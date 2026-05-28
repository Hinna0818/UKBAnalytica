#!/usr/bin/env python3
"""Verify PMID/DOI references in a Markdown research brief.

This script intentionally uses only Python's standard library so it can run in
minimal RAP or local environments. It extracts explicit PMID and DOI mentions
from a Markdown file, optionally validates them through PubMed E-utilities and
Crossref, and writes a machine-readable JSON audit.

Examples
--------
Parse only, without network:
    python3 verify_references.py report.md --offline

Validate with PubMed/Crossref:
    python3 verify_references.py report.md --output reference_audit.json

Validate and export lightweight reference-manager files:
    python3 verify_references.py report.md --bibtex references.bib --ris references.ris
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable


PMID_RE = re.compile(r"\bPMID\s*[:=]?\s*(\d{6,9})\b", re.IGNORECASE)
DOI_RE = re.compile(r"\b(10\.\d{4,9}/[^\s\]\)>,;]+)", re.IGNORECASE)
URL_RE = re.compile(r"https?://[^\s\]\)>]+", re.IGNORECASE)


@dataclass
class ReferenceRecord:
    kind: str
    identifier: str
    status: str
    title: str | None = None
    year: str | None = None
    source: str | None = None
    url: str | None = None
    error: str | None = None


def unique_keep_order(values: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for value in values:
        cleaned = value.strip().rstrip(".,;:)")
        if cleaned and cleaned.lower() not in seen:
            seen.add(cleaned.lower())
            out.append(cleaned)
    return out


def extract_references(text: str) -> tuple[list[str], list[str], list[str]]:
    pmids = unique_keep_order(PMID_RE.findall(text))
    dois = unique_keep_order(DOI_RE.findall(text))
    urls = unique_keep_order(URL_RE.findall(text))
    return pmids, dois, urls


def fetch_json(url: str, timeout: int = 20) -> dict:
    req = urllib.request.Request(
        url,
        headers={
            "User-Agent": "UKBAnalytica-reference-verifier/1.0 "
            "(mailto:maintainer@example.com)"
        },
    )
    with urllib.request.urlopen(req, timeout=timeout) as response:
        payload = response.read().decode("utf-8")
    return json.loads(payload)


def verify_pmid(pmid: str, timeout: int) -> ReferenceRecord:
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?"
        + urllib.parse.urlencode({"db": "pubmed", "id": pmid, "retmode": "json"})
    )
    try:
        data = fetch_json(url, timeout=timeout)
        result = data.get("result", {})
        item = result.get(pmid)
        if not item:
            return ReferenceRecord("PMID", pmid, "not_found", url=url)
        return ReferenceRecord(
            kind="PMID",
            identifier=pmid,
            status="verified",
            title=item.get("title"),
            year=(item.get("pubdate") or "")[:4] or None,
            source=item.get("source"),
            url=f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
        )
    except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
        return ReferenceRecord("PMID", pmid, "unverified_network_error", url=url, error=str(exc))


def verify_doi(doi: str, timeout: int) -> ReferenceRecord:
    clean = doi.strip().rstrip(".,;:)")
    url = "https://api.crossref.org/works/" + urllib.parse.quote(clean, safe="")
    try:
        data = fetch_json(url, timeout=timeout)
        msg = data.get("message", {})
        titles = msg.get("title") or []
        years = (
            msg.get("published-print", {}).get("date-parts")
            or msg.get("published-online", {}).get("date-parts")
            or msg.get("issued", {}).get("date-parts")
            or []
        )
        year = str(years[0][0]) if years and years[0] else None
        return ReferenceRecord(
            kind="DOI",
            identifier=clean,
            status="verified",
            title=titles[0] if titles else None,
            year=year,
            source=(msg.get("container-title") or [None])[0],
            url="https://doi.org/" + clean,
        )
    except urllib.error.HTTPError as exc:
        status = "not_found" if exc.code == 404 else "unverified_http_error"
        return ReferenceRecord("DOI", clean, status, url=url, error=str(exc))
    except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
        return ReferenceRecord("DOI", clean, "unverified_network_error", url=url, error=str(exc))


def audit_markdown(path: Path, offline: bool, timeout: int, sleep: float) -> dict:
    text = path.read_text(encoding="utf-8")
    pmids, dois, urls = extract_references(text)

    records: list[ReferenceRecord] = []
    if offline:
        records.extend(ReferenceRecord("PMID", x, "not_checked_offline", url=f"https://pubmed.ncbi.nlm.nih.gov/{x}/") for x in pmids)
        records.extend(ReferenceRecord("DOI", x, "not_checked_offline", url=f"https://doi.org/{x}") for x in dois)
    else:
        for pmid in pmids:
            records.append(verify_pmid(pmid, timeout=timeout))
            time.sleep(sleep)
        for doi in dois:
            records.append(verify_doi(doi, timeout=timeout))
            time.sleep(sleep)

    records.extend(ReferenceRecord("URL", x, "extracted_only", url=x) for x in urls)
    return {
        "input": str(path),
        "offline": offline,
        "counts": {
            "pmid": len(pmids),
            "doi": len(dois),
            "url": len(urls),
            "records": len(records),
        },
        "records": [asdict(record) for record in records],
    }


def clean_key(record: dict, index: int) -> str:
    base = record.get("identifier") or f"ref{index}"
    base = re.sub(r"^https?://", "", base, flags=re.IGNORECASE)
    base = re.sub(r"[^A-Za-z0-9]+", "_", base).strip("_")
    if not base:
        base = f"ref{index}"
    return f"ref_{base[:48]}_{index}"


def bibtex_escape(value: str | None) -> str:
    if not value:
        return ""
    return value.replace("\\", "\\\\").replace("{", "\\{").replace("}", "\\}")


def records_to_bibtex(records: list[dict]) -> str:
    entries: list[str] = []
    for i, record in enumerate(records, start=1):
        if record.get("kind") not in {"PMID", "DOI"}:
            continue
        key = clean_key(record, i)
        title = bibtex_escape(record.get("title") or record.get("identifier"))
        year = bibtex_escape(record.get("year"))
        journal = bibtex_escape(record.get("source"))
        doi = record.get("identifier") if record.get("kind") == "DOI" else ""
        url = record.get("url") or ""
        fields = [
            f"  title = {{{title}}}",
        ]
        if journal:
            fields.append(f"  journal = {{{journal}}}")
        if year:
            fields.append(f"  year = {{{year}}}")
        if doi:
            fields.append(f"  doi = {{{bibtex_escape(doi)}}}")
        if url:
            fields.append(f"  url = {{{bibtex_escape(url)}}}")
        fields.append(f"  note = {{{record.get('kind')} {bibtex_escape(record.get('identifier'))}; status: {record.get('status')}}}")
        entries.append("@article{" + key + ",\n" + ",\n".join(fields) + "\n}")
    return "\n\n".join(entries) + ("\n" if entries else "")


def records_to_ris(records: list[dict]) -> str:
    lines: list[str] = []
    for record in records:
        if record.get("kind") not in {"PMID", "DOI"}:
            continue
        lines.append("TY  - JOUR")
        if record.get("title"):
            lines.append(f"TI  - {record['title']}")
        if record.get("source"):
            lines.append(f"JO  - {record['source']}")
        if record.get("year"):
            lines.append(f"PY  - {record['year']}")
        if record.get("kind") == "DOI":
            lines.append(f"DO  - {record['identifier']}")
        if record.get("kind") == "PMID":
            lines.append(f"AN  - PMID:{record['identifier']}")
        if record.get("url"):
            lines.append(f"UR  - {record['url']}")
        lines.append(f"N1  - verification_status: {record.get('status')}")
        lines.append("ER  -")
        lines.append("")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("markdown", type=Path, help="Markdown file to audit")
    parser.add_argument("--output", type=Path, default=None, help="JSON output path")
    parser.add_argument("--offline", action="store_true", help="Only parse references; do not call PubMed/Crossref")
    parser.add_argument("--timeout", type=int, default=20, help="HTTP timeout in seconds")
    parser.add_argument("--sleep", type=float, default=0.34, help="Delay between API calls")
    parser.add_argument("--bibtex", type=Path, default=None, help="Optional BibTeX output path")
    parser.add_argument("--ris", type=Path, default=None, help="Optional RIS output path")
    args = parser.parse_args(argv)

    if not args.markdown.exists():
        parser.error(f"Input file not found: {args.markdown}")

    audit = audit_markdown(args.markdown, offline=args.offline, timeout=args.timeout, sleep=args.sleep)
    output = json.dumps(audit, ensure_ascii=False, indent=2)
    if args.output:
        args.output.write_text(output + "\n", encoding="utf-8")
    else:
        print(output)

    records = audit["records"]
    if args.bibtex:
        args.bibtex.write_text(records_to_bibtex(records), encoding="utf-8")
    if args.ris:
        args.ris.write_text(records_to_ris(records), encoding="utf-8")

    bad = [
        rec for rec in audit["records"]
        if rec["kind"] in {"PMID", "DOI"} and rec["status"] in {"not_found"}
    ]
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
