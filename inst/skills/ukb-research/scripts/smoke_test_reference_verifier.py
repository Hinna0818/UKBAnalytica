#!/usr/bin/env python3
"""Offline smoke test for the UKB research reference verifier."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path


def main() -> int:
    script = Path(__file__).with_name("verify_references.py")
    sample = (
        "# Mock report\n\n"
        "A cited PubMed record (PMID: 36138150) and DOI "
        "10.1038/s41591-022-01980-3 support this sentence. "
        "A source URL is https://www.ukbiobank.ac.uk/.\n"
    )

    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        md = tmpdir / "mock_report.md"
        audit_json = tmpdir / "audit.json"
        bib = tmpdir / "refs.bib"
        ris = tmpdir / "refs.ris"
        md.write_text(sample, encoding="utf-8")

        cmd = [
            sys.executable,
            str(script),
            str(md),
            "--offline",
            "--output",
            str(audit_json),
            "--bibtex",
            str(bib),
            "--ris",
            str(ris),
        ]
        result = subprocess.run(cmd, check=False, text=True, capture_output=True)
        if result.returncode != 0:
            print(result.stdout)
            print(result.stderr, file=sys.stderr)
            return result.returncode

        audit = json.loads(audit_json.read_text(encoding="utf-8"))
        assert audit["counts"]["pmid"] == 1
        assert audit["counts"]["doi"] == 1
        assert audit["counts"]["url"] >= 1
        assert "36138150" in bib.read_text(encoding="utf-8")
        assert "10.1038/s41591-022-01980-3" in ris.read_text(encoding="utf-8")

    print("UKB research reference verifier smoke test passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
