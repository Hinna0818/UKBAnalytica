#!/usr/bin/env python3
"""Validate one or more Codex skill folders.

This is intentionally lightweight: it catches structural mistakes before
running heavier package checks.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


FRONTMATTER_RE = re.compile(r"^---\n(.*?)\n---\n", re.DOTALL)
NAME_RE = re.compile(r"^[a-z0-9-]{1,64}$")


def parse_frontmatter(text: str) -> dict[str, str]:
    match = FRONTMATTER_RE.match(text)
    if not match:
        raise ValueError("missing YAML frontmatter delimited by ---")
    data: dict[str, str] = {}
    current_key: str | None = None
    for raw_line in match.group(1).splitlines():
        if raw_line.startswith((" ", "\t")) and current_key:
            data[current_key] = (data[current_key] + " " + raw_line.strip()).strip()
            continue
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if ":" not in line:
            raise ValueError(f"invalid frontmatter line: {raw_line!r}")
        key, value = line.split(":", 1)
        current_key = key.strip()
        value = value.strip().strip('"').strip("'")
        data[current_key] = "" if value == ">" else value
    return data


def validate_skill_dir(path: Path) -> list[str]:
    errors: list[str] = []
    skill_file = path / "SKILL.md"
    if not skill_file.exists():
        return errors

    text = skill_file.read_text(encoding="utf-8")
    try:
        meta = parse_frontmatter(text)
    except ValueError as exc:
        return [f"{skill_file}: {exc}"]

    allowed = {"name", "description"}
    extra = set(meta) - allowed
    if extra:
        errors.append(f"{skill_file}: unexpected frontmatter keys: {sorted(extra)}")

    name = meta.get("name", "")
    description = meta.get("description", "")
    if not name:
        errors.append(f"{skill_file}: missing frontmatter name")
    elif not NAME_RE.match(name):
        errors.append(f"{skill_file}: invalid skill name {name!r}")
    elif path.name != name:
        errors.append(f"{skill_file}: folder name {path.name!r} does not match name {name!r}")

    if not description or len(description) < 40:
        errors.append(f"{skill_file}: description is missing or too short")
    if "TODO" in text or "[TODO" in text:
        errors.append(f"{skill_file}: contains TODO placeholder text")

    openai_yaml = path / "agents" / "openai.yaml"
    if openai_yaml.exists() and name:
        ui = openai_yaml.read_text(encoding="utf-8")
        if f"${name}" not in ui:
            errors.append(f"{openai_yaml}: default_prompt should mention ${name}")

    return errors


def find_skill_dirs(root: Path) -> list[Path]:
    if (root / "SKILL.md").exists():
        return [root]
    return sorted(p.parent for p in root.rglob("SKILL.md"))


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate Codex skill folders.")
    parser.add_argument("path", type=Path, help="Skill folder or root containing skills")
    args = parser.parse_args()

    root = args.path.resolve()
    if not root.exists():
        print(f"Path not found: {root}", file=sys.stderr)
        return 2

    skill_dirs = find_skill_dirs(root)
    if not skill_dirs:
        print(f"No SKILL.md files found under {root}", file=sys.stderr)
        return 2

    all_errors: list[str] = []
    for skill_dir in skill_dirs:
        all_errors.extend(validate_skill_dir(skill_dir))

    if all_errors:
        print("Skill validation failed:", file=sys.stderr)
        for err in all_errors:
            print(f"- {err}", file=sys.stderr)
        return 1

    print(f"Validated {len(skill_dirs)} skill(s).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
