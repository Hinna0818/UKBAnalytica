# Install / loading recipes

Copy-paste blocks for each agent runtime that might consume this skill pack.
The pack lives at `skills/UKBAnalytica_skills/` in the repository root.

---

## 1. Resolve the pack path (any runtime)

```r
# From R — using the repo root directly:
pack_dir <- file.path(getwd(), "skills/UKBAnalytica_skills")
list.files(pack_dir, pattern = "^ukbsci")
#> [1] "ukbsci-baseline" "ukbsci-cohort" "ukbsci-imputation" ...

# Or with here::here() for scripts that may run from a subdirectory:
pack_dir <- here::here("skills/UKBAnalytica_skills")
```

```bash
# From shell (repo root):
ls skills/UKBAnalytica_skills/ | grep '^ukbsci'
```

---

## 2. Claude Code skills

Claude Code auto-discovers skills under `~/.claude/skills/` (user-scope) or
`<repo>/.claude/skills/` (workspace-scope).

### Option A — symlink the pack directory (recommended)

```bash
# Run from the repo root
PACK_DIR="$(pwd)/skills/UKBAnalytica_skills"
mkdir -p ~/.claude/skills
ln -snf "$PACK_DIR" ~/.claude/skills/UKBAnalytica_skills
```

Claude Code will see each skill as `~/.claude/skills/UKBAnalytica_skills/ukbsci-*/SKILL.md`.

### Option B — workspace-scope (.claude/settings.json)

Add to `.claude/settings.json`:

```json
{
  "skillDirs": ["skills/UKBAnalytica_skills"]
}
```

Verify in a Claude Code session that triggers such as
"build UKB survival dataset" or "extract UKB fields on RAP" resolve to
`ukbsci-cohort` and `ukbsci-rap-extract`.

---

## 3. OpenAI Assistants / Responses API

Each `SKILL.md` is a self-contained system-prompt fragment. Two common patterns:

### Pattern A — fold all skill bodies into `instructions`

```python
from pathlib import Path
from openai import OpenAI

PACK = Path("skills/UKBAnalytica_skills")
skills = "\n\n---\n\n".join(
    (p / "SKILL.md").read_text()
    for p in sorted(PACK.iterdir())
    if (p / "SKILL.md").exists()
)

client = OpenAI()
asst = client.beta.assistants.create(
    name="UKBAnalytica RAP Agent",
    model="gpt-4o",
    instructions=(
        "You are a UK Biobank RAP analysis assistant. "
        "Follow the skills below. Never export UK Biobank participant-level "
        "raw data, direct identifiers, or row-level source tables that can "
        "be linked back to participants.\n\n" + skills
    ),
)
```

### Pattern B — attach `references/*.md` via file-search (scales better)

Keep only the short `SKILL.md` headers in `instructions`; upload all
`references/*.md` to a Vector Store so the agent retrieves details on demand.

```python
files = []
for skill_dir in PACK.iterdir():
    refs = (skill_dir / "references")
    if refs.exists():
        for ref in refs.glob("*.md"):
            files.append(client.files.create(
                file=open(ref, "rb"), purpose="assistants"))
vs = client.beta.vector_stores.create(
    name="UKBAnalytica refs",
    file_ids=[f.id for f in files])
```

---

## 4. LangChain / LlamaIndex / generic tool-using agents

Treat each skill directory as one tool:

```python
from pathlib import Path
import yaml

def load_skill(skill_dir: Path) -> dict:
    txt = (skill_dir / "SKILL.md").read_text()
    parts = txt.split("---\n", 2)
    meta = yaml.safe_load(parts[1])
    body = parts[2] if len(parts) > 2 else ""
    refs_dir = skill_dir / "references"
    return {
        "name": meta["name"],
        "description": meta["description"],
        "body": body,
        "references": {
            p.stem: p.read_text()
            for p in refs_dir.glob("*.md")
        } if refs_dir.exists() else {},
    }

PACK = Path("skills/UKBAnalytica_skills")
tools = [load_skill(d) for d in sorted(PACK.iterdir()) if (d / "SKILL.md").exists()]
```

Bind each element to an LLM as a "knowledge tool" whose description is
`tools[i]["description"]` and whose call returns `body + references`.

---

## 5. MANIFEST.json (machine-readable index)

The pack ships a manifest you can ingest directly without parsing YAML:

```python
import json
from pathlib import Path

m = json.load(open("skills/UKBAnalytica_skills/MANIFEST.json"))
for skill in m["skills"]:
    print(skill["name"], "->", skill["path"])
    # path is relative to skills/UKBAnalytica_skills/
```

```r
library(jsonlite)
m <- fromJSON("skills/UKBAnalytica_skills/MANIFEST.json")
print(m$skills[, c("name", "phase")])
```

---

## 6. Sanity check after install

```bash
# 14 skill directories present
ls skills/UKBAnalytica_skills/ | grep '^ukbsci' | wc -l   # expected: 14

# All SKILL.md have valid YAML frontmatter
python3 - <<'PY'
import re, pathlib
for sk in sorted(pathlib.Path("skills/UKBAnalytica_skills").iterdir()):
    sm = sk / "SKILL.md"
    if not sm.exists(): continue
    txt = sm.read_text()
    m = re.match(r"^---\n(.*?)\n---\n", txt, re.DOTALL)
    name = re.search(r"^name:\s*(.+)$", m.group(1), re.MULTILINE).group(1) if m else "MISSING"
    print(f"{sk.name:32} name={name}")
PY
```
