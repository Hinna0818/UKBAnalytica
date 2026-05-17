---
name: skill-creator
description: Create, update, review, and validate Codex/OpenAI agent skills in this repository. Use when the user asks to add a new skill, revise an existing SKILL.md, design skill references/scripts/assets, audit trigger descriptions, add smoke tests, package skills under inst/skills, or align skills with UKBAnalytica/RAP privacy and reproducibility standards.
---

# Skill Creator

## Purpose

Use this skill to create and maintain repo-bundled Codex skills. Keep every
skill concise, triggerable, testable, and safe for UK Biobank RAP workflows.

## Workflow

1. Inspect the existing skill tree before editing:

```bash
find inst/skills -maxdepth 3 -type f | sort
```

2. Decide the target:

- New standalone skill: `inst/skills/<skill-name>/SKILL.md`.
- New UKBAnalytica workflow skill: `inst/skills/UKBAnalytica_skills/<ukbsci-name>/SKILL.md`.
- Existing skill update: patch only the relevant skill and directly referenced resources.

3. Write `SKILL.md` first. The frontmatter must contain only:

```yaml
---
name: lower-hyphen-name
description: Clear trigger description with what the skill does and when to use it.
---
```

4. Keep the body procedural and short. Put detailed examples, disease catalogs,
function lists, or long checklists in one-hop `references/` files.

5. Add resources only when they materially improve reuse:

- `references/`: detailed guidance loaded only when needed.
- `scripts/`: deterministic checks or generation utilities.
- `assets/`: templates or media used in outputs.
- `agents/openai.yaml`: UI metadata. Its `default_prompt` must mention `$<skill-name>`.

6. Validate after edits:

```bash
python3 inst/skills/skill-creator/scripts/validate_skill.py inst/skills/<skill-name>
```

For all repo skills:

```bash
python3 inst/skills/skill-creator/scripts/validate_skill.py inst/skills
```

## UKBAnalytica Standards

- Preserve RAP privacy boundaries: never instruct users to export UK Biobank
  individual-level raw data, `eid`, exact dates, raw RAP fields, or
  re-identifiable row-level tables.
- Allow de-identified figures, aggregate tables, model metrics, coefficients,
  feature-level tables, and bin/strata summaries when no identifying fields
  accompany them.
- Prefer package functions over ad hoc code in examples.
- Remind users to install or load the latest UKBAnalytica before high-level
  workflow examples.
- Add minimal smoke tests when a skill depends on package functions whose
  interfaces may change.

## Reference Files

Read only the needed reference:

- `references/skill-structure.md`: SKILL.md, resource layout, and UI metadata.
- `references/ukbanalytica-skill-standards.md`: UKBAnalytica-specific privacy,
  function, and analysis workflow standards.
- `references/evaluation-and-smoke-tests.md`: validation and smoke-test design.

## Completion Checklist

- `SKILL.md` has valid frontmatter and no placeholder text.
- `description` contains the trigger conditions.
- Long details are moved into direct `references/` files.
- `agents/openai.yaml`, if present, matches the skill and mentions `$skill-name`.
- Validation script passes.
- For package-facing skills, a small smoke test exists or the reason for
  omitting one is clear.
