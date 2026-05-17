# Skill Structure Reference

## Required Layout

Every skill is a directory with a required `SKILL.md` file:

```text
skill-name/
├── SKILL.md
├── agents/openai.yaml      # optional UI metadata
├── references/             # optional long guidance
├── scripts/                # optional deterministic utilities
└── assets/                 # optional templates/media
```

## SKILL.md Rules

- Use lowercase hyphenated names.
- Folder name and frontmatter `name` should match.
- Frontmatter should include only `name` and `description`.
- `description` is the trigger surface. Include what the skill does and when to
  use it there, not only in the body.
- Keep instructions imperative and action-oriented.
- Keep body content concise; move large examples or catalogs into references.
- Avoid README, changelog, installation guide, or other auxiliary files unless
  they are directly consumed by the skill.

## Progressive Disclosure

Use three levels:

1. Metadata: always available through `name` and `description`.
2. `SKILL.md`: loaded when the skill triggers.
3. References/scripts/assets: loaded or executed only when needed.

Prefer one-hop references directly linked from `SKILL.md`. Avoid deep chains of
references that make future agents search blindly.

## agents/openai.yaml

Use this only for UI metadata. Keep it short:

```yaml
interface:
  display_name: "Human Name"
  short_description: "Short scanning description."
  default_prompt: "Use $skill-name to ..."
```

The default prompt should mention `$skill-name`.
