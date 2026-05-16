# ukbsci-workflow

Agent skill that **plans and orchestrates** an end-to-end UK Biobank RAP study
using the `UKBAnalytica` R package. It does not do the math itself — it owns
the project plan, the directory layout, and the routing of each phase to the
right specialist skill.

## What this skill does

1. Writes the Phase 0 plan document and pauses for user approval.
2. Routes each subsequent phase to the correct `ukbsci-*` specialist skill.
3. Enforces a canonical project layout under `/mnt/project/<area>/`.
4. Pauses at known decision points (eight in total) so the user can pivot.
5. At the end, lists actual deliverables and writes the wrap-up README.

## What this skill does NOT do

- Run extraction, cohort building, regression, ML, or plotting itself — those
  are specialist skills.
- Re-derive UKB phenotypes, statistical methods, or plotting templates from
  scratch.
- Permit any deviation from RAP-resident storage.

## Compatibility

Same as the rest of the pack: provider-agnostic Markdown + YAML
(`name`/`description`). Works with Claude Code skills, OpenAI Assistants
instructions, and generic tool-using agents.

## File map

```
ukbsci-workflow/
├── SKILL.md
├── README.md
├── references/
│   ├── plan-template.md         ← Phase 0 fillable skeleton
│   ├── wrap-up-template.md      ← Phase 8 fillable skeleton
│   └── handoff-checklist.md     ← what to pass to each sub-skill
└── evals/
    └── evals.json
```
