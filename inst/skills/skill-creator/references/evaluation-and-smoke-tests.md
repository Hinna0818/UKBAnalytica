# Evaluation and Smoke Tests

## Structural Validation

Run the bundled validator after creating or editing a skill:

```bash
python3 inst/skills/skill-creator/scripts/validate_skill.py inst/skills/<skill-name>
```

For all bundled skills:

```bash
python3 inst/skills/skill-creator/scripts/validate_skill.py inst/skills
```

## Minimal Smoke Tests

Add a smoke test when a skill demonstrates code that may drift with package
interfaces. A good smoke test should:

- use simulated or toy data only;
- avoid UK Biobank individual-level data entirely, including de-identified
  row-level examples;
- run quickly;
- assert object classes, required columns, and that exported outputs are
  aggregate-only with no `eid`, exact dates, row indices, raw RAP fields, or
  row-level predictions / SHAP values;
- fail loudly when an interface changes.

## Suggested Locations

- Skill-specific checks: `<skill>/evals/`.
- Shared UKBAnalytica checks:
  `inst/skills/UKBAnalytica_skills/evals/smoke_test_core.R`.

## Forward Testing

For complex skills, test with a realistic prompt that asks an agent to use the
skill, not to review it. Pass only the skill path and task request, not the
expected answer.
