# jenner-check

Automated pipeline that selects real-world SAS repositories, generates
Jenner compatibility tests from their code, runs the tests ourselves to
self-validate, and submits a pull request to each repo with the test
bundle and a SAS-side runner that the repo owner can execute to confirm
(or contradict) our compatibility claims.

## Sub-pipelines

1. **Prevalence** (already in `sas_catalog/`) — discover SAS repos on
   GitHub, fetch per-repo metadata, decide which are PR-eligible.
2. **Runner** (`jenner_check/runner/`) — the `run_jenner_check.sas` master
   script we ship into every PR, plus reference test bundles that exercise
   its contract. See `runner/jenner-check/README.md` for the
   owner-facing view.
3. **Test generator** (planned) — given a `(repo, script_path)` pick,
   produce a `jenner-check/tNNN_*/` bundle with `script.sas`,
   `input/`, `expected/`, `validate.sas`, `meta.json`.
4. **Pre-PR self-check** (planned) — run the generated bundle through
   Jenner locally. If pass rate or numerical fidelity isn't good enough
   we DO NOT open a PR. If we would be embarrassed we stay quiet.
5. **PR service** (planned, adapted from `mseep/github_pr/`) — fork →
   branch → `add_test_directory()` → create PR with one body containing
   the intro, "how to run," and optional Jenner-Compatible badge
   markdown. One PR per repo, multiple test bundles per PR.
6. **PR status checker** (planned, from `mseep/pr_status_checker/`) —
   poll GH for merged/closed, scrape PR comments for attached
   `jenner_check_report.csv`, ingest into `jenner_check_result`.
7. **Diff-report ingestion** (planned) — endpoint + scrape path that
   records per-test outcomes against a `(pull_request_id, test_name)` key.

## Guiding principles (from discussion 2026-04-17)

- **One PR per repo, not per test.** Owners tolerate one PR before they
  get angry. Bundle everything we want to test into a single PR.
- **Validators beat listing diffs.** `validate.sas` checks numeric
  equivalence with tolerances. Listing diffs are for fallback.
- **Badge is static, not dynamic.** A single SVG at
  `jenneranalytics.com/badges/jenner-compatible.svg` linking to the
  homepage. The quality gate is implicit — we only PR if we'd stand
  behind the claim.
- **Cron, not Prefect.** Every phase's state lives in the DB with a
  resume predicate (`WHERE …_at IS NULL`). Row-level `FOR UPDATE SKIP
  LOCKED` on claims; stuck-row reset job every 15 min as a safety net.
- **do_not_contact is sacred.** Any owner reply containing
  `no-more-prs`, or an issue titled `jenner-check: opt out`, flips the
  flag — no future automation touches that repo.

## Not in scope

- Running Jenner for the owner. They have SAS; we ship SAS they can run.
- Dynamic per-repo compatibility pages — the repo's own history of
  merged/unmerged PRs + attached reports is the record.
