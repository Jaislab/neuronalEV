# Independent review follow-up

The review found robustness and maintenance gaps without identifying a corrupted numerical result. The following changes preserve the existing quantitative model and add explicit failures and regression coverage.

## Addressed findings

| Finding | Result |
| --- | --- |
| Empty vendor TEXT values could merge optional keywords | Parsed keywords containing the TEXT delimiter are now rejected before normalization. Empty `CELLS`, `PnV` and `SPILLOVER` fixtures followed by `CYT` fail explicitly with two delimiter styles. Escaped delimiters inside ordinary values remain supported. This is strict rejection of ambiguous input, not a guessed vendor-file repair. |
| Scatter settings absent from acquisition signature | Signature comparison now includes FITC/PE/APC and every axis used by bead, singlet or other pre-gates. Valid FCS fixtures with changed FSC-A or SSC-A voltages fail setup matching. |
| Unknown control types accepted | Readiness checks every control against the documented enumeration, including unused controls. A typo such as `assay blank` produces an explicit error. |
| Inconsistent donor IDs discovered late | Same-sample acquisitions mapped to different donor IDs fail metadata readiness before event processing or output creation. |
| Disabled normalization could defer a bad floor until processing | Every supplied floor is validated during configuration loading. Disabled normalization still accepts `null`; strings, booleans, negative and nonfinite floors fail early. |
| Calibration could create output before validating references | All setup references, compensation policies and FITC/PE transforms are checked before creating output. Missing/malformed references and unknown compensation leave no partial directory. Unexpected errors after writing starts retain `failure.json`. |
| Readability and unused imports | Removed unused demo imports, replaced the JSON-constant lambda with a named function, expanded long gating dictionaries and derived the inspect isotype count from panel roles. Gating formatting preserved the AST. |
| Documentation and package metadata | README commands use the repository root. Distribution metadata includes the Markdown README. A pinned development linter and GitHub Actions workflow were added. |
| Missing end-to-end coverage | CLI tests now exercise a technical replicate through analyze and report, plus demo/report success. Three acquisitions from two donors retain 117 measurement rows and distinct matrix rows; no technical-replicate averaging is introduced. |

## Verification

- **119 tests passed locally on Python 3.12**, with no failures, errors or skips.
- Ruff checks passed using the pinned development version.
- A fresh 24-donor synthetic CLI run completed with 936 measurements. All raw, background-corrected and normalized matrices were exactly equal to the pre-review synthetic baseline.
- A wheel built successfully; its Markdown README metadata and annotation resources were checked.
- CI is configured for Python 3.12 and 3.13. These local results do not establish the outcome of a remote Actions run.

The biological annotation work is documented separately in [marker-annotations.md](marker-annotations.md). It covers all 37 capture markers and two controls and preserves the quantitative values; its groups do not establish EV cellular origin or lineage fractions.

## Deliberately left open

- **Project license:** none has been chosen. A license should reflect the repository owner's decision rather than an assumed permission grant.
- **Repeated gating computation:** computation may still repeat for summaries, identity audits and plots. At the planned event counts this is a performance opportunity, not a correctness defect; no cache/refactor was introduced during this hardening pass.
- **Real FACSDiva compatibility and assay validation:** no experimental FCS files have been provided. In particular, ambiguous empty TEXT values are rejected explicitly; support for a concrete vendor variant should be added and validated against an actual export rather than inferred from synthetic fixtures.
