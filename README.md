# neuronalEV — MACSPlex EV Neuro FCS analysis

A reproducible command-line pipeline for the planned **24 independent human donors**, with **ExoQuick-enriched plasma**, **BD LSRFortessa** acquisition, and combined **CD9/CD63/CD81** detection.

FITC/PE identify 39 capture-bead populations. APC supplies one median signal per population: **37 marker measurements and two isotype controls per acquisition**. Bead counts are measurement counts, not additional biological replicates. No patient FCS files have been supplied; real-file compatibility, gate identity and assay performance remain unvalidated.

The input types are sample FCS files, a matched assay-blank FCS and the kit's setup-bead FCS. Gate definitions are derived and checked from that setup file; no separate provider gate/workspace file is required. See the [Figure 2/Table 1 identity guide](docs/bead-identification.md) for the complete 39-population layout and interpretation of 5,000 total beads per sample. Setup reference processing checks observed FITC/PE ordering, catching inconsistent placement before sample APC summaries.

The [biological annotation guide](docs/marker-annotations.md) now defines neuron-associated, glia/progenitor-associated and shared neural/progenitor views, alongside immune/myeloid, vascular/stromal, broad receptor, lysosomal and tetraspanin views. These organize individual marker signals with overlapping tags and source-specific caveats; they do not assign EV cell of origin or calculate lineage fractions.

## Use the installed environment

From the repository root:

```bash
.venv/bin/facs-neuro --help
.venv/bin/facs-neuro demo --output results/synthetic-demo
```

`demo` generates 24 explicitly synthetic FCS3.1 files, a setup reference and assay blank, analyzes all 39 populations, and creates `results/synthetic-demo/results/report.html`. Its known gates and thresholds apply **only to the generated data**. The standalone HTML embeds heatmaps and gating diagnostics. Keep its run folder for the linked numerical CSV files. Output directories must be new; nothing is silently overwritten.

## Install or reproduce

Python 3.11+ is required; the current environment was exercised with Python 3.12. Run:

```bash
bash scripts/setup.sh
```

The script reuses the current runtime's exact locked versions when available, installs the lock with normal TLS verification when needed, and installs the local CLI without network build isolation. A standard independent installation is:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.lock
.venv/bin/python -m pip install --no-deps --no-build-isolation -e .
```

No GPU, database, web server or live background process is needed. The code uses NumPy, pandas and Matplotlib. The strict FCS reader and explicit polygon gates in this first version do not require FlowKit/FlowIO downloads; see supported formats below.

## Prepare the real experiment

1. Put original event-level FCS files under ignored `data/raw/` or another retained data directory. Keep originals unchanged.
2. Fill `metadata/samples.csv`. The 24 donor/sample IDs and known preparation metadata are already present; no group assignments or file paths are invented. File paths are relative to the CSV containing them, not the shell working directory.
3. Fill `metadata/controls.csv` with actual setup-bead and matched no-EV assay-blank exports. A blank includes capture beads and detection cocktail processed without EV input. Isotype beads and setup beads cannot substitute for this blank. Budget **more than 24 assay reactions** for 24 donors plus blanks and any replicates.
4. Edit `config/analysis.json` for the actual detector names and acquisition setup. Its channel names and asinh cofactor are placeholders requiring review. Document compensation input state and policy before plotting or analysis.
5. Inspect files and metadata, then generate reference overviews:

```bash
.venv/bin/facs-neuro inspect --config config/analysis.json --output results/inventory.json
.venv/bin/facs-neuro calibrate --config config/analysis.json --output results/gate-review
```

Create the parent `results/` directory before writing an inspection JSON file if needed. `inspect` works on the incomplete template and reports what is missing; exit code 2 is expected until inputs are complete. It does not certify gates. `calibrate` supplies unreviewed overview PNGs and **empty** population polygons. It does not invent gates, align clusters, or assign identities automatically.

6. Establish scatter bead and singlet gates, then all 39 FITC/PE population polygons against setup beads and the manufacturer identity map. `scripts/identify_setup_beads.py` derives an unreviewed template from a setup-bead export (see below). Confirm the polygons with an expert-reviewed pilot including low/high signals and challenging populations. Configure a justified minimum bead count and optional normalization floor. Only then set `reviewed: true`, a nonempty `gate_version`, and `review_note`. See [configuration details](docs/configuration.md).
7. Run the reviewed configuration:

```bash
.venv/bin/facs-neuro analyze --config config/analysis.json --output results/run-001
.venv/bin/facs-neuro report --run results/run-001
```

The second command refreshes the HTML for a completed run. Failed runs retain `failure.json` and incomplete provenance; they do not produce a completed report. Analysis exit code 3 means processing completed but measurement QC failed; review the report/CSV before interpreting values. Exit code 0 means software processing completed, not biological or assay validation. Exit code 2 indicates invalid/incomplete input or a processing error.

## Derive population gates from a setup-bead file

```bash
.venv/bin/python scripts/identify_setup_beads.py data/raw/20261029_Beads_5000.fcs \
    --setup-id fortessa_20261029 --output results/identify-20261029
```

The script gates singlet beads on scatter, finds the density peaks of FITC versus PE in the chosen gating transform, matches them to the manufacturer's Figure 2 staircase and names them from Table 1. It writes `setup-template.json` (a setup object with `reviewed: false`), `populations.csv` with observed medians and counts, `identified.png` and `summary.json`, and exits with code 2 when the 39 populations are not resolved or their Figure 2 order is violated. It never reads APC for identity and never marks anything reviewed.

The supplied LSRFortessa setup-bead export of 29 October 2026 resolved all 39 populations with consistent Figure 2 order; `config/analysis.json` carries that unreviewed template as setup `fortessa_20261029`, with the evidence in [`docs/setup-bead-identification-20261029.md`](docs/setup-bead-identification-20261029.md) and `config/gates/`. Review the overlay, set `min_events_per_population`, `gate_version` and `review_note`, then set `reviewed: true`.

## Outputs and interpretation

- `measurements_long.csv`: per-acquisition/per-bead counts, raw median APC, matched blank median, signed correction, normalization and explicit QC states.
- `raw_median_apc.csv`, `background_corrected.csv`, `normalized.csv`: 37-marker matrices, indexed by acquisition ID. Raw medians are retained for audit even when QC-invalid; consult validity flags. Invalid corrected/normalized measurements are empty cells, never biological zeros.
- `control_measurements.csv`, `isotype_controls.csv`, `qc.csv`, `qc.json`: measurement/control evidence, counts through pre-gates, missing/ambiguous populations and saturation.
- `bead-identity-reference.png`, `bead-identity-reference.csv`, `bead_identity_audit.csv`: schematic Figure 2/Table 1 guide and observed setup-bead ordering/count evidence. The reference contains ordinal positions, not fluorescence coordinates; consistent order does not itself certify correct identity.
- `sample_metadata.csv`, original input manifests/config, `provenance.json`: donor/acquisition mapping, hashes, versions and processing status.
- `marker_annotations.csv`, `marker_annotation_sources.json`, `marker_views/`: versioned interpretive annotations and grouped individual-marker matrices, with annotation hashes recorded in provenance. Groups do not alter gating, numerical correction, normalization or detection rules.
- `report.html`, `diagnostics/`: self-contained descriptive report and individual diagnostic PNGs.

In the default assay-blank mode, correction is **sample median APC minus its matched bead-specific blank median** on the linear analysis scale. Negative results are retained and labeled `not_detected`; positive corrected values are labeled `unclassified_positive`, not automatically detected. No validated LOD/LOQ is inferred from one blank. Isotypes diagnose nonspecific signal and are not a second automatic subtraction.

An explicitly selected `background.mode: "setup_baseline"` subtracts a setup-bead reference when a valid assay blank is unavailable. It records the rationale, reference and unmeasured detection-reagent background. Its states are below/at/above setup baseline, with no detection or absence claim; the same limitation remains after optional normalization. Known sample failures can be removed through `sample_exclusions` with their reasons retained in `excluded_samples.csv` and provenance. Original manifests remain unchanged, and excluded samples receive no invented zero measurements. See [configuration details](docs/configuration.md).

Optional normalization divides by the arithmetic mean of corrected CD9/CD63/CD81 medians from a fixed declared set. Every chosen normalizer must be valid and positive and its mean above the configured floor. No per-sample fallback is used. ExoQuick does not deliberately affinity-block those markers, but it can co-enrich other plasma material. Normalization is not an EV-purity correction; signals do not establish neuronal origin, absolute EV counts or single-EV co-expression.

The current release keeps technical acquisitions separate and produces descriptive results. It does **not** combine technical replicates automatically, test donor groups, implement batch correction or run PCA. Group comparisons require the still-unspecified scientific contrast and design. High assay-blank/isotype signal and acquisition drift require pilot-defined review; no manufacturer-unsupported numerical thresholds are invented here.

## Supported FCS and gating scope

- FCS3.0/3.1, single-dataset list mode, linear float32/float64 or byte-aligned unsigned integer parameters, common endianness, unique `$PnN` names.
- Nonzero `$PnE`, nonunit signal/scatter `$PnG`, supplemental TEXT, multiple datasets, malformed offsets and nonfinite event values are explicitly rejected. A narrow FACSDiva floating-point `Time` exception permits a positive gain matching `$TIMESTEP`; native Time ticks are preserved, and seconds equal native Time multiplied once by `$TIMESTEP`. No fluorescence values are rescaled by that exception. Request a compatible export or extend/validate unsupported variants separately; do not strip metadata to bypass a rejection. Real FACSDiva files still require comparison against acquisition software.
- Explicit compensation policies: `none` with rationale, `embedded` matrix, or `explicit` matrix. Already-compensated inputs cannot receive another matrix. An embedded matrix does not establish whether compensation was already applied.
- Reviewed polygon gates in linear or asinh coordinates. Pre-gates use conventional Fortessa `FSC-*`, `SSC-*`, and `Time` names, with actual area-vs-height/width singlet geometry. No gating on APC. Every manufacturer bead ID must have one population definition; overlapping events are excluded and affected populations flagged.
- Saturation is checked on original APC before computational compensation using a conservative `$PnR - 1` boundary or a documented explicit threshold. Already-compensated inputs need a `saturation_review_note` attesting review of original acquisition QC; otherwise measurements are invalid because original clipping cannot be established.

The initial [analysis plan](docs/analysis-plan.md) documents the scientific rationale. The panel identity map was transcribed from the uploaded Miltenyi manual, document **140-006-716.02**, Table 1 (printed p.10). Acquisition guidance is referenced at p.19; correction/normalization at pp.20–21. Obtain that separate acquisition application note before accepting real gates.

## Verification

```bash
MPLCONFIGDIR=/tmp/facs-neuro-mpl XDG_CACHE_HOME=/tmp/facs-neuro-cache .venv/bin/python -m unittest discover -s tests -v
```

The tests exercise FCS encoding/offsets/scaling rejection, compensation, gated population identity and ambiguity, signed correction, blank matching, normalization failures, donor/acquisition identity, report safety, and the FCS-to-CSV/report workflow. The synthetic demo verifies software behavior; it does not validate real sample preparation, detector settings, reference gates or biological findings.

Development linting and CI use `requirements-dev.lock`:

```bash
.venv/bin/python -m pip install -r requirements-dev.lock
.venv/bin/ruff check .
```

GitHub Actions is configured to run lint, the test suite, a synthetic CLI smoke test and a wheel build on Python 3.12 and 3.13. A local pass does not establish the result of a remote Actions run. No project license has been selected; adding one requires the repository owner's choice.
