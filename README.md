# neuronalEV — MACSPlex EV Neuro FCS analysis

A reproducible command-line pipeline for the planned **24 independent human donors**, with **ExoQuick-enriched plasma**, **BD LSRFortessa** acquisition, and combined **CD9/CD63/CD81** detection.

FITC/PE identify 39 capture-bead populations. APC supplies one median signal per population: **37 marker measurements and two isotype controls per acquisition**. Bead counts are measurement counts, not additional biological replicates. No patient FCS files have been supplied; real-file compatibility, gate identity and assay performance remain unvalidated.

## Use the installed environment

From `/workspace/neuronalEV`:

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

6. Establish scatter bead and singlet gates, then all 39 FITC/PE population polygons against setup beads and the manufacturer identity map. Confirm them with an expert-reviewed pilot including low/high signals and challenging populations. Configure a justified minimum bead count and optional normalization floor. Only then set `reviewed: true`, a nonempty `gate_version`, and `review_note`. See [configuration details](docs/configuration.md).
7. Run the reviewed configuration:

```bash
.venv/bin/facs-neuro analyze --config config/analysis.json --output results/run-001
.venv/bin/facs-neuro report --run results/run-001
```

The second command refreshes the HTML for a completed run. Failed runs retain `failure.json` and incomplete provenance; they do not produce a completed report. Analysis exit code 3 means processing completed but measurement QC failed; review the report/CSV before interpreting values. Exit code 0 means software processing completed, not biological or assay validation. Exit code 2 indicates invalid/incomplete input or a processing error.

## Outputs and interpretation

- `measurements_long.csv`: per-acquisition/per-bead counts, raw median APC, matched blank median, signed correction, normalization and explicit QC states.
- `raw_median_apc.csv`, `background_corrected.csv`, `normalized.csv`: 37-marker matrices, indexed by acquisition ID. Raw medians are retained for audit even when QC-invalid; consult validity flags. Invalid corrected/normalized measurements are empty cells, never biological zeros.
- `control_measurements.csv`, `isotype_controls.csv`, `qc.csv`, `qc.json`: measurement/control evidence, counts through pre-gates, missing/ambiguous populations and saturation.
- `sample_metadata.csv`, original input manifests/config, `provenance.json`: donor/acquisition mapping, hashes, versions and processing status.
- `report.html`, `diagnostics/`: self-contained descriptive report and individual diagnostic PNGs.

Correction is **sample median APC minus its matched bead-specific blank median** on the linear analysis scale. Negative results are retained and labeled `not_detected`; positive corrected values are labeled `unclassified_positive`, not automatically detected. No validated LOD/LOQ is inferred from one blank. Isotypes diagnose nonspecific signal and are not a second automatic subtraction.

Optional normalization divides by the arithmetic mean of corrected CD9/CD63/CD81 medians from a fixed declared set. Every chosen normalizer must be valid and positive and its mean above the configured floor. No per-sample fallback is used. ExoQuick does not deliberately affinity-block those markers, but it can co-enrich other plasma material. Normalization is not an EV-purity correction; signals do not establish neuronal origin, absolute EV counts or single-EV co-expression.

The current release keeps technical acquisitions separate and produces descriptive results. It does **not** combine technical replicates automatically, test donor groups, implement batch correction or run PCA. Group comparisons require the still-unspecified scientific contrast and design. High assay-blank/isotype signal and acquisition drift require pilot-defined review; no manufacturer-unsupported numerical thresholds are invented here.

## Supported FCS and gating scope

- FCS3.0/3.1, single-dataset list mode, linear float32/float64 or byte-aligned unsigned integer parameters, common endianness, unique `$PnN` names.
- Nonzero `$PnE`, nonunit `$PnG`, supplemental TEXT, multiple datasets, malformed offsets and nonfinite event values are explicitly rejected. Request a compatible export or extend/validate the reader separately; do not strip metadata to bypass a rejection. Real FACSDiva files still require comparison against acquisition software.
- Explicit compensation policies: `none` with rationale, `embedded` matrix, or `explicit` matrix. Already-compensated inputs cannot receive another matrix. An embedded matrix does not establish whether compensation was already applied.
- Reviewed polygon gates in linear or asinh coordinates. Pre-gates use conventional Fortessa `FSC-*`, `SSC-*`, and `Time` names, with actual area-vs-height/width singlet geometry. No gating on APC. Every manufacturer bead ID must have one population definition; overlapping events are excluded and affected populations flagged.
- Saturation is checked on original APC before computational compensation using a conservative `$PnR - 1` boundary or a documented explicit threshold. Already-compensated inputs need a `saturation_review_note` attesting review of original acquisition QC; otherwise measurements are invalid because original clipping cannot be established.

The initial [analysis plan](docs/analysis-plan.md) documents the scientific rationale. The panel identity map was transcribed from the uploaded Miltenyi manual, document **140-006-716.02**, Table 1 (printed p.10). Acquisition guidance is referenced at p.19; correction/normalization at pp.20–21. Obtain that separate acquisition application note before accepting real gates.

## Verification

```bash
MPLCONFIGDIR=/tmp/facs-neuro-mpl XDG_CACHE_HOME=/tmp/facs-neuro-cache .venv/bin/python -m unittest discover -s tests -v
```

The tests exercise FCS encoding/offsets/scaling rejection, compensation, gated population identity and ambiguity, signed correction, blank matching, normalization failures, donor/acquisition identity, report safety, and the FCS-to-CSV/report workflow. The synthetic demo verifies software behavior; it does not validate real sample preparation, detector settings, reference gates or biological findings.
