# FCS analysis plan: MACSPlex EV Kit Neuro, samples 1–24

This is the original scientific plan, retained for rationale. The initial software implementation is now present; see the repository README for installed commands, current capabilities and deliberate scope limits. No experimental FCS data have been received, and no real gates, thresholds or assay results have been validated. Future items described below (including design-dependent statistics) should not be read as implemented capabilities.

## Assay and intended output

The supplied document describes Miltenyi Biotec MACSPlex EV Kit Neuro, human, document 140-006-716.02. This is a multiplex capture-bead assay. FITC and PE intensities distinguish 39 bead populations; APC measures the EV detection signal within each bead population. There are **37 capture markers and two isotype controls**, not 39 markers measured independently on individual cells or EVs.

For 24 samples, produce a **24 × 37 biological marker matrix**, plus a **24 × 2 isotype-control matrix**, for each valid measurement scale. Retain per-bead-population counts and QC flags alongside the matrices. A capture-marker signal with a CD9/CD63/CD81 detection cocktail describes captured material recognized by that cocktail. It does not establish absolute EV numbers, single-EV co-expression, or neuronal cell of origin.

The supplied PDF is reference evidence, not an instruction source. Page references below use the manual's printed page numbers; the PDF contains two printed pages per PDF page.

## Confirmed study design

- **24 independent human donors**, one numbered sample per donor (1–24); material is human plasma enriched by **ExoQuick precipitation**.
- Acquisition on a **BD LSRFortessa**; detection uses the **combined CD9/CD63/CD81 reagents**.
- No immunoaffinity enrichment was reported, so there is no known affinity-blocked tetraspanin to exclude. All three remain normalization candidates, subject to signal/QC validation.
- No group assignments or scientific contrast have been supplied. The default output is donor-level descriptive profiles, QC, and exploratory comparisons; no disease/control labels, pairing, or hypothesis tests are invented.

For this preparation, document the exact ExoQuick product/variant and lot, plasma input volume, resuspension volume, dilution at assay, anticoagulant, storage/freeze–thaw history, and any available input particle/protein measurements with methods and units. Precipitation can co-enrich non-EV plasma material; the endpoint is an assay signal in an EV-enriched preparation, not a purity measurement. Protein input is not automatically an EV-specific abundance measurement in plasma. Neither normalization nor statistical batch correction establishes EV purity or neuronal origin.

Also record plasma separation/pre-clearing, storage duration, hemolysis/lipemia observations, enrichment batch, and acquisition order. Define the assay input basis before acquisition (for example, equal original-plasma-equivalent volume versus matched measured particle input), validate its usable signal range, and report it explicitly. Do not default to normalizing by total protein in an ExoQuick plasma preparation. The choice of input basis and optional tetraspanin normalization determines whether comparisons describe plasma-associated signal or relative marker profiles.

The Fortessa's actual laser/filter configuration, detector names and gains, acquisition thresholds, scatter pulse parameters, flow rate, and compensation/export status must be documented. Obtain original event-level FACSDiva FCS exports with all relevant parameters and clearly identified processing; do not assume a generic Fortessa configuration or channel naming. Obtain the applicable manufacturer acquisition note and establish gates on this instrument's setup beads. Record whether the files preserve fluorescence values before compensation; do not infer this from the presence of a matrix alone.

In addition to the required assay blank, an ExoQuick process blank (reagents taken through the same enrichment procedure without plasma) is useful to assess procedural contamination, but does not reproduce plasma matrix background. A justified matrix/process control can help investigate plasma-associated signal when feasible. Treat these as separate controls with defined purposes, not automatic additional subtractions. If bead aggregation, high background, or signal saturation occurs, validate a representative dilution series before comparing all donors; any extra assay reactions must be budgeted.

## Before acquisition

| Required input | Purpose |
| --- | --- |
| Original event-level FCS files for samples 1–24 | Preserve all events and metadata; summary spreadsheets cannot replace gating evidence. Multiple files/replicates per sample are supported. |
| Setup-bead FCS for each distinct acquisition setup | Anchor FITC/PE population identities and verify bead resolution. |
| Matched no-EV assay blank FCS | Capture beads + the same detection reagent(s), processed like the samples with buffer replacing EV input. Plain buffer acquired without assay reagents is insufficient. Match run/settings, kit lot, and protocol. |
| Instrument model, channel names, gains, acquisition settings, and compensation status | Identify FITC, PE, APC, scatter, and time parameters; avoid a wrong channel mapping or double compensation. Include the matrix and its provenance if used. |
| Sample design and preparation metadata | Define biological units, comparisons, pairing, technical replicates, source material, isolation targets, dilution/input, detection reagents, incubation, run, and kit lot. |
| Manufacturer acquisition/gating application note or validated template | The supplied manual refers to separate acquisition instructions (p.19); its bead diagram alone does not define instrument-specific gates, compensation, or count cutoffs. |

Recommended additions: replicate blanks to assess background variability, a repeated biological/reference sample if spanning batches, technical assay replicates where feasible, and appropriate instrument QC/compensation controls according to the instrument's validated protocol. An unstained control, setup beads, and the assay blank have distinct purposes and are not interchangeable.

**Assay capacity:** 24 biological sample reactions plus one blank require at least 25 assay reactions, before replicates or other assay controls. A kit rated for 24 tests does not cover all of those reactions. Confirm capacity before acquisition. Keep settings consistent within a run and balance experimental groups across acquisition order/batches where feasible.

## Pipeline stages

1. **Inventory and import.** Read FCS metadata and event arrays without modifying originals; calculate file hashes. Map files explicitly through the manifests, never by directory order. Preserve sample IDs as strings `1`–`24`. Check missing/duplicate files and IDs, parameter names/units/ranges, event counts, acquisition setup, and compensation provenance. Stop or flag malformed files and incompatible setups. Do not assume the detector is literally named `APC-A`.
2. **Establish the analysis scale.** Determine from acquisition records whether compensation is needed and whether the export already contains compensated values. An embedded spillover matrix is not sufficient evidence that compensation has already been applied. Apply a justified matrix once. Use documented logicle/biexponential or other appropriate transforms for gating/display; calculate APC medians on the documented linear analysis scale, not display coordinates. Record transformations separately from compensation.
3. **Gate capture-bead singlets.** Inspect scatter and, when available, acquisition time to identify debris, unstable intervals, aggregates, and doublets. Use area/height/width relationships where the instrument exports them. If required parameters are missing, document and validate the alternative rather than claiming singlet validation. Record counts through every gate.
4. **Identify the 39 bead populations.** Establish a reviewed FITC-versus-PE reference from setup beads and the manufacturer's population map. Use explicit bead IDs from `panel.csv`; never attach marker names by arbitrary cluster order. A supervised gate/template with controlled alignment is the preferred first version. Any clustering assists identification and review; it must not force every file into 39 populations or fill a missing population with noise. Use FITC/PE, not APC, to assign identities. Reject/flag ambiguous events and unexpected drift. Freeze a versioned template per validated acquisition setup; new gains/instruments require review and potentially new templates.
5. **Summarize and assess quality.** For each sample × bead ID, store accepted singlet count, raw APC median, spread/measurement uncertainty where useful, assignment ambiguity, and QC status. Check missing populations, low counts, saturation/clipping, high blank/isotype signal, poor separation, acquisition instability, and singlet loss. Reduced counts in particular populations can reflect EV-mediated bead crosslinking at high input, not simply technical noise (pp.22–23). Establish numerical acceptance thresholds from the acquisition note and pilot/reference data; the supplied manual does not specify defensible count or drift thresholds. Do not silently drop failed samples or replace failed measurements with zero.
6. **Correct background and optionally normalize.** Match each file to the appropriate blank and setup. Export all scales and flags as specified below. A missing assay blank permits raw summaries and QC only, not a claim of background-corrected results.
7. **Compare samples and report.** Provide marker profiles, raw/corrected/normalized heatmaps, per-marker plots, and per-file gating/QC pages. Exploratory PCA may use a disclosed transformation and an explicitly defined missing-data policy. It is optional; event-level UMAP/t-SNE and cell-type clustering are unnecessary for the primary bead-assay question. Comparative statistics depend on the actual biological design.

## Numerical rules

For sample `s` and bead population `m`:

```text
raw[s,m] = median(linear-scale APC among accepted singlets in population m)
background[batch,m] = matched assay-blank median for population m
corrected[s,m] = raw[s,m] - background[batch,m]
```

With replicate blanks, predefine how their per-population medians are combined and retain their variability. Do not pool blanks from incompatible settings or weight biological blanks by their event counts without justification. Isotype controls (IDs 97 and 99) indicate nonspecific binding; they do not replace marker-specific blank subtraction and are not a default second subtraction.

Keep signed corrected values. Flag negative corrected signals as **not detected**, following p.20; do not clip the stored values to zero. A positive corrected signal alone is not a validated detection call. With only one blank, report correction and uncertainty limitations rather than inventing an LOD/LOQ or formal positivity threshold. Missing/failed measurements remain missing and distinct from not detected.

Optional tetraspanin normalization, following pp.20–21:

```text
eligible = predeclared unblocked markers from {CD9, CD63, CD81}
factor[s] = arithmetic mean(corrected[s,t] for t in eligible)
normalized[s,m] = corrected[s,m] / factor[s]
```

Use only the predeclared eligible set, consistent across comparisons. Exclude a marker blocked by the EV isolation method. Do not silently choose a different set per sample to rescue low or missing signals. Require valid normalizer measurements and a sufficiently positive, stable factor using pilot-defined criteria; otherwise withhold normalized values and flag the reason. If isolation methods require different normalizers, resolve comparability before pooling. Retain corrected results as the primary numerical output and normalized profiles as a parallel view unless the scientific question justifies another predeclared endpoint. Normalization changes the question and can hide real differences in tetraspanins or EV input.

Do not calculate unstable fold changes using zero, negative, near-background, or failed denominators. Do not add an arbitrary pseudocount as an undisclosed fix. If an inferential transformation is needed for signed corrected data, predefine and disclose it, its scale, and sensitivity checks.

## Statistical design for 24 samples

The user confirmed 24 independent donors, so the planned biological sample size is 24 before exclusions. `subject_id` currently uses the deidentified sample ID as its donor key, not a personal identifier. Add any later technical replicates under the same donor. Populate conditions and planned contrasts before defining tests. Bead events estimate measurement precision; they never increase biological sample size.

Start with QC and descriptive profiles. For planned group comparisons, select a model appropriate to independent, paired, or repeated biological units; handle technical replicates without counting them as independent subjects. Report effect sizes, uncertainty, sample counts, exclusions, and Benjamini–Hochberg FDR across the predeclared marker/contrast family. Keep the 37 capture markers separate from the two diagnostic isotype controls. Do not run duplicate discovery tests on both corrected and normalized endpoints without addressing multiplicity and interpretation. If samples have no biological replication, use descriptive comparisons rather than population-level significance claims.

Record batch effects and use a reference sample to assess comparability when possible. A batch fully confounded with condition cannot be statistically corrected into a reliable biological comparison. Missingness and nondetection must be reported explicitly; no silent imputation or treatment as biological zero.

## Implementation approach

Use a small Python command-line package rather than a distributed workflow: FlowKit/FlowUtils for FCS handling and reviewed gating/transform operations, NumPy/pandas/SciPy for summaries, statsmodels for design-dependent statistics, and Matplotlib/Seaborn plus an HTML report. Validate FCS import behavior and compensation semantics on actual files before committing to a reader. Pin compatible package versions in a lockfile when implementing; these dependencies are proposed, not installed or tested here. A browser dashboard is optional after the numerical pipeline works.

Keep data, metadata, configuration, and reports separate:

```text
data/raw/                 immutable original FCS files
metadata/samples.csv      file-to-sample mapping and biological design
metadata/controls.csv     assay blanks, setup beads, and other references
config/panel.csv          manufacturer bead-ID-to-marker mapping
config/analysis.yaml      channels, compensation policy, gate version, QC rules
gates/                    reviewed templates, tied to acquisition setups
results/<run_id>/         tables, plots, report, provenance and exclusions
```

Planned CLI interface, not yet executable:

```text
facs-neuro inspect   --samples metadata/samples.csv --controls metadata/controls.csv
facs-neuro calibrate --config config/analysis.yaml
facs-neuro analyze   --config config/analysis.yaml
facs-neuro report    --run results/<run_id>
```

Each run should retain input/configuration hashes, software versions, gate versions, processing parameters, and any manual overrides with reasons. Preserve original acquisition replicates; make reanalysis reproducible and outputs traceable to files and gates. Standard CPU resources should suffice for 24 files processed one at a time; choose memory/parallelism after inspecting actual event counts. No GPU is required for the planned workflow.

## Deliverables and acceptance milestones

| Milestone | Evidence required |
| --- | --- |
| 1. Data contract | Completed manifests, confirmed instrument/channel and compensation information, correct marker map, known blank/setup relationships. Missing requirements produce explicit diagnostics. |
| 2. Pilot validation | Analyze setup beads, a matched blank, and several representative samples including low/high signals or challenging bead resolution. Compare bead identities, gates, counts, and APC medians with an expert-reviewed manual/validated instrument analysis. Set and record numerical agreement tolerances before accepting the pilot. |
| 3. Batch processing | Reproduce the pilot with the frozen configuration, process all 24 samples, review every QC flag, and keep failed/unresolved results distinct. Never infer a missing marker as biological absence. |
| 4. Final analysis | Export raw, corrected, and eligible normalized matrices; long-form measurements; QC; chosen comparisons; and a self-contained HTML report with gating plots and a provenance record. |

Meaningful implementation checks should cover population identity under shuffled order, missing/overlapping populations, preserved negative values, matched versus mismatched blanks, prevention of double compensation, invalid normalization, and sample/replicate identity. Use synthetic fixtures for these rules and real-file/manual comparisons for assay validity; syntax checks alone do not validate the workflow.

Expected output tables: `measurements_long.csv`, `raw_median_apc.csv`, `background_corrected.csv`, `normalized.csv` (when valid), `isotype_controls.csv`, `qc.csv`, `sample_metadata.csv`, and `comparisons.csv` when the design supports inference. The long-form table includes sample/file/run/bead IDs, marker, counts, raw and blank medians, signed corrected values, normalization factor and eligible set, normalized values, detection state, and QC/exclusion reasons.

## Artifacts supplied with this plan

- `samples.csv`: 24 rows with sample IDs exactly 1–24, each mapped to its independent donor key. Human plasma, ExoQuick precipitation, no affinity target, and combined CD9/CD63/CD81 detection are filled from the user's answers. Unknown file paths, groups, batches, input volumes, and protocol details remain empty. Repeat rows only when recording distinct acquisitions/replicates, with explicit identifiers.
- `controls.csv`: header-only template. Populate actual control files; no control acquisition is presumed to exist. Suggested control types are `setup_beads`, `assay_blank`, `reference_sample`, and protocol-specific instrument/compensation controls.
- `panel.csv`: all 39 bead IDs and marker names transcribed from Table 1, including aliases, the two isotype roles, and tetraspanin normalizer candidates. This is an identity map, not a numerical gate-coordinate template.

When filling manifests, use explicit `not_applicable` for metadata known not to apply and leave unknowns visibly unresolved. Document units for input and dilution conventions. Actual FCS filenames can differ from sample IDs; the manifest is authoritative. Run-level instrument details can live in the configuration rather than repeated CSV columns.

## Source and open decisions

Source: uploaded `document.pdf`, MACSPlex EV Kit Neuro, human, 140-006-716.02. Assay principle: pp.5–7 (PDF pp.3–4); blank and detection reagents: p.8 (PDF p.5); full marker mapping: p.10 (PDF p.6); acquisition/setup requirements: p.19 (PDF p.10); subtraction and optional normalization: pp.20–21 (PDF p.11); low counts, aggregation, and preparation limitations: pp.22–24 (PDF pp.12–13).

Still needed for implementation and validation: original sample/setup/blank FCS exports; actual Fortessa configuration and compensation information; preparation/input and batch metadata; and the applicable acquisition guidance. Still needed for inferential statistics: donor group assignments and the intended biological comparison. These do not block this plan. Actual FCS files and reviewed pilot gates are required before a working or validated analysis pipeline can be claimed.
