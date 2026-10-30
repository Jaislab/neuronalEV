# Configuration and acquisition contract

`config/analysis.json` contains paths relative to itself. FCS paths inside a manifest are relative to that manifest. Unknown real-experiment values are intentionally not filled.

Sample FCS, assay-blank FCS and the kit's setup-bead FCS are the required input types. Gate definitions are established from the setup-bead events during analysis; you do not need a separate provider gate file. The [Figure 2/Table 1 guide](bead-identification.md) explains the 39 assignments and the setup-reference ordering audit. Relative order is checked using observed FITC/PE medians; schematic ordinal levels are never used as numeric gate coordinates.

## Manifest keys

Each acquired sample row requires:

- `sample_id`: original biological sample label (1–24 here).
- `acquisition_id`: unique event-file acquisition label; defaults to sample ID when omitted. Give repeated acquisitions distinct IDs while preserving `sample_id`/`subject_id`.
- `subject_id`: deidentified donor key. The provided study has 24 independent donor keys.
- `fcs_path`, `batch_id`, `setup_id`, `instrument_settings_id`, `kit_lot`, `detection_reagent`, `incubation_protocol`.
- `blank_id`: exact `control_id` for the applicable `assay_blank`.
- `setup_control_id`: exact `control_id` for the applicable `setup_beads` reference file. This differs from `setup_id`, which selects the configuration object.

Controls require unique `control_id`, `control_type`, `fcs_path`, batch, setup, instrument-settings ID and kit lot. Assay blanks additionally require the same detection reagent and incubation protocol as their sample. Setup/blank references cannot cross incompatible batches, settings or kit lots. No implicit pooling of replicate blanks occurs: each sample selects one file. Preserve replicate blanks as separate controls; pooling/background-uncertainty models need an explicit later implementation.

Accepted control types are `setup_beads`, `assay_blank`, `reference_sample`, `process_blank`, `compensation_control`, `unstained_control`, and `instrument_qc`. Unknown spellings are rejected even for unused controls. All files in the analysis manifest are processed with the configured bead-gating setup; separate non-assay instrument controls can instead be retained outside that manifest as acquisition evidence.

Readiness checks reject inconsistent donor IDs for repeated acquisitions of the same sample. A supplied normalization floor must be numeric, finite and nonnegative even when normalization is disabled; `null` remains valid for a disabled normalization template. Invalid or missing calibration-reference files and unresolved compensation settings are checked before a calibration output directory is created. Unexpected failures after writing begins retain `failure.json`.

Keep groups empty until known. Record ExoQuick variant/lot, original plasma and resuspension volumes, assay dilution/input and units, anticoagulant, enrichment batch, pre-clearing, storage/freeze–thaw and other preanalytical notes. No automatic protein/volume normalization is performed. `affinity_target` is `not_applicable` for the stated ExoQuick-only preparation; any future affinity-blocked tetraspanin must be excluded from the fixed normalizer set before comparisons.

Computed-field names such as `bead_id`, `marker`, `raw_median_apc`, `measurement_valid` and `qc_flags` are reserved and cannot appear in manifests. Keep marker aliases in the panel's `alias` column. All 39 canonical bead/marker identities and isotype roles are checked; the pipeline cannot be silently repurposed for another panel.

## Setup object

The `setups` keys are setup IDs used in manifests. Each setup specifies distinct detector names:

```json
"channels": {"fitc": "FITC-A", "pe": "PE-A", "apc": "APC-A"}
```

These names are examples. The instrument model alone does not determine installed lasers/filters, voltages or names. Within a setup, the pipeline compares available range, voltage, amplification/gain, filter/laser, instrument identity and embedded spillover metadata. When metadata are absent, the supplied setup/settings IDs and acquisition review remain essential; the software cannot reconstruct missing settings.

The comparison covers FITC, PE and APC plus every scatter/time parameter used by a pre-gate, including pulse geometry. Changes to FSC/SSC settings therefore require a separately reviewed setup just as fluorescence-setting changes do. The strict reader rejects ambiguous TEXT keywords caused by unsupported empty vendor values instead of silently losing optional voltage or compensation metadata. Escaped delimiters inside ordinary TEXT values remain supported.

Compensation examples (the matrix below illustrates the **format**, not an instrument calibration):

```json
"compensation": {
  "input_state": "uncompensated",
  "mode": "none",
  "reason": "Replace with the validated acquisition-specific rationale"
}
```

Or use `mode: "embedded"` with acquisition-confirmed `input_state: "uncompensated"`. For `mode: "explicit"`, supply ordered `channels` and a square spillover `matrix` with unit diagonal. The convention is **measured = true @ matrix**; inversion is performed with the stated channel order. Scatter/time channels are not compensation targets. Do not substitute an arbitrary identity matrix to avoid establishing compensation requirements. Already-compensated exports use `input_state: "compensated"`, `mode: "none"`, rationale, and reviewed original acquisition saturation evidence through `saturation_review_note`.

Each gate uses vertices in its declared coordinate space:

```json
{
  "name": "beads",
  "x": "FSC-A",
  "y": "SSC-A",
  "x_transform": {"kind": "linear"},
  "y_transform": {"kind": "linear"},
  "vertices": []
}
```

Empty vertices intentionally cannot analyze data. Draw real polygons from reviewed controls, or start from the unreviewed template that `scripts/identify_setup_beads.py` derives from a setup-bead export (its `derived_from` entry records the source file hash and is informational). Add a distinct `singlets` pre-gate using FSC-A vs FSC-H/FSC-W or the analogous SSC pulse pair. Gates run sequentially in listed order; additional time/scatter polygons are allowed. Pre-gate axis names must be conventional FSC/SSC/Time. Missing pulse geometry needs another validated acquisition/export, rather than a fabricated singlet gate.

`bead_space` specifies the transforms for population gates (`x_transform` for FITC, `y_transform` for PE). Supported forms are `{"kind":"linear"}` and `{"kind":"asinh","cofactor":150}`; **150 is an illustrative value, not a validated setting**. Draw all population coordinates after the chosen transform. Population objects have `bead_id` and polygon `vertices`; IDs determine biological labels regardless of list order. Do not transfer gate coordinates across different units, transforms, gains or instruments without validating them.

Set `min_events_per_population` to a positive integer justified by acquisition guidance/pilot precision. No real-data default is supplied. `apc_saturation_threshold` may replace the conservative `$PnR - 1` upper boundary when the acquisition review establishes a different finite positive limit. Any ambiguous or clipped accepted event flags its affected population in this conservative first version. Raw values remain auditable; invalid rows cannot enter correction/normalization. Population centroid drift and high blank/isotype signal need manual review; the current release has no automatic alignment or validated cutoff for them.

Complete the review with `reviewed: true`, a version string and a substantive `review_note`, only after setup-bead recognition and pilot agreement. A Boolean review flag records a human decision; it cannot itself establish assay validity. All experimental reports remain labeled as requiring real-file assay validation.

## Normalization

The default is disabled. To enable it, declare a fixed unique subset of CD9/CD63/CD81 and an explicit `minimum_factor`. Every chosen corrected normalizer must be valid and positive; their arithmetic mean must exceed the floor. Determine that floor from pilot/background behavior. Do not copy the synthetic demonstration's value. Invalid normalizers suppress normalized output and append reasons; raw/corrected results remain available where valid.

In the default assay-blank mode, negative corrected target-marker values remain signed, including after valid normalization. They are labeled `not_detected`, not removed. Neither positive correction nor successful normalization is a validated detection limit or an EV-purity assessment. The setup-baseline fallback below uses baseline-relative states instead.

## Explicit setup-bead baseline fallback

The default background mode remains a matched no-EV assay blank. If the assay blank is known to be invalid and a setup-bead reference is deliberately selected as a limited fallback, declare the choice explicitly:

```json
"background": {
  "mode": "setup_baseline",
  "rationale": "The assay blank was processed incorrectly; use the specified setup-bead file as a limited baseline. Detection-reagent background is unmeasured."
}
```

The existing sample `blank_id` column identifies the selected reference in either mode. In `setup_baseline` mode it must reference a control of type `setup_beads`; no detection cocktail is falsely attributed to that control. Match batch, setup, instrument settings and kit lot. A single physical control row may serve both as `setup_control_id` and the baseline reference without duplicating its FCS path. Separate setup controls can be used for identification and baseline when explicitly selected and compatible.

This subtracts the bead/setup baseline only. It does **not** correct nonspecific detection-reagent binding or replace a fully processed no-EV assay blank. Every result retains `background_mode`, `background_reference_id`, `background_limitation` and a `setup_baseline_only` QC flag. Signed values and missingness remain unchanged. Baseline states are `below_setup_baseline`, `at_setup_baseline`, and `above_setup_baseline`; none is a validated detection/absence call. Optional normalization does not remove this limitation. Count, ambiguity and saturation requirements remain unchanged; a sparse baseline cannot be rescued by silently relaxing QC.

## Audited sample exclusions

Known preparation failures can be excluded explicitly without converting missing assay evidence into biological zeros:

```json
"sample_exclusions": {
  "10": "Capture beads lost during preparation; exclude biological interpretation."
}
```

All acquisitions of an excluded sample are omitted from numerical analysis. The original full input manifest is retained, while `excluded_samples.csv` and provenance preserve the sample IDs and exact reasons. Excluded samples do not require usable FCS files and do not acquire invented zero measurements. Unknown sample IDs, empty reasons and an exclusion leaving no active samples are errors. `inspect` distinguishes input, active and excluded acquisitions. The default empty exclusion map preserves every sample.

## Native FACSDiva Time values

The reader supports the observed floating-point FACSDiva convention in which a channel named `Time` has nonunit `$PnG` equal to a positive `$TIMESTEP`. It preserves the stored native Time values and metadata; a Time pre-gate must therefore use native ticks. Multiply by `$TIMESTEP` once for seconds in separate diagnostics. Signal/scatter channels still require unit gain, and missing or mismatched timestep metadata is rejected. This exception never changes FITC, PE or APC measurements.
