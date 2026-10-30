# Analysis run of 29 October 2026, setup-bead baseline

This records the first analysis of real donor files. The run had **no working assay blank**, so the setup-bead acquisition `20261029_Beads_5000.fcs` serves as the per-bead baseline (`background.mode: setup_baseline`). Setup beads carry no detection reagent. Baseline subtraction therefore removes bead autofluorescence but **not** nonspecific binding of the detection reagent, and no detection or nondetection calls are made.

Reproduce with the FCS files in `data/raw/` (git-ignored):

```bash
.venv/bin/facs-neuro analyze --config config/analysis.json --output results/run-20261029-setup-baseline
```

## Inputs

| Item | Value |
| --- | --- |
| Samples | 1 to 24 except 10; one acquisition each, 13:15 to 13:52 |
| Excluded | Sample 10: no FCS file in the export archive |
| Reference | Setup beads, acquired 13:54 after the samples |
| Instrument settings | Identical in all 24 files: FSC 150, SSC 120, FITC 280, PE 480, APC 520 V; FSC threshold 5000; identity spill matrix |
| Gates | `fortessa-20261029-v1`, derived from the setup beads; minimum 50 events per population (provisional) |
| Not recorded | Kit lot and incubation protocol; entered as `unrecorded_single_*_run_20261029` placeholders |

## Results

- **Identity:** the Figure 2 order audit is consistent for all 39 populations in every sample, with no ambiguous events. The setup-bead gates transfer to the samples.
- **Measurements:** 819 of 897 are valid. The run status is `completed_with_qc_failures`.
- **Low bead counts:** samples 5, 22 and 20 lose 30, 20 and 9 populations below 50 events. These samples have many FSC-W doublets, consistent with EV-mediated bead crosslinking at high input (manual p.22); the manual recommends repeating with diluted EVs. Population 53 GLAST is under 50 events in ten further samples.
- **APC saturation:** samples 18 and 24 each have three populations with clipped APC events (18: CD29, CD9, CD36; 24: CD29, CD47, CD49f). These are invalid.
- **High nonspecific signal:** in samples 3, 8, 11, 4 and 9 both isotype controls read 5,000 to 20,000 above baseline, close to the tetraspanin signal. Their marker values are inflated by nonspecific binding that a real assay blank would have measured. Compare them with other samples only with that caveat, or consider excluding them.
- **Strongest signals overall:** CD36, CD9, CD81, CD63 and CD49f.

## Recommendations

1. Include a buffer-only assay blank (beads plus detection reagent, no EVs) in the next run. It is the only control that measures reagent background.
2. Repeat samples 5, 20 and 22 diluted.
3. Record the kit lot and incubation protocol in `metadata/` so the placeholders can be replaced.
4. Decide on a validated `min_events_per_population` and on whether to exclude the high-isotype samples.
