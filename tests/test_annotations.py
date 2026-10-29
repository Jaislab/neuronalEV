"""Curated association metadata must not alter quantitative assay results."""

import csv
import io
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal

from facs_neuro import annotations
from facs_neuro.annotations import (
    ANNOTATION_COLUMNS, ANNOTATION_VERSION, GROUP_LABELS,
    annotate_measurements, load_annotations, load_sources,
)
from facs_neuro.config import EXPECTED_MARKERS


def measurement_fixture():
    frame = load_annotations()[["bead_id", "marker", "role"]].copy()
    frame["sample_id"] = "1"
    frame["n_events"] = np.arange(39, dtype=np.int64)
    frame["raw_median_apc"] = np.arange(39, dtype=float) * 0.125
    frame["corrected_apc"] = np.arange(39, dtype=float) - 17.5
    frame.loc[1, "corrected_apc"] = np.nan
    frame["measurement_valid"] = frame["n_events"] > 0
    frame["qc_flags"] = ""
    return frame


def csv_text(frame):
    handle = io.StringIO()
    writer = csv.DictWriter(handle, fieldnames=frame.columns.tolist())
    writer.writeheader()
    writer.writerows(frame.to_dict("records"))
    return handle.getvalue()


class AnnotationTests(unittest.TestCase):
    def test_all_39_canonical_identities_are_present_with_control_roles(self):
        table = load_annotations()
        self.assertEqual(len(table), 39)
        self.assertEqual(table["bead_id"].nunique(), 39)
        self.assertEqual(table.set_index("bead_id")["marker"].to_dict(), EXPECTED_MARKERS)
        self.assertEqual(int(table["role"].eq("capture_marker").sum()), 37)
        controls = table[table["role"] == "isotype_control"]
        self.assertEqual(set(controls["bead_id"]), {"97", "99"})
        self.assertTrue(controls["biological_group"].eq("isotype_control").all())
        self.assertTrue(controls["association_tags"].eq("assay_control;isotype_background").all())
        self.assertTrue(table["annotation_version"].eq(ANNOTATION_VERSION).all())

    def test_primary_display_groups_are_explicit_and_do_not_overlap(self):
        table = load_annotations()
        expected = {
            "neuron_associated": {"32", "33", "75", "96"},
            "glial_progenitor_associated": {"24", "52", "53", "63", "65"},
            "shared_neural_progenitor": {"46", "98"},
            "immune_myeloid": {"34", "43", "45", "55", "77", "87", "89"},
            "vascular_stromal_adhesion": {"23", "35", "54", "57", "64", "76", "78", "79", "85", "86"},
            "broad_receptor_signaling": {"22", "42", "44", "56", "74"},
            "lysosomal_trafficking": {"88"},
            "tetraspanin_reference": {"66", "67", "68"},
            "isotype_control": {"97", "99"},
        }
        self.assertEqual(set(table["biological_group"]), set(GROUP_LABELS))
        for group, beads in expected.items():
            with self.subTest(group=group):
                self.assertEqual(set(table.loc[table["biological_group"] == group, "bead_id"]), beads)

    def test_numbers_negative_values_nan_dtypes_and_order_are_unchanged(self):
        frame = measurement_fixture()
        original = frame.copy(deep=True)
        result = annotate_measurements(frame)
        assert_frame_equal(result[frame.columns], original)
        assert_frame_equal(frame, original)
        self.assertTrue(np.isnan(result.loc[1, "corrected_apc"]))
        self.assertLess(result.loc[0, "corrected_apc"], 0)
        self.assertEqual(result.columns.tolist(), frame.columns.tolist() + list(ANNOTATION_COLUMNS))
        self.assertFalse(any("score" in column or "fraction" in column for column in result))

    def test_shuffled_rows_and_duplicate_indexes_are_safe(self):
        original = measurement_fixture().iloc[[38, 0, 15, 0, 32, 11]].copy()
        original.index = pd.Index([7, 7, 2, 7, 2, 7], name="original_row")
        result = annotate_measurements(original)
        assert_frame_equal(result[original.columns], original)
        expected = ["isotype_control", "broad_receptor_signaling", "immune_myeloid", "broad_receptor_signaling", "immune_myeloid", "shared_neural_progenitor"]
        self.assertEqual(result["biological_group"].tolist(), expected)

    def test_integer_bead_ids_are_preserved(self):
        frame = measurement_fixture()
        frame["bead_id"] = frame["bead_id"].astype(np.int64)
        result = annotate_measurements(frame)
        assert_frame_equal(result[frame.columns], frame)

    def test_mismatched_id_marker_or_role_is_rejected(self):
        for field, value in [("bead_id", "999"), ("marker", "CD56"), ("role", "isotype_control")]:
            with self.subTest(field=field):
                frame = measurement_fixture()
                frame.loc[0, field] = value
                with self.assertRaisesRegex(ValueError, "canonical"):
                    annotate_measurements(frame)

    def test_conflicting_annotations_are_rejected_matching_ones_are_idempotent(self):
        frame = annotate_measurements(measurement_fixture())
        assert_frame_equal(annotate_measurements(frame), frame)
        frame.loc[0, "biological_group"] = "neuron_associated"
        original = frame.copy(deep=True)
        with self.assertRaisesRegex(ValueError, "Conflicting preexisting annotation"):
            annotate_measurements(frame)
        assert_frame_equal(frame, original)

    def test_overlapping_tags_preserve_important_alternative_associations(self):
        table = load_annotations().set_index("bead_id")
        required_tags = {
            "24": {"oligodendrocyte_progenitor_associated", "mesenchymal_associated"},
            "33": {"neural_associated", "nk_cell_associated"},
            "35": {"endothelial_associated", "platelet_associated", "leukocyte_associated"},
            "43": {"microglia_associated", "peripheral_myeloid"},
            "55": {"microglia_associated", "peripheral_myeloid"},
            "63": {"neural_glial_progenitor_associated", "neuronal_associated", "ganglioside_epitope"},
            "65": {"oligodendrocyte_progenitor_associated", "pericyte_associated"},
            "87": {"platelet_associated", "peripheral_myeloid", "microglia_associated"},
            "89": {"microglia_associated", "peripheral_myeloid"},
        }
        for bead, expected in required_tags.items():
            with self.subTest(bead=bead):
                self.assertTrue(expected <= set(table.at[bead, "association_tags"].split(";")))
        all_tags = set(";".join(table["association_tags"]).split(";"))
        self.assertNotIn("neuron_specific", all_tags)
        self.assertNotIn("microglia_specific", all_tags)

    def test_critical_interpretation_caveats_and_epitope_identity(self):
        table = load_annotations().set_index("bead_id")
        notes = table["interpretation_note"].str.lower()
        for bead, words in {
            "32": ("soluble", "plasma", "l1cam"),
            "33": ("nk", "non-neural"),
            "24": ("mesenchymal", "non-neural"),
            "65": ("pericyte", "mesenchymal"),
            "96": ("vesicular", "epitope", "accessibility"),
            "43": ("microglia", "peripheral"),
            "55": ("microglial", "peripheral"),
            "89": ("microglia", "peripheral"),
            "35": ("platelet", "leukocyte"),
            "87": ("platelet", "peripheral"),
        }.items():
            with self.subTest(bead=bead):
                for word in words:
                    self.assertIn(word, notes.at[bead])
        for bead in ("52", "63", "98"):
            self.assertIn("not a single gene", notes.at[bead])
            self.assertIn("epitope", table.at[bead, "gene_or_epitope"].lower())
        for bead in ("66", "67", "68"):
            self.assertIn("not a universal ev marker", notes.at[bead])

    def test_source_registry_resolves_references_and_correct_vglut2_mapping(self):
        registry = load_sources()
        self.assertEqual(registry["annotation_version"], ANNOTATION_VERSION)
        self.assertIn("not manufacturer-validated", registry["curation_notice"])
        sources = registry["sources"]
        self.assertIn("identity only", sources["kit_table1"]["evidence_scope"])
        self.assertEqual(sources["MISEV2023"]["doi"], "10.1002/jev2.12404")
        self.assertEqual(sources["L1CAM2021"]["doi"], "10.1038/s41592-021-01174-8")
        self.assertEqual(sources["GD2_non_neural_2022"]["doi"], "10.1111/nyas.14700")
        table = load_annotations().set_index("bead_id")
        self.assertEqual(table.at["96", "gene_or_epitope"], "SLC17A6")
        self.assertIn("uniprot:Q9P2U8", table.at["96", "annotation_source_ids"])
        self.assertNotIn("uniprot:Q9P2U7", sources)
        self.assertIn("GD2_non_neural_2022", table.at["98", "annotation_source_ids"])
        self.assertIn("Non-neural/tumor", sources["GD2_non_neural_2022"]["evidence_scope"])
        for references in table["annotation_source_ids"]:
            self.assertTrue(set(references.split(";")) <= set(sources))

    def test_corrupt_resource_fails_closed(self):
        original_reader = annotations._read_resource
        original = load_annotations()
        changes = [
            ("marker", "CD56"),
            ("role", "isotype_control"),
            ("biological_group", "neuron_associated"),
            ("annotation_version", "0.0.0"),
            ("annotation_source_ids", "does_not_exist"),
            ("association_tags", "neuron_specific"),
        ]
        for field, value in changes:
            with self.subTest(field=field):
                modified = original.copy()
                modified.loc[0, field] = value
                payload = csv_text(modified)
                with patch.object(annotations, "_read_resource", side_effect=lambda name: payload if name == "marker_annotations.csv" else original_reader(name)):
                    with self.assertRaises(ValueError):
                        load_annotations()
        for modified in (original.iloc[:-1], pd.concat([original, original.iloc[[0]]])):
            payload = csv_text(modified)
            with patch.object(annotations, "_read_resource", side_effect=lambda name: payload if name == "marker_annotations.csv" else original_reader(name)):
                with self.assertRaisesRegex(ValueError, "all 39"):
                    load_annotations()

    def test_empty_measurement_table_is_supported_without_inventing_rows(self):
        frame = measurement_fixture().iloc[:0]
        result = annotate_measurements(frame)
        self.assertEqual(len(result), 0)
        assert_frame_equal(result[frame.columns], frame)


if __name__ == "__main__":
    unittest.main()
