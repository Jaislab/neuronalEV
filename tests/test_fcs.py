"""Independent binary fixtures exercise the reader without clinical files."""

from pathlib import Path
import hashlib
import struct
import tempfile
import unittest

import numpy as np

from facs_neuro.fcs import FCSFormatError, read_fcs


def encode_fcs(events, *, datatype="F", endian="<", widths=None,
               changes=None, text_offsets=False, version=b"FCS3.1", delimiter="|"):
    """Construct fixtures from the FCS HEADER/TEXT/DATA specification directly."""
    events = np.asarray(events)
    if events.ndim != 2:
        raise ValueError("Fixture needs an explicit event-by-parameter array")
    n, p = events.shape
    widths = widths or [32 if datatype != "D" else 64] * p
    fields = {
        "$MODE": "L", "$DATATYPE": datatype, "$TOT": str(n), "$PAR": str(p),
        "$BYTEORD": "1,2,3,4" if endian == "<" else "4,3,2,1",
        "$NEXTDATA": "0",
    }
    for i, width in enumerate(widths, 1):
        fields.update({f"$P{i}N": f"channel_{i}", f"$P{i}B": str(width),
                       f"$P{i}E": "0,0", f"$P{i}G": "1", f"$P{i}R": "65536"})
    if text_offsets:
        fields.update({"$BEGINDATA": "0000000000", "$ENDDATA": "0000000000"})
    for key, value in (changes or {}).items():
        if value is None:
            fields.pop(key, None)
        else:
            fields[key] = str(value)
    def text_bytes():
        components = [part.replace(delimiter, delimiter * 2)
                      for pair in fields.items() for part in pair]
        return (delimiter + delimiter.join(components) + delimiter).encode("latin-1")

    if datatype in {"F", "D"}:
        format_character = "f" if datatype == "F" else "d"
        data = b"".join(struct.pack(endian + format_character, float(x)) for x in events.flat)
    else:
        characters = {8: "B", 16: "H", 32: "I", 64: "Q"}
        data = b"".join(struct.pack(endian + characters[widths[i]], int(row[i]))
                        for row in events for i in range(p))
    text = text_bytes()
    start = 58 + len(text)
    end = start + len(data) - 1
    if text_offsets:
        fields["$BEGINDATA"] = f"{start:010d}"
        fields["$ENDDATA"] = f"{end:010d}"
        text = text_bytes()
        assert start == 58 + len(text)
    header = version + b"    " + b"".join(
        f"{value:8d}".encode("ascii")
        for value in (58, 57 + len(text), 0 if text_offsets else start,
                      0 if text_offsets else end, 0, 0)
    )
    assert len(header) == 58
    return header + text + data


class FCSReaderTests(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tempdir.cleanup)
        self.path = Path(self.tempdir.name) / "fixture.fcs"

    def read(self, payload):
        self.path.write_bytes(payload)
        return read_fcs(self.path)

    def test_signed_float_values_channels_metadata_and_hash(self):
        values = [[-4.5, 10], [0, 65535], [7.25, -1.5]]
        payload = encode_fcs(values, changes={"$P1N": "FITC-A", "$P1S": "FITC",
                                             "$SPILLOVER": "2,FITC-A,channel_2,1,0,0,1"})
        actual = self.read(payload)
        np.testing.assert_array_equal(actual.events, values)
        self.assertEqual(actual.channels, ("FITC-A", "channel_2"))
        self.assertEqual(actual.ranges, {"FITC-A": 65536.0, "channel_2": 65536.0})
        self.assertIn("spillover", actual.metadata)
        self.assertEqual(actual.sha256, hashlib.sha256(payload).hexdigest())

    def test_big_endian_float_and_both_double_orders(self):
        values = [[-11.5, 32], [15.125, 9]]
        for datatype, endian in [("F", ">"), ("D", "<"), ("D", ">")]:
            with self.subTest(datatype=datatype, endian=endian):
                actual = self.read(encode_fcs(values, datatype=datatype, endian=endian))
                np.testing.assert_array_equal(actual.events, values)

    def test_header_zero_offsets_use_text(self):
        actual = self.read(encode_fcs([[1, 2], [3, 4]], text_offsets=True))
        np.testing.assert_array_equal(actual.events, [[1, 2], [3, 4]])

    def test_fcs30_and_escaped_text_delimiters(self):
        actual = self.read(encode_fcs([[1]], version=b"FCS3.0",
                                     changes={"$P1N": "APC|A", "COMMENT": "one|two|"}))
        self.assertEqual(actual.channels, ("APC|A",))
        self.assertEqual(actual.metadata["comment"], "one|two|")

    def test_empty_optional_vendor_values_do_not_silently_merge_keywords(self):
        for keyword in ["$CELLS", "$P1V", "$SPILLOVER"]:
            for delimiter in ["|", "!"]:
                with self.subTest(keyword=keyword, delimiter=delimiter):
                    # This otherwise valid FCS file has all required metadata.
                    # The adjacent separators in e.g. |$CELLS||$CYT| must not
                    # consume CYT into a bogus 'cells|$cyt' keyword silently.
                    payload = encode_fcs(
                        [[1, 2]], changes={keyword: "", "$CYT": "LSRFortessa"},
                        delimiter=delimiter,
                    )
                    with self.assertRaisesRegex(FCSFormatError, "empty vendor values"):
                        self.read(payload)

    def test_escaped_delimiters_remain_supported_in_values(self):
        for delimiter in ["|", "!"]:
            with self.subTest(delimiter=delimiter):
                expected = f"human{delimiter}plasma{delimiter}"
                actual = self.read(encode_fcs(
                    [[1]], changes={"$CELLS": expected, "$CYT": "LSRFortessa"},
                    delimiter=delimiter,
                ))
                self.assertEqual(actual.metadata["cells"], expected)
                self.assertEqual(actual.metadata["cyt"], "LSRFortessa")

    def test_mixed_integer_widths_and_byte_orders(self):
        for endian in ["<", ">"]:
            with self.subTest(endian=endian):
                values = [[8, 512, 131072], [0, 40000, 2**32 - 1]]
                actual = self.read(encode_fcs(values, datatype="I", endian=endian,
                                             widths=[8, 16, 32],
                                             changes={"$P1R": 256, "$P3R": 2**32}))
                np.testing.assert_array_equal(actual.events, values)

    def test_integer64_rejects_loss_of_precision(self):
        with self.assertRaisesRegex(FCSFormatError, "exact integer range"):
            self.read(encode_fcs(np.array([[2**53 + 1]], dtype=np.uint64),
                                 datatype="I", widths=[64]))

    def test_missing_optional_gain_means_unity(self):
        actual = self.read(encode_fcs([[3]], changes={"$P1G": None}))
        np.testing.assert_array_equal(actual.events, [[3]])

    def test_unsupported_scaling_is_not_silently_interpreted(self):
        for changes, message in [({"$P1E": "4,1"}, "amplification"),
                                 ({"$P1E": None}, "P1E"),
                                 ({"$P1G": "2"}, "gain"),
                                 ({"$P1G": "nan"}, "finite")]:
            with self.subTest(changes=changes):
                with self.assertRaisesRegex(FCSFormatError, message):
                    self.read(encode_fcs([[100]], changes=changes))

    def test_missing_duplicate_channels_and_invalid_ranges(self):
        for changes, message in [({"$P1N": None}, "P1N"),
                                 ({"$P2N": "channel_1"}, "unique"),
                                 ({"$P1R": 0}, "positive"),
                                 ({"$P1B": 64}, "inconsistent")]:
            with self.subTest(changes=changes):
                with self.assertRaisesRegex(FCSFormatError, message):
                    self.read(encode_fcs([[1, 2]], changes=changes))

    def test_unsupported_structure(self):
        for changes, message in [({"$NEXTDATA": 500}, "Multi-dataset"),
                                 ({"$BEGINSTEXT": 500}, "Supplemental TEXT"),
                                 ({"$MODE": "C"}, "event-list"),
                                 ({"$BYTEORD": "2,4,1,3"}, "BYTEORD")]:
            with self.subTest(changes=changes):
                with self.assertRaisesRegex(FCSFormatError, message):
                    self.read(encode_fcs([[1]], changes=changes))
        with self.assertRaisesRegex(FCSFormatError, "FCS3.0 and FCS3.1"):
            self.read(encode_fcs([[1]], version=b"FCS2.0"))

    def test_truncated_and_inconsistent_data_are_rejected(self):
        payload = encode_fcs([[1, 2], [3, 4]])
        with self.assertRaisesRegex(FCSFormatError, "truncated DATA"):
            self.read(payload[:-1])
        with self.assertRaisesRegex(FCSFormatError, "byte count"):
            self.read(encode_fcs([[1, 2]], changes={"$TOT": 2}))
        payload = bytearray(encode_fcs([[1]], text_offsets=True))
        payload[26:34] = b"     100"
        with self.assertRaisesRegex(FCSFormatError, "offsets disagree"):
            self.read(payload)

    def test_nonfinite_values_are_rejected(self):
        for value in [float("nan"), float("inf"), -float("inf")]:
            with self.subTest(value=value):
                with self.assertRaisesRegex(FCSFormatError, "nonfinite"):
                    self.read(encode_fcs([[value]]))

    def test_invalid_text_and_offset_fields_are_rejected(self):
        payload = bytearray(encode_fcs([[1]]))
        payload[10:18] = b"   hello"
        with self.assertRaisesRegex(FCSFormatError, "integer"):
            self.read(payload)
        payload = bytearray(encode_fcs([[1]]))
        text_end = int(payload[18:26])
        payload[text_end] = ord("!")
        with self.assertRaisesRegex(FCSFormatError, "end with its delimiter"):
            self.read(payload)

    def test_empty_event_export(self):
        actual = self.read(encode_fcs(np.empty((0, 2))))
        self.assertEqual(actual.events.shape, (0, 2))


if __name__ == "__main__":
    unittest.main()
