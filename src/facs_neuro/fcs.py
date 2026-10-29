"""Strict reader for the linear, event-list FCS subset used by this pipeline.

Supported: one FCS3.0/3.1 dataset; float32 (F), float64 (D), and byte-aligned
unsigned integer (I) parameters of 8/16/32/64 bits; common little/big endian
encodings; escaped TEXT delimiters; DATA offsets from HEADER or TEXT. Returned
events are float64 on the file's *uncompensated or already compensated* linear
measurement scale. Compensation status cannot be inferred by this reader.

Nonzero amplification decades ($PnE), nonunit gain ($PnG), supplemental TEXT,
multi-dataset files, nonfinite events, and inconsistent offsets fail explicitly.
This intentionally avoids guessing physical scaling, repairing malformed files,
or silently applying a compensation matrix. Actual instrument exports still need
comparison with the acquisition software before assay use. No display transform
is applied; negative floating-point values are preserved.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
from pathlib import Path
import math

import numpy as np


class FCSFormatError(ValueError):
    """The file is malformed or outside the reader's supported FCS subset."""


@dataclass(frozen=True)
class FCSData:
    events: np.ndarray
    channels: tuple[str, ...]
    metadata: dict[str, str]
    ranges: dict[str, float]
    sha256: str


def _integer(value: str, name: str) -> int:
    try:
        result = int(value.strip())
    except (TypeError, ValueError) as exc:
        raise FCSFormatError(f"{name} must be an integer, got {value!r}") from exc
    return result


def _number(value: str, name: str) -> float:
    try:
        result = float(value.strip())
    except (TypeError, ValueError) as exc:
        raise FCSFormatError(f"{name} must be a finite number") from exc
    if not math.isfinite(result):
        raise FCSFormatError(f"{name} must be a finite number")
    return result


def _parse_text(segment: bytes) -> dict[str, str]:
    """Tokenize FCS TEXT without losing escaped delimiters inside values."""
    if len(segment) < 3:
        raise FCSFormatError("TEXT segment is too short")
    text = segment.decode("latin-1")
    delimiter = text[0]
    if text[-1] != delimiter:
        raise FCSFormatError("TEXT segment must end with its delimiter")
    tokens: list[str] = []
    token: list[str] = []
    i = 1
    while i < len(text):
        if text[i] != delimiter:
            token.append(text[i])
            i += 1
        elif i + 1 < len(text) and text[i + 1] == delimiter:
            token.append(delimiter)
            i += 2
        else:
            tokens.append("".join(token))
            token = []
            i += 1
    if token or len(tokens) % 2:
        raise FCSFormatError("TEXT contains incomplete keyword/value pairs")
    metadata: dict[str, str] = {}
    for key, value in zip(tokens[::2], tokens[1::2]):
        key = key.strip().lstrip("$").lower()
        if not key:
            raise FCSFormatError("TEXT contains an empty keyword")
        if key in metadata:
            raise FCSFormatError(f"Duplicate TEXT keyword: {key}")
        metadata[key] = value
    return metadata


def _required(metadata: dict[str, str], name: str) -> str:
    if name not in metadata:
        raise FCSFormatError(f"Required FCS keyword ${name.upper()} is missing")
    return metadata[name]


def _offset(header: bytes, start: int, name: str) -> int:
    raw = header[start : start + 8].decode("ascii").strip()
    value = _integer(raw or "0", name)
    if value < 0:
        raise FCSFormatError(f"{name} cannot be negative")
    return value


def _bounds(start: int, end: int, size: int, name: str) -> None:
    if start < 58 or end < start or end >= size:
        raise FCSFormatError(f"Invalid or truncated {name} segment offsets: {start}..{end}")


def _overlap(a: tuple[int, int], b: tuple[int, int]) -> bool:
    return a[0] <= b[1] and b[0] <= a[1]


def read_fcs(path: str | Path) -> FCSData:
    """Read a supported FCS file, raising FCSFormatError for unsafe assumptions.

    Metadata keys are lowercase with a leading ``$`` removed. Channel names are
    $PnN values, not optional stain labels ($PnS). Integer counts are converted
    exactly to float64 or rejected if they exceed its exact integer range.
    """
    raw = Path(path).read_bytes()
    if len(raw) < 58:
        raise FCSFormatError("File is too short for an FCS header")
    header = raw[:58]
    if header[:6] not in {b"FCS3.0", b"FCS3.1"}:
        raise FCSFormatError("Only FCS3.0 and FCS3.1 files are supported")
    try:
        text_start = _offset(header, 10, "HEADER TEXT start")
        text_end = _offset(header, 18, "HEADER TEXT end")
        header_data_start = _offset(header, 26, "HEADER DATA start")
        header_data_end = _offset(header, 34, "HEADER DATA end")
        analysis_start = _offset(header, 42, "HEADER ANALYSIS start")
        analysis_end = _offset(header, 50, "HEADER ANALYSIS end")
    except UnicodeDecodeError as exc:
        raise FCSFormatError("FCS header offsets must be ASCII integers") from exc
    _bounds(text_start, text_end, len(raw), "TEXT")
    metadata = _parse_text(raw[text_start : text_end + 1])
    if _integer(metadata.get("nextdata", "0"), "$NEXTDATA") != 0:
        raise FCSFormatError("Multi-dataset FCS files are unsupported; export individual datasets")
    if any(_integer(metadata.get(k, "0"), f"${k.upper()}") != 0 for k in ("beginstext", "endstext")):
        raise FCSFormatError("Supplemental TEXT is unsupported; export a single primary TEXT segment")
    if _required(metadata, "mode").strip().upper() != "L":
        raise FCSFormatError("Only event-list FCS files ($MODE=L) are supported")

    n_parameters = _integer(_required(metadata, "par"), "$PAR")
    n_events = _integer(_required(metadata, "tot"), "$TOT")
    if not 1 <= n_parameters <= 4096 or n_events < 0:
        raise FCSFormatError("$PAR must be 1..4096 and $TOT must be nonnegative")
    datatype = _required(metadata, "datatype").strip().upper()
    if datatype not in {"F", "D", "I"}:
        raise FCSFormatError(f"Unsupported $DATATYPE={datatype!r}; expected F, D, or I")
    order = _required(metadata, "byteord").replace(" ", "")
    if order in {"1,2", "1,2,3,4", "1,2,3,4,5,6,7,8"}:
        endian = "<"
    elif order in {"2,1", "4,3,2,1", "8,7,6,5,4,3,2,1"}:
        endian = ">"
    else:
        raise FCSFormatError(f"Unsupported $BYTEORD={order!r}")

    channels: list[str] = []
    widths: list[int] = []
    ranges: dict[str, float] = {}
    for index in range(1, n_parameters + 1):
        prefix = f"p{index}"
        name = _required(metadata, prefix + "n").strip()
        if not name or name in ranges:
            raise FCSFormatError("$PnN channel names must be nonempty and unique")
        width = _integer(_required(metadata, prefix + "b"), f"$P{index}B")
        if datatype == "F" and width != 32 or datatype == "D" and width != 64:
            raise FCSFormatError(f"$P{index}B is inconsistent with $DATATYPE={datatype}")
        if datatype == "I" and width not in {8, 16, 32, 64}:
            raise FCSFormatError("Integer FCS parameters must be byte-aligned 8/16/32/64-bit values")
        amplification = _required(metadata, prefix + "e").split(",")
        if len(amplification) != 2:
            raise FCSFormatError(f"$P{index}E must contain two numbers")
        decades, zero = (_number(v, f"$P{index}E") for v in amplification)
        if decades != 0 or zero != 0:
            raise FCSFormatError(f"Nonlinear $P{index}E amplification is unsupported; export linear values")
        gain = _number(metadata.get(prefix + "g", "1"), f"$P{index}G")
        if gain != 1:
            raise FCSFormatError(f"Nonunit $P{index}G gain is unsupported; export calibrated linear values")
        value_range = _number(_required(metadata, prefix + "r"), f"$P{index}R")
        if value_range <= 0:
            raise FCSFormatError(f"$P{index}R must be positive")
        if datatype == "I" and value_range > 2**width:
            raise FCSFormatError(f"$P{index}R exceeds its integer bit width")
        channels.append(name)
        widths.append(width)
        ranges[name] = value_range

    text_data_start = _integer(metadata.get("begindata", "0"), "$BEGINDATA")
    text_data_end = _integer(metadata.get("enddata", "0"), "$ENDDATA")
    if min(text_data_start, text_data_end) < 0:
        raise FCSFormatError("DATA segment offsets cannot be negative")
    if header_data_start and text_data_start and header_data_start != text_data_start:
        raise FCSFormatError("HEADER and TEXT DATA start offsets disagree")
    if header_data_end and text_data_end and header_data_end != text_data_end:
        raise FCSFormatError("HEADER and TEXT DATA end offsets disagree")
    data_start = header_data_start or text_data_start
    data_end = header_data_end or text_data_end
    expected_bytes = n_events * sum(widths) // 8
    if not expected_bytes:
        if (data_start or data_end) and data_end != data_start - 1:
            raise FCSFormatError("Zero-event FCS file must have an empty DATA segment")
        events = np.empty((0, n_parameters), dtype=np.float64)
    else:
        _bounds(data_start, data_end, len(raw), "DATA")
        if _overlap((text_start, text_end), (data_start, data_end)):
            raise FCSFormatError("TEXT and DATA segments overlap")
        if data_end - data_start + 1 != expected_bytes:
            raise FCSFormatError("DATA byte count disagrees with $TOT, $PAR, or $PnB")
        data = memoryview(raw)[data_start : data_end + 1]
        if datatype in {"F", "D"}:
            dtype = np.dtype(endian + ("f4" if datatype == "F" else "f8"))
            events = np.frombuffer(data, dtype=dtype).reshape(n_events, n_parameters).astype(np.float64)
        else:
            dtype = np.dtype([(f"p{i}", endian + f"u{width // 8}") for i, width in enumerate(widths)])
            records = np.frombuffer(data, dtype=dtype, count=n_events)
            events = np.empty((n_events, n_parameters), dtype=np.float64)
            for i in range(n_parameters):
                column = records[f"p{i}"]
                if widths[i] == 64 and np.any(column > 2**53):
                    raise FCSFormatError("64-bit integer counts exceed float64's exact integer range")
                events[:, i] = column
        if not np.isfinite(events).all():
            raise FCSFormatError("DATA contains nonfinite event values")
    if analysis_start or analysis_end:
        _bounds(analysis_start, analysis_end, len(raw), "ANALYSIS")
        if _overlap((analysis_start, analysis_end), (text_start, text_end)) or (
            expected_bytes and _overlap((analysis_start, analysis_end), (data_start, data_end))
        ):
            raise FCSFormatError("ANALYSIS overlaps another FCS segment")
    return FCSData(
        events=events,
        channels=tuple(channels),
        metadata=metadata,
        ranges=ranges,
        sha256=hashlib.sha256(raw).hexdigest(),
    )
