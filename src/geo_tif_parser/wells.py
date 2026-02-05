"""Well data parsing and Petrel well format output."""

import csv
import re
import shlex
from dataclasses import dataclass, field
from pathlib import Path
from typing import TextIO

from pyproj import Transformer

# US survey foot to meter conversion (exact definition)
US_SURVEY_FT_TO_M = 1200 / 3937  # ≈ 0.3048006096

# Column-to-surface name mapping for the Columbia River Basalt dataset
TOP_COLUMNS: dict[str, str] = {
    "OVBD_TOP": "Overburden",
    "SDMB_TOP": "Saddle_Mountains_Basalt",
    "MBTN_TOP": "Mabton_Interbed",
    "WNB_TOP": "Wanapum_Basalt",
    "VNTG_TOP": "Vantage_Interbed",
    "GRB_TOP": "Grande_Ronde_Basalt",
    "PreM_TOP": "Pre_Miocene",
}

# Surface name for ground/KB elevation in Petrel tops export
SURFACE_ELEVATION_NAME = "Surface"


@dataclass
class WellRecord:
    """A single well with location and formation top picks."""

    name: str
    x: float
    y: float
    surface_elevation: float | None
    tops: dict[str, float] = field(default_factory=dict)


def detect_file_format(path: Path) -> str:
    """Detect whether file is columnar TSV or Petrel tops export format.

    Args:
        path: Path to well data file.

    Returns:
        "columnar" for tab-delimited one-well-per-row format,
        "petrel_tops" for Petrel well tops export format.
    """
    with open(path) as f:
        first_line = f.readline()
        if first_line.startswith("# Petrel well tops"):
            return "petrel_tops"
        # Check for VERSION header
        for line in [first_line] + [f.readline() for _ in range(5)]:
            if line.strip().startswith("VERSION"):
                return "petrel_tops"
    return "columnar"


def _parse_petrel_tops_line(line: str, header_cols: list[str]) -> dict[str, str]:
    """Parse a Petrel tops data line with mixed quoted/unquoted fields.

    Args:
        line: Data line to parse.
        header_cols: List of column names from header.

    Returns:
        Dictionary mapping column names to values.
    """
    # Use shlex to handle quoted strings
    try:
        tokens = shlex.split(line)
    except ValueError:
        # Fallback: simple whitespace split
        tokens = line.split()

    result = {}
    for i, col in enumerate(header_cols):
        if i < len(tokens):
            result[col] = tokens[i]
        else:
            result[col] = ""
    return result


def read_petrel_tops(path: Path, nodata: float = -999) -> list[WellRecord]:
    """Read Petrel well tops export file (VERSION 2).

    Format: One row per well top, with X, Y, Z, Surface, Well columns.
    Groups rows by Well name to build WellRecord objects.

    Args:
        path: Path to Petrel tops export file.
        nodata: Value indicating missing data (default: -999).

    Returns:
        List of WellRecord objects.
    """
    # Accumulator: well_name -> {x, y, surface_elev, tops}
    well_data: dict[str, dict] = {}

    with open(path) as f:
        header_cols: list[str] = []
        in_header = False

        for line in f:
            line = line.strip()
            if not line:
                continue

            # Skip comment lines at start
            if line.startswith("#"):
                continue

            # Skip VERSION line
            if line.startswith("VERSION"):
                continue

            # Header section
            if line == "BEGIN HEADER":
                in_header = True
                continue
            if line == "END HEADER":
                in_header = False
                continue
            if in_header:
                header_cols.append(line)
                continue

            # Data line
            if not header_cols:
                continue

            row = _parse_petrel_tops_line(line, header_cols)

            # Extract key fields
            try:
                x = float(row.get("X", ""))
                y = float(row.get("Y", ""))
                z = float(row.get("Z", ""))
            except ValueError:
                continue

            well_name = row.get("Well", "").strip()
            surface_name = row.get("Surface", "").strip()

            if not well_name:
                continue

            # Skip nodata values
            if z == nodata:
                continue

            # Initialize well entry
            if well_name not in well_data:
                well_data[well_name] = {
                    "x": x,
                    "y": y,
                    "surface_elev": None,
                    "tops": {},
                }

            # Determine if this is surface elevation or a formation top
            if surface_name == SURFACE_ELEVATION_NAME:
                well_data[well_name]["surface_elev"] = z
            else:
                # Normalize surface name (replace spaces with underscores)
                normalized_name = surface_name.replace(" ", "_")
                well_data[well_name]["tops"][normalized_name] = z

    # Convert to WellRecord list
    wells: list[WellRecord] = []
    for name, data in well_data.items():
        wells.append(
            WellRecord(
                name=name.replace(" ", "_"),  # Normalize well name too
                x=data["x"],
                y=data["y"],
                surface_elevation=data["surface_elev"],
                tops=data["tops"],
            )
        )

    return wells


def read_well_data(path: Path, nodata: float = -9999) -> list[WellRecord]:
    """Read tab-delimited well data file (columnar format).

    Expects columns: SEQ_NUM, SOURCE, SOURCE_ID, X, Y, SURF_ELEV,
    plus formation top columns matching TOP_COLUMNS keys.

    Args:
        path: Path to tab-delimited well data file.
        nodata: Value indicating missing data (default: -9999).

    Returns:
        List of WellRecord objects.
    """
    wells: list[WellRecord] = []

    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            # Build well name from SOURCE and SEQ_NUM
            source = row.get("SOURCE", "").strip()
            seq_num = row.get("SEQ_NUM", "").strip()
            name = f"{source}_{seq_num}" if source else seq_num

            # Parse coordinates
            try:
                x = float(row["X"])
                y = float(row["Y"])
            except (KeyError, ValueError):
                continue

            # Parse surface elevation
            surf_elev: float | None = None
            try:
                val = float(row["SURF_ELEV"])
                if val != nodata:
                    surf_elev = val
            except (KeyError, ValueError):
                pass

            # Parse formation tops (elevation values)
            tops: dict[str, float] = {}
            for col, surface_name in TOP_COLUMNS.items():
                try:
                    val = float(row[col])
                    if val != nodata:
                        tops[surface_name] = val
                except (KeyError, ValueError):
                    continue

            wells.append(
                WellRecord(
                    name=name,
                    x=x,
                    y=y,
                    surface_elevation=surf_elev,
                    tops=tops,
                )
            )

    return wells


def convert_elevations_to_meters(wells: list[WellRecord]) -> None:
    """Convert all elevation values from US survey feet to meters in-place.

    Args:
        wells: List of WellRecord objects to convert.
    """
    for w in wells:
        if w.surface_elevation is not None:
            w.surface_elevation *= US_SURVEY_FT_TO_M
        w.tops = {k: v * US_SURVEY_FT_TO_M for k, v in w.tops.items()}


def transform_wells(
    wells: list[WellRecord], src_crs: str, dst_crs: str
) -> None:
    """Batch-transform well coordinates from src_crs to dst_crs in-place.

    Args:
        wells: List of WellRecord objects to transform.
        src_crs: Source CRS as EPSG code string (e.g., "2927").
        dst_crs: Destination CRS as EPSG code string (e.g., "32611").
    """
    transformer = Transformer.from_crs(
        f"EPSG:{src_crs}", f"EPSG:{dst_crs}", always_xy=True
    )
    xs = [w.x for w in wells]
    ys = [w.y for w in wells]
    xs_out, ys_out = transformer.transform(xs, ys)

    for well, x, y in zip(wells, xs_out, ys_out):
        well.x = x
        well.y = y


def write_petrel_headers(
    wells: list[WellRecord], output: TextIO, decimals: int = 2
) -> int:
    """Write Petrel-format well header file (tab-delimited).

    Format: WELL  X  Y  ELEV

    Args:
        wells: List of WellRecord objects.
        output: File-like object to write to.
        decimals: Number of decimal places for coordinates.

    Returns:
        Number of wells written.
    """
    output.write("WELL\tX\tY\tELEV\n")
    count = 0
    for w in wells:
        elev = f"{w.surface_elevation:.{decimals}f}" if w.surface_elevation is not None else ""
        output.write(f"{w.name}\t{w.x:.{decimals}f}\t{w.y:.{decimals}f}\t{elev}\n")
        count += 1
    return count


def write_petrel_tops(
    wells: list[WellRecord], output: TextIO, decimals: int = 1
) -> int:
    """Write Petrel-format well tops file (tab-delimited).

    Format: WELL  SURFACE  ELEV  TVD

    TVD is computed as surface_elevation - top_elevation (depth below
    ground surface). If surface elevation is missing, TVD is left blank.

    Args:
        wells: List of WellRecord objects.
        output: File-like object to write to.
        decimals: Number of decimal places for values.

    Returns:
        Number of top records written.
    """
    output.write("WELL\tSURFACE\tELEV\tTVD\n")
    count = 0
    for w in wells:
        for surface_name, elev in w.tops.items():
            if w.surface_elevation is not None:
                tvd = f"{w.surface_elevation - elev:.{decimals}f}"
            else:
                tvd = ""
            output.write(
                f"{w.name}\t{surface_name}\t{elev:.{decimals}f}\t{tvd}\n"
            )
            count += 1
    return count


def convert_wells(
    input_path: Path,
    output_dir: Path,
    source_crs: str = "2927",
    output_crs: str = "32611",
    nodata: float = -9999,
    decimals: int = 2,
    feet_to_meters: bool = True,
    input_format: str | None = None,
) -> dict:
    """Read well data, transform coordinates, and write Petrel output files.

    Args:
        input_path: Path to well data file.
        output_dir: Directory for output files.
        source_crs: Source EPSG code (default: 2927 = NAD83 / Washington South).
        output_crs: Output EPSG code (default: 32611 = WGS 84 / UTM zone 11N).
        nodata: Nodata sentinel value (default: -9999 for columnar, -999 for petrel_tops).
        decimals: Decimal places for coordinates.
        feet_to_meters: Convert elevation values from US survey feet to meters.
        input_format: "columnar", "petrel_tops", or None to auto-detect.

    Returns:
        Dictionary with conversion statistics.
    """
    # Auto-detect format if not specified
    if input_format is None:
        input_format = detect_file_format(input_path)

    if input_format == "petrel_tops":
        # Petrel tops export typically uses -999 as nodata
        petrel_nodata = nodata if nodata != -9999 else -999
        wells = read_petrel_tops(input_path, nodata=petrel_nodata)
    else:
        wells = read_well_data(input_path, nodata=nodata)

    transform_wells(wells, source_crs, output_crs)

    if feet_to_meters:
        convert_elevations_to_meters(wells)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    headers_path = output_dir / "well_headers.txt"
    tops_path = output_dir / "well_tops.txt"

    with open(headers_path, "w") as f:
        well_count = write_petrel_headers(wells, f, decimals=decimals)

    with open(tops_path, "w") as f:
        tops_count = write_petrel_tops(wells, f, decimals=max(1, decimals))

    # Compute per-surface counts
    surface_counts: dict[str, int] = {}
    for w in wells:
        for surface_name in w.tops:
            surface_counts[surface_name] = surface_counts.get(surface_name, 0) + 1

    # Coordinate ranges
    xs = [w.x for w in wells]
    ys = [w.y for w in wells]

    return {
        "input_file": str(input_path),
        "input_format": input_format,
        "headers_file": str(headers_path),
        "tops_file": str(tops_path),
        "well_count": well_count,
        "tops_count": tops_count,
        "surface_counts": surface_counts,
        "x_min": min(xs) if xs else 0.0,
        "x_max": max(xs) if xs else 0.0,
        "y_min": min(ys) if ys else 0.0,
        "y_max": max(ys) if ys else 0.0,
        "feet_to_meters": feet_to_meters,
    }
