"""Core converter functionality for GeoTIFF to ASCII formats."""

from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import TextIO

import numpy as np
import rasterio
from rasterio.crs import CRS


class OutputFormat(str, Enum):
    """Supported ASCII output formats."""

    PETREL_POINTS = "petrel_points"  # X Y Z point format
    PETREL_GRID = "petrel_grid"  # CPS-3 grid format (ZMAP compatible)
    ESRI_ASCII = "esri_ascii"  # Standard ESRI ASCII Grid


@dataclass
class GeoTiffInfo:
    """Metadata extracted from a GeoTIFF file."""

    path: Path
    width: int
    height: int
    crs: CRS | None
    transform: rasterio.Affine
    nodata: float | None
    dtype: str
    bounds: tuple[float, float, float, float]  # left, bottom, right, top
    pixel_size: tuple[float, float]  # x, y (note: y is typically negative)

    @property
    def origin(self) -> tuple[float, float]:
        """Return the upper-left corner coordinates."""
        return (self.transform.c, self.transform.f)

    @property
    def cell_size_x(self) -> float:
        """Return the cell size in X direction."""
        return abs(self.transform.a)

    @property
    def cell_size_y(self) -> float:
        """Return the cell size in Y direction."""
        return abs(self.transform.e)


def get_geotiff_info(tif_path: Path) -> GeoTiffInfo:
    """Extract metadata from a GeoTIFF file."""
    with rasterio.open(tif_path) as src:
        return GeoTiffInfo(
            path=tif_path,
            width=src.width,
            height=src.height,
            crs=src.crs,
            transform=src.transform,
            nodata=src.nodata,
            dtype=str(src.dtypes[0]),
            bounds=src.bounds,
            pixel_size=(src.transform.a, src.transform.e),
        )


def read_geotiff_data(tif_path: Path) -> tuple[np.ndarray, GeoTiffInfo]:
    """Read GeoTIFF data and metadata.

    Returns:
        Tuple of (data array, GeoTiffInfo)
    """
    info = get_geotiff_info(tif_path)
    with rasterio.open(tif_path) as src:
        data = src.read(1)  # Read first band
    return data, info


def write_petrel_points(
    data: np.ndarray,
    info: GeoTiffInfo,
    output: TextIO,
    nodata_value: float | None = None,
    skip_nodata: bool = True,
    decimals: int = 3,
) -> int:
    """Write data as Petrel-compatible XYZ points.

    Format: X Y Z (space-delimited, one point per line)

    Args:
        data: 2D numpy array of Z values
        info: GeoTiffInfo with coordinate information
        output: File-like object to write to
        nodata_value: Value to use for nodata in output (None = skip nodata)
        skip_nodata: If True, skip nodata points entirely
        decimals: Number of decimal places for coordinates

    Returns:
        Number of points written
    """
    count = 0
    nodata = info.nodata

    # Get transform parameters
    x_origin = info.transform.c
    y_origin = info.transform.f
    cell_x = info.transform.a
    cell_y = info.transform.e  # Negative for north-up images

    for row in range(info.height):
        for col in range(info.width):
            z = data[row, col]

            # Handle nodata
            if nodata is not None and (z == nodata or np.isnan(z) or np.isinf(z)):
                if skip_nodata:
                    continue
                elif nodata_value is not None:
                    z = nodata_value
                else:
                    continue

            # Calculate center of cell coordinates
            x = x_origin + (col + 0.5) * cell_x
            y = y_origin + (row + 0.5) * cell_y

            output.write(f"{x:.{decimals}f} {y:.{decimals}f} {z:.{decimals}f}\n")
            count += 1

    return count


def write_petrel_grid(
    data: np.ndarray,
    info: GeoTiffInfo,
    output: TextIO,
    nodata_value: float = 0.1e31,
    decimals: int = 3,
) -> None:
    """Write data as Petrel CPS-3 grid format.

    This format is compatible with Petrel's surface import.
    Uses the Petrel-native CPS-3 format with FSASCI header.

    Header format:
        FSASCI 0 1 "COMPUTED" 0 <nodata>
        FSATTR 0 0
        FSLIMI <xmin> <xmax> <ymin> <ymax> <zmin> <zmax>
        FSNROW <nrows> <nrows>
        FSXINC <xinc> <xinc>
        -><label>
        <data>

    Args:
        data: 2D numpy array of Z values
        info: GeoTiffInfo with coordinate information
        output: File-like object to write to
        nodata_value: Value to use for nodata cells (default: 0.1E+31)
        decimals: Number of decimal places for values
    """
    # Make a copy to avoid modifying original
    data = data.copy()

    # Replace nodata values
    if info.nodata is not None:
        mask = (data == info.nodata) | np.isnan(data) | np.isinf(data)
        data = np.where(mask, nodata_value, data)

    # Also handle any remaining NaN/Inf
    mask = np.isnan(data) | np.isinf(data)
    data = np.where(mask, nodata_value, data)

    # Calculate Z min/max from valid data
    valid_data = data[data != nodata_value]
    if len(valid_data) > 0:
        zmin = float(np.min(valid_data))
        zmax = float(np.max(valid_data))
    else:
        zmin = 0.0
        zmax = 0.0

    # Get bounds
    xmin, ymin, xmax, ymax = info.bounds

    # Petrel CPS-3 Header
    # FSASCI: format=0, type=1, name="COMPUTED", unknown=0, nodata
    output.write(f'FSASCI 0 1 "COMPUTED" 0 0.1E+31\n')

    # FSATTR: attribute flags
    output.write("FSATTR 0 0\n")

    # FSLIMI: xmin xmax ymin ymax zmin zmax
    output.write(f"FSLIMI {xmin:.{decimals}f} {xmax:.{decimals}f} {ymin:.{decimals}f} {ymax:.{decimals}f} {zmin:.{decimals}f} {zmax:.{decimals}f}\n")

    # FSNROW: number of columns (ncols) and rows (nrows)
    output.write(f"FSNROW {info.width} {info.height}\n")

    # FSXINC: X increment and Y increment
    output.write(f"FSXINC {info.cell_size_x:.{decimals}f} {info.cell_size_y:.{decimals}f}\n")

    # Label line (arrow indicates start of data)
    output.write(f"->Converted from {info.path.name}\n")

    # Write data (row by row, from top to bottom)
    # Petrel expects 5 values per line
    values_per_line = 5

    for row in range(info.height):
        row_data = data[row, :]
        for i in range(0, len(row_data), values_per_line):
            chunk = row_data[i : i + values_per_line]
            # Format with appropriate precision, using scientific notation for nodata
            formatted = []
            for v in chunk:
                if v == nodata_value or v >= 1e30:
                    formatted.append("0.1E+31")
                else:
                    formatted.append(f"{v:.{decimals}f}")
            line = "  ".join(formatted)
            output.write(line + "\n")


def write_esri_ascii(
    data: np.ndarray,
    info: GeoTiffInfo,
    output: TextIO,
    nodata_value: float = -9999,
    decimals: int = 3,
) -> None:
    """Write data as ESRI ASCII Grid format.

    This is a widely compatible format that Petrel can also import.

    Args:
        data: 2D numpy array of Z values
        info: GeoTiffInfo with coordinate information
        output: File-like object to write to
        nodata_value: Value to use for nodata cells
        decimals: Number of decimal places for values
    """
    # Replace nodata values
    if info.nodata is not None:
        mask = (data == info.nodata) | np.isnan(data) | np.isinf(data)
        data = np.where(mask, nodata_value, data)

    # ESRI ASCII Grid header
    # Note: xllcorner/yllcorner are lower-left corner coordinates
    xll = info.bounds[0]
    yll = info.bounds[1]

    # Check if cells are square
    if abs(info.cell_size_x - info.cell_size_y) < 1e-6:
        output.write(f"ncols         {info.width}\n")
        output.write(f"nrows         {info.height}\n")
        output.write(f"xllcorner     {xll:.{decimals}f}\n")
        output.write(f"yllcorner     {yll:.{decimals}f}\n")
        output.write(f"cellsize      {info.cell_size_x:.{decimals}f}\n")
        output.write(f"NODATA_value  {nodata_value}\n")
    else:
        # Use dx/dy for non-square cells (ArcGIS extension)
        output.write(f"ncols         {info.width}\n")
        output.write(f"nrows         {info.height}\n")
        output.write(f"xllcorner     {xll:.{decimals}f}\n")
        output.write(f"yllcorner     {yll:.{decimals}f}\n")
        output.write(f"dx            {info.cell_size_x:.{decimals}f}\n")
        output.write(f"dy            {info.cell_size_y:.{decimals}f}\n")
        output.write(f"NODATA_value  {nodata_value}\n")

    # Write data (row by row, from top to bottom)
    for row in range(info.height):
        row_data = data[row, :]
        line = " ".join(f"{v:.{decimals}f}" for v in row_data)
        output.write(line + "\n")


def convert_geotiff(
    input_path: Path,
    output_path: Path,
    output_format: OutputFormat = OutputFormat.PETREL_POINTS,
    nodata_value: float | None = None,
    skip_nodata: bool = True,
    decimals: int = 3,
) -> dict:
    """Convert a GeoTIFF file to ASCII format.

    Args:
        input_path: Path to input GeoTIFF
        output_path: Path to output ASCII file
        output_format: Output format type
        nodata_value: Custom nodata value for output
        skip_nodata: Skip nodata values (points format only)
        decimals: Decimal precision for output

    Returns:
        Dictionary with conversion statistics
    """
    data, info = read_geotiff_data(input_path)

    stats = {
        "input_file": str(input_path),
        "output_file": str(output_path),
        "format": output_format.value,
        "width": info.width,
        "height": info.height,
        "total_cells": info.width * info.height,
    }

    with open(output_path, "w") as f:
        if output_format == OutputFormat.PETREL_POINTS:
            default_nodata = None if skip_nodata else -999.25
            count = write_petrel_points(
                data,
                info,
                f,
                nodata_value=nodata_value or default_nodata,
                skip_nodata=skip_nodata,
                decimals=decimals,
            )
            stats["points_written"] = count
        elif output_format == OutputFormat.PETREL_GRID:
            write_petrel_grid(
                data,
                info,
                f,
                nodata_value=nodata_value or -999.25,
                decimals=decimals,
            )
        elif output_format == OutputFormat.ESRI_ASCII:
            write_esri_ascii(
                data,
                info,
                f,
                nodata_value=nodata_value or -9999,
                decimals=decimals,
            )

    return stats
