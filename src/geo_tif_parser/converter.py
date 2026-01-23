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
    nodata_value: float = -999.25,
    decimals: int = 3,
) -> None:
    """Write data as Petrel CPS-3 grid format.

    This format is compatible with Petrel's surface import.
    The CPS-3 format is a ZMAP-compatible grid format.

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

    # CPS-3 Header
    # @GRID HEADER
    output.write("@GRID HEADER\n")

    # Grid description line
    output.write(f"!  Grid converted from: {info.path.name}\n")
    if info.crs:
        output.write(f"!  CRS: {info.crs.to_string()}\n")

    # Format specification
    # FSASCI (Field Size ASCII), number of values per line, nodata value
    values_per_line = 6
    output.write(f"@ FSASCI,  {values_per_line}, {nodata_value}, , 4, 1\n")

    # Grid dimensions and bounds
    # rows, cols, xmin, xmax, ymin, ymax
    xmin, ymin, xmax, ymax = info.bounds
    output.write(f"@ {info.height}, {info.width}, {xmin:.{decimals}f}, {xmax:.{decimals}f}, {ymin:.{decimals}f}, {ymax:.{decimals}f}\n")

    # Cell sizes
    output.write(f"@ {info.cell_size_x:.{decimals}f}, {info.cell_size_y:.{decimals}f}\n")

    # End header marker
    output.write("@\n")

    # Write data (row by row, from top to bottom)
    # Petrel expects values in column-major order for CPS-3
    for row in range(info.height):
        row_data = data[row, :]
        # Write values_per_line values per line
        for i in range(0, len(row_data), values_per_line):
            chunk = row_data[i : i + values_per_line]
            line = " ".join(f"{v:.{decimals}f}" for v in chunk)
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
