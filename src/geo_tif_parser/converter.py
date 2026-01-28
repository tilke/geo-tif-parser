"""Core converter functionality for GeoTIFF to ASCII formats."""

from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import TextIO

import numpy as np
import rasterio
from rasterio.crs import CRS
from rasterio.warp import (
    Resampling,
    calculate_default_transform,
    reproject,
    transform,
    transform_bounds,
)


# Mapping of resampling method names to Resampling enum values
RESAMPLING_METHODS = {
    "nearest": Resampling.nearest,
    "bilinear": Resampling.bilinear,
    "cubic": Resampling.cubic,
    "lanczos": Resampling.lanczos,
}


def parse_crs(crs_string: str) -> CRS:
    """Parse CRS from EPSG code (26910), EPSG:26910, or proj4 string.

    Args:
        crs_string: CRS specification string

    Returns:
        Parsed CRS object

    Raises:
        ValueError: If CRS cannot be parsed
    """
    crs_string = crs_string.strip()

    # Try numeric EPSG code first (e.g., "26910")
    if crs_string.isdigit():
        return CRS.from_epsg(int(crs_string))

    # Try EPSG:XXXXX format (e.g., "EPSG:26910")
    if crs_string.upper().startswith("EPSG:"):
        try:
            epsg_code = int(crs_string.split(":")[1])
            return CRS.from_epsg(epsg_code)
        except (ValueError, IndexError) as e:
            raise ValueError(f"Invalid EPSG format: {crs_string}") from e

    # Try as generic CRS string (WKT, proj4, etc.)
    try:
        return CRS.from_string(crs_string)
    except Exception as e:
        raise ValueError(f"Cannot parse CRS: {crs_string}") from e


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


def reproject_raster(
    data: np.ndarray,
    info: GeoTiffInfo,
    dst_crs: CRS,
    resampling: Resampling = Resampling.bilinear,
) -> tuple[np.ndarray, GeoTiffInfo]:
    """Reproject raster data to a new CRS.

    Args:
        data: 2D numpy array of values
        info: GeoTiffInfo with source coordinate information
        dst_crs: Destination CRS
        resampling: Resampling method (default: bilinear)

    Returns:
        Tuple of (reprojected data array, updated GeoTiffInfo)
    """
    # Calculate the transform and dimensions for the destination CRS
    dst_transform, dst_width, dst_height = calculate_default_transform(
        info.crs,
        dst_crs,
        info.width,
        info.height,
        *info.bounds,
    )

    # Create destination array
    dst_data = np.zeros((dst_height, dst_width), dtype=data.dtype)

    # Perform the reprojection
    reproject(
        source=data,
        destination=dst_data,
        src_transform=info.transform,
        src_crs=info.crs,
        dst_transform=dst_transform,
        dst_crs=dst_crs,
        resampling=resampling,
        src_nodata=info.nodata,
        dst_nodata=info.nodata,
    )

    # Calculate new bounds from the transform and dimensions
    dst_bounds = (
        dst_transform.c,  # left (xmin)
        dst_transform.f + dst_height * dst_transform.e,  # bottom (ymin)
        dst_transform.c + dst_width * dst_transform.a,  # right (xmax)
        dst_transform.f,  # top (ymax)
    )

    # Create updated GeoTiffInfo
    new_info = GeoTiffInfo(
        path=info.path,
        width=dst_width,
        height=dst_height,
        crs=dst_crs,
        transform=dst_transform,
        nodata=info.nodata,
        dtype=info.dtype,
        bounds=dst_bounds,
        pixel_size=(dst_transform.a, dst_transform.e),
    )

    return dst_data, new_info


def write_petrel_points(
    data: np.ndarray,
    info: GeoTiffInfo,
    output: TextIO,
    nodata_value: float | None = None,
    skip_nodata: bool = True,
    decimals: int = 3,
    dst_crs: CRS | None = None,
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
        dst_crs: Optional destination CRS for coordinate transformation

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

    # Collect all valid points first for batch transformation
    xs = []
    ys = []
    zs = []

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

            xs.append(x)
            ys.append(y)
            zs.append(z)

    # Transform coordinates if dst_crs is specified
    if dst_crs is not None and info.crs is not None:
        xs_out, ys_out = transform(info.crs, dst_crs, xs, ys)
    else:
        xs_out, ys_out = xs, ys

    # Write points
    for x, y, z in zip(xs_out, ys_out, zs):
        output.write(f"{x:.{decimals}f} {y:.{decimals}f} {z:.{decimals}f}\n")
        count += 1

    return count


def write_petrel_grid(
    data: np.ndarray,
    info: GeoTiffInfo,
    output: TextIO,
    nodata_value: float = 0.1e31,
    decimals: int = 3,
    flip_z: bool = False,
    dst_crs: CRS | None = None,
) -> None:
    """Write data as ZMAP+ ASCII grid format.

    This format is compatible with Petrel's CPS-3 grid surface import.
    Uses the standard ZMAP+ header with @ delimiters.

    Header format:
        ! <comment>
        @<name>, GRID, <nodes_per_line>
        <field_width>, <null_value>, , <decimal_places>, <start_col>
        <nrows>, <ncols>, <xmin>, <xmax>, <ymin>, <ymax>
        0.0, 0.0, 0.0
        @
        <data: column-major, N→S per column, columns W→E>

    Args:
        data: 2D numpy array of Z values
        info: GeoTiffInfo with coordinate information
        output: File-like object to write to
        nodata_value: Value to use for nodata cells (default: 0.1E+31)
        decimals: Number of decimal places for values
        flip_z: If True, indicates Z values were negated (for label comment)
        dst_crs: Optional destination CRS for coordinate transformation (bounds-only)
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

    # Get bounds and transform if dst_crs specified
    xmin, ymin, xmax, ymax = info.bounds
    cell_size_x = info.cell_size_x
    cell_size_y = info.cell_size_y

    if dst_crs is not None and info.crs is not None:
        # Transform bounds to destination CRS
        xmin, ymin, xmax, ymax = transform_bounds(
            info.crs, dst_crs, xmin, ymin, xmax, ymax
        )
        # Recalculate cell sizes for transformed bounds
        cell_size_x = (xmax - xmin) / info.width
        cell_size_y = (ymax - ymin) / info.height

    # ZMAP+ Header
    nodes_per_line = 5
    field_width = 15

    # Comment line
    if flip_z:
        output.write(f"! Converted from {info.path.name} | Z: Depth (positive down)\n")
    else:
        output.write(f"! Converted from {info.path.name} | Z: Elevation (positive up)\n")

    # Grid definition header (between @ markers)
    output.write(f"@{info.path.stem}, GRID, {nodes_per_line}\n")
    output.write(f"{field_width}, 0.1000000E+31, , {decimals}, 1\n")
    output.write(
        f"{info.height}, {info.width}, "
        f"{xmin:.{decimals}f}, {xmax:.{decimals}f}, "
        f"{ymin:.{decimals}f}, {ymax:.{decimals}f}\n"
    )
    output.write("0.0, 0.0, 0.0\n")
    output.write("@\n")

    # Write data column by column (ZMAP+ column-major order)
    # Each column: N→S (top to bottom of raster), columns: W→E (left to right)
    for col in range(info.width):
        col_data = data[:, col]  # top to bottom (N→S)
        for i in range(0, len(col_data), nodes_per_line):
            chunk = col_data[i : i + nodes_per_line]
            formatted = []
            for v in chunk:
                if v == nodata_value or v >= 1e30:
                    formatted.append(f"{'0.1000000E+31':>{field_width}}")
                else:
                    formatted.append(f"{v:{field_width}.{decimals}f}")
            output.write("".join(formatted) + "\n")


def write_esri_ascii(
    data: np.ndarray,
    info: GeoTiffInfo,
    output: TextIO,
    nodata_value: float = -9999,
    decimals: int = 3,
    dst_crs: CRS | None = None,
) -> None:
    """Write data as ESRI ASCII Grid format.

    This is a widely compatible format that Petrel can also import.

    Args:
        data: 2D numpy array of Z values
        info: GeoTiffInfo with coordinate information
        output: File-like object to write to
        nodata_value: Value to use for nodata cells
        decimals: Number of decimal places for values
        dst_crs: Optional destination CRS for coordinate transformation (bounds-only)
    """
    # Replace nodata values
    if info.nodata is not None:
        mask = (data == info.nodata) | np.isnan(data) | np.isinf(data)
        data = np.where(mask, nodata_value, data)

    # ESRI ASCII Grid header
    # Note: xllcorner/yllcorner are lower-left corner coordinates
    xll = info.bounds[0]
    yll = info.bounds[1]
    xur = info.bounds[2]
    yur = info.bounds[3]
    cell_size_x = info.cell_size_x
    cell_size_y = info.cell_size_y

    # Transform bounds if dst_crs specified
    if dst_crs is not None and info.crs is not None:
        xll, yll, xur, yur = transform_bounds(
            info.crs, dst_crs, xll, yll, xur, yur
        )
        # Recalculate cell sizes for transformed bounds
        cell_size_x = (xur - xll) / info.width
        cell_size_y = (yur - yll) / info.height

    # Check if cells are square
    if abs(cell_size_x - cell_size_y) < 1e-6:
        output.write(f"ncols         {info.width}\n")
        output.write(f"nrows         {info.height}\n")
        output.write(f"xllcorner     {xll:.{decimals}f}\n")
        output.write(f"yllcorner     {yll:.{decimals}f}\n")
        output.write(f"cellsize      {cell_size_x:.{decimals}f}\n")
        output.write(f"NODATA_value  {nodata_value}\n")
    else:
        # Use dx/dy for non-square cells (ArcGIS extension)
        output.write(f"ncols         {info.width}\n")
        output.write(f"nrows         {info.height}\n")
        output.write(f"xllcorner     {xll:.{decimals}f}\n")
        output.write(f"yllcorner     {yll:.{decimals}f}\n")
        output.write(f"dx            {cell_size_x:.{decimals}f}\n")
        output.write(f"dy            {cell_size_y:.{decimals}f}\n")
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
    flip_z: bool = False,
    output_crs: str | None = None,
    reproject_mode: str = "bounds-only",
    resampling: str = "bilinear",
) -> dict:
    """Convert a GeoTIFF file to ASCII format.

    Args:
        input_path: Path to input GeoTIFF
        output_path: Path to output ASCII file
        output_format: Output format type
        nodata_value: Custom nodata value for output
        skip_nodata: Skip nodata values (points format only)
        decimals: Decimal precision for output
        flip_z: If True, negate Z values (convert elevation to depth or vice versa)
        output_crs: Output CRS (e.g., "26910", "EPSG:26910"). If None, use source CRS
        reproject_mode: "bounds-only" (transform coordinates only) or "full" (reproject raster)
        resampling: Resampling method for full mode ("nearest", "bilinear", "cubic", "lanczos")

    Returns:
        Dictionary with conversion statistics
    """
    data, info = read_geotiff_data(input_path)

    # Parse and validate output CRS
    dst_crs = None
    if output_crs is not None:
        if info.crs is None:
            raise ValueError(
                f"Cannot transform coordinates: source file '{input_path.name}' has no CRS defined"
            )
        dst_crs = parse_crs(output_crs)

        # Check if source and destination CRS are the same
        if info.crs == dst_crs:
            dst_crs = None  # No transformation needed

    # Handle full reprojection mode
    if dst_crs is not None and reproject_mode == "full":
        resampling_method = RESAMPLING_METHODS.get(resampling, Resampling.bilinear)
        data, info = reproject_raster(data, info, dst_crs, resampling_method)
        # After full reprojection, no need for bounds-only transform
        dst_crs = None

    # Flip Z values if requested (negate non-nodata values)
    if flip_z:
        if info.nodata is not None:
            mask = (data != info.nodata) & ~np.isnan(data) & ~np.isinf(data)
            data = np.where(mask, -data, data)
        else:
            mask = ~np.isnan(data) & ~np.isinf(data)
            data = np.where(mask, -data, data)

    stats = {
        "input_file": str(input_path),
        "output_file": str(output_path),
        "format": output_format.value,
        "width": info.width,
        "height": info.height,
        "total_cells": info.width * info.height,
    }

    # Include CRS info in stats
    if output_crs is not None:
        stats["output_crs"] = output_crs
        stats["reproject_mode"] = reproject_mode

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
                dst_crs=dst_crs,
            )
            stats["points_written"] = count
        elif output_format == OutputFormat.PETREL_GRID:
            write_petrel_grid(
                data,
                info,
                f,
                nodata_value=nodata_value or 0.1e31,
                decimals=decimals,
                flip_z=flip_z,
                dst_crs=dst_crs,
            )
        elif output_format == OutputFormat.ESRI_ASCII:
            write_esri_ascii(
                data,
                info,
                f,
                nodata_value=nodata_value or -9999,
                decimals=decimals,
                dst_crs=dst_crs,
            )

    return stats
