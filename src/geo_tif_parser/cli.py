"""Command-line interface for GeoTIFF to ASCII converter."""

from pathlib import Path
from typing import Annotated, Optional

import typer
from rich.console import Console
from rich.table import Table

import numpy as np
import rasterio

from rasterio.warp import transform_bounds

from .converter import (
    OutputFormat,
    convert_geotiff,
    get_geotiff_info,
    parse_crs,
)
from .wells import convert_wells

app = typer.Typer(
    name="geotif",
    help="Convert GeoTIFF files to Petrel-compatible ASCII formats.",
    add_completion=False,
)
console = Console()


@app.command()
def info(
    input_file: Annotated[
        Path,
        typer.Argument(
            help="Path to GeoTIFF file",
            exists=True,
            dir_okay=False,
            readable=True,
        ),
    ],
) -> None:
    """Display information about a GeoTIFF file."""
    try:
        tif_info = get_geotiff_info(input_file)

        table = Table(title=f"GeoTIFF Info: {input_file.name}")
        table.add_column("Property", style="cyan")
        table.add_column("Value", style="green")

        table.add_row("File", str(tif_info.path))
        table.add_row("Dimensions", f"{tif_info.width} x {tif_info.height}")
        table.add_row("Total Cells", f"{tif_info.width * tif_info.height:,}")
        table.add_row("Data Type", tif_info.dtype)
        table.add_row("NoData Value", str(tif_info.nodata))
        table.add_row("Cell Size X", f"{tif_info.cell_size_x:.3f}")
        table.add_row("Cell Size Y", f"{tif_info.cell_size_y:.3f}")
        table.add_row("Origin (UL)", f"({tif_info.origin[0]:.3f}, {tif_info.origin[1]:.3f})")

        bounds = tif_info.bounds
        table.add_row("Bounds (L,B,R,T)", f"({bounds[0]:.3f}, {bounds[1]:.3f}, {bounds[2]:.3f}, {bounds[3]:.3f})")

        if tif_info.crs:
            table.add_row("CRS", tif_info.crs.to_string())

        console.print(table)

    except Exception as e:
        console.print(f"[red]Error reading GeoTIFF: {e}[/red]")
        raise typer.Exit(1)


@app.command()
def header(
    input_file: Annotated[
        Path,
        typer.Argument(
            help="Path to GeoTIFF file",
            exists=True,
            dir_okay=False,
            readable=True,
        ),
    ],
    output_crs: Annotated[
        Optional[str],
        typer.Option(
            "--output-crs",
            "-c",
            help="Preview transformed bounds in this CRS (e.g., 26910, EPSG:26910)",
        ),
    ] = None,
) -> None:
    """Display detailed header information including CRS for Petrel import."""
    try:
        tif_info = get_geotiff_info(input_file)

        # Read data to get Z min/max
        with rasterio.open(input_file) as src:
            data = src.read(1)

        # Calculate Z stats from valid data
        nodata = tif_info.nodata
        if nodata is not None:
            valid_mask = (data != nodata) & ~np.isnan(data) & ~np.isinf(data)
            valid_data = data[valid_mask]
        else:
            valid_data = data[~np.isnan(data) & ~np.isinf(data)]

        if len(valid_data) > 0:
            zmin = float(np.min(valid_data))
            zmax = float(np.max(valid_data))
            zmean = float(np.mean(valid_data))
        else:
            zmin = zmax = zmean = 0.0

        # CRS Information
        console.print("\n[bold cyan]Coordinate Reference System[/bold cyan]")
        console.print("─" * 60)

        if tif_info.crs:
            crs = tif_info.crs

            # Get CRS name
            if crs.is_projected:
                console.print(f"[green]Type:[/green] Projected CRS")
            elif crs.is_geographic:
                console.print(f"[green]Type:[/green] Geographic CRS")

            # Try to get readable name
            if hasattr(crs, 'name') and crs.name:
                console.print(f"[green]Name:[/green] {crs.name}")

            # EPSG code if available
            epsg = crs.to_epsg()
            if epsg:
                console.print(f"[green]EPSG:[/green] {epsg}")

            # Units
            if crs.linear_units:
                console.print(f"[green]Linear Units:[/green] {crs.linear_units}")

            # WKT (truncated)
            try:
                wkt = crs.to_wkt()
            except Exception:
                wkt = str(crs)
            console.print(f"\n[green]WKT:[/green]")
            # Show first portion
            if len(wkt) > 500:
                console.print(f"  {wkt[:500]}...")
            else:
                console.print(f"  {wkt}")
        else:
            console.print("[yellow]No CRS defined[/yellow]")

        # Grid Information
        console.print("\n[bold cyan]Grid Geometry[/bold cyan]")
        console.print("─" * 60)

        bounds = tif_info.bounds
        table = Table(show_header=False, box=None)
        table.add_column("Property", style="green", width=20)
        table.add_column("Value")

        table.add_row("Columns (X)", str(tif_info.width))
        table.add_row("Rows (Y)", str(tif_info.height))
        table.add_row("Cell Size X", f"{tif_info.cell_size_x:.6f}")
        table.add_row("Cell Size Y", f"{tif_info.cell_size_y:.6f}")
        table.add_row("X Min", f"{bounds[0]:.6f}")
        table.add_row("X Max", f"{bounds[2]:.6f}")
        table.add_row("Y Min", f"{bounds[1]:.6f}")
        table.add_row("Y Max", f"{bounds[3]:.6f}")
        table.add_row("Z Min", f"{zmin:.3f}")
        table.add_row("Z Max", f"{zmax:.3f}")
        table.add_row("Z Mean", f"{zmean:.3f}")
        table.add_row("NoData Value", str(tif_info.nodata))
        table.add_row("Valid Cells", f"{len(valid_data):,} of {tif_info.width * tif_info.height:,}")

        console.print(table)

        # Show transformed bounds if output_crs specified
        if output_crs is not None:
            if tif_info.crs is None:
                console.print("\n[yellow]Cannot transform bounds: source file has no CRS defined[/yellow]")
            else:
                try:
                    dst_crs = parse_crs(output_crs)
                    xmin, ymin, xmax, ymax = transform_bounds(
                        tif_info.crs, dst_crs, bounds[0], bounds[1], bounds[2], bounds[3]
                    )
                    cell_size_x = (xmax - xmin) / tif_info.width
                    cell_size_y = (ymax - ymin) / tif_info.height

                    console.print(f"\n[bold cyan]Transformed Bounds ({output_crs})[/bold cyan]")
                    console.print("─" * 60)

                    # Get destination CRS info
                    dst_epsg = dst_crs.to_epsg()
                    dst_name = dst_crs.name if hasattr(dst_crs, 'name') else None
                    if dst_name:
                        console.print(f"[green]Target CRS:[/green] {dst_name}")
                    if dst_epsg:
                        console.print(f"[green]Target EPSG:[/green] {dst_epsg}")
                    if dst_crs.linear_units:
                        console.print(f"[green]Target Units:[/green] {dst_crs.linear_units}")

                    trans_table = Table(show_header=False, box=None)
                    trans_table.add_column("Property", style="green", width=20)
                    trans_table.add_column("Value")
                    trans_table.add_row("X Min", f"{xmin:.6f}")
                    trans_table.add_row("X Max", f"{xmax:.6f}")
                    trans_table.add_row("Y Min", f"{ymin:.6f}")
                    trans_table.add_row("Y Max", f"{ymax:.6f}")
                    trans_table.add_row("Cell Size X", f"{cell_size_x:.6f}")
                    trans_table.add_row("Cell Size Y", f"{cell_size_y:.6f}")
                    console.print(trans_table)
                except Exception as e:
                    console.print(f"\n[red]Error transforming bounds: {e}[/red]")

        # ZMAP+ Header Preview
        console.print("\n[bold cyan]ZMAP+ Header Preview[/bold cyan]")
        console.print("─" * 60)

        # Use transformed bounds for preview if output_crs specified
        preview_bounds = bounds
        preview_cell_x = tif_info.cell_size_x
        preview_cell_y = tif_info.cell_size_y

        if output_crs is not None and tif_info.crs is not None:
            try:
                dst_crs = parse_crs(output_crs)
                preview_bounds = transform_bounds(
                    tif_info.crs, dst_crs, bounds[0], bounds[1], bounds[2], bounds[3]
                )
                preview_cell_x = (preview_bounds[2] - preview_bounds[0]) / tif_info.width
                preview_cell_y = (preview_bounds[3] - preview_bounds[1]) / tif_info.height
                console.print(f"[dim](Using transformed bounds for CRS {output_crs})[/dim]")
            except Exception:
                pass  # Fall back to original bounds

        console.print(f"! Converted from {input_file.name}")
        console.print(f"@{input_file.stem}, GRID, 5")
        console.print(f"15, 0.1000000E+31, , 3, 1")
        console.print(
            f"{tif_info.height}, {tif_info.width}, "
            f"{preview_bounds[0]:.3f}, {preview_bounds[2]:.3f}, "
            f"{preview_bounds[1]:.3f}, {preview_bounds[3]:.3f}"
        )
        console.print(f"0.0, 0.0, 0.0")
        console.print(f"@")
        console.print()

    except Exception as e:
        console.print(f"[red]Error reading GeoTIFF: {e}[/red]")
        raise typer.Exit(1)


@app.command()
def convert(
    input_file: Annotated[
        Path,
        typer.Argument(
            help="Path to GeoTIFF file",
            exists=True,
            dir_okay=False,
            readable=True,
        ),
    ],
    output_file: Annotated[
        Optional[Path],
        typer.Option(
            "-o",
            "--output",
            help="Output file path. Default: input name with new extension",
        ),
    ] = None,
    format: Annotated[
        OutputFormat,
        typer.Option(
            "-f",
            "--format",
            help="Output format",
            case_sensitive=False,
        ),
    ] = OutputFormat.PETREL_POINTS,
    nodata: Annotated[
        Optional[float],
        typer.Option(
            "--nodata",
            help="NoData value for output file",
        ),
    ] = None,
    skip_nodata: Annotated[
        bool,
        typer.Option(
            "--skip-nodata/--include-nodata",
            help="Skip or include nodata values (points format only)",
        ),
    ] = True,
    decimals: Annotated[
        int,
        typer.Option(
            "-d",
            "--decimals",
            help="Number of decimal places",
            min=0,
            max=10,
        ),
    ] = 3,
    flip_z: Annotated[
        bool,
        typer.Option(
            "--flip-z",
            help="Negate Z values (convert elevation to depth or vice versa)",
        ),
    ] = False,
    output_crs: Annotated[
        Optional[str],
        typer.Option(
            "--output-crs",
            "-c",
            help="Output CRS for coordinate transformation (e.g., 26910, EPSG:26910)",
        ),
    ] = None,
    reproject_mode: Annotated[
        str,
        typer.Option(
            "--reproject-mode",
            help="Reprojection mode: 'bounds-only' (fast, default) or 'full' (resamples raster)",
        ),
    ] = "bounds-only",
    resampling: Annotated[
        str,
        typer.Option(
            "--resampling",
            help="Resampling method for full mode: nearest, bilinear (default), cubic, lanczos",
        ),
    ] = "bilinear",
) -> None:
    """Convert a GeoTIFF file to ASCII format for Petrel import."""
    # Validate reproject_mode
    if reproject_mode not in ("bounds-only", "full"):
        console.print(f"[red]Invalid reproject-mode: {reproject_mode}. Use 'bounds-only' or 'full'.[/red]")
        raise typer.Exit(1)

    # Validate resampling method
    valid_resampling = ("nearest", "bilinear", "cubic", "lanczos")
    if resampling not in valid_resampling:
        console.print(f"[red]Invalid resampling method: {resampling}. Use one of: {', '.join(valid_resampling)}[/red]")
        raise typer.Exit(1)

    # Determine output path
    if output_file is None:
        suffix_map = {
            OutputFormat.PETREL_POINTS: ".xyz",
            OutputFormat.PETREL_GRID: ".cps3",
            OutputFormat.ESRI_ASCII: ".asc",
        }
        output_file = input_file.with_suffix(suffix_map[format])

    console.print(f"[cyan]Converting:[/cyan] {input_file.name}")
    console.print(f"[cyan]Format:[/cyan] {format.value}")
    if flip_z:
        console.print(f"[cyan]Flip Z:[/cyan] Yes (negating Z values)")
    if output_crs:
        console.print(f"[cyan]Output CRS:[/cyan] {output_crs}")
        console.print(f"[cyan]Reproject Mode:[/cyan] {reproject_mode}")
        if reproject_mode == "full":
            console.print(f"[cyan]Resampling:[/cyan] {resampling}")
    console.print(f"[cyan]Output:[/cyan] {output_file}")

    try:
        with console.status("[bold green]Processing..."):
            stats = convert_geotiff(
                input_path=input_file,
                output_path=output_file,
                output_format=format,
                nodata_value=nodata,
                skip_nodata=skip_nodata,
                decimals=decimals,
                flip_z=flip_z,
                output_crs=output_crs,
                reproject_mode=reproject_mode,
                resampling=resampling,
            )

        console.print("[green]✓ Conversion complete![/green]")

        # Show statistics
        if "points_written" in stats:
            console.print(f"  Points written: {stats['points_written']:,}")
        console.print(f"  Grid size: {stats['width']} x {stats['height']} ({stats['total_cells']:,} cells)")
        if "output_crs" in stats:
            console.print(f"  CRS transformed: {stats['output_crs']} ({stats['reproject_mode']} mode)")

    except Exception as e:
        console.print(f"[red]Error during conversion: {e}[/red]")
        raise typer.Exit(1)


@app.command()
def batch(
    input_dir: Annotated[
        Path,
        typer.Argument(
            help="Directory containing GeoTIFF files",
            exists=True,
            file_okay=False,
            readable=True,
        ),
    ],
    output_dir: Annotated[
        Optional[Path],
        typer.Option(
            "-o",
            "--output",
            help="Output directory. Default: same as input",
        ),
    ] = None,
    format: Annotated[
        OutputFormat,
        typer.Option(
            "-f",
            "--format",
            help="Output format",
            case_sensitive=False,
        ),
    ] = OutputFormat.PETREL_POINTS,
    pattern: Annotated[
        str,
        typer.Option(
            "-p",
            "--pattern",
            help="Glob pattern for input files",
        ),
    ] = "*.tif",
    nodata: Annotated[
        Optional[float],
        typer.Option(
            "--nodata",
            help="NoData value for output files",
        ),
    ] = None,
    skip_nodata: Annotated[
        bool,
        typer.Option(
            "--skip-nodata/--include-nodata",
            help="Skip or include nodata values (points format only)",
        ),
    ] = True,
    decimals: Annotated[
        int,
        typer.Option(
            "-d",
            "--decimals",
            help="Number of decimal places",
            min=0,
            max=10,
        ),
    ] = 3,
    flip_z: Annotated[
        bool,
        typer.Option(
            "--flip-z",
            help="Negate Z values (convert elevation to depth or vice versa)",
        ),
    ] = False,
    output_crs: Annotated[
        Optional[str],
        typer.Option(
            "--output-crs",
            "-c",
            help="Output CRS for coordinate transformation (e.g., 26910, EPSG:26910)",
        ),
    ] = None,
    reproject_mode: Annotated[
        str,
        typer.Option(
            "--reproject-mode",
            help="Reprojection mode: 'bounds-only' (fast, default) or 'full' (resamples raster)",
        ),
    ] = "bounds-only",
    resampling: Annotated[
        str,
        typer.Option(
            "--resampling",
            help="Resampling method for full mode: nearest, bilinear (default), cubic, lanczos",
        ),
    ] = "bilinear",
) -> None:
    """Batch convert multiple GeoTIFF files in a directory."""
    # Validate reproject_mode
    if reproject_mode not in ("bounds-only", "full"):
        console.print(f"[red]Invalid reproject-mode: {reproject_mode}. Use 'bounds-only' or 'full'.[/red]")
        raise typer.Exit(1)

    # Validate resampling method
    valid_resampling = ("nearest", "bilinear", "cubic", "lanczos")
    if resampling not in valid_resampling:
        console.print(f"[red]Invalid resampling method: {resampling}. Use one of: {', '.join(valid_resampling)}[/red]")
        raise typer.Exit(1)

    # Set output directory
    if output_dir is None:
        output_dir = input_dir
    else:
        output_dir.mkdir(parents=True, exist_ok=True)

    # Find input files
    input_files = sorted(input_dir.glob(pattern))

    if not input_files:
        console.print(f"[yellow]No files matching '{pattern}' found in {input_dir}[/yellow]")
        raise typer.Exit(0)

    console.print(f"[cyan]Found {len(input_files)} files to convert[/cyan]")
    if flip_z:
        console.print(f"[cyan]Flip Z:[/cyan] Yes (negating Z values)")
    if output_crs:
        console.print(f"[cyan]Output CRS:[/cyan] {output_crs} ({reproject_mode} mode)")

    suffix_map = {
        OutputFormat.PETREL_POINTS: ".xyz",
        OutputFormat.PETREL_GRID: ".cps3",
        OutputFormat.ESRI_ASCII: ".asc",
    }

    success_count = 0
    error_count = 0

    for input_file in input_files:
        output_file = output_dir / (input_file.stem + suffix_map[format])

        try:
            console.print(f"  [cyan]{input_file.name}[/cyan] → {output_file.name}")
            convert_geotiff(
                input_path=input_file,
                output_path=output_file,
                output_format=format,
                flip_z=flip_z,
                nodata_value=nodata,
                skip_nodata=skip_nodata,
                decimals=decimals,
                output_crs=output_crs,
                reproject_mode=reproject_mode,
                resampling=resampling,
            )
            success_count += 1
        except Exception as e:
            console.print(f"    [red]Error: {e}[/red]")
            error_count += 1

    console.print()
    console.print(f"[green]✓ Converted: {success_count}[/green]")
    if error_count:
        console.print(f"[red]✗ Errors: {error_count}[/red]")


@app.command()
def wells(
    input_file: Annotated[
        Path,
        typer.Argument(
            help="Path to well data file (auto-detects format)",
            exists=True,
            dir_okay=False,
            readable=True,
        ),
    ],
    output_dir: Annotated[
        Optional[Path],
        typer.Option(
            "-o",
            "--output",
            help="Output directory. Default: same directory as input file",
        ),
    ] = None,
    input_format: Annotated[
        Optional[str],
        typer.Option(
            "--input-format",
            "-i",
            help="Input format: 'columnar' (TSV one-well-per-row) or 'petrel_tops' (Petrel export). Auto-detects if not specified.",
        ),
    ] = None,
    source_crs: Annotated[
        str,
        typer.Option(
            "--source-crs",
            "-s",
            help="Source EPSG code (e.g., 2927 for NAD83 / Washington South)",
        ),
    ] = "2927",
    output_crs: Annotated[
        str,
        typer.Option(
            "--output-crs",
            "-c",
            help="Output EPSG code (e.g., 32611 for WGS 84 / UTM zone 11N)",
        ),
    ] = "32611",
    nodata: Annotated[
        float,
        typer.Option(
            "--nodata",
            help="NoData sentinel value in input file (default: -9999 for columnar, -999 for petrel_tops)",
        ),
    ] = -9999,
    decimals: Annotated[
        int,
        typer.Option(
            "-d",
            "--decimals",
            help="Number of decimal places for coordinates",
            min=0,
            max=10,
        ),
    ] = 2,
    feet_to_meters: Annotated[
        bool,
        typer.Option(
            "--feet-to-meters/--no-feet-to-meters",
            help="Convert elevation values from US survey feet to meters",
        ),
    ] = True,
) -> None:
    """Parse well data, transform coordinates, and write Petrel well header/tops files."""
    if output_dir is None:
        output_dir = input_file.parent

    # Validate input format if specified
    if input_format is not None and input_format not in ("columnar", "petrel_tops"):
        console.print(f"[red]Invalid input format: {input_format}. Use 'columnar' or 'petrel_tops'.[/red]")
        raise typer.Exit(1)

    console.print(f"[cyan]Input:[/cyan] {input_file}")
    console.print(f"[cyan]Output dir:[/cyan] {output_dir}")
    console.print(f"[cyan]CRS:[/cyan] EPSG:{source_crs} → EPSG:{output_crs}")
    if feet_to_meters:
        console.print(f"[cyan]Z units:[/cyan] US survey ft → m")

    try:
        with console.status("[bold green]Processing wells..."):
            stats = convert_wells(
                input_path=input_file,
                output_dir=output_dir,
                source_crs=source_crs,
                output_crs=output_crs,
                nodata=nodata,
                decimals=decimals,
                feet_to_meters=feet_to_meters,
                input_format=input_format,
            )

        console.print("[green]Done.[/green]")

        # Summary table
        z_unit = "m" if feet_to_meters else "ft"
        table = Table(title="Well Data Summary")
        table.add_column("Metric", style="cyan")
        table.add_column("Value", style="green")

        table.add_row("Input format", stats.get("input_format", "columnar"))
        table.add_row("Wells", f"{stats['well_count']:,}")
        table.add_row("Top records", f"{stats['tops_count']:,}")
        table.add_row("X range (m)", f"{stats['x_min']:.2f} – {stats['x_max']:.2f}")
        table.add_row("Y range (m)", f"{stats['y_min']:.2f} – {stats['y_max']:.2f}")
        table.add_row("Z units", z_unit)

        console.print(table)

        # Per-surface breakdown
        surface_counts = stats.get("surface_counts", {})
        if surface_counts:
            st = Table(title="Tops per Surface")
            st.add_column("Surface", style="cyan")
            st.add_column("Count", style="green", justify="right")
            for name, count in sorted(surface_counts.items(), key=lambda x: -x[1]):
                st.add_row(name, f"{count:,}")
            console.print(st)

        console.print(f"\n[dim]Headers:[/dim] {stats['headers_file']}")
        console.print(f"[dim]Tops:[/dim]    {stats['tops_file']}")

    except Exception as e:
        console.print(f"[red]Error processing wells: {e}[/red]")
        raise typer.Exit(1)


if __name__ == "__main__":
    app()
