"""Command-line interface for GeoTIFF to ASCII converter."""

from pathlib import Path
from typing import Annotated, Optional

import typer
from rich.console import Console
from rich.table import Table

from .converter import (
    OutputFormat,
    convert_geotiff,
    get_geotiff_info,
)

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
) -> None:
    """Convert a GeoTIFF file to ASCII format for Petrel import."""
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
            )

        console.print("[green]✓ Conversion complete![/green]")

        # Show statistics
        if "points_written" in stats:
            console.print(f"  Points written: {stats['points_written']:,}")
        console.print(f"  Grid size: {stats['width']} x {stats['height']} ({stats['total_cells']:,} cells)")

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
) -> None:
    """Batch convert multiple GeoTIFF files in a directory."""
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
                nodata_value=nodata,
                skip_nodata=skip_nodata,
                decimals=decimals,
            )
            success_count += 1
        except Exception as e:
            console.print(f"    [red]Error: {e}[/red]")
            error_count += 1

    console.print()
    console.print(f"[green]✓ Converted: {success_count}[/green]")
    if error_count:
        console.print(f"[red]✗ Errors: {error_count}[/red]")


if __name__ == "__main__":
    app()
