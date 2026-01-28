# CLAUDE.md

Project guidance for Claude Code when working with geo-tif-parser.

## Project Overview

**geo-tif-parser** is a Python CLI tool that converts GeoTIFF raster files to ASCII formats compatible with Schlumberger Petrel for surface import. Primary use case is geological surface data (tops, thicknesses) for reservoir modeling.

## Architecture

```
geo_tif_parser/
├── src/geo_tif_parser/
│   ├── __init__.py      # Package version
│   ├── cli.py           # Typer CLI interface (info, convert, batch commands)
│   └── converter.py     # Core conversion logic and format writers
├── pyproject.toml       # Project metadata and dependencies
└── README.md
```

### Key Components

- **converter.py**: Core module with `GeoTiffInfo` dataclass, format writers (`write_petrel_points`, `write_petrel_grid`, `write_esri_ascii`), and main `convert_geotiff` function
- **cli.py**: Three commands (`info`, `convert`, `batch`) using Typer with Rich formatting

## Dependencies

- **rasterio**: GeoTIFF reading (uses GDAL under the hood)
- **numpy**: Array operations
- **typer**: CLI framework
- **rich**: Terminal formatting

## Output Formats

| Format | Extension | Petrel Import As |
|--------|-----------|------------------|
| `petrel_points` | `.xyz` | Points → Surface |
| `petrel_grid` | `.cps3` | CPS-3 Grid Surface |
| `esri_ascii` | `.asc` | ESRI ASCII Grid |

## Development Commands

```bash
# Activate virtual environment
source .venv/bin/activate

# Install in editable mode
pip install -e .

# Run CLI
geotif --help
geotif info <file.tif>
geotif convert <file.tif> -f petrel_grid
geotif batch <directory> -f petrel_points
```

## Testing Data

Test GeoTIFFs are located at:
`/Volumes/SSD-2TB/Year/2026/TrapRock/ColumbiaRiverPlateauGeomodel/`

These are Columbia River Basalt group surfaces:
- `grtop_f.tif`, `grthk_f.tif` - Grande Ronde top/thickness
- `wntop_f.tif`, `wnthk_f.tif` - Wanapum top/thickness
- `smtop_f.tif`, `smthk_f.tif` - Saddle Mountains top/thickness
- `obtop_f.tif`, `obthk_f.tif` - Overburden top/thickness
- `pmtop_f.tif` - Pre-Miocene top
- `geomap_f.tif`, `geomapbr_f.tif` - Geologic maps

Coordinate system: NAD83 / Washington South (Lambert Conformal Conic), units in US survey feet.

## Code Conventions

- Type hints throughout
- Dataclasses for structured data
- Enum for output format options
- Rich console output for user feedback
- NoData handling: skip by default for points, replace with 0.1E+31 for grids (Petrel standard)

## Common Tasks

### Adding a new output format

1. Add enum value to `OutputFormat` in `converter.py`
2. Create `write_<format>()` function following existing patterns
3. Add case to `convert_geotiff()` function
4. Update suffix map in `cli.py` commands
5. Update README format table

### Modifying CLI options

All CLI code is in `cli.py`. Uses Typer's `Annotated` syntax for option definitions with full type hints.

## ZMAP+ Grid Format Reference

The `petrel_grid` format uses standard ZMAP+ ASCII grid format (imported as CPS-3 grid in Petrel):

```
! <comment>                                         # Comment line
@<name>, GRID, <nodes_per_line>                     # Grid identifier
<field_width>, <null_value>, , <decimals>, <start>  # Format spec
<nrows>, <ncols>, <xmin>, <xmax>, <ymin>, <ymax>   # Grid dimensions & bounds
0.0, 0.0, 0.0                                      # Rotation (unused)
@                                                   # End of header
<values>                                            # Column-major, N→S per column, W→E
```

Data is written column-by-column (column-major order), starting from the upper-left (NW)
corner, going down each column (N→S), then advancing columns left to right (W→E).
