# geo-tif-parser

A Python CLI tool for converting GeoTIFF raster files to ASCII formats compatible with Schlumberger Petrel and other geological modeling software.

## Features

- Convert GeoTIFF files to multiple ASCII formats
- Batch processing of entire directories
- Preserves coordinate reference system information
- Configurable NoData handling
- Adjustable decimal precision
- Rich terminal output with progress feedback

## Installation

### Prerequisites

- Python 3.10 or higher
- GDAL libraries (required by rasterio)

On macOS with Homebrew:
```bash
brew install gdal
```

On Ubuntu/Debian:
```bash
sudo apt-get install gdal-bin libgdal-dev
```

### Install from source

```bash
git clone https://github.com/tilke/geo-tif-parser.git
cd geo-tif-parser
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
```

## Usage

### Display GeoTIFF Information

View metadata about a GeoTIFF file including dimensions, coordinate system, cell size, and bounds:

```bash
geotif info input.tif
```

Example output:
```
                    GeoTIFF Info: grtop_f.tif
┏━━━━━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
│ Property         │ Value                                   │
┡━━━━━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┩
│ Dimensions       │ 3151 x 2958                             │
│ Total Cells      │ 9,320,658                               │
│ Data Type        │ float32                                 │
│ Cell Size X      │ 500.000                                 │
│ Cell Size Y      │ 500.000                                 │
│ CRS              │ NAD83 / Washington South                │
└──────────────────┴─────────────────────────────────────────┘
```

### Convert Single File

Convert a GeoTIFF to ASCII format:

```bash
# Default: XYZ points format (.xyz)
geotif convert input.tif

# Specify output format
geotif convert input.tif --format petrel_grid
geotif convert input.tif -f esri_ascii

# Custom output path
geotif convert input.tif -o /path/to/output.xyz

# Adjust decimal precision
geotif convert input.tif --decimals 2

# Include NoData values instead of skipping them
geotif convert input.tif --include-nodata --nodata -9999
```

### Batch Convert Multiple Files

Convert all GeoTIFF files in a directory:

```bash
# Convert all .tif files in directory
geotif batch /path/to/geotiffs/

# Specify output directory and format
geotif batch /path/to/geotiffs/ -o /path/to/output/ -f petrel_grid

# Use custom file pattern
geotif batch /path/to/geotiffs/ --pattern "surface_*.tif"
```

## Output Formats

### Petrel Points (`petrel_points`)

Space-delimited X Y Z points, one per line. This is the default format.

- **Extension**: `.xyz`
- **Petrel Import**: Import as Points, then create surface from points
- **NoData Handling**: Skips NoData cells by default

```
1886381.621 1095389.234 2567.263
1886381.621 1094889.234 2593.142
1886881.621 1094889.234 2556.512
```

### Petrel Grid (`petrel_grid`)

Petrel-native CPS-3 format preserving the full grid structure.

- **Extension**: `.cps3`
- **Petrel Import**: Import directly as CPS-3 Grid Surface
- **NoData Handling**: Replaced with 0.1E+31 (Petrel standard)

```
FSASCI 0 1 "COMPUTED" 0 0.1E+31
FSATTR 0 0
FSLIMI 1313631.621 2889131.621 -353360.766 1125639.234 -2293.685 9063.338
FSNROW 3151 2958
FSXINC 500.000 500.000
->Converted from surface.tif
2567.263  2593.142  2556.512  2654.270  2579.251
```

### ESRI ASCII Grid (`esri_ascii`)

Standard ESRI ASCII raster format, widely compatible.

- **Extension**: `.asc`
- **Petrel Import**: Import as ESRI ASCII Grid
- **NoData Handling**: Replaced with -9999 (configurable)

```
ncols         3151
nrows         2958
xllcorner     1313631.621
yllcorner     -353360.766
cellsize      500.000
NODATA_value  -9999
2567.263 2593.142 2556.512 ...
```

## Command Reference

### `geotif info`

```
Usage: geotif info [OPTIONS] INPUT_FILE

  Display information about a GeoTIFF file.

Arguments:
  INPUT_FILE  Path to GeoTIFF file [required]
```

### `geotif convert`

```
Usage: geotif convert [OPTIONS] INPUT_FILE

  Convert a GeoTIFF file to ASCII format for Petrel import.

Arguments:
  INPUT_FILE                      Path to GeoTIFF file [required]

Options:
  -o, --output PATH               Output file path (default: input name with new extension)
  -f, --format [petrel_points|petrel_grid|esri_ascii]
                                  Output format [default: petrel_points]
  --nodata FLOAT                  NoData value for output file
  --skip-nodata / --include-nodata
                                  Skip or include nodata values [default: skip-nodata]
  -d, --decimals INTEGER          Number of decimal places (0-10) [default: 3]
```

### `geotif batch`

```
Usage: geotif batch [OPTIONS] INPUT_DIR

  Batch convert multiple GeoTIFF files in a directory.

Arguments:
  INPUT_DIR                       Directory containing GeoTIFF files [required]

Options:
  -o, --output PATH               Output directory (default: same as input)
  -f, --format [petrel_points|petrel_grid|esri_ascii]
                                  Output format [default: petrel_points]
  -p, --pattern TEXT              Glob pattern for input files [default: *.tif]
  --nodata FLOAT                  NoData value for output files
  --skip-nodata / --include-nodata
                                  Skip or include nodata values [default: skip-nodata]
  -d, --decimals INTEGER          Number of decimal places (0-10) [default: 3]
```

## Importing into Petrel

### XYZ Points

1. In Petrel, go to **Import** > **Points**
2. Select the `.xyz` file
3. Set delimiter to **Space**
4. Map columns: Column 1 → X, Column 2 → Y, Column 3 → Z
5. Set the coordinate reference system to match your data
6. Create surface: Right-click points → **Create Surface**

### CPS-3 Grid

1. In Petrel, go to **Import** > **Surface (ZMAP/CPS-3)**
2. Select the `.cps3` file
3. Set the coordinate reference system
4. The surface imports directly with full grid structure

### ESRI ASCII Grid

1. In Petrel, go to **Import** > **Surface (ESRI ASCII)**
2. Select the `.asc` file
3. Set the coordinate reference system

## Development

```bash
# Clone and setup
git clone https://github.com/tilke/geo-tif-parser.git
cd geo-tif-parser
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"

# Run tests
pytest

# Run CLI in development
geotif --help
```

## License

MIT License - see LICENSE file for details.

## Author

Peter Tilke <tilke@alum.mit.edu>
