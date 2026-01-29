# geo-tif-parser

A Python CLI tool for converting GeoTIFF raster files to ASCII formats compatible with Schlumberger Petrel and other geological modeling software.

## Features

- Convert GeoTIFF files to multiple ASCII formats
- Batch processing of entire directories
- **Well data parsing**: Read tab-delimited well/stratigraphic data, transform coordinates, and export Petrel-format well headers and tops
- **CRS transformation**: Transform coordinates to a different CRS (e.g., UTM)
- **Unit conversion**: Convert elevation values from US survey feet to meters
- Two reprojection modes: fast bounds-only or full raster reprojection
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

# Flip Z values (convert elevation to depth, positive down)
geotif convert input.tif -f petrel_grid --flip-z

# Transform coordinates to UTM Zone 10N (bounds-only mode, fast)
geotif convert input.tif -f petrel_grid --output-crs 26910

# Full reprojection with bilinear resampling
geotif convert input.tif -f petrel_grid --output-crs 26910 --reproject-mode full

# Full reprojection with cubic resampling
geotif convert input.tif -f petrel_points --output-crs EPSG:26910 --reproject-mode full --resampling cubic
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

# Batch convert with Z flipped to depth
geotif batch /path/to/geotiffs/ -f petrel_grid --flip-z

# Batch convert with CRS transformation
geotif batch /path/to/geotiffs/ -f petrel_grid --output-crs 26910
```

### Parse Well Data

Parse tab-delimited well/stratigraphic point data, transform coordinates, convert elevations, and write Petrel-format well header and well tops files:

```bash
# Convert well data from NAD83 Washington South (ft) to UTM 11N (m)
geotif wells wells.txt -s 2927 -c 32611

# Specify output directory
geotif wells wells.txt -s 2927 -c 32611 -o /path/to/output/

# Keep elevations in feet (no vertical unit conversion)
geotif wells wells.txt -s 2927 -c 32611 --no-feet-to-meters

# Adjust coordinate decimal precision
geotif wells wells.txt -s 2927 -c 32611 --decimals 3
```

This produces two files:
- **`well_headers.txt`** — tab-delimited `WELL  X  Y  ELEV` (surface elevation)
- **`well_tops.txt`** — tab-delimited `WELL  SURFACE  ELEV  TVD` (formation top elevation and true vertical depth from surface)

### Display Detailed Header Information

View detailed CRS and grid information, optionally with transformed bounds:

```bash
# Show detailed header info
geotif header input.tif

# Preview transformed bounds in a different CRS
geotif header input.tif --output-crs 26910
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

Standard ZMAP+ ASCII grid format, importable as CPS-3 Grid Surface in Petrel.

- **Extension**: `.cps3`
- **Petrel Import**: Import directly as CPS-3 Grid Surface
- **NoData Handling**: Replaced with 0.1E+31 (Petrel standard)

```
! Converted from surface.tif | Z: Elevation (positive up)
@surface, GRID, 5
15, 0.1000000E+31, , 3, 1
2958, 3151, 1313631.621, 2889131.621, -353360.766, 1125639.234
0.0, 0.0, 0.0
@
   2567.263   2593.142   2556.512   2654.270   2579.251
```

The comment line indicates whether Z values represent elevation (positive up) or depth (positive down, when using `--flip-z`). Data is written column-major (N→S per column, columns W→E).

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

### Well Headers (`well_headers.txt`)

Tab-delimited well location file with surface elevation.

```
WELL	X	Y	ELEV
NWIS_1	102471.49	5049515.40	382.83
NWIS_2	334771.85	5038653.22	592.53
```

### Well Tops (`well_tops.txt`)

Tab-delimited formation top picks with elevation and true vertical depth from surface.

- **TVD** = surface elevation − top elevation (depth below ground surface)

```
WELL	SURFACE	ELEV	TVD
NWIS_3	Grande_Ronde_Basalt	588.3	0.6
NWIS_4	Wanapum_Basalt	66.8	103.6
```

## CRS Transformation

Transform output coordinates to a different Coordinate Reference System (CRS) using the `--output-crs` option. This is useful when your GeoTIFF is in one CRS (e.g., NAD83 State Plane) but you need to import into Petrel using a different CRS (e.g., UTM).

### CRS Specification

The `--output-crs` option accepts:
- **EPSG code as number**: `26910` (UTM Zone 10N)
- **EPSG code with prefix**: `EPSG:26910`
- **Proj4 string**: `+proj=utm +zone=10 +datum=NAD83`

### Reprojection Modes

| Mode | Description | Best For |
|------|-------------|----------|
| `bounds-only` (default) | Transform bounds/coordinates only, preserves original pixel data | Fast, minor projection changes |
| `full` | Reproject entire raster with resampling | Major projection changes, accurate results |

**bounds-only** mode transforms coordinate values but keeps the original pixel data unchanged. This is fast and works well when the source and destination CRS are similar (e.g., both projected coordinate systems).

**full** mode reprojects the entire raster, resampling pixel values to the new grid. This is more accurate for major projection changes but slower. Available resampling methods: `nearest`, `bilinear` (default), `cubic`, `lanczos`.

### Examples

```bash
# Preview transformed bounds before converting
geotif header surface.tif --output-crs 26910

# Convert State Plane (US Survey Feet) to UTM (meters)
geotif convert surface.tif -f petrel_grid --output-crs 26910

# Full reprojection for maximum accuracy
geotif convert surface.tif -f petrel_grid --output-crs 26910 --reproject-mode full --resampling cubic

# Batch convert directory with CRS transformation
geotif batch /data/surfaces/ -f petrel_grid --output-crs EPSG:26910
```

## Command Reference

### `geotif info`

```
Usage: geotif info [OPTIONS] INPUT_FILE

  Display information about a GeoTIFF file.

Arguments:
  INPUT_FILE  Path to GeoTIFF file [required]
```

### `geotif header`

```
Usage: geotif header [OPTIONS] INPUT_FILE

  Display detailed header information including CRS for Petrel import.

Arguments:
  INPUT_FILE                      Path to GeoTIFF file [required]

Options:
  -c, --output-crs TEXT           Preview transformed bounds in this CRS (e.g., 26910, EPSG:26910)
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
  --flip-z                        Negate Z values (convert elevation to depth or vice versa)
  -c, --output-crs TEXT           Output CRS for coordinate transformation (e.g., 26910, EPSG:26910)
  --reproject-mode TEXT           Reprojection mode: 'bounds-only' (fast, default) or 'full'
  --resampling TEXT               Resampling method for full mode: nearest, bilinear, cubic, lanczos
```

### `geotif wells`

```
Usage: geotif wells [OPTIONS] INPUT_FILE

  Parse well data, transform coordinates, and write Petrel well header/tops files.

Arguments:
  INPUT_FILE                      Path to tab-delimited well data file [required]

Options:
  -o, --output PATH               Output directory (default: same as input file)
  -s, --source-crs TEXT           Source EPSG code [default: 2927]
  -c, --output-crs TEXT           Output EPSG code [default: 32611]
  --nodata FLOAT                  NoData sentinel value in input file [default: -9999]
  -d, --decimals INTEGER          Number of decimal places (0-10) [default: 2]
  --feet-to-meters / --no-feet-to-meters
                                  Convert elevations from US survey feet to meters [default: feet-to-meters]
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
  --flip-z                        Negate Z values (convert elevation to depth or vice versa)
  -c, --output-crs TEXT           Output CRS for coordinate transformation (e.g., 26910, EPSG:26910)
  --reproject-mode TEXT           Reprojection mode: 'bounds-only' (fast, default) or 'full'
  --resampling TEXT               Resampling method for full mode: nearest, bilinear, cubic, lanczos
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

### Well Headers and Tops

1. In Petrel, go to **Import** > **Well Heads**
2. Select `well_headers.txt`
3. Map columns: WELL → Well name, X → X, Y → Y, ELEV → Elevation
4. Set the coordinate reference system to match the output CRS
5. Import well tops: **Import** > **Well Tops**
6. Select `well_tops.txt`
7. Map columns: WELL → Well name, SURFACE → Surface, TVD → TVD

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
