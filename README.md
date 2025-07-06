# scbamop

> Last updated: 2025-07-06

A high-performance single-cell BAM operations toolkit with UMI-based deduplication and cell barcode splitting, powered by vendored `htslib` and [`uthash`](https://troydhanson.github.io/uthash/) under MIT license.

## Features

- **Cell barcode-based BAM splitting**: Subset BAM files by cell barcodes in parallel
- **UMI-based deduplication**: Memory-efficient 3-pass algorithm for removing PCR duplicates
- **Multiple platform support**: 10X Genomics v2/v3, sci-RNA-seq3, and custom configurations
- **Automatic label sanitization**: Secure handling of metadata file labels
- **High performance**: Optimized for large single-cell datasets

## Installation

### Prerequisites

**Build Tools:**
- CMake (≥ 3.18)
- C compiler with C99 support (gcc, clang)
- Git (for submodule management)

**Runtime Libraries:**
- zlib development headers
- bzip2 development headers  
- liblzma development headers
- libcurl development headers

**Installation on Debian/Ubuntu:**
```bash
apt-get update
apt-get install build-essential cmake git zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev
```

**Installation on RHEL/CentOS/Fedora:**
```bash
yum install gcc make cmake git zlib-devel bzip2-devel xz-devel libcurl-devel
# or dnf install ... on newer systems
```

### Building from Source

```bash
# Clone the repository
git clone https://github.com/chenyenchung/scbamop.git
cd scbamop

# Initialize submodules (includes vendored htslib)
git submodule update --init --recursive

# Build
mkdir build && cd build
cmake ..
make
```

The compiled binary `scbamop` will be available in the `build` directory.

### Debug Builds

For development and debugging, additional sanitizer options are available:

```bash
# Default debug build (UndefinedBehaviorSanitizer)
cmake ..
make

# Memory debugging (AddressSanitizer + UBSan)
cmake -DENABLE_ASAN=ON ..
make

# Thread safety testing (ThreadSanitizer only)
cmake -DENABLE_UBSAN=OFF -DENABLE_TSAN=ON ..
make

# Release build (no sanitizers)
cmake -DCMAKE_BUILD_TYPE=Release -DENABLE_UBSAN=OFF ..
make
```

## Usage

### Basic Command Structure

```bash
scbamop split -f input.bam -m metadata.csv [options]
```

### Required Arguments

- `-f, --file`: Input BAM file path
- `-m, --meta`: Metadata CSV file (two columns: barcode, label)

### Common Options

- `-o, --output`: Output directory (default: current directory)
- `-d, --dedup`: Enable UMI-based deduplication
- `-q, --mapq`: Minimum MAPQ threshold (default: 0)
- `-v, --verbose`: Verbosity level (0-5, default: 2)
- `-h, --help`: Show help message

### Platform-Specific Options

- `-p, --platform`: Pre-configured platform settings
  - `10Xv2`: 10X Genomics v2 chemistry
  - `10Xv3`: 10X Genomics v3 chemistry  
  - `sciRNAseq3`: sci-RNA-seq3 pipeline
- `-b, --cbc-location`: Cell barcode tag/field (default: CB)
- `-u, --umi-location`: UMI tag/field (default: UB)

### Example Usage

**Basic splitting without deduplication:**
```bash
scbamop split -f sample.bam -m metadata.csv -o output/
```

**10X Genomics data with UMI deduplication:**
```bash
scbamop split -f sample.bam -m metadata.csv -d -p 10Xv3 -q 30 -v 3
```

**sci-RNA-seq3 data (barcodes in read names):**
```bash
scbamop split -f sample.bam -m metadata.csv -p sciRNAseq3 -d
```

**Custom barcode/UMI locations:**
```bash
scbamop split -f sample.bam -m metadata.csv -b CR -u UR -d
```

## Metadata File Format

The metadata file must be a two-column CSV with headers:

```csv
barcode,label
AAACCCAAGAAACACT,CD4_T_cells
AAACCCAAGAAACCAT,B_cells
AAACCCAAGAAACCCA,NK_cells
```

## Security Features

### Automatic Label Sanitization

Output labels from metadata files are automatically sanitized to prevent security vulnerabilities:

- **Path traversal sequences** (`..`) → `__`
- **Directory separators** (`/`, `\`) → `_`
- **Hidden file prefixes** (`.`) → `_`
- **Special characters** → `_`

**Example sanitization:**
```csv
barcode,label
AAACCCAAGAAACACT,CD4/CD8 T-cells        # → CD4_CD8 T-cells
AAACCCAAGAAACCAT,../../../etc/passwd    # → ___________etc_passwd
AAACCCAAGAAACCCA,.hidden_file           # → _hidden_file
```

## UMI Deduplication

The tool uses a memory-efficient 3-pass algorithm:

1. **Pass 1**: Extract read information used for deduplication (CB, UMI, coordinates, and MAPQ)
2. **Pass 2**: In-memory duplicate marking
3. **Pass 3**: Write deduplicated reads to output files

When deduplication is enabled (`-d`), reads with identical cell barcode + UMI + genomic coordinates are considered duplicates. The primary mapping with the highest MAPQ is retained.
Memory usage scales with the number of unique molecules when deduplication is enabled

## License

MIT License - see LICENSE file for details.

## Support

For issues, questions, or feature requests, please open an issue on the GitHub repository.
