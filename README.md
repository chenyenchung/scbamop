# scbamsplit

> Last updated: 2023-09-07

An open-source and light-weight tool to subset BAM files based on read tags in parallel
powered by `htslib` and [`uthash`](https://troydhanson.github.io/uthash/) under MIT license.

## Installation

Before installing `scbamsplit`, you will need to install [`samtools`](http://www.htslib.org/) so you get
`htslib` in your environment. You also need [`CMake`](https://cmake.org/) to build the tool.

Once you have both installed, you should run `samtools --version` in the terminal and expect to see:

```
samtools 1.17
Using htslib 1.17
Copyright (C) 2023 Genome Research Ltd.
......
```

Also try `cmake --version`, and you should see this in your terminal:

```
cmake version 3.26.0

CMake suite maintained and supported by Kitware (kitware.com/cmake).
```

If you have both tools ready, get the source code of `scbamsplit` by:

```
git clone https://github.com/chenyenchung/scbamsplit.git
```

or from [this link](https://github.com/chenyenchung/scbamsplit/archive/refs/heads/main.zip).

If you clone the repository, please run the following command in
your terminal:

```bash
cd scbamsplit
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release .. && make && cmake --install .
```
This would install the app to `/usr/local/bin`. If you want to install it elsewhere,
replace the last line above with:
```
# Replace the square bracket and enclosed content with your path
cmake -DCMAKE_BUILD_TYPE=Relase -DCMAKE_INSTALL_PREFIX=[Your preferred path] ..
make && cmake --install .
```


If you download from the above link, please decompress it and run the following command in
your terminal:

```bash
cd scbamsplit-main
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release .. && make && cmake --install .
```
## Usage

After installation, you should be able to see the following usage decomentation if
by typing `scbamsplit -h` in your terminal.

```
Program: scbamsplit (Parallel BAM file subset by CBC/UMI)
Version: v0.3.1 (Dependent on htslib v1.17)

Usage: scbamsplit -f path -m path
Options:

    Generic:
        [-o path] [-q MAPQ] [-d] [-r read name length] [-M memory usage (in GB)] [-n] [-v (verbosity)] [-h]
    CBC/UMI related:
        [-p platform] [-b CBC tag/field] [-L CBC length] [-u UMI tag/field] [-l UMI length]

    -f/--file: the path for input bam file
    -m/--meta: the path for input metadata an unquoted two-column csv with column names)
    -o/--output: the path to export bam files to default: ./)
    -q/--mapq: Minimal MAPQ threshold for output default: 0)
    -p/--platform: Pre-fill locations and lengths for CBC and UMI (Supported platform: 10Xv2, 10Xv3, sciRNAseq3
         (e.g., for 10Xv3, both are stored as read tags. CBC is 16 mers tagged CB, while UMI is 12mers tagged UB).
    -d/--dedup: Remove duplicated reads with the same cell barcode/UMI combination
    -b/--cbc-location: If CBC is a read tag, provide the name (e.g., CB); if it is in the read name,
         provide the field number (e.g., 3) (default: CB)
    -L/--cbc-length: The length of the barcode you want to filter against (default: 20)
    -u/--umi-location: If UMI is a read tag, provide the name (e.g., UB); if it is in the read name,
        provide the field number (e.g., 3) (default: CB)
    -l/--umi-length: The length of the UMI default: 20)
    -r/--rn-length: The length of the read name (default: 70)
    -M/--mem: The estimated maximum amount of memory to use (In GB, default: 4)
    -@/--threads: Setting the number of threads to use (default: 1)
    -n/--dry-run: Only print out parameters
    -v/--verbose: Set verbosity level (1 - 5) (default: 2, 3 if -v provided without a value)
    -h/--help: Show this documentation
```

`scbamsplit` requires:

1. A BAM file with reads that are marked by read tags (-f/--file)
2. A comma-separated file (.csv) (-m/--meta) in which the first column is the value of the tag to filtered
while the second is the subset identity (e.g., cluster, sample...). `scbamsplit` will generate a BAM file
for each identity and export reads that contains a tag that belongs to this identity in the file.

By default, the read tag used to query the provided metadata is `CB` (the read tag that contains
corrected cell barcode in `cellranger`-aligned BAMs). While if deduplication is set (`-d`), UMI
is retrieved from `UB`.

To adapted for various platform, you can use `-p`. Currently, we have 10X Genomics v2 and
v3 chemistry (invoke with `-p 10Xv2` and `-p 10Xv3` respectively), and sci-RNAseq3 pipeline
that stores cell barcode and UMI information in the read name (`-p scirnaseq3`).

It is important to make sure the retrieved cell barcode and UMI information is compatible
with the metadata provided. To use a customized setting, you can use `-b` to specify where
cell barcode is: `-b CR` or any two letter tag will be interpreted as a read tag;
`-b 3` or any other number will be interpreted as the n-th field delimited by commas in
the read name (e.g, if the read name is AAA,BBB,CCC,DDD, `-b 3` will retrieve CCC as
the barcode). `-u` follows the same convention but for UMI.

If the CBC/UMI lengths parsed are not correct, `-L [number]` (CBC) and `-l [number]` (UMI)
can be used to specify how long a CBC/UMI should be parsed. 

It is advisable to add `-n`/`--dry-run` to the command, so when you execute the command, nothing but
**printing the information for the run** will happen. If you are happy with what you see, remove
`-n`/`--dry-run` to really run the tool.

`scbamsplit` is pretty easy on memory usage if you don't need deduplication. On the otherhand,
since deduplication involves read sorting in memory, you might want to use `-M [number]` to
restrict the maximum amount of memory `scbamsplit` uses (e.g., `-M 4` will restrict memory
usage under 4GB).

There are some extra functionalities that are optional:

### MAPQ filtering

If you want to set an lower limit of MAPQ for reads to be exported, try `-q`
(e.g., if you want only reads with MAPQ>=30, try `-q 30` or `--mapq 30`).

### UMI-based deduplication

Some sequencing techniques involves adding a unique molecule index (UMI) to
each fragments of interest. In such cases, one might want to keep one read
from each unique fragment once. To enable this function, please use
`-d` (or `--dedup`) flag to turn on deduplication.

When deduplication is on, all reads that contain CBC and UMI will be sorted by CBC-UMI
combinations, and the primary mapping record with the highest MAPQ will be kept along with
all secondary mappings (if present) (New in v0.2.0; in previous versions the first read
with the same CBC-UMI combo will be kept, which will be confounded by genomic location and
does not consider MAPQ and primary/secondary mapping info).

Please note that this function assumes UMIs and filter tags are corrected and contains
no PCR errors.
