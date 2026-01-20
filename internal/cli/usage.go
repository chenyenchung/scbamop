package cli

import (
	"fmt"
	"io"
)

func PrintGlobalUsage(writer io.Writer) {
	fmt.Fprintln(writer, "Program: scbamop (Single-cell BAM operations toolkit)")
	fmt.Fprintln(writer, "Version: v0.5.0 (subcommand structure)")
	fmt.Fprintln(writer, "")
	fmt.Fprintln(writer, "Usage: scbamop <command> [options]")
	fmt.Fprintln(writer, "")
	fmt.Fprintln(writer, "Commands:")
	fmt.Fprintln(writer, "  split    Split BAM file by cell barcodes with optional deduplication")
	fmt.Fprintln(writer, "")
	fmt.Fprintln(writer, "Use 'scbamop <command> --help' for command-specific help")
	fmt.Fprintln(writer, "")
}

func PrintSplitUsage(writer io.Writer) {
	fmt.Fprintln(writer, "Usage: scbamop split -f FILE -m FILE [options]")
	fmt.Fprintln(writer, "")
	fmt.Fprintln(writer, "Split BAM file by cell barcodes with optional UMI-based deduplication")
	fmt.Fprintln(writer, "")
	fmt.Fprintln(writer, "Required arguments:")
	fmt.Fprintln(writer, "  -f, --file FILE        Input BAM file path")
	fmt.Fprintln(writer, "  -m, --meta FILE        Metadata file with cell barcode assignments")
	fmt.Fprintln(writer, "")
	fmt.Fprintln(writer, "Optional arguments:")
	fmt.Fprintln(writer, "  -o, --output DIR       Output directory prefix (default: ./)")
	fmt.Fprintln(writer, "  -q, --mapq INT         MAPQ threshold (default: 0)")
	fmt.Fprintln(writer, "  -d, --dedup            Enable UMI-based deduplication")
	fmt.Fprintln(writer, "  -b, --cbc-location STR Cell barcode tag name or field number (default: CB)")
	fmt.Fprintln(writer, "  -u, --umi-location STR UMI tag name or field number (default: UB, use 0 to ignore)")
	fmt.Fprintln(writer, "  -v, --verbose [INT]    Verbosity level: -v (INFO), -v 5 or --verbose=5 (DEBUG)")
	fmt.Fprintln(writer, "  --atac                 Ignore UMI checks (disables deduplication)")
	fmt.Fprintln(writer, "  -h, --help             Show this help message")
	fmt.Fprintln(writer, "")
}
