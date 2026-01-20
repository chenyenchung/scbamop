package cli

import (
	"strings"
	"testing"

	"scbamop/internal/config"
	"scbamop/internal/logging"
)

func TestParseSplitArgsDefaults(t *testing.T) {
	splitConfig, parseErr := ParseSplitArgs([]string{"-f", "input.bam", "-m", "meta.csv"})
	if parseErr != nil {
		t.Fatalf("expected no error, got %v", parseErr)
	}
	if splitConfig.InputPath != "input.bam" {
		t.Fatalf("expected input path to be set")
	}
	if splitConfig.MetaPath != "meta.csv" {
		t.Fatalf("expected meta path to be set")
	}
	if splitConfig.OutputPrefix != "./" {
		t.Fatalf("expected default output prefix './', got %q", splitConfig.OutputPrefix)
	}
	if splitConfig.MapQThreshold != 0 {
		t.Fatalf("expected default mapq 0, got %d", splitConfig.MapQThreshold)
	}
	if splitConfig.Deduplicate {
		t.Fatalf("expected dedup false by default")
	}
	if splitConfig.DryRun {
		t.Fatalf("expected dry run false by default")
	}
	if splitConfig.LogLevel != logging.Warning {
		t.Fatalf("expected warning log level by default")
	}
	if splitConfig.CellBarcode.TagName != "CB" || splitConfig.CellBarcode.Length != 20 {
		t.Fatalf("unexpected default cell barcode meta: %+v", splitConfig.CellBarcode)
	}
	if splitConfig.UMI.TagName != "UB" || splitConfig.UMI.Length != 20 {
		t.Fatalf("unexpected default UMI meta: %+v", splitConfig.UMI)
	}
}

func TestParseSplitArgsOutputPrefix(t *testing.T) {
	splitConfig, parseErr := ParseSplitArgs([]string{"-f", "input.bam", "-m", "meta.csv", "-o", "out"})
	if parseErr != nil {
		t.Fatalf("expected no error, got %v", parseErr)
	}
	if splitConfig.OutputPrefix != "out/" {
		t.Fatalf("expected output prefix to end with '/', got %q", splitConfig.OutputPrefix)
	}
}

func TestParseSplitArgsVerboseLevel(t *testing.T) {
	splitConfig, parseErr := ParseSplitArgs([]string{"-f", "input.bam", "-m", "meta.csv", "-v", "5"})
	if parseErr != nil {
		t.Fatalf("expected no error, got %v", parseErr)
	}
	if splitConfig.LogLevel != logging.Debug {
		t.Fatalf("expected debug level, got %v", splitConfig.LogLevel)
	}
	if !splitConfig.Verbose {
		t.Fatalf("expected verbose flag to be true")
	}

	splitConfig, parseErr = ParseSplitArgs([]string{"-f", "input.bam", "-m", "meta.csv", "--verbose=2"})
	if parseErr != nil {
		t.Fatalf("expected no error, got %v", parseErr)
	}
	if splitConfig.LogLevel != logging.Level(2) {
		t.Fatalf("expected level 2, got %v", splitConfig.LogLevel)
	}
}

func TestParseSplitArgsTagField(t *testing.T) {
	splitConfig, parseErr := ParseSplitArgs([]string{"-f", "input.bam", "-m", "meta.csv", "-b", "2", "-u", "3"})
	if parseErr != nil {
		t.Fatalf("expected no error, got %v", parseErr)
	}
	if splitConfig.CellBarcode.Location != config.TagLocationReadName || splitConfig.CellBarcode.Field != 2 {
		t.Fatalf("unexpected cell barcode location: %+v", splitConfig.CellBarcode)
	}
	if splitConfig.UMI.Location != config.TagLocationReadName || splitConfig.UMI.Field != 3 {
		t.Fatalf("unexpected UMI location: %+v", splitConfig.UMI)
	}
	if splitConfig.UMIIgnored {
		t.Fatalf("expected UMI to be enabled")
	}
}

func TestParseSplitArgsAtac(t *testing.T) {
	splitConfig, parseErr := ParseSplitArgs([]string{"-f", "input.bam", "-m", "meta.csv", "--atac"})
	if parseErr != nil {
		t.Fatalf("expected no error, got %v", parseErr)
	}
	if !splitConfig.UMIIgnored {
		t.Fatalf("expected UMI to be ignored")
	}

	splitConfig, parseErr = ParseSplitArgs([]string{"-f", "input.bam", "-m", "meta.csv", "-u", "0"})
	if parseErr != nil {
		t.Fatalf("expected no error, got %v", parseErr)
	}
	if !splitConfig.UMIIgnored {
		t.Fatalf("expected UMI to be ignored for -u 0")
	}
}

func TestParseSplitArgsMissingRequired(t *testing.T) {
	_, parseErr := ParseSplitArgs([]string{"-f", "input.bam"})
	if parseErr == nil {
		t.Fatalf("expected missing required arguments error")
	}
	if !strings.Contains(parseErr.Error(), "missing required arguments") {
		t.Fatalf("unexpected error: %v", parseErr)
	}
}

func TestParseSplitArgsInvalidValues(t *testing.T) {
	cases := []struct {
		name string
		args []string
	}{
		{name: "mapq", args: []string{"-f", "input.bam", "-m", "meta.csv", "-q", "-1"}},
		{name: "cb-length", args: []string{"-f", "input.bam", "-m", "meta.csv", "-L", "0"}},
		{name: "umi-length", args: []string{"-f", "input.bam", "-m", "meta.csv", "-l", "0"}},
		{name: "tag", args: []string{"-f", "input.bam", "-m", "meta.csv", "-b", "TOO"}},
		{name: "verbose", args: []string{"-f", "input.bam", "-m", "meta.csv", "-v", "nope"}},
		{name: "atac-dedup", args: []string{"-f", "input.bam", "-m", "meta.csv", "--atac", "-d"}},
		{name: "unknown", args: []string{"-f", "input.bam", "-m", "meta.csv", "--unknown"}},
	}
	for _, testCase := range cases {
		_, parseErr := ParseSplitArgs(testCase.args)
		if parseErr == nil {
			t.Fatalf("expected error for %s", testCase.name)
		}
	}
}

func TestParseSplitArgsVerboseInline(t *testing.T) {
	splitConfig, parseErr := ParseSplitArgs([]string{"-f", "input.bam", "-m", "meta.csv", "-v5"})
	if parseErr != nil {
		t.Fatalf("expected no error, got %v", parseErr)
	}
	if splitConfig.LogLevel != logging.Debug {
		t.Fatalf("expected debug level, got %v", splitConfig.LogLevel)
	}
	if !splitConfig.Verbose {
		t.Fatalf("expected verbose flag to be true")
	}
}
