package app

import (
	"fmt"
	"io"
	"os"
	"path/filepath"
	"testing"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"

	"scbamop/internal/config"
	"scbamop/internal/logging"
)

func TestCheckFileReadableMissing(t *testing.T) {
	if err := checkFileReadable("", "Input BAM"); err == nil {
		t.Fatalf("expected error for empty path")
	}
	if err := checkFileReadable("/path/does/not/exist", "Input BAM"); err == nil {
		t.Fatalf("expected error for missing file")
	}
}

func TestCheckFileReadableDirectory(t *testing.T) {
	tempDir := t.TempDir()
	if err := checkFileReadable(tempDir, "Input BAM"); err == nil {
		t.Fatalf("expected directory error")
	}
}

func TestEnsureOutputDirectoryWithFile(t *testing.T) {
	tempDir := t.TempDir()
	filePath := filepath.Join(tempDir, "output")
	if err := os.WriteFile(filePath, []byte("data"), 0644); err != nil {
		t.Fatalf("failed to write temp file: %v", err)
	}
	logger := logging.Logger{Level: logging.Error}
	if err := ensureOutputDirectory(filePath, logger); err == nil {
		t.Fatalf("expected error for file output path")
	}
}

func TestEnsureOutputDirectoryCreates(t *testing.T) {
	tempDir := t.TempDir()
	outputDir := filepath.Join(tempDir, "new-output")
	outputPrefix := outputDir + string(os.PathSeparator)
	logger := logging.Logger{Level: logging.Error}
	if err := ensureOutputDirectory(outputPrefix, logger); err != nil {
		t.Fatalf("expected output directory created, got %v", err)
	}
	if info, err := os.Stat(outputDir); err != nil || !info.IsDir() {
		t.Fatalf("expected output directory to exist")
	}
}

func TestRunSplitUMIIgnoredSkipsUMICheck(t *testing.T) {
	tempDir := t.TempDir()
	inputPath := filepath.Join(tempDir, "input.bam")
	metadataPath := filepath.Join(tempDir, "metadata.csv")
	outputDir := filepath.Join(tempDir, "output")
	outputPrefix := outputDir + string(os.PathSeparator)
	barcode := "CELL"
	label := "cluster1"

	content := fmt.Sprintf("barcode,label\n%s,%s\n", barcode, label)
	if err := os.WriteFile(metadataPath, []byte(content), 0644); err != nil {
		t.Fatalf("failed to write metadata: %v", err)
	}

	header := buildTestHeader(t)
	if err := writeTestBAM(inputPath, header, barcode); err != nil {
		t.Fatalf("failed to write input BAM: %v", err)
	}

	splitConfig := config.SplitConfig{
		InputPath:     inputPath,
		MetaPath:      metadataPath,
		OutputPrefix:  outputPrefix,
		MapQThreshold: 0,
		Deduplicate:   false,
		DryRun:        false,
		Verbose:       false,
		LogLevel:      logging.Error,
		CellBarcode: config.TagMeta{
			Location: config.TagLocationReadTag,
			TagName:  "CB",
			Length:   len(barcode),
		},
		UMI: config.TagMeta{
			Location: config.TagLocationReadTag,
			TagName:  "UB",
			Length:   10,
		},
		UMIIgnored: true,
	}

	if err := RunSplit(splitConfig); err != nil {
		t.Fatalf("expected split to succeed, got %v", err)
	}

	outputPath := filepath.Join(outputDir, label+".bam")
	outputFile, err := os.Open(outputPath)
	if err != nil {
		t.Fatalf("failed to open output BAM: %v", err)
	}
	defer outputFile.Close()

	reader, err := bam.NewReader(outputFile, 0)
	if err != nil {
		t.Fatalf("failed to read output BAM: %v", err)
	}
	defer reader.Close()

	if _, err := reader.Read(); err != nil {
		t.Fatalf("expected one record, got %v", err)
	}
	if _, err := reader.Read(); err != io.EOF {
		t.Fatalf("expected EOF after one record, got %v", err)
	}
}

func buildTestHeader(t *testing.T) *sam.Header {
	header, err := sam.NewHeader(nil, nil)
	if err != nil {
		t.Fatalf("failed to create header: %v", err)
	}
	header.Version = "1.6"
	header.SortOrder = sam.Unsorted
	return header
}

func writeTestBAM(path string, header *sam.Header, barcode string) error {
	file, err := os.Create(path)
	if err != nil {
		return err
	}
	defer file.Close()

	writer, err := bam.NewWriter(file, header, 0)
	if err != nil {
		return err
	}
	defer writer.Close()

	aux, err := sam.NewAux(sam.NewTag("CB"), barcode)
	if err != nil {
		return err
	}

	record, err := sam.NewRecord("read1", nil, nil, -1, -1, 0, 0, nil, []byte("A"), nil, []sam.Aux{aux})
	if err != nil {
		return err
	}
	record.Flags = sam.Unmapped | sam.MateUnmapped

	return writer.Write(record)
}
