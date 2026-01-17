package metadata

import (
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/biogo/hts/sam"

	"scbamop/internal/logging"
)

func TestLoadMetadataInvalidFieldCount(t *testing.T) {
	tempDir := t.TempDir()
	metadataPath := filepath.Join(tempDir, "metadata.csv")
	content := "barcode,label\nCB1,label1,extra\n"
	if err := os.WriteFile(metadataPath, []byte(content), 0644); err != nil {
		t.Fatalf("failed to write metadata: %v", err)
	}
	outputPrefix := filepath.Join(tempDir, "output") + string(os.PathSeparator)
	header := buildHeader(t)
	_, err := LoadMetadata(metadataPath, outputPrefix, header, 3, logging.Logger{Level: logging.Error})
	if err == nil {
		t.Fatalf("expected metadata field count error")
	}
}

func TestLoadMetadataBarcodeLengthMismatch(t *testing.T) {
	tempDir := t.TempDir()
	metadataPath := filepath.Join(tempDir, "metadata.csv")
	barcode := strings.Repeat("A", 4)
	content := "barcode,label\n" + barcode + ",label1\n"
	if err := os.WriteFile(metadataPath, []byte(content), 0644); err != nil {
		t.Fatalf("failed to write metadata: %v", err)
	}
	outputPrefix := filepath.Join(tempDir, "output") + string(os.PathSeparator)
	header := buildHeader(t)
	_, err := LoadMetadata(metadataPath, outputPrefix, header, 3, logging.Logger{Level: logging.Error})
	if err == nil {
		t.Fatalf("expected barcode length mismatch error")
	}
}

func TestLoadMetadataLabelTooLong(t *testing.T) {
	tempDir := t.TempDir()
	metadataPath := filepath.Join(tempDir, "metadata.csv")
	label := strings.Repeat("B", 64)
	content := "barcode,label\nCB1," + label + "\n"
	if err := os.WriteFile(metadataPath, []byte(content), 0644); err != nil {
		t.Fatalf("failed to write metadata: %v", err)
	}
	outputPrefix := filepath.Join(tempDir, "output") + string(os.PathSeparator)
	header := buildHeader(t)
	_, err := LoadMetadata(metadataPath, outputPrefix, header, 3, logging.Logger{Level: logging.Error})
	if err == nil {
		t.Fatalf("expected label length error")
	}
}

func TestLoadMetadataSanitizeLabel(t *testing.T) {
	tempDir := t.TempDir()
	metadataPath := filepath.Join(tempDir, "metadata.csv")
	content := "barcode,label\nCB1,../bad/label\n"
	if err := os.WriteFile(metadataPath, []byte(content), 0644); err != nil {
		t.Fatalf("failed to write metadata: %v", err)
	}
	outputPrefix := filepath.Join(tempDir, "output") + string(os.PathSeparator)
	if err := os.MkdirAll(outputPrefix, 0755); err != nil {
		t.Fatalf("failed to create output dir: %v", err)
	}
	header := buildHeader(t)
	directMap, err := LoadMetadata(metadataPath, outputPrefix, header, 3, logging.Logger{Level: logging.Error})
	if err != nil {
		t.Fatalf("expected metadata load to succeed: %v", err)
	}
	defer directMap.Close()

	if _, ok := directMap.LabelToWriter["_._bad_label"]; !ok {
		t.Fatalf("expected sanitized label key to exist")
	}
	outputPath := filepath.Join(outputPrefix, "_._bad_label.bam")
	if _, err := os.Stat(outputPath); err != nil {
		t.Fatalf("expected output BAM to exist: %v", err)
	}
}

func buildHeader(t *testing.T) *sam.Header {
	header, err := sam.NewHeader(nil, nil)
	if err != nil {
		t.Fatalf("failed to create header: %v", err)
	}
	header.Version = "1.6"
	header.SortOrder = sam.Unsorted
	return header
}
