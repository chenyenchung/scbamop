package app

import (
	"os"
	"path/filepath"
	"testing"

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
