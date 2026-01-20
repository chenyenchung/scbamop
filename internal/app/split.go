package app

import (
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strings"

	"github.com/biogo/hts/bam"

	"scbamop/internal/config"
	"scbamop/internal/dedup"
	"scbamop/internal/logging"
	"scbamop/internal/metadata"
	"scbamop/internal/tag"
)

func RunSplit(splitConfig config.SplitConfig) error {
	logger := logging.Logger{Level: splitConfig.LogLevel}

	if splitConfig.Verbose || splitConfig.DryRun {
		printRunConfig(splitConfig)
	}

	if splitConfig.UMIIgnored && splitConfig.Deduplicate {
		return fmt.Errorf("UMI is required for deduplication")
	}

	if err := checkFileReadable(splitConfig.InputPath, "Input BAM"); err != nil {
		logger.Logf(logging.Error, err.Error())
		return err
	}
	if err := checkFileReadable(splitConfig.MetaPath, "Metadata"); err != nil {
		logger.Logf(logging.Error, err.Error())
		return err
	}

	if splitConfig.DryRun {
		fmt.Fprintln(os.Stderr, "Dry run completed successfully.")
		return nil
	}

	if err := ensureOutputDirectory(splitConfig.OutputPrefix, logger); err != nil {
		logger.Logf(logging.Error, err.Error())
		return err
	}

	inputFile, openErr := os.Open(splitConfig.InputPath)
	if openErr != nil {
		logger.Logf(logging.Error, "Failed to open BAM file: %v", openErr)
		return openErr
	}

	reader, readerErr := bam.NewReader(inputFile, 0)
	if readerErr != nil {
		_ = inputFile.Close()
		logger.Logf(logging.Error, "Failed to read BAM header: %v", readerErr)
		return readerErr
	}

	header := reader.Header()
	directMap, mapErr := metadata.LoadMetadata(splitConfig.MetaPath, splitConfig.OutputPrefix, header, splitConfig.CellBarcode.Length, logger)
	if mapErr != nil {
		_ = reader.Close()
		_ = inputFile.Close()
		return mapErr
	}
	defer directMap.Close()

	if splitConfig.Deduplicate {
		_ = reader.Close()
		_ = inputFile.Close()
		return dedup.Deduplicate(splitConfig.InputPath, header, directMap, splitConfig.CellBarcode, splitConfig.UMI, splitConfig.MapQThreshold, logger)
	}

	defer inputFile.Close()
	return splitWithoutDedup(reader, directMap, splitConfig, logger)
}

func splitWithoutDedup(reader *bam.Reader, directMap *metadata.DirectMap, splitConfig config.SplitConfig, logger logging.Logger) error {
	defer reader.Close()

	for {
		record, readErr := reader.Read()
		if readErr == io.EOF {
			break
		}
		if readErr != nil {
			return readErr
		}

		cellBarcode, tagErr := tag.Extract(record, splitConfig.CellBarcode)
		if tagErr != nil {
			continue
		}
		umiValue := ""
		if !splitConfig.UMIIgnored {
			var tagErr error
			umiValue, tagErr = tag.Extract(record, splitConfig.UMI)
			if tagErr != nil {
				continue
			}
		}
		if int(record.MapQ) < splitConfig.MapQThreshold {
			continue
		}
		if !directMap.HasBarcode(cellBarcode) {
			continue
		}
		if lengthErr := tag.ValidateLength(cellBarcode, splitConfig.CellBarcode.Length, "cell barcode"); lengthErr != nil {
			return fmt.Errorf("read %s: %w", record.Name, lengthErr)
		}
		if !splitConfig.UMIIgnored {
			if lengthErr := tag.ValidateLength(umiValue, splitConfig.UMI.Length, "UMI"); lengthErr != nil {
				return fmt.Errorf("read %s: %w", record.Name, lengthErr)
			}
		}

		writeErr := directMap.WriteRecord(cellBarcode, record)
		if writeErr != nil {
			logger.Logf(logging.Error, "Failed to write read: %v", writeErr)
			return writeErr
		}
	}

	return nil
}

func printRunConfig(splitConfig config.SplitConfig) {
	fmt.Fprintln(os.Stderr, "- Run configuration:")
	fmt.Fprintf(os.Stderr, "\tInput BAM: %s\n", splitConfig.InputPath)
	fmt.Fprintf(os.Stderr, "\tMetadata: %s\n", splitConfig.MetaPath)
	fmt.Fprintf(os.Stderr, "\tMAPQ threshold: %d\n", splitConfig.MapQThreshold)
	fmt.Fprintf(os.Stderr, "\tOutput prefix: %s\n", splitConfig.OutputPrefix)
	printTagMeta("Cell barcode", splitConfig.CellBarcode)
	printTagMeta("UMI", splitConfig.UMI)
	fmt.Fprintf(os.Stderr, "\tDeduplication: %s\n\n", boolLabel(splitConfig.Deduplicate))
}

func printTagMeta(header string, tagMeta config.TagMeta) {
	location := "Read tag"
	if tagMeta.Location == config.TagLocationReadName {
		location = "Read name"
	}
	fmt.Fprintf(os.Stderr, "\t%s:\n", header)
	fmt.Fprintf(os.Stderr, "\t\tLocation: %s\n", location)
	if tagMeta.Location == config.TagLocationReadName {
		fmt.Fprintf(os.Stderr, "\t\tSeparator: %s\n", tagMeta.Separator)
		fmt.Fprintf(os.Stderr, "\t\tField number: %d\n", tagMeta.Field)
	} else {
		fmt.Fprintf(os.Stderr, "\t\tTag name: %s\n", tagMeta.TagName)
	}
	fmt.Fprintf(os.Stderr, "\t\tTag length: %d\n\n", tagMeta.Length)
}

func boolLabel(value bool) string {
	if value {
		return "enabled"
	}
	return "disabled"
}

func checkFileReadable(path string, label string) error {
	if path == "" {
		return fmt.Errorf("%s file path is empty", label)
	}
	fileInfo, statErr := os.Stat(path)
	if statErr != nil {
		return fmt.Errorf("%s file not found: %s", label, path)
	}
	if fileInfo.IsDir() {
		return fmt.Errorf("%s path is a directory: %s", label, path)
	}
	fileHandle, openErr := os.Open(path)
	if openErr != nil {
		return fmt.Errorf("%s file not readable: %s", label, path)
	}
	_ = fileHandle.Close()
	return nil
}

func ensureOutputDirectory(outputPrefix string, logger logging.Logger) error {
	outputDir := strings.TrimSuffix(outputPrefix, "/")
	if outputDir == "" {
		outputDir = "."
	}

	directoryInfo, statErr := os.Stat(outputDir)
	if statErr == nil {
		if !directoryInfo.IsDir() {
			return fmt.Errorf("output path is not a directory: %s", outputDir)
		}
		if !isWritableDirectory(outputDir) {
			return fmt.Errorf("output directory not writable: %s", outputDir)
		}
		logger.Logf(logging.Warning, "Output directory already exists: %s", outputPrefix)
		return nil
	}
	if !os.IsNotExist(statErr) {
		return fmt.Errorf("failed to stat output directory: %s", outputDir)
	}

	parentDir := filepath.Dir(outputDir)
	if parentDir == "." {
		parentDir = filepath.Clean(".")
	}
	if !isWritableDirectory(parentDir) {
		return fmt.Errorf("parent directory not writable for output creation: %s", parentDir)
	}

	if createErr := os.Mkdir(outputDir, 0755); createErr != nil {
		return fmt.Errorf("failed to create directory: %s", outputDir)
	}
	logger.Logf(logging.Info, "Created output directory: %s", outputPrefix)
	return nil
}

func isWritableDirectory(path string) bool {
	tempFile, createErr := os.CreateTemp(path, ".scbamop")
	if createErr != nil {
		return false
	}
	tempName := tempFile.Name()
	_ = tempFile.Close()
	_ = os.Remove(tempName)
	return true
}
