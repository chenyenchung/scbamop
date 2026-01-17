package metadata

import (
	"bufio"
	"fmt"
	"os"
	"strings"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"

	"scbamop/internal/logging"
)

const (
	maxLabelLength      = 64
	maxOutputPathLength = 512
)

type DirectMap struct {
	BarcodeToWriter map[string]*bam.Writer
	LabelToWriter   map[string]*bam.Writer
}

func LoadMetadata(path string, outputPrefix string, header *sam.Header, cellBarcodeLength int, logger logging.Logger) (*DirectMap, error) {
	fileHandle, openErr := os.Open(path)
	if openErr != nil {
		return nil, fmt.Errorf("cannot open file (%s): %w", path, openErr)
	}
	defer fileHandle.Close()

	if header == nil {
		return nil, fmt.Errorf("missing BAM header for metadata parsing")
	}

	directMap := &DirectMap{
		BarcodeToWriter: make(map[string]*bam.Writer),
		LabelToWriter:   make(map[string]*bam.Writer),
	}
	warnedLabels := make(map[string]bool)

	scanner := bufio.NewScanner(fileHandle)
	firstLine := true
	for scanner.Scan() {
		line := scanner.Text()
		if firstLine {
			firstLine = false
			continue
		}
		fields := strings.Split(line, ",")
		if len(fields) != 2 {
			closeErr := directMap.Close()
			if closeErr != nil {
				logger.Logf(logging.Warning, "Failed to close output files after metadata error: %v", closeErr)
			}
			return nil, fmt.Errorf("there are %d fields in the metadata but only 2 are expected", len(fields))
		}
		barcode := fields[0]
		label := fields[1]
		if cellBarcodeLength > 0 && len(barcode) != cellBarcodeLength {
			closeErr := directMap.Close()
			if closeErr != nil {
				logger.Logf(logging.Warning, "Failed to close output files after metadata error: %v", closeErr)
			}
			return nil, fmt.Errorf("cell barcode length %d does not match expected %d: %s", len(barcode), cellBarcodeLength, barcode)
		}
		if len(label) >= maxLabelLength {
			closeErr := directMap.Close()
			if closeErr != nil {
				logger.Logf(logging.Warning, "Failed to close output files after metadata error: %v", closeErr)
			}
			return nil, fmt.Errorf("label too long (max %d chars): %s", maxLabelLength-1, label)
		}

		sanitizedLabel, modified := sanitizeLabel(label)
		if modified && !warnedLabels[label] {
			logger.Logf(logging.Warning, "Sanitized label: '%s' -> '%s'", label, sanitizedLabel)
			warnedLabels[label] = true
		}

		writer, ok := directMap.LabelToWriter[sanitizedLabel]
		if !ok {
			outputPath := fmt.Sprintf("%s%s.bam", outputPrefix, sanitizedLabel)
			if len(outputPrefix)+len(sanitizedLabel)+len(".bam") >= maxOutputPathLength {
				closeErr := directMap.Close()
				if closeErr != nil {
					logger.Logf(logging.Warning, "Failed to close output files after metadata error: %v", closeErr)
				}
				return nil, fmt.Errorf("output path too long for label: %s", sanitizedLabel)
			}
			outputFile, createErr := os.Create(outputPath)
			if createErr != nil {
				closeErr := directMap.Close()
				if closeErr != nil {
					logger.Logf(logging.Warning, "Failed to close output files after metadata error: %v", closeErr)
				}
				return nil, fmt.Errorf("failed to create output file: %s", outputPath)
			}
			writer, createErr = bam.NewWriter(outputFile, header, 0)
			if createErr != nil {
				_ = outputFile.Close()
				closeErr := directMap.Close()
				if closeErr != nil {
					logger.Logf(logging.Warning, "Failed to close output files after metadata error: %v", closeErr)
				}
				return nil, fmt.Errorf("failed to write header to: %s", outputPath)
			}
			directMap.LabelToWriter[sanitizedLabel] = writer
			logger.Logf(logging.Info, "Created output file: %s", outputPath)
		}

		directMap.BarcodeToWriter[barcode] = writer
	}

	if scanErr := scanner.Err(); scanErr != nil {
		closeErr := directMap.Close()
		if closeErr != nil {
			logger.Logf(logging.Warning, "Failed to close output files after metadata error: %v", closeErr)
		}
		return nil, scanErr
	}

	return directMap, nil
}

func (directMap *DirectMap) Close() error {
	if directMap == nil {
		return nil
	}
	var firstErr error
	for _, writer := range directMap.LabelToWriter {
		if writer == nil {
			continue
		}
		closeErr := writer.Close()
		if closeErr != nil && firstErr == nil {
			firstErr = closeErr
		}
	}
	return firstErr
}

func (directMap *DirectMap) HasBarcode(barcode string) bool {
	if directMap == nil {
		return false
	}
	_, ok := directMap.BarcodeToWriter[barcode]
	return ok
}

func (directMap *DirectMap) WriteRecord(barcode string, record *sam.Record) error {
	if directMap == nil {
		return nil
	}
	writer, ok := directMap.BarcodeToWriter[barcode]
	if !ok || writer == nil {
		return nil
	}
	return writer.Write(record)
}

func sanitizeLabel(label string) (string, bool) {
	if label == "" {
		return label, false
	}
	replaced := []rune(label)
	modified := false
	for index, value := range replaced {
		if value == '/' || value == '\\' || value == '~' {
			replaced[index] = '_'
			modified = true
			continue
		}
		if !isAllowedLabelRune(value) {
			replaced[index] = '_'
			modified = true
		}
	}

	if len(replaced) > 0 && replaced[0] == '.' {
		replaced[0] = '_'
		modified = true
	}

	sanitized := string(replaced)
	if strings.Contains(sanitized, "..") {
		sanitized = strings.ReplaceAll(sanitized, "..", "__")
		modified = true
	}

	return sanitized, modified
}

func isAllowedLabelRune(value rune) bool {
	if value == '_' || value == '-' || value == ' ' || value == '.' {
		return true
	}
	if value >= '0' && value <= '9' {
		return true
	}
	if value >= 'A' && value <= 'Z' {
		return true
	}
	if value >= 'a' && value <= 'z' {
		return true
	}
	return false
}
