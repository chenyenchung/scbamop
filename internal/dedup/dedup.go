package dedup

import (
	"fmt"
	"io"
	"os"
	"sort"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"

	"scbamop/internal/config"
	"scbamop/internal/logging"
	"scbamop/internal/metadata"
	"scbamop/internal/tag"
)

type readDecision struct {
	ReadIndex   uint64
	CellBarcode string
	UMI         string
	ReferenceID int
	Coordinate  int
	Strand      uint8
	MapQ        uint8
	Keep        bool
}

func Deduplicate(inputPath string, header *sam.Header, directMap *metadata.DirectMap, cellMeta config.TagMeta, umiMeta config.TagMeta, mapqThreshold int, logger logging.Logger) error {
	if directMap == nil {
		return fmt.Errorf("missing metadata mapping")
	}

	decisions, readCount, extractErr := buildDecisionBitset(inputPath, header, directMap, cellMeta, umiMeta, mapqThreshold, logger)
	if extractErr != nil {
		return extractErr
	}

	writeErr := writeDeduplicated(inputPath, header, directMap, cellMeta, &decisions, logger, readCount)
	if writeErr != nil {
		return writeErr
	}

	return nil
}

func writeDeduplicated(inputPath string, header *sam.Header, directMap *metadata.DirectMap, cellMeta config.TagMeta, decisions *decisionBitset, logger logging.Logger, totalReads uint64) error {
	fileHandle, openErr := os.Open(inputPath)
	if openErr != nil {
		return fmt.Errorf("failed to open BAM file for pass 3: %w", openErr)
	}
	defer fileHandle.Close()

	reader, readerErr := bam.NewReader(fileHandle, 0)
	if readerErr != nil {
		return fmt.Errorf("failed to read BAM header for pass 3: %w", readerErr)
	}
	defer reader.Close()

	if header == nil {
		header = reader.Header()
	}

	logger.Logf(logging.Info, "Pass 3: Writing deduplicated reads to output files")

	var readIndex uint64
	readsWritten := 0
	readsSkipped := 0

	for {
		record, readErr := reader.Read()
		if readErr == io.EOF {
			break
		}
		if readErr != nil {
			return readErr
		}

		shouldWrite := decisions != nil && decisions.Has(readIndex)
		if shouldWrite {
			cellBarcode, tagErr := tag.Extract(record, cellMeta)
			if tagErr == nil {
				writeErr := directMap.WriteRecord(cellBarcode, record)
				if writeErr != nil {
					logger.Logf(logging.Error, "Failed to write read: %v", writeErr)
					readsSkipped++
				} else {
					readsWritten++
				}
			} else {
				logger.Logf(logging.Error, "Failed to extract cell barcode for output: %v", tagErr)
				readsSkipped++
			}
		} else {
			readsSkipped++
		}

		readIndex++
	}

	logger.Logf(logging.Info, "Pass 3 complete: %d reads written, %d reads skipped", readsWritten, readsSkipped)
	if totalReads > readIndex {
		logger.Logf(logging.Warning, "Pass 3 processed fewer reads (%d) than expected (%d)", readIndex, totalReads)
	}

	return nil
}

func markDuplicates(decisions []readDecision) {
	if len(decisions) == 0 {
		return
	}

	sort.Slice(decisions, func(leftIndex int, rightIndex int) bool {
		left := decisions[leftIndex]
		right := decisions[rightIndex]
		if left.CellBarcode != right.CellBarcode {
			return left.CellBarcode < right.CellBarcode
		}
		if left.ReferenceID != right.ReferenceID {
			return left.ReferenceID < right.ReferenceID
		}
		if left.Coordinate != right.Coordinate {
			return left.Coordinate < right.Coordinate
		}
		if left.Strand != right.Strand {
			return left.Strand < right.Strand
		}
		if left.UMI != right.UMI {
			return left.UMI < right.UMI
		}
		if left.MapQ != right.MapQ {
			return left.MapQ > right.MapQ
		}
		return left.ReadIndex < right.ReadIndex
	})

	for decisionIndex := 1; decisionIndex < len(decisions); decisionIndex++ {
		previous := decisions[decisionIndex-1]
		current := &decisions[decisionIndex]
		if previous.CellBarcode == current.CellBarcode &&
			previous.ReferenceID == current.ReferenceID &&
			previous.Coordinate == current.Coordinate &&
			previous.Strand == current.Strand &&
			previous.UMI == current.UMI {
			current.Keep = false
		}
	}

	sort.Slice(decisions, func(leftIndex int, rightIndex int) bool {
		return decisions[leftIndex].ReadIndex < decisions[rightIndex].ReadIndex
	})
}
