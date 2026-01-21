package dedup

import (
	"bufio"
	"bytes"
	"encoding/binary"
	"fmt"
	"io"
	"os"
	"sort"
	"strconv"

	"github.com/biogo/hts/bam"
	"github.com/biogo/hts/sam"

	"scbamop/internal/config"
	"scbamop/internal/logging"
	"scbamop/internal/metadata"
	"scbamop/internal/tag"
)

const (
	recordHeaderSize          = 22
	defaultDedupChunkBytes    = int64(4 * 1024 * 1024 * 1024)
	dedupChunkBytesEnv        = "SCBAMOP_DEDUP_CHUNK_BYTES"
	initialChunkCapacityBytes = 4 * 1024 * 1024
)

type decisionBitset struct {
	bits []uint64
}

func newDecisionBitset(size uint64) decisionBitset {
	if size == 0 {
		return decisionBitset{}
	}
	return decisionBitset{bits: make([]uint64, (size+63)/64)}
}

func (bitset decisionBitset) Set(index uint64) {
	if len(bitset.bits) == 0 {
		return
	}
	slot := index / 64
	if slot >= uint64(len(bitset.bits)) {
		return
	}
	bitset.bits[slot] |= 1 << (index % 64)
}

func (bitset decisionBitset) Has(index uint64) bool {
	if len(bitset.bits) == 0 {
		return false
	}
	slot := index / 64
	if slot >= uint64(len(bitset.bits)) {
		return false
	}
	return bitset.bits[slot]&(1<<(index%64)) != 0
}

func buildDecisionBitset(inputPath string, header *sam.Header, directMap *metadata.DirectMap, cellMeta config.TagMeta, umiMeta config.TagMeta, mapqThreshold int, logger logging.Logger) (decisionBitset, uint64, error) {
	fileHandle, openErr := os.Open(inputPath)
	if openErr != nil {
		return decisionBitset{}, 0, fmt.Errorf("failed to open BAM file: %w", openErr)
	}
	defer fileHandle.Close()

	reader, readerErr := bam.NewReader(fileHandle, 0)
	if readerErr != nil {
		return decisionBitset{}, 0, fmt.Errorf("failed to read BAM header: %w", readerErr)
	}
	defer reader.Close()

	if header == nil {
		header = reader.Header()
	}

	recordSize, sizeErr := decisionRecordSize(cellMeta.Length, umiMeta.Length)
	if sizeErr != nil {
		return decisionBitset{}, 0, sizeErr
	}
	chunkLimit, limitErr := dedupChunkLimit()
	if limitErr != nil {
		return decisionBitset{}, 0, limitErr
	}
	if int64(recordSize) > chunkLimit {
		return decisionBitset{}, 0, fmt.Errorf("dedup chunk size (%d bytes) is smaller than decision record size (%d bytes)", chunkLimit, recordSize)
	}
	maxInt := int64(^uint(0) >> 1)
	if chunkLimit > maxInt {
		return decisionBitset{}, 0, fmt.Errorf("dedup chunk size (%d bytes) exceeds supported memory limit", chunkLimit)
	}
	maxChunkBytes := int(chunkLimit)
	data := make([]byte, 0, minInt(maxChunkBytes, initialChunkCapacityBytes))
	var chunkPaths []string
	defer func() {
		for _, path := range chunkPaths {
			_ = os.Remove(path)
		}
	}()

	var readIndex uint64
	var keptCount uint64

	logger.Logf(logging.Info, "Pass 1: Extracting read information")

	for {
		record, readErr := reader.Read()
		if readErr == io.EOF {
			break
		}
		if readErr != nil {
			return decisionBitset{}, readIndex, readErr
		}

		cellBarcode, tagErr := tag.Extract(record, cellMeta)
		if tagErr != nil {
			readIndex++
			continue
		}
		umiValue, tagErr := tag.Extract(record, umiMeta)
		if tagErr != nil {
			readIndex++
			continue
		}
		if int(record.MapQ) < mapqThreshold {
			readIndex++
			continue
		}
		if record.Flags&(sam.Secondary|sam.Supplementary) != 0 {
			readIndex++
			continue
		}
		if !directMap.HasBarcode(cellBarcode) {
			readIndex++
			continue
		}
		if lengthErr := tag.ValidateLength(cellBarcode, cellMeta.Length, "cell barcode"); lengthErr != nil {
			return decisionBitset{}, readIndex, fmt.Errorf("read %s: %w", record.Name, lengthErr)
		}
		if lengthErr := tag.ValidateLength(umiValue, umiMeta.Length, "UMI"); lengthErr != nil {
			return decisionBitset{}, readIndex, fmt.Errorf("read %s: %w", record.Name, lengthErr)
		}

		if len(data)+recordSize > maxChunkBytes {
			path, spillErr := spillDecisionChunk(data, recordSize, cellMeta.Length, umiMeta.Length)
			if spillErr != nil {
				return decisionBitset{}, readIndex, spillErr
			}
			chunkPaths = append(chunkPaths, path)
			data = data[:0]
		}

		data = ensureDecisionCapacity(data, recordSize, maxChunkBytes)
		recordStart := len(data)
		data = append(data, make([]byte, recordSize)...)
		recordBytes := data[recordStart : recordStart+recordSize]

		referenceID := int32(-1)
		if record.Ref != nil {
			referenceID = int32(record.Ref.ID())
		}
		strand := uint8(0)
		if record.Flags&sam.Reverse != 0 {
			strand = 1
		}

		encodeDecisionRecord(recordBytes, readIndex, referenceID, int64(record.Pos), strand, record.MapQ, cellBarcode, umiValue, cellMeta.Length, umiMeta.Length)
		keptCount++
		readIndex++
	}

	if len(data) > 0 {
		path, spillErr := spillDecisionChunk(data, recordSize, cellMeta.Length, umiMeta.Length)
		if spillErr != nil {
			return decisionBitset{}, readIndex, spillErr
		}
		chunkPaths = append(chunkPaths, path)
	}

	logger.Logf(logging.Info, "Pass 1 complete: %d reads processed, %d kept for deduplication", readIndex, keptCount)

	decisions := newDecisionBitset(readIndex)
	if len(chunkPaths) == 0 {
		return decisions, readIndex, nil
	}

	logger.Logf(logging.Info, "Pass 2: Deduplicating %d decision chunk(s)", len(chunkPaths))
	mergeErr := mergeDecisionChunks(chunkPaths, recordSize, cellMeta.Length, umiMeta.Length, &decisions)
	if mergeErr != nil {
		return decisions, readIndex, mergeErr
	}

	return decisions, readIndex, nil
}

func decisionRecordSize(cellLength int, umiLength int) (int, error) {
	if cellLength <= 0 {
		return 0, fmt.Errorf("cell barcode length must be larger than 0")
	}
	if umiLength <= 0 {
		return 0, fmt.Errorf("UMI length must be larger than 0")
	}
	return recordHeaderSize + cellLength + umiLength, nil
}

func encodeDecisionRecord(buffer []byte, readIndex uint64, referenceID int32, position int64, strand uint8, mapq uint8, cellBarcode string, umiValue string, cellLength int, umiLength int) {
	binary.LittleEndian.PutUint64(buffer[0:8], readIndex)
	binary.LittleEndian.PutUint32(buffer[8:12], uint32(referenceID))
	binary.LittleEndian.PutUint64(buffer[12:20], uint64(position))
	buffer[20] = strand
	buffer[21] = mapq
	cellStart := recordHeaderSize
	copy(buffer[cellStart:cellStart+cellLength], cellBarcode)
	umiStart := cellStart + cellLength
	copy(buffer[umiStart:umiStart+umiLength], umiValue)
}

func dedupChunkLimit() (int64, error) {
	override := os.Getenv(dedupChunkBytesEnv)
	if override == "" {
		return defaultDedupChunkBytes, nil
	}
	parsed, err := strconv.ParseInt(override, 10, 64)
	if err != nil || parsed <= 0 {
		return 0, fmt.Errorf("invalid %s: %s", dedupChunkBytesEnv, override)
	}
	return parsed, nil
}

func spillDecisionChunk(data []byte, recordSize int, cellLength int, umiLength int) (string, error) {
	if len(data) == 0 {
		return "", nil
	}
	sorter := &recordSorter{
		data:       data,
		recordSize: recordSize,
		cellLength: cellLength,
		umiLength:  umiLength,
	}
	sort.Sort(sorter)

	file, err := os.CreateTemp("", "scbamop-dedup-")
	if err != nil {
		return "", fmt.Errorf("failed to create dedup temp file: %w", err)
	}
	if err := writeAll(file, data); err != nil {
		_ = file.Close()
		_ = os.Remove(file.Name())
		return "", err
	}
	if closeErr := file.Close(); closeErr != nil {
		_ = os.Remove(file.Name())
		return "", closeErr
	}
	return file.Name(), nil
}

func writeAll(writer io.Writer, data []byte) error {
	written := 0
	for written < len(data) {
		count, err := writer.Write(data[written:])
		if err != nil {
			return err
		}
		written += count
	}
	return nil
}

type recordSorter struct {
	data       []byte
	recordSize int
	cellLength int
	umiLength  int
	swapBuf    []byte
}

func (sorter *recordSorter) Len() int {
	return len(sorter.data) / sorter.recordSize
}

func (sorter *recordSorter) Less(leftIndex int, rightIndex int) bool {
	left := recordAt(sorter.data, sorter.recordSize, leftIndex)
	right := recordAt(sorter.data, sorter.recordSize, rightIndex)
	return lessRecord(left, right, sorter.cellLength, sorter.umiLength)
}

func (sorter *recordSorter) Swap(leftIndex int, rightIndex int) {
	if leftIndex == rightIndex {
		return
	}
	leftStart := leftIndex * sorter.recordSize
	rightStart := rightIndex * sorter.recordSize
	left := sorter.data[leftStart : leftStart+sorter.recordSize]
	right := sorter.data[rightStart : rightStart+sorter.recordSize]
	if len(sorter.swapBuf) < sorter.recordSize {
		sorter.swapBuf = make([]byte, sorter.recordSize)
	}
	copy(sorter.swapBuf, left)
	copy(left, right)
	copy(right, sorter.swapBuf)
}

func mergeDecisionChunks(paths []string, recordSize int, cellLength int, umiLength int, decisions *decisionBitset) error {
	if len(paths) == 0 {
		return nil
	}

	heap := &recordHeap{cellLength: cellLength, umiLength: umiLength}
	readers := make([]*chunkReader, 0, len(paths))
	defer func() {
		for _, reader := range readers {
			if reader != nil {
				reader.closeAndRemove()
			}
		}
	}()
	for _, path := range paths {
		reader, err := newChunkReader(path, recordSize)
		if err != nil {
			if err == io.EOF {
				_ = os.Remove(path)
				continue
			}
			return err
		}
		heap.items = append(heap.items, reader)
		readers = append(readers, reader)
	}
	heap.Init()

	if heap.Len() == 0 {
		return nil
	}

	prevCell := make([]byte, cellLength)
	prevUMI := make([]byte, umiLength)
	var prevRefID int32
	var prevPos int64
	var prevStrand uint8
	hasPrev := false

	for heap.Len() > 0 {
		reader := heap.Pop()
		record := reader.record

		refID := decodeRefID(record)
		pos := decodePosition(record)
		strand := record[20]
		cellBytes := recordCell(record, cellLength)
		umiBytes := recordUMI(record, cellLength, umiLength)

		if !hasPrev || !sameDedupKey(cellBytes, umiBytes, refID, pos, strand, prevCell, prevUMI, prevRefID, prevPos, prevStrand) {
			decisions.Set(decodeReadIndex(record))
			copy(prevCell, cellBytes)
			copy(prevUMI, umiBytes)
			prevRefID = refID
			prevPos = pos
			prevStrand = strand
			hasPrev = true
		}

		if err := reader.advance(); err != nil {
			if err == io.EOF {
				reader.closeAndRemove()
				continue
			}
			reader.closeAndRemove()
			return err
		}
		heap.Push(reader)
	}

	return nil
}

type recordHeap struct {
	items      []*chunkReader
	cellLength int
	umiLength  int
}

func (heap *recordHeap) Init() {
	sort.Slice(heap.items, func(left int, right int) bool {
		return lessRecord(heap.items[left].record, heap.items[right].record, heap.cellLength, heap.umiLength)
	})
}

func (heap *recordHeap) Len() int {
	return len(heap.items)
}

func (heap *recordHeap) Less(left int, right int) bool {
	return lessRecord(heap.items[left].record, heap.items[right].record, heap.cellLength, heap.umiLength)
}

func (heap *recordHeap) Swap(left int, right int) {
	heap.items[left], heap.items[right] = heap.items[right], heap.items[left]
}

func (heap *recordHeap) Push(reader *chunkReader) {
	heap.items = append(heap.items, reader)
	heap.up(len(heap.items) - 1)
}

func (heap *recordHeap) Pop() *chunkReader {
	if len(heap.items) == 0 {
		return nil
	}
	lastIndex := len(heap.items) - 1
	heap.Swap(0, lastIndex)
	item := heap.items[lastIndex]
	heap.items = heap.items[:lastIndex]
	heap.down(0)
	return item
}

func (heap *recordHeap) up(index int) {
	for {
		parent := (index - 1) / 2
		if index == 0 || !heap.Less(index, parent) {
			break
		}
		heap.Swap(index, parent)
		index = parent
	}
}

func (heap *recordHeap) down(index int) {
	for {
		left := 2*index + 1
		right := left + 1
		smallest := index
		if left < len(heap.items) && heap.Less(left, smallest) {
			smallest = left
		}
		if right < len(heap.items) && heap.Less(right, smallest) {
			smallest = right
		}
		if smallest == index {
			break
		}
		heap.Swap(index, smallest)
		index = smallest
	}
}

type chunkReader struct {
	path       string
	file       *os.File
	reader     *bufio.Reader
	record     []byte
	recordSize int
}

func newChunkReader(path string, recordSize int) (*chunkReader, error) {
	file, err := os.Open(path)
	if err != nil {
		return nil, fmt.Errorf("failed to open dedup temp file: %w", err)
	}
	chunk := &chunkReader{
		path:       path,
		file:       file,
		reader:     bufio.NewReader(file),
		record:     make([]byte, recordSize),
		recordSize: recordSize,
	}
	if err := chunk.advance(); err != nil {
		if err == io.EOF {
			chunk.closeAndRemove()
			return nil, io.EOF
		}
		chunk.closeAndRemove()
		return nil, err
	}
	return chunk, nil
}

func (chunk *chunkReader) advance() error {
	if chunk.reader == nil {
		return io.EOF
	}
	_, err := io.ReadFull(chunk.reader, chunk.record)
	if err == io.EOF {
		return io.EOF
	}
	if err == io.ErrUnexpectedEOF {
		return fmt.Errorf("unexpected end of chunk file %s", chunk.path)
	}
	if err != nil {
		return err
	}
	return nil
}

func (chunk *chunkReader) closeAndRemove() {
	if chunk.file != nil {
		_ = chunk.file.Close()
		chunk.file = nil
	}
	if chunk.path != "" {
		_ = os.Remove(chunk.path)
		chunk.path = ""
	}
}

func recordAt(data []byte, recordSize int, index int) []byte {
	start := index * recordSize
	return data[start : start+recordSize]
}

func decodeReadIndex(record []byte) uint64 {
	return binary.LittleEndian.Uint64(record[0:8])
}

func decodeRefID(record []byte) int32 {
	return int32(binary.LittleEndian.Uint32(record[8:12]))
}

func decodePosition(record []byte) int64 {
	return int64(binary.LittleEndian.Uint64(record[12:20]))
}

func recordCell(record []byte, cellLength int) []byte {
	return record[recordHeaderSize : recordHeaderSize+cellLength]
}

func recordUMI(record []byte, cellLength int, umiLength int) []byte {
	start := recordHeaderSize + cellLength
	return record[start : start+umiLength]
}

func lessRecord(left []byte, right []byte, cellLength int, umiLength int) bool {
	leftCell := recordCell(left, cellLength)
	rightCell := recordCell(right, cellLength)
	if cmp := bytes.Compare(leftCell, rightCell); cmp != 0 {
		return cmp < 0
	}

	leftRefID := decodeRefID(left)
	rightRefID := decodeRefID(right)
	if leftRefID != rightRefID {
		return leftRefID < rightRefID
	}

	leftPos := decodePosition(left)
	rightPos := decodePosition(right)
	if leftPos != rightPos {
		return leftPos < rightPos
	}

	leftStrand := left[20]
	rightStrand := right[20]
	if leftStrand != rightStrand {
		return leftStrand < rightStrand
	}

	leftUMI := recordUMI(left, cellLength, umiLength)
	rightUMI := recordUMI(right, cellLength, umiLength)
	if cmp := bytes.Compare(leftUMI, rightUMI); cmp != 0 {
		return cmp < 0
	}

	leftMapq := left[21]
	rightMapq := right[21]
	if leftMapq != rightMapq {
		return leftMapq > rightMapq
	}

	return decodeReadIndex(left) < decodeReadIndex(right)
}

func sameDedupKey(cellBytes []byte, umiBytes []byte, refID int32, pos int64, strand uint8, prevCell []byte, prevUMI []byte, prevRefID int32, prevPos int64, prevStrand uint8) bool {
	if refID != prevRefID || pos != prevPos || strand != prevStrand {
		return false
	}
	if !bytes.Equal(cellBytes, prevCell) {
		return false
	}
	if !bytes.Equal(umiBytes, prevUMI) {
		return false
	}
	return true
}

func minInt(left int, right int) int {
	if left < right {
		return left
	}
	return right
}

func ensureDecisionCapacity(data []byte, recordSize int, maxChunkBytes int) []byte {
	needed := len(data) + recordSize
	if needed <= cap(data) {
		return data
	}
	newCap := cap(data) * 2
	if newCap < needed {
		newCap = needed
	}
	if newCap > maxChunkBytes {
		newCap = maxChunkBytes
	}
	if newCap < needed {
		newCap = needed
	}
	resized := make([]byte, len(data), newCap)
	copy(resized, data)
	return resized
}
