//
// 3-Pass Algorithm for UMI-based Deduplication (Clean implementation)
//
// Efficient in-memory deduplication without temporary files

#ifndef SCBAMSPLIT_DEDUP_3PASS_H
#define SCBAMSPLIT_DEDUP_3PASS_H

// Standard library includes
#include <stdint.h>
#include <stdbool.h>

// External library includes
#include "htslib/sam.h"

// Project includes
#include "utils.h"
#include "hash.h"

// Core data structure for read decisions
typedef struct {
    uint64_t read_idx;      // Position in original BAM (0-based)
    char cb[32];            // Cell barcode 
    char ub[32];            // UMI
    int32_t coord;          // Genomic position
    uint8_t strand;         // 0 for +, 1 for -
    uint8_t mapq;           // Mapping quality
    bool keep;              // Set in pass 2
} read_decision_t;

// Container for region-based processing
typedef struct {
    read_decision_t *decisions;     // Array of decisions
    uint64_t capacity;              // Allocated size
    uint64_t count;                 // Current number of decisions
} region_decisions_t;

// Context for deduplication operations
typedef struct {
    region_decisions_t *region;     // Current region being processed
    cb2fp *direct_map;              // Direct cell barcode to file pointer mapping
    tag_meta_t *cb_meta;            // Cell barcode metadata
    tag_meta_t *ub_meta;            // UMI metadata
    int16_t mapq_threshold;         // MAPQ threshold
} dedup_context_t;

// Comparison functions for qsort
int compare_by_molecule(const void *a, const void *b);
int compare_by_read_idx(const void *a, const void *b);

// Core functions
uint64_t estimate_capacity_from_file_size(const char* bampath);
region_decisions_t *create_region_decisions(uint64_t initial_capacity);
void destroy_region_decisions(region_decisions_t *region);

// 3-pass algorithm functions
int extract_region_decisions(samFile *fp, sam_hdr_t *header, 
                           region_decisions_t *region, 
                           dedup_context_t *ctx);

void mark_duplicates_in_region(region_decisions_t *region);

int write_deduplicated_region(samFile *fp, sam_hdr_t *header,
                            region_decisions_t *region,
                            cb2fp *direct_map,
                            tag_meta_t *cb_meta);

// Main deduplication function
int dedup_3pass(const char *bampath, sam_hdr_t *header, 
               cb2fp *direct_map,
               tag_meta_t *cb_meta, tag_meta_t *ub_meta,
               int16_t mapq_threshold);

#endif //SCBAMSPLIT_DEDUP_3PASS_H