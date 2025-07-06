//
// 3-Pass Algorithm for UMI-based Deduplication (Clean implementation)
//

#include "dedup_3pass.h"
#include "sort.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <sys/stat.h>
#include "htslib/bgzf.h"
#include "htslib/hts.h"
#include "htslib/tbx.h"

// Comparison function for sorting by molecule (CB, coord, strand, UB, MAPQ desc)
int compare_by_molecule(const void *a, const void *b) {
    const read_decision_t *read_a = (const read_decision_t *)a;
    const read_decision_t *read_b = (const read_decision_t *)b;
    
    // 1. Compare cell barcode
    int cmp = strcmp(read_a->cb, read_b->cb);
    if (cmp != 0) return cmp;
    
    // 2. Compare genomic coordinate
    if (read_a->coord != read_b->coord) {
        return (read_a->coord < read_b->coord) ? -1 : 1;
    }
    
    // 3. Compare strand
    if (read_a->strand != read_b->strand) {
        return read_a->strand - read_b->strand;
    }
    
    // 4. Compare UMI
    cmp = strcmp(read_a->ub, read_b->ub);
    if (cmp != 0) return cmp;
    
    // 5. Compare MAPQ (higher quality first - descending order)
    if (read_a->mapq != read_b->mapq) {
        return read_b->mapq - read_a->mapq;
    }
    
    // 6. Tie-breaker: read index (for stable sorting)
    if (read_a->read_idx != read_b->read_idx) {
        return (read_a->read_idx < read_b->read_idx) ? -1 : 1;
    }
    
    return 0;
}

// Comparison function for sorting by read index (to restore original order)
int compare_by_read_idx(const void *a, const void *b) {
    const read_decision_t *read_a = (const read_decision_t *)a;
    const read_decision_t *read_b = (const read_decision_t *)b;
    
    if (read_a->read_idx != read_b->read_idx) {
        return (read_a->read_idx < read_b->read_idx) ? -1 : 1;
    }
    
    return 0;
}

// Estimate initial capacity based on BAM file size
uint64_t estimate_capacity_from_file_size(const char* bampath) {
    struct stat st;
    if (stat(bampath, &st) != 0) {
        log_msg("Cannot stat BAM file, using default capacity", WARNING);
        return 1000000;
    }
    
    // Heuristic: ~20 bytes per read in compressed BAM (more realistic for modern data)
    uint64_t estimated_reads = st.st_size / 20;
    
    // Add 50% buffer for safety (fewer reallocations)
    estimated_reads = (estimated_reads * 3) / 2;
    
    // Clamp to reasonable bounds
    const uint64_t MIN_CAPACITY = 10000;
    const uint64_t MAX_CAPACITY = 200000000;
    
    if (estimated_reads < MIN_CAPACITY) return MIN_CAPACITY;
    if (estimated_reads > MAX_CAPACITY) return MAX_CAPACITY;
    
    log_msg("Estimated %llu reads from file size %lld bytes (%.1f MB initial allocation)", 
            DEBUG, estimated_reads, (long long)st.st_size, 
            (estimated_reads * sizeof(read_decision_t)) / (1024.0 * 1024.0));
    
    return estimated_reads;
}

// Create a new region decisions container
region_decisions_t *create_region_decisions(uint64_t initial_capacity) {
    region_decisions_t *region = malloc(sizeof(region_decisions_t));
    if (!region) {
        log_msg("Failed to allocate region_decisions_t", ERROR);
        return NULL;
    }
    
    region->decisions = malloc(initial_capacity * sizeof(read_decision_t));
    if (!region->decisions) {
        log_msg("Failed to allocate decisions array", ERROR);
        free(region);
        return NULL;
    }
    
    region->capacity = initial_capacity;
    region->count = 0;
    
    return region;
}

// Destroy region decisions container
void destroy_region_decisions(region_decisions_t *region) {
    if (region) {
        free(region->decisions);
        free(region);
    }
}

// Pass 1: Extract minimal information from BAM file
int extract_region_decisions(samFile *fp, sam_hdr_t *header, 
                           region_decisions_t *region, 
                           dedup_context_t *ctx) {
    bam1_t *read = bam_init1();
    if (!read) {
        log_msg("Failed to initialize BAM read", ERROR);
        return -1;
    }
    
    uint64_t read_idx = 0;
    int read_stat;
    
    log_msg("Pass 1: Extracting read information", INFO);
    
    while ((read_stat = sam_read1(fp, header, read)) >= 0) {
        // Check if we need to expand the array
        if (region->count >= region->capacity) {
            // Check for overflow before doubling capacity
            if (region->capacity > UINT64_MAX / 2) {
                log_msg("Cannot expand decisions array: capacity would overflow", ERROR);
                bam_destroy1(read);
                return -1;
            }
            uint64_t new_capacity = region->capacity * 2;
            
            // Additional check for multiplication overflow with sizeof
            if (new_capacity > UINT64_MAX / sizeof(read_decision_t)) {
                log_msg("Cannot expand decisions array: allocation size would overflow", ERROR);
                bam_destroy1(read);
                return -1;
            }
            
            read_decision_t *new_decisions = realloc(region->decisions, 
                                                   new_capacity * sizeof(read_decision_t));
            if (!new_decisions) {
                log_msg("Failed to expand decisions array", ERROR);
                bam_destroy1(read);
                return -1;
            }
            region->decisions = new_decisions;
            region->capacity = new_capacity;
            log_msg("Expanded decisions array to %llu entries", DEBUG, new_capacity);
        }
        
        read_decision_t *decision = &region->decisions[region->count];
        decision->read_idx = read_idx;
        
        // Extract cell barcode
        char cb_temp[CB_LENGTH];
        int8_t cb_stat = get_CB(read, ctx->cb_meta, cb_temp);
        if (cb_stat != 0) {
            // Skip reads without valid CB
            read_idx++;
            continue;
        }
        strncpy(decision->cb, cb_temp, sizeof(decision->cb) - 1);
        decision->cb[sizeof(decision->cb) - 1] = '\0';
        
        // Extract UMI
        char ub_temp[UB_LENGTH];
        int8_t ub_stat = get_UB(read, ctx->ub_meta, ub_temp);
        if (ub_stat != 0) {
            // Skip reads without valid UB
            read_idx++;
            continue;
        }
        strncpy(decision->ub, ub_temp, sizeof(decision->ub) - 1);
        decision->ub[sizeof(decision->ub) - 1] = '\0';
        
        // Check MAPQ threshold
        decision->mapq = read->core.qual;
        if (decision->mapq < ctx->mapq_threshold) {
            // Skip reads below MAPQ threshold
            read_idx++;
            continue;
        }
        
        // Check for secondary alignment (0x100 flag)
        if (read->core.flag & 0x100) {
            // Skip secondary alignments
            read_idx++;
            continue;
        }
        
        // Extract genomic coordinate and strand
        decision->coord = read->core.pos;
        decision->strand = bam_is_rev(read) ? 1 : 0;
        
        // Check if cell barcode exists in metadata (skip if not found)
        cb2fp *cluster_entry;
        HASH_FIND_STR(ctx->direct_map, decision->cb, cluster_entry);
        if (!cluster_entry) {
            // Skip reads not in any cluster
            read_idx++;
            continue;
        }
        
        // Initialize as keep=true, will be updated in Pass 2
        decision->keep = true;
        
        region->count++;
        read_idx++;
    }
    
    log_msg("Pass 1 complete: %llu reads processed, %llu kept for deduplication", 
            INFO, read_idx, region->count);
    
    bam_destroy1(read);
    return (read_stat == -1) ? 0 : -1;  // -1 is normal EOF
}

// Pass 2: Mark duplicates in memory
void mark_duplicates_in_region(region_decisions_t *region) {
    if (region->count == 0) {
        log_msg("No reads to deduplicate", INFO);
        return;
    }
    
    log_msg("Pass 2: Sorting %llu reads by molecule", INFO, region->count);
    
    // Sort by molecule (CB, coord, strand, UB, MAPQ desc)
    qsort(region->decisions, region->count, sizeof(read_decision_t), compare_by_molecule);
    
    log_msg("Pass 2: Marking duplicates", INFO);
    
    // Mark duplicates: keep only the highest MAPQ read per molecule
    uint64_t duplicates_marked = 0;
    
    for (uint64_t i = 0; i < region->count; i++) {
        if (i == 0) {
            // First read is always kept
            continue;
        }
        
        read_decision_t *prev = &region->decisions[i - 1];
        read_decision_t *curr = &region->decisions[i];
        
        // Check if this read is from the same molecule as the previous one
        bool same_molecule = (
            strcmp(prev->cb, curr->cb) == 0 &&
            prev->coord == curr->coord &&
            prev->strand == curr->strand &&
            strcmp(prev->ub, curr->ub) == 0
        );
        
        if (same_molecule) {
            // This is a duplicate - mark for removal
            curr->keep = false;
            duplicates_marked++;
        }
    }
    
    log_msg("Pass 2: Marked %llu duplicates for removal", INFO, duplicates_marked);
    
    // Sort back by read index to restore original order
    log_msg("Pass 2: Restoring original read order", DEBUG);
    qsort(region->decisions, region->count, sizeof(read_decision_t), compare_by_read_idx);
    
    log_msg("Pass 2 complete: %llu reads to keep, %llu duplicates to discard", 
            INFO, region->count - duplicates_marked, duplicates_marked);
}

// Pass 3: Write deduplicated reads to output files
int write_deduplicated_region(samFile *fp, sam_hdr_t *header,
                            region_decisions_t *region,
                            cb2fp *direct_map,
                            tag_meta_t *cb_meta) {
    
    // No need to read header - we've already seeked to data start position
    
    bam1_t *read = bam_init1();
    if (!read) {
        log_msg("Failed to initialize BAM read for Pass 3", ERROR);
        return -1;
    }
    
    uint64_t read_idx = 0;
    uint64_t decision_idx = 0;
    uint64_t reads_written = 0;
    uint64_t reads_skipped = 0;
    int read_stat;
    
    log_msg("Pass 3: Writing deduplicated reads to output files", INFO);
    log_msg("Pass 3: Total decisions: %llu", DEBUG, region->count);
    
    while ((read_stat = sam_read1(fp, header, read)) >= 0) {
        // Find the corresponding decision
        while (decision_idx < region->count && 
               region->decisions[decision_idx].read_idx < read_idx) {
            decision_idx++;
        }
        
        bool should_write = false;
        
        if (decision_idx < region->count && 
            region->decisions[decision_idx].read_idx == read_idx) {
            // We have a decision for this read
            should_write = region->decisions[decision_idx].keep;
        }
        
        if (should_write) {
            // Extract cell barcode and use read_dump like non-deduplication path
            char this_CB[CB_LENGTH];
            int8_t cb_stat = get_CB(read, cb_meta, this_CB);
            if (cb_stat == 0) {
                // Use read_dump exactly like non-deduplication path
                int8_t rdump_stat = read_dump(direct_map, this_CB, header, read);
                if (rdump_stat == 0) {
                    reads_written++;
                } else {
                    log_msg("Failed to write read using read_dump", ERROR);
                    reads_skipped++;
                }
            } else {
                log_msg("Failed to extract cell barcode for output", ERROR);
                reads_skipped++;
            }
        } else {
            reads_skipped++;
        }
        
        read_idx++;
    }
    
    log_msg("Pass 3 complete: %llu reads written, %llu reads skipped", 
            INFO, reads_written, reads_skipped);
    
    bam_destroy1(read);
    return (read_stat == -1) ? 0 : -1;  // -1 is normal EOF
}

// Main 3-pass deduplication function
int dedup_3pass(const char *bampath, sam_hdr_t *header, 
               cb2fp *direct_map,
               tag_meta_t *cb_meta, tag_meta_t *ub_meta,
               int16_t mapq_threshold) {
    
    log_msg("Starting 3-pass deduplication algorithm", INFO);
    
    // Open BAM file for Pass 1
    samFile *fp = sam_open(bampath, "r");
    if (!fp) {
        log_msg("Failed to open BAM file for Pass 1", ERROR);
        return -1;
    }
    
    // Skip the header for Pass 1 and save position
    sam_hdr_t *temp_header = sam_hdr_read(fp);
    if (!temp_header) {
        log_msg("Failed to read header in Pass 1", ERROR);
        sam_close(fp);
        return -1;
    }
    sam_hdr_destroy(temp_header);
    
    // Save the position after header for later seeking
    int64_t data_start_offset = bgzf_tell(hts_get_bgzfp(fp));
    if (data_start_offset < 0) {
        log_msg("Failed to get file position after header", ERROR);
        sam_close(fp);
        return -1;
    }
    
    // Create region decisions container with estimated capacity
    uint64_t initial_capacity = estimate_capacity_from_file_size(bampath);
    region_decisions_t *region = create_region_decisions(initial_capacity);
    if (!region) {
        sam_close(fp);
        return -1;
    }
    
    // Create deduplication context
    dedup_context_t ctx = {
        .region = region,
        .direct_map = direct_map,
        .cb_meta = cb_meta,
        .ub_meta = ub_meta,
        .mapq_threshold = mapq_threshold
    };
    
    // Pass 1: Extract minimal information
    if (extract_region_decisions(fp, header, region, &ctx) != 0) {
        log_msg("Pass 1 failed", ERROR);
        destroy_region_decisions(region);
        sam_close(fp);
        return -1;
    }
    
    // Pass 2: Mark duplicates in memory (no file I/O)
    mark_duplicates_in_region(region);
    
    // Pass 3: Seek back to data start and write deduplicated reads
    if (bgzf_seek(hts_get_bgzfp(fp), data_start_offset, SEEK_SET) < 0) {
        log_msg("Failed to seek back to data start for Pass 3", ERROR);
        destroy_region_decisions(region);
        sam_close(fp);
        return -1;
    }
    
    if (write_deduplicated_region(fp, header, region, direct_map, cb_meta) != 0) {
        log_msg("Pass 3 failed", ERROR);
        destroy_region_decisions(region);
        sam_close(fp);
        return -1;
    }
    
    sam_close(fp);
    
    log_msg("3-pass deduplication completed successfully", INFO);
    
    destroy_region_decisions(region);
    return 0;
}