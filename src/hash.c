//
// Created by Yen-Chung Chen on 2/23/23.
//
#include "hash.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <ctype.h>
#include "utils.h"

// Count unique labels in metadata file for pre-allocation
uint32_t count_unique_labels(const char *path) {
    FILE* fp = fopen(path, "r");
    if (!fp) return 0;
    
    // Simple hash table for counting unique labels
    typedef struct {
        char label[64];
        UT_hash_handle hh;
    } label_count_t;
    label_count_t *labels = NULL;
    
    char line[MAX_LINE_LENGTH];
    bool first_line = true;
    uint32_t unique_count = 0;
    
    while (fgets(line, MAX_LINE_LENGTH, fp) != NULL) {
        if (first_line) {
            first_line = false;
            continue;
        }
        
        line[strcspn(line, "\n")] = 0;
        char *tokens = strtok(line, ",");
        
        // Skip to second field (label)
        if (tokens) tokens = strtok(NULL, ",");
        if (!tokens) continue;
        
        // Check if we've seen this label before
        label_count_t *existing;
        HASH_FIND_STR(labels, tokens, existing);
        if (!existing) {
            label_count_t *new_label = calloc(1, sizeof(label_count_t));
            if (new_label) {
                strncpy(new_label->label, tokens, sizeof(new_label->label) - 1);
                new_label->label[sizeof(new_label->label) - 1] = '\0';
                HASH_ADD_STR(labels, label, new_label);
                unique_count++;
            }
        }
    }
    
    // Cleanup
    label_count_t *current, *tmp;
    HASH_ITER(hh, labels, current, tmp) {
        HASH_DEL(labels, current);
        free(current);
    }
    
    fclose(fp);
    log_msg("Pre-counted %u unique labels", DEBUG, unique_count);
    return unique_count;
}

cb2fp* hash_readtag_direct(char *path, const char *prefix, sam_hdr_t *header) {
    // Initialize all resources to NULL for cleanup
    FILE* meta_fp = NULL;
    cb2fp *direct_map = NULL;
    cb2fp *direct_entry = NULL;
    int ret = 0;  // 0 for success, -1 for error
    
    // Temporary hash table to track unique labels and their file pointers
    typedef struct {
        char label[64];                   /* consistent with cb2fp label size */
        samFile* fp;
        UT_hash_handle hh;
    } label_to_fp_t;
    label_to_fp_t *label_fps = NULL;
    
    // Hash table to track which labels we've already warned about
    typedef struct {
        char label[MAX_LINE_LENGTH];
        UT_hash_handle hh;
    } warned_label_t;
    warned_label_t *warned_labels = NULL;
    
    // Open metadata file
    meta_fp = fopen(path, "r");
    if (meta_fp == NULL) {
        // Exit and print error message if the file does not exist
        log_msg("Cannot open file (%s)", ERROR, path);
        return NULL;
    }
    
    // Pre-count unique labels for better memory allocation
    uint32_t estimated_labels = count_unique_labels(path);
    log_msg("Estimated %u unique labels for hash table pre-allocation", DEBUG, estimated_labels);

    // Read every line
    char meta_line[MAX_LINE_LENGTH]; // Temporary variable to store each line
    bool first_line = true;
    char* tokens; // Temporary variable for iterating tokens
    uint32_t field_num = 0; // Counting numbers to examine if expected field numbers are present
    char trt[MAX_LINE_LENGTH]; // Temporary variable for read tag content
    char tlabel[MAX_LINE_LENGTH]; // Temporary variable for corresponding label content

    while (fgets(meta_line, MAX_LINE_LENGTH, meta_fp) != NULL) {
        // Assuming header and skip it
        if (first_line) {
            first_line = false;
            continue;
        }

        // Strip the linebreak
        meta_line[strcspn(meta_line, "\n")] = 0;

        // Tokenize by comma
        tokens = strtok(meta_line, ",");

        // Reset field number
        field_num = 0;

        // Iterate through tokens
        while (tokens != NULL) {
            switch(field_num) {
                case 0:
                    // Expect the first field to be read tags
                    strncpy(trt, tokens, MAX_LINE_LENGTH - 1);
                    trt[MAX_LINE_LENGTH - 1] = '\0';
                    break;
                case 1:
                    // Expect the second to be labels
                    strncpy(tlabel, tokens, MAX_LINE_LENGTH - 1);
                    tlabel[MAX_LINE_LENGTH - 1] = '\0';
                    break;
                default:
                    log_msg("There are %d fields in the metadata but only 2 are expected",
                            ERROR, field_num);
                    ret = -1;
                    goto cleanup;
            }
            field_num++;
            tokens = strtok(NULL, ",");
        }

        // Deal with metadata that has < 2 fields
        if (field_num < 2) {
            log_msg("There is only %d field in the metadata (expecting 2)", ERROR, field_num);
            ret = -1;
            goto cleanup;
        }

        // Sanitize label by replacing invalid characters with underscores
        char original_label[MAX_LINE_LENGTH];
        strncpy(original_label, tlabel, MAX_LINE_LENGTH - 1);
        original_label[MAX_LINE_LENGTH - 1] = '\0';
        int label_modified = 0;
        
        // Replace path traversal sequences and invalid characters
        for (char *p = tlabel; *p; p++) {
            if (*p == '/' || *p == '\\' || *p == '~') {
                *p = '_';
                label_modified = 1;
            } else if (!isalnum(*p) && *p != '_' && *p != '-' && *p != ' ' && *p != '.') {
                *p = '_';
                label_modified = 1;
            }
        }
        
        // Handle leading dots (hidden files)
        if (tlabel[0] == '.') {
            tlabel[0] = '_';
            label_modified = 1;
        }
        
        // Handle ".." sequences
        char *dot_dot = strstr(tlabel, "..");
        while (dot_dot) {
            dot_dot[0] = '_';
            dot_dot[1] = '_';
            label_modified = 1;
            dot_dot = strstr(dot_dot + 2, "..");
        }
        
        // Log if label was sanitized (only once per unique original label)
        if (label_modified) {
            warned_label_t *existing_warning;
            HASH_FIND_STR(warned_labels, original_label, existing_warning);
            
            if (!existing_warning) {
                // First time seeing this label, add warning
                log_msg("Sanitized label: '%s' -> '%s'", WARNING, original_label, tlabel);
                
                // Track that we've warned about this label
                warned_label_t *new_warning = calloc(1, sizeof(warned_label_t));
                if (new_warning) {
                    strncpy(new_warning->label, original_label, sizeof(new_warning->label) - 1);
                    new_warning->label[sizeof(new_warning->label) - 1] = '\0';
                    HASH_ADD_STR(warned_labels, label, new_warning);
                }
            }
        }

        // Check if we already have a file pointer for this label
        label_to_fp_t *existing_label;
        HASH_FIND_STR(label_fps, tlabel, existing_label);
        
        samFile* output_fp = NULL;
        
        if (existing_label == NULL) {
            // Create new output file for this label
            char output_path[512];
            size_t prefix_len = strlen(prefix);
            size_t label_len = strlen(tlabel);
            
            // Check if the combined path would exceed buffer size
            if (prefix_len + label_len + 5 >= sizeof(output_path)) {  // 5 = ".bam" + null terminator
                log_msg("Output path too long for label: %s", ERROR, tlabel);
                ret = -1;
                goto cleanup;
            }
            
            snprintf(output_path, sizeof(output_path), "%s%s.bam", prefix, tlabel);
            
            output_fp = sam_open(output_path, "wb");
            if (!output_fp) {
                log_msg("Failed to create output file: %s", ERROR, output_path);
                ret = -1;
                goto cleanup;
            }
            
            // Write header to the new file
            if (sam_hdr_write(output_fp, header) < 0) {
                log_msg("Failed to write header to: %s", ERROR, output_path);
                sam_close(output_fp);
                ret = -1;
                goto cleanup;
            }
            
            // Add to label tracking hash table
            label_to_fp_t *new_label = calloc(1, sizeof(label_to_fp_t));
            if (!new_label) {
                log_msg("Failed to allocate memory for label tracking", ERROR);
                sam_close(output_fp);
                ret = -1;
                goto cleanup;
            }
            
            strncpy(new_label->label, tlabel, sizeof(new_label->label) - 1);
            new_label->label[sizeof(new_label->label) - 1] = '\0';
            new_label->fp = output_fp;
            HASH_ADD_STR(label_fps, label, new_label);
            
            log_msg("Created output file: %s", INFO, output_path);
        } else {
            // Use existing file pointer
            output_fp = existing_label->fp;
        }

        // Create direct mapping entry: cell_barcode -> file_pointer
        direct_entry = calloc(1, sizeof(cb2fp));
        if (!direct_entry) {
            log_msg("Failed to allocate memory for direct mapping entry", ERROR);
            ret = -1;
            goto cleanup;
        }

        // Assign the cell barcode, label, and file pointer with bounds checking
        if (strlen(trt) >= sizeof(direct_entry->cb)) {
            log_msg("Cell barcode too long (max %zu chars): %s", ERROR, sizeof(direct_entry->cb) - 1, trt);
            free(direct_entry);
            direct_entry = NULL;
            ret = -1;
            goto cleanup;
        }
        if (strlen(tlabel) >= sizeof(direct_entry->label)) {
            log_msg("Label too long (max %zu chars): %s", ERROR, sizeof(direct_entry->label) - 1, tlabel);
            free(direct_entry);
            direct_entry = NULL;
            ret = -1;
            goto cleanup;
        }
        
        strncpy(direct_entry->cb, trt, sizeof(direct_entry->cb) - 1);
        direct_entry->cb[sizeof(direct_entry->cb) - 1] = '\0';
        strncpy(direct_entry->label, tlabel, sizeof(direct_entry->label) - 1);
        direct_entry->label[sizeof(direct_entry->label) - 1] = '\0';
        direct_entry->fp = output_fp;

        HASH_ADD_STR(direct_map, cb, direct_entry);
        direct_entry = NULL;  // Successfully added, don't free in cleanup
    }
    
cleanup:
    // Always close metadata file
    if (meta_fp) {
        fclose(meta_fp);
    }
    
    // Free any unadded direct_entry
    if (direct_entry) {
        free(direct_entry);
    }
    
    // Clean up label_fps hash table
    if (label_fps) {
        label_to_fp_t *current_label, *tmp_label;
        HASH_ITER(hh, label_fps, current_label, tmp_label) {
            HASH_DEL(label_fps, current_label);
            if (ret != 0 && current_label->fp) {
                // On error, close all file handles
                sam_close(current_label->fp);
            }
            free(current_label);
        }
    }
    
    // Clean up warned_labels hash table
    if (warned_labels) {
        warned_label_t *current_warning, *tmp_warning;
        HASH_ITER(hh, warned_labels, current_warning, tmp_warning) {
            HASH_DEL(warned_labels, current_warning);
            free(current_warning);
        }
    }
    
    // On error, clean up partial direct_map
    if (ret != 0 && direct_map) {
        cb2fp *entry, *tmp;
        HASH_ITER(hh, direct_map, entry, tmp) {
            HASH_DEL(direct_map, entry);
            // Note: Don't close file pointers here as they're shared with label_fps
            // and already closed above
            free(entry);
        }
        direct_map = NULL;
    }
    
    return direct_map;
}
