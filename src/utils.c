//
// Utility functions (Clean implementation)
//

#include "utils.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <errno.h>
#include <ctype.h>
#include <sys/stat.h>

void show_global_usage() {
    fprintf(stderr, "Program: scbamop (Single-cell BAM operations toolkit)\n");
    fprintf(stderr, "Version: v0.5.0 (subcommand structure)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: scbamop <command> [options]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Commands:\n");
    fprintf(stderr, "  split    Split BAM file by cell barcodes with optional deduplication\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Use 'scbamop <command> --help' for command-specific help\n");
    fprintf(stderr, "\n");
}

void show_split_usage() {
    fprintf(stderr, "Usage: scbamop split -f FILE -m FILE [options]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Split BAM file by cell barcodes with optional UMI-based deduplication\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Required arguments:\n");
    fprintf(stderr, "  -f, --file FILE        Input BAM file path\n");
    fprintf(stderr, "  -m, --meta FILE        Metadata file with cell barcode assignments\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Optional arguments:\n");
    fprintf(stderr, "  -o, --output DIR       Output directory prefix (default: ./)\n");
    fprintf(stderr, "  -q, --mapq INT         MAPQ threshold (default: 0)\n");
    fprintf(stderr, "  -d, --dedup            Enable UMI-based deduplication\n");
    fprintf(stderr, "  -b, --cbc-location STR Cell barcode tag name or field number (default: CB)\n");
    fprintf(stderr, "  -u, --umi-location STR UMI tag name or field number (default: UB)\n");
    fprintf(stderr, "  -v, --verbose [INT]    Verbosity level: -v (INFO), -v 5 or --verbose=5 (DEBUG)\n");
    fprintf(stderr, "  -h, --help             Show this help message\n");
    fprintf(stderr, "\n");
}

int8_t read_dump(cb2fp *direct_map, char *this_CB, 
                sam_hdr_t *header, bam1_t *read) {
    cb2fp *entry;
    
    HASH_FIND_STR(direct_map, this_CB, entry);
    if (entry == NULL) {
        // Cell barcode not found in metadata, skip
        return 0;
    }
    
    int write_stat = sam_write1(entry->fp, header, read);
    if (write_stat < 0) {
        log_msg("Failed to write read to output file", ERROR);
        return 1;
    }
    
    return 0;
}

void log_message(char* log_path, log_level_t out_level, char* message, log_level_t level, ...) {
    if (level > out_level) return;
    
    va_list args;
    va_start(args, level);
    
    // Print level prefix
    char* level_names[] = {"ERROR", "WARNING", "", "INFO", "", "DEBUG"};
    if (level < 6) {
        fprintf(stderr, "[%s] ", level_names[level]);
    }
    
    // Print formatted message
    vfprintf(stderr, message, args);
    fprintf(stderr, "\n");
    
    va_end(args);
}

int create_directory(char* pathname) {
    struct stat st = {0};
    
    if (stat(pathname, &st) == -1) {
        if (mkdir(pathname, 0755) == -1) {
            log_msg("Failed to create directory: %s", ERROR, pathname);
            return -1;
        }
        log_msg("Created output directory: %s", INFO, pathname);
    } else {
        log_msg("Output directory already exists: %s", WARNING, pathname);
        // Ask user for confirmation or just proceed
        return 0;
    }
    return 0;
}

tag_meta_t *initialize_tag_meta() {
    tag_meta_t *tag_meta = calloc(1, sizeof(tag_meta_t));
    if (!tag_meta) return NULL;
    
    tag_meta->tag_name = calloc(3, sizeof(char));
    tag_meta->sep = calloc(2, sizeof(char));
    
    if (!tag_meta->tag_name || !tag_meta->sep) {
        free(tag_meta->tag_name);
        free(tag_meta->sep);
        free(tag_meta);
        return NULL;
    }
    
    // Set default values
    tag_meta->location = READ_TAG;
    strcpy(tag_meta->tag_name, "CB");
    strcpy(tag_meta->sep, ",");
    tag_meta->field = 1;
    tag_meta->length = 21;
    
    return tag_meta;
}

void destroy_tag_meta(tag_meta_t *tag_meta) {
    if (tag_meta) {
        free(tag_meta->sep);
        free(tag_meta->tag_name);
        free(tag_meta);
    }
}

void set_CB(tag_meta_t *tag_meta, char *platform) {
    if (!platform) return;
    
    // Convert to lowercase for comparison
    char *platform_lower = strdup(platform);
    if (!platform_lower) {
        log_msg("Failed to allocate memory for platform string", ERROR);
        return;
    }
    for (int i = 0; platform_lower[i]; i++) {
        platform_lower[i] = tolower(platform_lower[i]);
    }
    
    if (strcmp("10xv2", platform_lower) == 0) {
        tag_meta->length = 18 + 1;
    } else if (strcmp("scirnaseq3", platform_lower) == 0) {
        tag_meta->location = READ_NAME;
        tag_meta->length = 20 + 1;
        tag_meta->field = 1;
    } else {
        // Default (10Xv3)
        tag_meta->length = 18 + 1;
    }
    
    free(platform_lower);
}

void set_UB(tag_meta_t *tag_meta, char *platform) {
    if (!platform) return;
    
    // Convert to lowercase for comparison
    char *platform_lower = strdup(platform);
    if (!platform_lower) {
        log_msg("Failed to allocate memory for platform string", ERROR);
        return;
    }
    for (int i = 0; platform_lower[i]; i++) {
        platform_lower[i] = tolower(platform_lower[i]);
    }
    
    if (strcmp("10xv2", platform_lower) == 0) {
        strcpy(tag_meta->tag_name, "UB");
        tag_meta->length = 10 + 1;
    } else if (strcmp("scirnaseq3", platform_lower) == 0) {
        tag_meta->location = READ_NAME;
        tag_meta->length = 8 + 1;
        tag_meta->field = 2;
    } else {
        // Default (10Xv3)
        strcpy(tag_meta->tag_name, "UB");
        tag_meta->length = 12 + 1;
    }
    
    free(platform_lower);
}

void print_tag_meta(tag_meta_t *tag_meta, const char *header) {
    char* location_names[] = {"Read tag", "Read name"};
    
    if (header) {
        fprintf(stderr, "\t%s:\n", header);
    } else {
        fprintf(stderr, "Tag Information\n");
    }
    
    fprintf(stderr, "\t\tLocation: %s\n", location_names[tag_meta->location]);
    if (tag_meta->location == READ_NAME) {
        fprintf(stderr, "\t\tSeparator: %s\n", tag_meta->sep);
        fprintf(stderr, "\t\tField number: %d\n", tag_meta->field);
    } else {
        fprintf(stderr, "\t\tTag name: %s\n", tag_meta->tag_name);
    }
    fprintf(stderr, "\t\tTag length: %d\n\n", tag_meta->length - 1);
}