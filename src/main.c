#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <stdbool.h>
#include <unistd.h>
#include <errno.h>
#include <limits.h>
#include <ctype.h>

#ifndef PATH_MAX
#define PATH_MAX 4096
#endif
#include <htslib/sam.h>
#include "uthash.h"
#include "hash.h"
#include "utils.h"
#include "sort.h"
#include "dedup_3pass.h"

// Global variables
char *OUT_PATH = "";
char *LEVEL_FLAG[] = {"ERROR", "WARNING", "", "INFO", "", "DEBUG"};
log_level_t OUT_LEVEL = WARNING;
int64_t CB_LENGTH = 21;
int64_t UB_LENGTH = 21;


// Command function for split subcommand
int cmd_split(int argc, char *argv[]) {
    int32_t opt;
    int64_t mapq_thres = 0;
    int64_t out_level_raw = 0;
    bool dedup = false, dryrun = false, verbose = false;
    char *bampath = NULL;
    char *metapath = NULL;
    char *oprefix = NULL;
    tag_meta_t *cb_meta = initialize_tag_meta();
    tag_meta_t *ub_meta = initialize_tag_meta();
    
    if (!cb_meta || !ub_meta) {
        log_msg("Failed to allocate tag metadata structures", ERROR);
        destroy_tag_meta(cb_meta);
        destroy_tag_meta(ub_meta);
        return 1;
    }
    
    strcpy(ub_meta->tag_name, "UB");
    int64_t cb_field = 0;
    int64_t ub_field = 0;
    int32_t return_val = 0;

    // Command line options
    static struct option cl_opts[] = {
        {"file", required_argument, NULL, 'f'},
        {"meta", required_argument, NULL, 'm'},
        {"output", required_argument, NULL, 'o'},
        {"mapq", required_argument, NULL, 'q'},
        {"platform", required_argument, NULL, 'p'},
        {"dedup", no_argument, NULL, 'd'},
        {"cbc-location", required_argument, NULL, 'b'},
        {"cbc-length", required_argument, NULL, 'L'},
        {"umi-location", required_argument, NULL, 'u'},
        {"umi-length", required_argument, NULL, 'l'},
        {"dry-run", no_argument, NULL, 'n'},
        {"verbose", optional_argument, NULL, 'v'},
        {"help", no_argument, NULL, 'h'}
    };

    while ((opt = getopt_long(argc, argv, ":f:m:o:q:p:db:L:u:l:nv::h", cl_opts, NULL)) != -1) {
        switch (opt) {
            case 'f':
                bampath = optarg;
                break;
            case 'm':
                metapath = optarg;
                break;
            case 'o':
                oprefix = optarg;
                break;
            case 'q':
                {
                    char *endptr;
                    errno = 0;
                    mapq_thres = strtol(optarg, &endptr, 10);
                    if (errno == ERANGE || *endptr != '\0' || mapq_thres < 0) {
                        log_msg("Invalid MAPQ threshold: %s", ERROR, optarg);
                        goto error_out_and_free;
                    }
                }
                break;
            case 'p':
                set_CB(cb_meta, optarg);
                set_UB(ub_meta, optarg);
                CB_LENGTH = cb_meta->length;
                UB_LENGTH = ub_meta->length;
                if (CB_LENGTH >= 32 || UB_LENGTH >= 32) {
                    log_msg("Platform barcode/UMI lengths exceed maximum supported size (32)", ERROR);
                    goto error_out_and_free;
                }
                break;
            case 'd':
                dedup = true;
                break;
            case 'b':
                {
                    char *endptr;
                    errno = 0;
                    cb_field = strtol(optarg, &endptr, 10);
                    if (errno == ERANGE) {
                        log_msg("Invalid cell barcode field (overflow): %s", ERROR, optarg);
                        goto error_out_and_free;
                    } else if (*endptr != '\0' && cb_field == 0) {
                        // Not a number, treat as tag name
                        cb_field = 0;
                    } else if (*endptr != '\0' || cb_field < 0) {
                        log_msg("Invalid cell barcode field: %s", ERROR, optarg);
                        goto error_out_and_free;
                    }
                }
                if (cb_field == 0) {
                    // Validate tag name length
                    if (strlen(optarg) > 2) {
                        log_msg("Tag name too long (max 2 chars): %s", ERROR, optarg);
                        goto error_out_and_free;
                    }
                    strncpy(cb_meta->tag_name, optarg, 2);
                    cb_meta->tag_name[2] = '\0';
                } else {
                    cb_meta->location = READ_NAME;
                    cb_meta->field = cb_field;
                }
                break;
            case 'L':
                {
                    char *endptr;
                    errno = 0;
                    long tmp = strtol(optarg, &endptr, 10);
                    if (errno == ERANGE || *endptr != '\0' || tmp < 0) {
                        log_msg("Invalid cell barcode length: %s", ERROR, optarg);
                        goto error_out_and_free;
                    }
                    CB_LENGTH = tmp + 1;
                }
                if (CB_LENGTH > 0 && CB_LENGTH < 32) {
                    cb_meta->length = CB_LENGTH;
                } else if (CB_LENGTH >= 32) {
                    log_msg("Cell barcode length must be less than 32", ERROR);
                    goto error_out_and_free;
                } else {
                    log_msg("Cell barcode length must be larger than 0", ERROR);
                    goto error_out_and_free;
                }
                break;
            case 'u':
                {
                    char *endptr;
                    errno = 0;
                    ub_field = strtol(optarg, &endptr, 10);
                    if (errno == ERANGE) {
                        log_msg("Invalid UMI field (overflow): %s", ERROR, optarg);
                        goto error_out_and_free;
                    } else if (*endptr != '\0' && ub_field == 0) {
                        // Not a number, treat as tag name
                        ub_field = 0;
                    } else if (*endptr != '\0' || ub_field < 0) {
                        log_msg("Invalid UMI field: %s", ERROR, optarg);
                        goto error_out_and_free;
                    }
                }
                if (ub_field == 0) {
                    // Validate tag name length
                    if (strlen(optarg) > 2) {
                        log_msg("Tag name too long (max 2 chars): %s", ERROR, optarg);
                        goto error_out_and_free;
                    }
                    strncpy(ub_meta->tag_name, optarg, 2);
                    ub_meta->tag_name[2] = '\0';
                } else {
                    ub_meta->location = READ_NAME;
                    ub_meta->field = ub_field;
                }
                break;
            case 'l':
                {
                    char *endptr;
                    errno = 0;
                    long tmp = strtol(optarg, &endptr, 10);
                    if (errno == ERANGE || *endptr != '\0' || tmp < 0) {
                        log_msg("Invalid UMI length: %s", ERROR, optarg);
                        goto error_out_and_free;
                    }
                    UB_LENGTH = tmp + 1;
                }
                if (UB_LENGTH > 0 && UB_LENGTH < 32) {
                    ub_meta->length = UB_LENGTH;
                } else if (UB_LENGTH >= 32) {
                    log_msg("UMI length must be less than 32", ERROR);
                    goto error_out_and_free;
                } else {
                    log_msg("UMI length must be larger than 0", ERROR);
                    goto error_out_and_free;
                }
                break;
            case 'n':
                dryrun = true;
                break;
            case 'v':
                {
                    char *level_arg = NULL;
                    
                    if (optarg) {
                        // --verbose=5 case
                        level_arg = optarg;
                    } else if (optind < argc && argv[optind][0] != '-') {
                        // Check if next argument looks like a verbosity level (digit)
                        char *next_arg = argv[optind];
                        if (strlen(next_arg) == 1 && isdigit(next_arg[0])) {
                            // -v 5 case - consume the next argument
                            level_arg = next_arg;
                            optind++;
                        }
                        // If next arg doesn't look like a level, treat as -v with no argument
                    }
                    
                    if (level_arg) {
                        // Parse the verbosity level
                        char *endptr;
                        errno = 0;
                        out_level_raw = strtol(level_arg, &endptr, 10);
                        if (errno == ERANGE || *endptr != '\0') {
                            log_msg("Invalid verbosity level: %s", ERROR, level_arg);
                            goto error_out_and_free;
                        }
                        if (out_level_raw == 0) {
                            OUT_LEVEL = INFO;
                        } else if (out_level_raw > 0 && out_level_raw < 5) {
                            OUT_LEVEL = out_level_raw;
                        } else {
                            OUT_LEVEL = DEBUG;
                        }
                    } else {
                        // -v with no argument - default to INFO level
                        OUT_LEVEL = INFO;
                    }
                    verbose = true;
                }
                break;
            case 'h':
                show_split_usage();
                return 0;
            case ':':
                log_msg("Option -%c requires an argument", ERROR, optopt);
                goto error_out_and_free;
            default:
                log_msg("Unknown option", ERROR);
                error_out_and_free:
                    goto cleanup;
        }
    }

    // Check required arguments
    if (bampath == NULL || metapath == NULL) {
        log_msg("Error: Missing required arguments (-f and -m)", ERROR);
        show_split_usage();
        return_val = 1;
        goto cleanup;
    }

    // Set default output prefix
    if (oprefix == NULL) {
        oprefix = "./";
    }

    // Add trailing slash if needed (using stack allocation)
    char oprefix_buffer[PATH_MAX];
    size_t oplen = strlen(oprefix);
    if (oplen > 0 && oprefix[oplen - 1] != '/') {
        if (oplen + 2 >= sizeof(oprefix_buffer)) {
            log_msg("Output prefix path too long (max %zu chars)", ERROR, sizeof(oprefix_buffer) - 2);
            goto cleanup;
        }
        strcpy(oprefix_buffer, oprefix);
        strcat(oprefix_buffer, "/");
        oprefix = oprefix_buffer;
    }

    if (verbose || dryrun) {
        fprintf(stderr, "- Run configuration:\n");
        fprintf(stderr, "\tInput BAM: %s\n", bampath);
        fprintf(stderr, "\tMetadata: %s\n", metapath);
        fprintf(stderr, "\tMAPQ threshold: %lld\n", (long long)mapq_thres);
        fprintf(stderr, "\tOutput prefix: %s\n", oprefix);
        print_tag_meta(cb_meta, "Cell barcode");
        print_tag_meta(ub_meta, "UMI");
        fprintf(stderr, "\tDeduplication: %s\n\n", dedup ? "enabled" : "disabled");
    }

    if (access(bampath, F_OK) != 0) {
        log_msg("Input BAM file not found: %s", ERROR, bampath);
        goto cleanup;
    }

    if (access(metapath, F_OK) != 0) {
        log_msg("Metadata file not found: %s", ERROR, metapath);
        goto cleanup;
    }

    if (access(metapath, R_OK) != 0) {
        log_msg("Metadata file not readable: %s", ERROR, metapath);
        goto cleanup;
    }

    if (dryrun) {
        fprintf(stderr, "Dry run completed successfully.\n");
        goto cleanup;
    }

    // Check output directory write permissions (using stack allocation)
    char output_parent[PATH_MAX];
    size_t oprefix_len = strlen(oprefix);
    if (oprefix_len >= sizeof(output_parent)) {
        log_msg("Output prefix path too long for validation", ERROR);
        goto cleanup;
    }
    strcpy(output_parent, oprefix);
    
    {
        // Remove trailing slash for parent directory check
        size_t len = strlen(output_parent);
        if (len > 1 && output_parent[len-1] == '/') {
            output_parent[len-1] = '\0';
        }
        
        // Check if output directory exists, or parent directory is writable
        if (access(output_parent, F_OK) == 0) {
            if (access(output_parent, W_OK) != 0) {
                log_msg("Output directory not writable: %s", ERROR, output_parent);
                goto cleanup;
            }
        } else {
            // Check parent directory for creation permission (using stack allocation)
            char parent_dir[PATH_MAX];
            if (strlen(output_parent) >= sizeof(parent_dir)) {
                log_msg("Parent directory path too long for validation", ERROR);
                goto cleanup;
            }
            strcpy(parent_dir, output_parent);
            
            {
                char *last_slash = strrchr(parent_dir, '/');
                if (last_slash) {
                    *last_slash = '\0';
                    if (access(parent_dir, W_OK) != 0) {
                        log_msg("Parent directory not writable for output creation: %s", ERROR, parent_dir);
                        goto cleanup;
                    }
                }
            }
        }
    }

    // Create output directory
    if (create_directory(oprefix) != 0) {
        goto cleanup;
    }

    // Open BAM file
    samFile *fp = sam_open(bampath, "r");
    if (!fp) {
        log_msg("Failed to open BAM file: %s", ERROR, bampath);
        goto cleanup;
    }

    // Read header
    sam_hdr_t *header = sam_hdr_read(fp);
    if (!header) {
        log_msg("Failed to read BAM header", ERROR);
        sam_close(fp);
        goto cleanup;
    }

    if (dedup) {
        sam_hdr_change_HD(header, "SO", "scbamsplit");
    }

    // Load metadata and create direct mapping
    cb2fp *direct_map = hash_readtag_direct(metapath, oprefix, header);
    if (!direct_map) {
        log_msg("Failed to load metadata and create output files from: %s", ERROR, metapath);
        sam_hdr_destroy(header);
        sam_close(fp);
        goto cleanup;
    }

    // Process reads
    if (!dedup) {
        // Simple splitting without deduplication
        bam1_t *read = bam_init1();
        char this_CB[CB_LENGTH];
        char this_UB[UB_LENGTH];
        int read_stat;

        while ((read_stat = sam_read1(fp, header, read)) >= 0) {
            int8_t cb_stat = get_CB(read, cb_meta, this_CB);
            int8_t ub_stat = get_UB(read, ub_meta, this_UB);
            int16_t mapq = read->core.qual;

            if (cb_stat != 0 || ub_stat != 0 || mapq < mapq_thres) {
                continue;
            }

            if (read_dump(direct_map, this_CB, header, read) != 0) {
                log_msg("Failed to write read", ERROR);
                break;
            }
        }

        bam_destroy1(read);
    } else {
        // Use 3-pass deduplication
        log_msg("Using 3-pass deduplication algorithm", INFO);
        
        int dedup_result = dedup_3pass(bampath, header, direct_map, cb_meta, ub_meta, mapq_thres);
        if (dedup_result != 0) {
            log_msg("3-pass deduplication failed", ERROR);
            return_val = 1;
        }
    }

    // Cleanup
    sam_close(fp);
    sam_hdr_destroy(header);

    // Close output files and free direct mapping hash table
    // First, collect unique file pointers to avoid double-closing
    samFile **unique_fps = NULL;
    int n_fps = 0;
    int fps_capacity = 10;
    unique_fps = malloc(fps_capacity * sizeof(samFile*));
    
    cb2fp *entry, *tmp;
    HASH_ITER(hh, direct_map, entry, tmp) {
        if (entry->fp) {
            // Check if we've already seen this file pointer
            int found = 0;
            for (int i = 0; i < n_fps; i++) {
                if (unique_fps[i] == entry->fp) {
                    found = 1;
                    break;
                }
            }
            
            // If not found, add to unique list
            if (!found) {
                if (n_fps >= fps_capacity) {
                    fps_capacity *= 2;
                    unique_fps = realloc(unique_fps, fps_capacity * sizeof(samFile*));
                }
                unique_fps[n_fps++] = entry->fp;
            }
        }
        HASH_DEL(direct_map, entry);
        free(entry);
    }
    
    // Now close each unique file pointer once
    for (int i = 0; i < n_fps; i++) {
        sam_close(unique_fps[i]);
    }
    free(unique_fps);

cleanup:
    destroy_tag_meta(cb_meta);
    destroy_tag_meta(ub_meta);
    return return_val;
}

// Main function with subcommand parsing
int main(int argc, char *argv[]) {
    if (argc < 2) {
        show_global_usage();
        return 1;
    }
    
    char *subcommand = argv[1];
    
    // Handle global help
    if (strcmp(subcommand, "--help") == 0 || strcmp(subcommand, "-h") == 0) {
        show_global_usage();
        return 0;
    }
    
    // Handle subcommands
    if (strcmp(subcommand, "split") == 0) {
        // Remove subcommand from argv and pass to cmd_split
        return cmd_split(argc - 1, &argv[1]);
    } else {
        fprintf(stderr, "Error: Unknown command '%s'\n", subcommand);
        fprintf(stderr, "\n");
        show_global_usage();
        return 1;
    }
}