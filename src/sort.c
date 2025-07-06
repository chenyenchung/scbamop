//
// BAM file processing functions (Clean implementation)
//

#include "sort.h"
#include <string.h>
#include <stdlib.h>

// Helper function to fetch tag from BAM record
static int8_t fetch_tag2(bam1_t *read, char* tag_ptr, tag_meta_t *info) {
    char* tag_content = (char *) bam_aux_get(read, info->tag_name);
    if (NULL == tag_content) {
        // Tag not found
        return -1;
    }
    strncpy(tag_ptr, tag_content + 1, info->length - 1);
    tag_ptr[info->length - 1] = '\0';
    return 0;
}

// Helper function to fetch tag from read name
static int8_t fetch_name(bam1_t *read, char* tag_ptr, tag_meta_t *info) {
    char* rn = bam_get_qname(read);

    if (rn[0] == info->sep[0]) return -1;

    // Use static buffer to avoid malloc/free for every read
    // Note: Static variables are NOT thread-safe, but we'll keep simple for now
    static char rncpy[512];
    size_t rn_len = strlen(rn);
    
    // Check if read name fits in our buffer
    if (rn_len >= sizeof(rncpy)) {
        return -1; // Read name too long
    }
    
    strncpy(rncpy, rn, sizeof(rncpy) - 1);
    rncpy[sizeof(rncpy) - 1] = '\0';
    char *token;

    uint8_t field_num = 0;
    int8_t return_val = 1;
    token = strtok(rncpy, info->sep);
    while (token != NULL) {
        field_num++;
        if (field_num == info->field) {
            strncpy(tag_ptr, token, info->length - 1);
            tag_ptr[info->length - 1] = '\0';
            return_val = 0;
            break;
        }
        token = strtok(NULL, info->sep);
    }
    return return_val;
}

int8_t get_CB(bam1_t *read, tag_meta_t* info, char* tag_ptr) {
    int8_t exit_code = 0;
    switch (info->location) {
        case READ_TAG:
            exit_code = fetch_tag2(read, tag_ptr, info);
            break;
        case READ_NAME:
            exit_code = fetch_name(read, tag_ptr, info);
            break;
        default:
            log_msg("Unknown location type to fetch cell barcode", ERROR);
            exit_code = 1;
    }
    return exit_code;
}

int8_t get_UB(bam1_t *read, tag_meta_t* info, char* tag_ptr) {
    int8_t exit_code = 0;
    switch (info->location) {
        case READ_TAG:
            exit_code = fetch_tag2(read, tag_ptr, info);
            break;
        case READ_NAME:
            exit_code = fetch_name(read, tag_ptr, info);
            break;
        default:
            log_msg("Unknown location type to fetch UMI", ERROR);
            exit_code = 1;
    }
    return exit_code;
}