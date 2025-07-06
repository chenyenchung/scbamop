//
// Utility functions and data structures (Clean implementation)
//

#ifndef SCBAMSPLIT_UTILS_H
#define SCBAMSPLIT_UTILS_H

// Standard library includes
#include <stdbool.h>
#include <stdint.h>

// External library includes
#include "htslib/sam.h"

// Project includes
#include "hash.h"

// Tag metadata structure
enum location {
    READ_TAG,
    READ_NAME
};

typedef struct {
    enum location location;
    char *tag_name;
    char *sep;
    uint8_t field;
    uint8_t length;
} tag_meta_t;

// Logging utilities
typedef enum {
    DEBUG = 5,
    ERROR = 0,
    WARNING = 1,
    INFO = 3
} log_level_t;

void log_message(char* message, log_level_t level, char* log_path, log_level_t OUT_LEVEL, ...);
#define log_msg(...) log_message(OUT_PATH, OUT_LEVEL, __VA_ARGS__)

// Essential functions
void show_global_usage();
void show_split_usage();
int create_directory(char* pathname);
tag_meta_t *initialize_tag_meta();
void destroy_tag_meta(tag_meta_t *tag_meta);
void set_CB(tag_meta_t *tag_meta, char *platform);
void set_UB(tag_meta_t *tag_meta, char *platform);
void print_tag_meta(tag_meta_t *tag_meta, const char *header);
int8_t read_dump(cb2fp *direct_map, char *this_CB, 
                sam_hdr_t *header, bam1_t *read);

// Global variables
extern log_level_t OUT_LEVEL;
extern char *OUT_PATH;
extern char* LEVEL_FLAG[6];
extern int64_t CB_LENGTH;
extern int64_t UB_LENGTH;

#endif //SCBAMSPLIT_UTILS_H