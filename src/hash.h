//
// Created by Yen-Chung Chen on 2/23/23.
//
// Hash table definitions for read-to-label and label-to-file mappings
// Dependencies: htslib, uthash (external), shared_const.h

#ifndef SCBAMSPLIT_HASH_H
#define SCBAMSPLIT_HASH_H

// External library includes
#include "htslib/sam.h"
#include "uthash.h"

// Project includes
#include "shared_const.h"

// Direct mapping: cell_barcode -> file_pointer
typedef struct {
    char cb[32];                          /* key: cell barcode (sufficient for all platforms) */
    char label[64];                       /* cluster label (reasonable for most labels) */
    samFile* fp;                          /* direct file pointer */
    UT_hash_handle hh;                    /* makes this structure hashable */
} cb2fp;

cb2fp* hash_readtag_direct(char *path, const char *prefix, sam_hdr_t *header);


#endif //SCBAMSPLIT_HASH_H
