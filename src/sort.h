//
// BAM file processing functions (Clean implementation)
//

#ifndef SCBAMSPLIT_SORT_H
#define SCBAMSPLIT_SORT_H

// External library includes
#include "htslib/sam.h"

// Project includes
#include "utils.h"

// Essential BAM tag extraction functions
int8_t get_CB(bam1_t *read, tag_meta_t* info, char* tag_ptr);
int8_t get_UB(bam1_t *read, tag_meta_t* info, char* tag_ptr);

#endif //SCBAMSPLIT_SORT_H