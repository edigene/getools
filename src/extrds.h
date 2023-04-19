#ifndef EXT_READS_H
#define EXT_READS_H

#include "htslib/sam.h"
#include <algorithm>
#include <assert.h>
#include "bamsorter.h"
#include "common.h"
#include "ksw4get.h"
#include <vector>

#ifndef KSW_FTYPE
#define KSW_FTYPE uint32_t
#endif

struct ext_opt_t{
    const char* inbam = NULL;
    const char* outbam = NULL;
    const char* contig = NULL;
    const char* cellbc = NULL;
    KSW_FTYPE imask = 0;
    KSW_FTYPE emask = 0;
    int mergs = -1;
    bool sbyc = false;
    

    ext_opt_t(){}

    ~ext_opt_t(){}

    bool valid(){
        bool valid = true;
        if(inbam == NULL){
            valid = false;
            fprintf(stderr, "input/output must be provided\n");
        }
        return valid;
    }

    void ext();
};

extern int extrds_main(int argc, char** argv);

#endif
