#ifndef BAMSORTER_H
#define BAMSORTER_H

#include <string.h>
#include "htslib/sam.h"

struct BamSortByCoord{
    inline bool operator()(const bam1_t* b1, const bam1_t* b2) const {
        if(b1 == b2) return false; // strick weak
        if(b1->core.tid >= 0){ // b1 is mapped
            if(b2->core.tid < 0) return true; // b2 is unmapped
            else if(b2->core.tid >  b1->core.tid) // both mapped
               return true;
            else if(b2->core.tid == b1->core.tid && b2->core.pos >  b1->core.pos)
                return true;
            else if(b2->core.tid == b1->core.tid && b2->core.pos == b1->core.pos && b2->core.mtid >  b1->core.mtid)
                return true;
            else if(b2->core.tid == b1->core.tid && b2->core.pos == b1->core.pos && b2->core.mtid == b1->core.mtid && b2->core.mpos >  b1->core.mpos)
                return true;
            else if(b2->core.tid == b1->core.tid && b2->core.pos == b1->core.pos && b2->core.mtid == b1->core.mtid && b2->core.mpos == b1->core.mpos){
                if(b2->core.isize > b1->core.isize) return true;
                else if(b2->core.isize == b1->core.isize && b2->data > b1->data) return true;
                else return false;
            }else return false;
        }else{ // b1 is unmapped
            if(b2->core.tid < 0){ // both are unmapped
                return b2->data > b1->data;
            }else return false;
        }
    }
};

struct BamSortByQname{
    inline bool operator()(const bam1_t* b1, const bam1_t* b2) const {
        if(b1 == b2) return false;
        int cmp = strcmp(bam_get_qname(b1), bam_get_qname(b2));
        if(cmp < 0) return true;
        else if(cmp == 0) return b1->core.flag & BAM_FREAD1;
        else return false;
    }
};

#endif
