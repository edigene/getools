#ifndef KSVAR_STAT_H
#define KSVAR_STAT_H

#ifndef KSW_FTYPE
#define KSW_FTYPE uint32_t
#endif

#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif

#include "htslib/kstring.h"
#include "htslib/sam.h"
#include "common.h"
#include <libgen.h>
#include <unistd.h>
#include <vector>

struct stat_opt_t{
     const char* inbam; // input bam
     const char* name; // contig name
     const char* outtsv; // output tsv of ref
     int rlen = 300; // ref len
     int buflen = 600; // buffer len
     int maxl = 0; // max used buffer len
     KSW_FTYPE imask = 0; // VC flag must met
     KSW_FTYPE emask = 0; // VC flag must not met

     stat_opt_t(){
         name = NULL;
         rlen = 300;
         buflen = 600;
         maxl = 0;
         imask = emask = 0;
     }

     ~stat_opt_t(){}

     void init(){
         buflen = rlen * 2;
     }

     void stat();
     void out_head(kstring_t* ks);
     void out_rec1d(kstring_t* ks, const char* name, const std::vector<int64_t>& res);
     void out_rec2d(kstring_t* ks, const char* name, int mh, const std::vector<std::vector<int64_t>>& res);
};

void stats_usage(stat_opt_t* opt, char* arg0);
int stats_main(int argc, char** argv);

#endif
