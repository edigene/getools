#ifndef VARCNT_H
#define VARCNT_H

#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <memory.h>
#include "htslib/kstring.h"

// variant quantitive number stat(count, length...)
struct var_num_t{
    uint64_t* vcdist_array = NULL; // variant count dist array
    int64_t vcdist_aloc = 0; // variant count dist allocated
    int64_t vcmax = 0; // variant count max item
    int64_t vcmin = 0; // variant count min item
    uint64_t vcmod = 0; // variant count mod item
    uint64_t vcmedian = 0; // variant count median item
    uint64_t vcquant[4] = {0}; // variant quantile

    var_num_t(int buflen){
        vcdist_aloc = buflen;
        vcdist_array = (uint64_t*)calloc(vcdist_aloc, sizeof(uint64_t));
    }

    ~var_num_t(){
        if(vcdist_array){
            free(vcdist_array);
            vcdist_array = NULL;
        }
    }

    // add/sub variant stat item c with count n
    void statvar(int c, int n, bool add){
        if(add){
            if(c >= vcdist_aloc){
                int ol = vcdist_aloc;
                vcdist_aloc = c + 1;
                vcdist_array = (uint64_t*)realloc(vcdist_array, sizeof(uint64_t)*vcdist_aloc);
                memset(vcdist_array + ol, 0, (vcdist_aloc-ol)*sizeof(uint64_t));
            }
            vcdist_array[c] += n;
        }else{
            assert(c < vcdist_aloc);
            vcdist_array[c] -= n;
        }
    }

    void summary(){
        uint64_t totitem = 0;
        vcmin = vcdist_aloc + 1;
        while(vcdist_aloc > 1 && vcdist_array[vcdist_aloc-1] == 0) --vcdist_aloc;
        for(int i = 0; i < vcdist_aloc; ++i) totitem += vcdist_array[i];
        uint64_t acc_sum = 0;
        uint64_t old_sum = acc_sum;
        uint64_t quat_arr[4] = {0};
        for(int i =  1; i <= 4; ++i) quat_arr[i-1] = totitem * (i * 0.25);
        for(int64_t i = 0; i < vcdist_aloc; ++i){
            if(vcdist_array[i]){
                if(i < vcmin) vcmin = i;
                if(i > vcmax) vcmax = i;
                if(vcdist_array[i] > vcdist_array[vcmod]) vcmod = i;
                old_sum = acc_sum;
                acc_sum += vcdist_array[i];
                for(int q = 0; q < 4; ++q){
                    if(quat_arr[q] >= old_sum && quat_arr[q] <= acc_sum){
                        vcquant[q] = i;
                    }
                }
            }
        }
    }

    void reportJSON(kstring_t* s, const char* dh, const char* dm, const char* name){
        ksprintf(s, "%s\"%sSummary\": {\n", dh, name);
        ksprintf(s, "%s%s\"%sMax\": %lld,\n", dh, dm, name, vcmax);
        ksprintf(s, "%s%s\"%sMin\": %lld,\n", dh, dm, name, vcmin);
        ksprintf(s, "%s%s\"%sMod\": %lld,\n", dh, dm, name, vcmod);
        ksprintf(s, "%s%s\"%sQuartile\": [", dh, dm, name);
        for(int i = 0; i < 4; ++i){
            ksprintf(s, "%lld,", vcquant[i]);
        }
        s->s[s->l-1] = ']';
        ksprintf(s, ",\n");
        ksprintf(s, "%s%s\"%sDist\": [", dh, dm, name);
        for(int64_t i = 0; i < vcdist_aloc; ++i){
            ksprintf(s, "%lld,", vcdist_array[i]);
        }
        if(s->s[s->l-1] == ',') s->s[s->l-1] = ']';
        else kputc(']', s);
        ksprintf(s, "\n%s},\n", dh);
    }
};

#endif
