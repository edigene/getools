#include "report.h"

void Reporter::reportJSON(SplitResult* sr, QcStat* qs1, QcStat* qs2, PairMerger* merge){
    kstring_t *ks = (kstring_t*)calloc(1, sizeof(kstring_t));
    const char* d1 = "    ";
    const char* d2 = "        ";
    // json beg
    ksprintf(ks, "{\n");
    if(qs1){
        // beg BasicQC
        ksprintf(ks, "%s\"BasicQC\": {\n", d1);
        qs1->reportJSON(ks, d2, d1, 1);
        if(qs2){
            ksprintf(ks, ",\n");
            qs2->reportJSON(ks, d2, d1, 2);
        }
        ksprintf(ks, "\n%s},\n", d1);
    }
    // end BasicQC
    // beg SplitResult
    if(sr){
        ksprintf(ks, "%s\"SplitResult\": {\n", d1);
        sr->reportJSON(ks, d2, d1);
        ksprintf(ks, "\n%s},\n", d1);
    }
    // end SplitResult
    // beg MergeStat
    if(merge){
        ksprintf(ks, "%s\"MergeResult\": {\n", d1);
        merge->reportJSON(ks, d2, d1);
        ksprintf(ks, "\n%s},\n", d1);
    }
    // end MergeStat
    // beg Software
    ksprintf(ks, "%s\"Software\": {\n", d1);
    mOpt->soft->reportJSON(ks, d1, d1);
    ksprintf(ks, "\n%s}", d1);
    // end Software
    // json end
    ksprintf(ks, "\n}\n"); 
    // write to file
    FILE* fout = fopen(mOpt->jsr.c_str(), "w");
    fprintf(fout, "%s", ks->s);
    fclose(fout);
    if(ks->l) free(ks->s);
}

void Reporter::reportTSV(SplitResult* sr, QcStat* qs1, QcStat* qs2, PairMerger* merge){
    kstring_t* ks = (kstring_t*)calloc(1, sizeof(kstring_t));
    kstring_t* ms = (kstring_t*)calloc(1, sizeof(kstring_t));
    // head
    if(qs1){ 
        qs1->tsvHead(ks); kputc('\t', ks);
        kputsn(ks->s, ks->l, ms);
    }
    if(sr){ 
        sr->tsvHead(ks, false); kputc('\t', ks);
        sr->tsvHead(ms, true); kputc('\t', ms);
    }
    if(merge){ 
        merge->tsvHead(ks); kputc('\t', ks); 
        merge->tsvHead(ms); kputc('\t', ms);
    }
    ks->s[ks->l-1] = '\n';
    ms->s[ms->l-1] = '\n';
    // body
    if(qs1){
        qs1->tsvBody(ks, qs2); kputc('\t', ks);
        qs1->tsvBody(ms, qs2); kputc('\t', ms);
    }
    if(sr){
        sr->tsvBody(ks, false); kputc('\t', ks);
        sr->tsvBody(ms, true); kputc('\t', ms);
    }
    if(merge){
        merge->tsvBody(ks); kputc('\t', ks); 
        merge->tsvBody(ms); kputc('\t', ms); 
    }
    ks->s[ks->l-1] = '\n';
    ms->s[ms->l-1] = '\n';
    // out 1
    FILE* fout = fopen(mOpt->tsr.c_str(), "w");
    fprintf(fout, "%s", ks->s);
    fclose(fout);
    // out 2
    fout = fopen(mOpt->tsm.c_str(), "w");
    fprintf(fout, "%s", ms->s);
    fclose(fout);
    // release mem
    if(ks->l) free(ks->s);
    if(ms->l) free(ms->s);
    free(ks); ks = NULL;
    free(ms); ms = NULL;
}

void Reporter::report(SplitResult* sr, QcStat* qs1, QcStat* qs2, PairMerger* merge){
    reportJSON(sr, qs1, qs2, merge);
    reportTSV(sr, qs1, qs2, merge);
}
