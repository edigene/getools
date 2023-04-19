#include "gedetect.h"
#include "paligner.h"

bam_hdr_t* PAligner::get_bam_hdr(){
    kstring_t* ks = (kstring_t*)calloc(1, sizeof(kstring_t));
    kputs("@SQ\tSN:", ks);
    kputs(name, ks);
    ksprintf(ks, "\tLN:%d\n", rlen);
    bam_hdr_t* h = sam_hdr_parse(ks->l, ks->s);
    h->l_text = ks->l;
    h->text = ks->s;
    free(ks);
    return h;
}

kswr_t* PAligner::align(char* qseq, int qlen){
    // init
    uint8_t* qints = kswge_seq2ints(qseq, qlen);
    // align, first round
    kswr_t *ret = kswge_semi_global(qlen, qints, rlen, refints, 5, score_mat, gapopen, gapext);
    if(ret){
        if(!ret->cigar){
            ret->smask = 0;
            free(qints);
            return ret;
        }
        kswge_mark_mismatch(ret, refints, rlen, qints, qlen);
    }else{
        ret = (kswr_t*)calloc(1, sizeof(kswr_t));
        ret->smask = 0;
        ret->cigar = NULL;
    }
    free(qints);
    return ret;
}

bam1_t* PAligner::ksw2bam(kswr_t* ret, char* qseq, int qlen, char* qname, int qnlen, int mask){
    bam1_t* b = bam_init1();
    b->core.tid = -1;
    b->core.pos = -1;
    b->core.qual = 255;
    b->core.flag = mask;
    b->core.n_cigar = 0;
    b->core.mtid = -1;
    b->core.isize = 0;
    int score = 0;
    if(ret->cigar){// aligned
        b->core.tid = rtid;
        b->core.pos = ret->tb;
        b->core.qual = 60*((double)ret->score/(double)((qlen - ret->lsc - ret->tsc) * match));
        b->core.n_cigar = ret->ncigar;
        score = ret->score;
    }else{
        b->core.flag |= BAM_FUNMAP;
        ret->nmm = qlen;
    }
    // allocate bam record memory
    b->core.l_qname = qnlen + 1;
    b->core.l_qseq = qlen;
    b->l_data = b->core.l_qname + (b->core.n_cigar << 2) + ((b->core.l_qseq + 1) >> 1) + (b->core.l_qseq);
    b->data = (uint8_t*)malloc(b->l_data*sizeof(uint8_t));
    memcpy(b->data, qname, qnlen + 1);
    if(ret->cigar) memcpy(b->data + b->core.l_qname, (uint8_t*)ret->cigar, ret->ncigar << 2);
    uint8_t* mbases = b->data + b->core.l_qname + (b->core.n_cigar << 2);
    for(int i = 0; i < qlen; ++i){
        uint8_t base = seq_nt16_table[(int)qseq[i]];
        mbases[i>>1] &= ~(0xF << ((~i & 1) << 2));
        mbases[i>>1] |= base << ((~i & 1) << 2);
    }
    uint8_t* quals = bam_get_qual(b);
    quals[0] = 0xff;
    // add NM, AS
    bam_aux_update_int(b, "NM", ret->nmm); // mismatch
    bam_aux_update_int(b, "AS", score); // alignment score
    return b;
}

void PAligner::add_var(bam1_t* b){
    bams.push_back(b);
}

void PAligner::stat_var(kswr_t* r1, kswr_t* r2){
    if(r1){
        if(r1->nmm >= rlen) ++nmmdist_se[rlen-1];
        else ++nmmdist_se[r1->nmm];
        if(r1->nins >= rlen) ++ninsdist_se[rlen-1];
        else ++ninsdist_se[r1->nins];
        if(r1->ndel >= rlen) ++ndeldist_se[rlen-1];
        else ++ndeldist_se[r1->ndel];
        if(r1->nvar >= rlen) ++nvardist_se[rlen-1];
        else ++nvardist_se[r1->nvar];
        if(r1->mm >= rlen) ++mmdist_se[rlen-1];
        else ++mmdist_se[r1->mm];
    }
    if(r2){
        if(r2->nmm >= rlen) ++nmmdist_se[rlen-1];
        else ++nmmdist_se[r2->nmm];
        if(r2->nins >= rlen) ++ninsdist_se[rlen-1];
        else ++ninsdist_se[r2->nins];
        if(r2->ndel >= rlen) ++ndeldist_se[rlen-1];
        else ++ndeldist_se[r2->ndel];
        if(r2->nvar >= rlen) ++nvardist_se[rlen-1];
        else ++nvardist_se[r2->nvar];
        if(r2->mm >= rlen) ++mmdist_se[rlen-1];
        else ++mmdist_se[r2->mm];
    }
    if(r1 && r2){
        if(r1->nmm+r2->nmm >= rlen) ++nmmdist_pe[rlen-1];
        else ++nmmdist_pe[r1->nmm+r2->nmm];
        if(r1->nins+r2->nins >= rlen) ++ninsdist_pe[rlen-1];
        else ++ninsdist_pe[r1->nins+r2->nins];
        if(r1->ndel+r2->ndel >= rlen) ++ndeldist_pe[rlen-1];
        else ++ndeldist_pe[r1->ndel+r2->ndel];
        if(r1->nvar+r2->nvar >= rlen) ++nvardist_pe[rlen-1];
        else ++nvardist_pe[r1->nvar+r2->nvar];
        if(r1->mm+r2->mm >= rlen) ++mmdist_pe[rlen-1];
        else ++mmdist_pe[r1->mm+r2->mm];
    }
}

void PAligner::reportJSON(kstring_t* s, const char* dh, const char* dm){
    summary();
    // summary
    ksprintf(s, "%s\"Summary\": {\n", dh);
    ksprintf(s, "%s%s\"TotalReads\": %lld,\n", dh, dm, mapped+unmapped);
    ksprintf(s, "%s%s\"MappedReads\": %lld,\n", dh, dm, mapped);
    ksprintf(s, "%s%s\"UnmappedReads\": %lld\n", dh, dm, unmapped);
    ksprintf(s, "%s},\n", dh);
    // mutdist
    ksprintf(s, "%s\"MutateDistSE\": {\n", dh);
    // nmmdist_se
    ksprintf(s, "%s%s\"ReadsMismatchCountDist\": [", dh, dm);
    for(int i = 0; i < maxdstl_se; ++i) ksprintf(s, "%lld,", nmmdist_se[i]);
    s->s[s->l-1] = ']';
    ksprintf(s, ",\n");
    // ninsdist_se
    ksprintf(s, "%s%s\"ReadsInsCountDist\": [", dh, dm);
    for(int i = 0; i < maxdstl_se; ++i) ksprintf(s, "%lld,", ninsdist_se[i]);
    s->s[s->l-1] = ']';
    ksprintf(s, ",\n");
    // ndeldist_se
    ksprintf(s, "%s%s\"ReadsDelCountDist\": [", dh, dm);
    for(int i = 0; i < maxdstl_se; ++i) ksprintf(s, "%lld,", ndeldist_se[i]);
    s->s[s->l-1] = ']';
    ksprintf(s, ",\n");
    // nvardist_se
    ksprintf(s, "%s%s\"ReadsVarCountDist\": [", dh, dm);
    for(int i = 0; i < maxdstl_se; ++i) ksprintf(s, "%lld,", nvardist_se[i]);
    s->s[s->l-1] = ']';
    ksprintf(s, ",\n");
    // mmdist_se
    ksprintf(s, "%s%s\"ReadsTotalDiffLengthDist\": [", dh, dm);
    for(int i = 0; i < maxdstl_se; ++i) ksprintf(s, "%lld,", mmdist_se[i]);
    s->s[s->l-1] = ']';
    ksprintf(s, "%s},\n", dh);

    ksprintf(s, "%s\"MutateDistPE\": {\n", dh);
    // nmmdist_pe
    ksprintf(s, "%s%s\"ReadPairMismatchCountDist\": [", dh, dm);
    for(int i = 0; i < maxdstl_se; ++i) ksprintf(s, "%lld,", nmmdist_pe[i]);
    s->s[s->l-1] = ']';
    ksprintf(s, ",\n");
    // ninsdist_pe
    ksprintf(s, "%s%s\"ReadPairInsCountDist\": [", dh, dm);
    for(int i = 0; i < maxdstl_se; ++i) ksprintf(s, "%lld,", ninsdist_pe[i]);
    s->s[s->l-1] = ']';
    ksprintf(s, ",\n");
    // ndeldist_pe
    ksprintf(s, "%s%s\"ReadPairDelCountDist\": [", dh, dm);
    for(int i = 0; i < maxdstl_se; ++i) ksprintf(s, "%lld,", ndeldist_pe[i]);
    s->s[s->l-1] = ']';
    ksprintf(s, ",\n");
    // nvardist_pe
    ksprintf(s, "%s%s\"ReadPairVarCountDist\": [", dh, dm);
    for(int i = 0; i < maxdstl_se; ++i) ksprintf(s, "%lld,", nvardist_pe[i]);
    s->s[s->l-1] = ']';
    ksprintf(s, ",\n");
    // mmdist_pe
    ksprintf(s, "%s%s\"ReadPairTotalDiffLengthDist\": [", dh, dm);
    for(int i = 0; i < maxdstl_se; ++i) ksprintf(s, "%lld,", mmdist_pe[i]);
    s->s[s->l-1] = ']';
    ksprintf(s, "%s}\n", dh);
}

void PAligner::tsvHead(kstring_t* s){
    ksprintf(s, "%s\t%s", "MappedReads", "UnmappedReads");
}

void PAligner::tsvBody(kstring_t* s){
    ksprintf(s, "%lld\t%lld", mapped, unmapped);
}

void PAligner::summary(){
    if(summarized) return;
    maxdstl_se = maxdstl_pe = 0;
    int idx = rlen-1;
    while(nmmdist_se[idx] == 0 && idx > 0) --idx;
    if(maxdstl_se < idx+1) maxdstl_se = idx+1;
    idx = rlen-1;
    while(ninsdist_se[idx] == 0 && idx > 0) --idx;
    if(maxdstl_se < idx+1) maxdstl_se = idx+1;
    idx = rlen-1;
    while(ndeldist_se[idx] == 0 && idx > 0) --idx;
    if(maxdstl_se < idx+1) maxdstl_se = idx+1;
    idx = rlen-1;
    while(nvardist_se[idx] == 0 && idx > 0) --idx;
    if(maxdstl_se < idx+1) maxdstl_se = idx+1;
    idx = rlen-1;
    while(mmdist_se[idx] == 0 && idx > 0) --idx;
    if(maxdstl_se < idx+1) maxdstl_se = idx+1;

    idx = rlen-1;
    while(nmmdist_pe[idx] == 0 && idx > 0) --idx;
    if(maxdstl_pe < idx+1) maxdstl_pe = idx+1;
    idx = rlen-1;
    while(ninsdist_pe[idx] == 0 && idx > 0) --idx;
    if(maxdstl_pe < idx+1) maxdstl_pe = idx+1;
    idx = rlen-1;
    while(ndeldist_pe[idx] == 0 && idx > 0) --idx;
    if(maxdstl_pe < idx+1) maxdstl_pe = idx+1;
    idx = rlen-1;
    while(nvardist_pe[idx] == 0 && idx > 0) --idx;
    if(maxdstl_pe < idx+1) maxdstl_pe = idx+1;
    idx = rlen-1;
    while(mmdist_pe[idx] == 0 && idx > 0) --idx;
    if(maxdstl_pe < idx+1) maxdstl_pe = idx+1;

    summarized = true;
}

void PAligner::tsv4DiffxHead(kstring_t* s, int maxl){
    ksprintf(s, "Ampliton\tItems");
    for(int i = 0; i <= maxl; ++i) ksprintf(s, "\t%d", i);
    kputc('\n', s);
}

void PAligner::tsv4DiffxBodySE(kstring_t* s){
    summary();
    ksprintf(s, "%s\tMismatch", name);
    for(int i = 0; i < maxdstl_se; ++i) ksprintf(s, "\t%lld", nmmdist_se[i]);
    kputc('\n', s);
    ksprintf(s, "%s\tInsertion", name);
    for(int i = 0; i < maxdstl_se; ++i) ksprintf(s, "\t%lld", ninsdist_se[i]);
    kputc('\n', s);
    ksprintf(s, "%s\tDeletion", name);
    for(int i = 0; i < maxdstl_se; ++i) ksprintf(s, "\t%lld", ndeldist_se[i]);
    kputc('\n', s);
    ksprintf(s, "%s\tVariants", name);
    for(int i = 0; i < maxdstl_se; ++i) ksprintf(s, "\t%lld", nvardist_se[i]);
    kputc('\n', s);
    ksprintf(s, "%s\tEditDistance", name);
    for(int i = 0; i < maxdstl_se; ++i) ksprintf(s, "\t%lld", mmdist_se[i]);
    kputc('\n', s);
}

void PAligner::tsv4DiffxBodyPE(kstring_t* s){
    summary();
    ksprintf(s, "%s\tMismatch", name);
    for(int i = 0; i < maxdstl_pe; ++i) ksprintf(s, "\t%lld", nmmdist_pe[i]);
    kputc('\n', s);
    ksprintf(s, "%s\tInsertion", name);
    for(int i = 0; i < maxdstl_pe; ++i) ksprintf(s, "\t%lld", ninsdist_pe[i]);
    kputc('\n', s);
    ksprintf(s, "%s\tDeletion", name);
    for(int i = 0; i < maxdstl_pe; ++i) ksprintf(s, "\t%lld", ndeldist_pe[i]);
    kputc('\n', s);
    ksprintf(s, "%s\tVariants", name);
    for(int i = 0; i < maxdstl_pe; ++i) ksprintf(s, "\t%lld", nvardist_pe[i]);
    kputc('\n', s);
    ksprintf(s, "%s\tEditDistance", name);
    for(int i = 0; i < maxdstl_pe; ++i) ksprintf(s, "\t%lld", mmdist_pe[i]);
    kputc('\n', s);
}
