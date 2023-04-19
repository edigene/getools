#include "gedetect.h"
#include "bamplp.h"

void GEDetector::bam2bcf(){
    if(!vsorted){
        std::sort(vars.begin(), vars.end(), BamSortByCoord());
        vsorted = 1;
    }
    // cal samples
    int n = 1;
    if(opt->ctrl.size()) n = 2;
    // init aux data
    mplp_aux_t** data = (mplp_aux_t**)calloc(n, sizeof(mplp_aux_t*));
    data[0] = (mplp_aux_t*)calloc(1, sizeof(mplp_aux_t));
    data[0]->bams = &vars;
    data[0]->idx = 0;
    data[0]->dropmask = dropmask;
    data[0]->iter = NULL;
    if(n == 2){
        data[1] = (mplp_aux_t*)calloc(1, sizeof(mplp_aux_t));
        data[1]->fp = sam_open(opt->ctrl.c_str(), "r");
        data[1]->hdr = sam_hdr_read(data[1]->fp);
        data[1]->sidx = sam_index_load(data[1]->fp, opt->ctrl.c_str());
        data[1]->iter = sam_itr_queryi(data[1]->sidx, sam_hdr_name2tid(data[1]->hdr, name), 0, rlen);
        data[1]->dropmask = dropmask;
    }
    // init pileup results
    mplp_pileup_t gplp;
    mplp_pileup_init(&gplp, n);
    // calling results
    call_ret_t** cr = (call_ret_t**)calloc(n, sizeof(call_ret_t*));
    for(int i = 0; i < n; ++i) cr[i] = new call_ret_t();
    cr[0]->totcnt = totcnt;
    if(n == 2) cr[1]->totcnt = ctrltt;
    // begin pileup and calling
    bam_mplp_t iter = bam_mplp_init(n, mplp_func, (void**)data);
    bam_mplp_set_maxcnt(iter, MAX(vars.size(), INT_MAX));
    bam_mplp_init_overlaps(iter);
    bam_mplp_constructor(iter, pileup_constructor);
    bam_mplp_destructor(iter, pileup_destructor);
    int tid = -1, pos = -1;
    while(bam_mplp_auto(iter, &tid, &pos, gplp.n_plp, gplp.plp) > 0){
        if(tid == rtid && pos >= 0 && pos < rlen){ // valid pos and tid
            for(int j = 0; j < n; ++j) plp_call(rtid, pos, gplp.plp[j], gplp.n_plp[j], cr[j], vvmask);
            // call BCF record
            if(n == 1) call2bcf1(cr[0], opt->bcfh, bcfs, ref);
            else if(n == 2) call2bcf2(cr[0], cr[1], opt->bcfh, bcfs, ref);
        }
    }
    // clean res
    for(int i = 0; i < n; ++i) mplp_aux_destroy(data[i]);
    free(data);
    mplp_pileup_destroy(&gplp);
    for(int i = 0; i < n; ++i) delete cr[i];
    free(cr);
    bam_mplp_destroy(iter);
}
