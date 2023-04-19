#include "processor.h"

void Processor::updateOneBCFHeader(int i){
    kstring_t tmps = {0, 0, 0};
    kstring_t tmpx = {0, 0, 0};
    // edit stats
    int cid = bcf_hdr_id2int(mBcfHeader, BCF_DT_CTG, mGEDetectors[i]->name);
    bcf_hrec_t* hrec = bcf_hdr_id2hrec(mBcfHeader, BCF_DT_CTG, BCF_HL_CTG, cid);
    tmps.l = 0;
    ksprintf(&tmps, "%lld", mGEDetectors[i]->allcnt);
    assert(bcf_hrec_add_key(hrec, "rawcnt", 6) >= 0);
    assert(bcf_hrec_set_val(hrec, hrec->nkeys-1, tmps.s, tmps.l, 0) >= 0);
    tmps.l = 0;
    ksprintf(&tmps, "%lld", mGEDetectors[i]->totcnt);
    assert(bcf_hrec_add_key(hrec, "totcnt", 6) >= 0);
    assert(bcf_hrec_set_val(hrec, hrec->nkeys-1, tmps.s, tmps.l, 0) >= 0);
    tmps.l = 0;
    ksprintf(&tmps, "%lf", mGEDetectors[i]->edieff);
    assert(bcf_hrec_add_key(hrec, "edieff", 6) >= 0);
    assert(bcf_hrec_set_val(hrec, hrec->nkeys-1, tmps.s, tmps.l, 0) >= 0);
    tmps.l = 0;
    ksprintf(&tmps, "%lf", mGEDetectors[i]->recpef);
    assert(bcf_hrec_add_key(hrec, "recpef", 6) >= 0);
    assert(bcf_hrec_set_val(hrec, hrec->nkeys-1, tmps.s, tmps.l, 0) >= 0);
    tmps.l = 0;
    ksprintf(&tmps, "%lf", mGEDetectors[i]->receff);
    assert(bcf_hrec_add_key(hrec, "receff", 6) >= 0);
    assert(bcf_hrec_set_val(hrec, hrec->nkeys-1, tmps.s, tmps.l, 0) >= 0);
    tmps.l = 0;
    ksprintf(&tmps, "%lf", mGEDetectors[i]->reeeff);
    assert(bcf_hrec_add_key(hrec, "reeeff", 6) >= 0);
    assert(bcf_hrec_set_val(hrec, hrec->nkeys-1, tmps.s, tmps.l, 0) >= 0);
    tmps.l = 0;
    ksprintf(&tmps, "%lf", mGEDetectors[i]->muteff);
    assert(bcf_hrec_add_key(hrec, "muteff", 6) >= 0);
    assert(bcf_hrec_set_val(hrec, hrec->nkeys-1, tmps.s, tmps.l, 0) >= 0);
    tmps.l = 0;
    ksprintf(&tmps, "%d", mGEDetectors[i]->donorseq ? 1 : 0);
    assert(bcf_hrec_add_key(hrec, "recstat", 7) >= 0);
    assert(bcf_hrec_set_val(hrec, hrec->nkeys-1, tmps.s, tmps.l, 0) >= 0);
    // single cell if needed
    if(mOptions->dosccal){
        for(auto& sr: mGEDetectors[i]->scretm){
            tmps.l = 0;
            ksprintf(&tmps, "%d-%d-%d-%d", 
                            sr.second->droped,
                            sr.second->edicnt,
                            sr.second->refcnt,
                            sr.second->othcnt);
            tmpx.l = 0;
            ksprintf(&tmpx, "%s%s", SC_ID_PREFIX, sr.first.c_str());
            assert(bcf_hrec_add_key(hrec, tmpx.s, tmpx.l) >= 0);
            assert(bcf_hrec_set_val(hrec, hrec->nkeys-1, tmps.s, tmps.l, 0) >= 0);
        }
    }
    // sync header
    if(bcf_hdr_sync(mBcfHeader) < 0){
        fprintf(stderr, "Error in sync bcf header after update\n");
        exit(EXIT_FAILURE);
    }
    // release res
    if(tmps.s) free(tmps.s);
    if(tmpx.s) free(tmpx.s);
}

void Processor::writeBCF(std::string outBcfPath){
    std::string out = mOutTotBcf;
    if(outBcfPath.size()) out = outBcfPath;
    // first round, update/write header
    for(size_t i = 0; i < mOptions->amplicons.size(); ++i) updateOneBCFHeader(i);
    htsFile* ofp = bcf_open(out.c_str(), "wb");
    assert(bcf_hdr_write(ofp, mBcfHeader) == 0);
    // second round, write rec
    for(size_t i = 0; i < mOptions->amplicons.size(); ++i){
        for(auto& b: mGEDetectors[i]->bcfs){
            assert(bcf_write1(ofp, mBcfHeader, b) >= 0);
        }
    }
    bcf_close(ofp);
    bcf_index_build(out.c_str(), 14);
}

void Processor::writeOneBCF(int i){
    // write the BCF header
    std::string outBCF = mOutBcfs[i];
    htsFile* bcf_fp = bcf_open(outBCF.c_str(), "wb");
    assert(bcf_hdr_write(bcf_fp, mBcfHeader) == 0);
    // write BCF record
    for(auto& b: mGEDetectors[i]->bcfs) assert(bcf_write1(bcf_fp, mBcfHeader, b) == 0);
    // close file
    bcf_close(bcf_fp);
}

void Processor::mergeBCF(std::string outBcfPath){
    std::string out = mOutTotBcf;
    if(outBcfPath.size()) out = outBcfPath;
    htsFile* fp = bcf_open(out.c_str(), "wb");
    assert(bcf_hdr_write(fp, mBcfHeader) >= 0);
    bcf1_t* rec = bcf_init1();
    for(auto& s: mOutBcfs){
        htsFile* tp = bcf_open(s.c_str(), "r");
        bcf_hdr_t* h = bcf_hdr_read(tp);
        rec->max_unpack = 0;
        while(bcf_read1(tp, h, rec) >= 0){
            assert(bcf_write1(fp, mBcfHeader, rec) >= 0);
        }
        bcf_close(tp);
        bcf_hdr_destroy(h);
    }
    bcf_destroy(rec);
    bcf_close(fp);
    bcf_index_build(out.c_str(), 14);
}
