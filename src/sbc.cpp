#include "sbc.h"

spl_res_t::spl_res_t(){
    memset(b1mm, 0, BARCODE_MATCH_STATUS_ARR_LEN*sizeof(int64_t));
    memset(b2mm, 0, BARCODE_MATCH_STATUS_ARR_LEN*sizeof(int64_t));
    memset(bbmm, 0, BARCODE_MATCH_STATUS_ARR_LEN*sizeof(int64_t));
    memset(r2drop, 0, BARCODE_MATCH_STATUS_ARR_LEN*sizeof(int64_t));
    memset(offc, 0, (MISSIONBIO_BARCODE_CONST_MAX_OFFSET+1)*sizeof(int64_t));
    rcmap = (int64_t**)malloc(MISSIONBIO_BARCODE_COUNT*sizeof(int64_t*));
    for(int i = 0; i < MISSIONBIO_BARCODE_COUNT; ++i){
        rcmap[i] = (int64_t*)calloc(MISSIONBIO_BARCODE_COUNT, sizeof(int64_t));
    }
}

spl_res_t::~spl_res_t(){
    if(rcmap){
        for(int i = 0; i < MISSIONBIO_BARCODE_COUNT; ++i){
            free(rcmap[i]);
            rcmap[i] = NULL;
        }
        free(rcmap);
        rcmap = NULL;
    }
}

void spl_res_t::res2jsn(kstring_t* s){
    ksprintf(s, "  \"CellBarcodeParseStat\":{\n");
    ksprintf(s, "    \"TotalReadToParse\": %lld,\n", totr);
    ksprintf(s, "    \"ReadWithCellBarcode\": %lld,\n", gotr);
    ksprintf(s, "    \"ReadWithoutCellBarcode\": %lld,\n", dropr);
    ksprintf(s, "    \"ReadWithCellBarcodeRate\": %lf\n  },\n", (double)gotr/(double)totr);
    ksprintf(s, "  \"ConstSeq1OffsetStat\":{\n");
    for(int i = 0; i <= MISSIONBIO_BARCODE_CONST_MAX_OFFSET; ++i){
        ksprintf(s, "    \"Offset%d\": %lld,\n", i, offc[i]);
    }
    s->l -= 2;
    ksprintf(s, "\n  },\n");
    ksprintf(s, "  \"CellBarcode1MismatchStat\":{\n");
    ksprintf(s, "    \"Mismatch0\": %lld,\n", b1mm[BARCODE_MISMATCH0_CNTARR_INDEX]);
    ksprintf(s, "    \"Mismatch1\": %lld,\n", b1mm[BARCODE_MISMATCH1_CNTARR_INDEX]);
    ksprintf(s, "    \"Mismatch2\": %lld,\n", b1mm[BARCODE_MISMATCH2_CNTARR_INDEX]);
    ksprintf(s, "    \"WithOneIns\": %lld,\n", b1mm[BARCODE_INS1_CNTARR_INDEX]);
    ksprintf(s, "    \"WithOneDel\": %lld\n  },\n", b1mm[BARCODE_DEL1_CNTARR_INDEX]);
    ksprintf(s, "  \"CellBarcode2MismatchStat\":{\n");
    ksprintf(s, "    \"Mismatch0\": %lld,\n", b2mm[BARCODE_MISMATCH0_CNTARR_INDEX]);
    ksprintf(s, "    \"Mismatch1\": %lld,\n", b2mm[BARCODE_MISMATCH1_CNTARR_INDEX]);
    ksprintf(s, "    \"Mismatch2\": %lld,\n", b2mm[BARCODE_MISMATCH2_CNTARR_INDEX]);
    ksprintf(s, "    \"WithOneIns\": %lld,\n", b2mm[BARCODE_INS1_CNTARR_INDEX]);
    ksprintf(s, "    \"WithOneDel\": %lld\n  },\n", b2mm[BARCODE_DEL1_CNTARR_INDEX]);
    ksprintf(s, "  \"WholeCellBarcodeMismatchStat\":{\n");
    ksprintf(s, "    \"Mismatch0\": %lld,\n", bbmm[BARCODE_MISMATCH0_CNTARR_INDEX]);
    ksprintf(s, "    \"Mismatch1\": %lld,\n", bbmm[BARCODE_MISMATCH1_CNTARR_INDEX]);
    ksprintf(s, "    \"Mismatch2\": %lld,\n", bbmm[BARCODE_MISMATCH2_CNTARR_INDEX]);
    ksprintf(s, "    \"Mismatch3\": %lld,\n", bbmm[BARCODE_MISMATCH3_CNTARR_INDEX]);
    ksprintf(s, "    \"Mismatch4\": %lld\n  },\n", bbmm[BARCODE_MISMATCH4_CNTARR_INDEX]);
    ksprintf(s, "  \"CellBarcodeParseFailReasonCount\":{\n");
    ksprintf(s, "    \"ConstSeq1NotFound\": %lld,\n", r2drop[BARCODE_DROP_BY_CONSTSEQ1_NOT_FOUND]);
    ksprintf(s, "    \"ConstSeq2NotFound\": %lld,\n", r2drop[BARCODE_DROP_BY_CONSTSEQ2_NOT_FOUND]);
    ksprintf(s, "    \"NineBpBarcodesMismatchTooMany\": %lld,\n", r2drop[BARCODE_DROP_BY_NINEBP_MM_TOO_MANY]);
    ksprintf(s, "    \"NineBpBarcodesSameMismatchHitTooMany\": %lld,\n", r2drop[BARCODE_DROP_BY_NINEBP_HIT_TOO_MANY]);
    ksprintf(s, "    \"IndelTooMany\": %lld,\n", r2drop[BARCODE_DROP_BY_INDEL_TOO_MANY]);
    ksprintf(s, "    \"SameIndelHitTooMany\": %lld,\n", r2drop[BARCODE_DROP_BY_INDEL_HIT_TOO_MANY]);
    ksprintf(s, "    \"ReadTooShort\": %lld\n  }", r2drop[BARCODE_DROP_BY_READ_TOO_SHORT]);
}

sbc_biom_t::sbc_biom_t(){
    bcidx = (int32_t**)malloc(MISSIONBIO_BARCODE_COUNT*sizeof(int64_t*));
    int accidx = 0;
    for(int i = 0; i < MISSIONBIO_BARCODE_COUNT; ++i){
        bcidx[i] = (int32_t*)calloc(MISSIONBIO_BARCODE_COUNT, sizeof(int64_t));
        for(int j = 0; j < MISSIONBIO_BARCODE_COUNT; ++j){
            bcidx[i][j] = ++accidx;
        }
    }
}

sbc_biom_t::~sbc_biom_t(){
    if(wr1){ delete wr1; wr1 = NULL; }
    if(wr2){ delete wr2; wr2 = NULL; }
    if(dw1){ delete dw1; dw1 = NULL; }
    if(dw2){ delete dw2; dw2 = NULL; }
    for(auto& r: br1){
        if(r && r->s){ free(r->s); free(r); }
    }
    for(auto& r: br2){
        if(r && r->s){ free(r->s); free(r); }
    }
    for(auto& r: dbr1){
        if(r && r->s){ free(r->s); free(r); }
    }
    for(auto& r: dbr2){
        if(r && r->s){ free(r->s); free(r); }
    }
    if(tpl) delete tpl;
    if(jks){
        if(jks->s) free(jks->s);
        free(jks);
        jks = NULL;
    }
    if(bcidx){
        for(int i = 0; i < MISSIONBIO_BARCODE_COUNT; ++i){
            free(bcidx[i]);
            bcidx[i] = NULL;
        }
        free(bcidx);
        bcidx = NULL;
    }
}

bool sbc_biom_t::valid4idx(){
    if(wlist == NULL || owidx == NULL){
        fprintf(stderr, "input cellbarcode list file and output index file must both be provided\n");
        return false;
    }
    return true;
}

bool sbc_biom_t::valid4spl(){
    if(inr1 == NULL){
        fprintf(stderr, "read1 input must be provided\n");
        return false;
    }
    if(thread <= 0){
        fprintf(stderr, "thread number must be positive number\n");
        return false;
    }
    if(maxp < 1){
        fprintf(stderr, "packs must be at least 1\n");
        return false;
    }
    if(maxr < 100000){
        fprintf(stderr, "reads in each packs should be at least 100000\n");
        return false;
    }
    if(wlist == NULL && iwidx == NULL){
        fprintf(stderr, "either the cellbarcode whitelist or the cellbarcode index file must be provided\n");
        return false;
    }
    return true;
}

void sbc_biom_t::init(){
    if(!util::exists(outdir)) util::makedir(outdir);
    outdir = util::abspath(outdir);
    outr1 = util::joinpath(outdir, outr1);
    outr2 = util::joinpath(outdir, outr2);
    outc1 = util::joinpath(outdir, outc1);
    outc2 = util::joinpath(outdir, outc2);
    outd1 = util::joinpath(outdir, outd1);
    outd2 = util::joinpath(outdir, outd2);
    outtsv = util::joinpath(outdir, outtsv);
    outjsn = util::joinpath(outdir, outjsn);
    jks = (kstring_t*)calloc(1, sizeof(kstring_t));
    ksprintf(jks, "{\n");
    // index
    if(iwidx){
        file2idx();
    }else if(wlist){
        util::makeListFromFileByLine(wlist, warray);
        if(warray.size() > MISSIONBIO_BARCODE_COUNT){
            fprintf(stderr, "ERROR: too many barcode provided[%ld], only %d is valid\n", warray.size(), MISSIONBIO_BARCODE_COUNT);
            exit(EXIT_FAILURE);
        }
        for(auto& e: warray){
            if(e.size() != MISSIONBIO_BARCODE_LENGTH){
                fprintf(stderr, "ERROR: invalid barcode %s [%ldbp], only %dbp is valid\n", e.c_str(),  e.size(), MISSIONBIO_BARCODE_LENGTH);
                exit(EXIT_FAILURE);
            }
        }
        bcs2mmidx();
        bcs2ididx();
    }
    // temp val
    cs1l = strlen(cs1s);
    cs2l = strlen(cs2s);
    sbeg1 = bc1l-1;
    sbeg2 = bc1l-1+cs1l+bc2l-1;
    srng1 = maxoff+cs1l+3;
    srng2 = maxoff+cs2l+6;
    wr1 = new Writer(outr1);
    dw1 = new Writer(outd1);
    if(inr2){
        wr2 = new Writer(outr2);
        dw2 = new Writer(outd2);
    }
    br1.resize(thread, NULL);
    for(int i = 0; i < thread; ++i) br1[i] = (kstring_t*)calloc(1, sizeof(kstring_t));
    br2.resize(thread, NULL);
    for(int i = 0; i < thread; ++i) br2[i] = (kstring_t*)calloc(1, sizeof(kstring_t));
    dbr1.resize(thread, NULL);
    for(int i = 0; i < thread; ++i) dbr1[i] = (kstring_t*)calloc(1, sizeof(kstring_t));
    dbr2.resize(thread, NULL);
    for(int i = 0; i < thread; ++i) dbr2[i] = (kstring_t*)calloc(1, sizeof(kstring_t));
    tpl = new ThreadPool(thread);
    minl = bc1l + bc2l + cs1l + cs2l;
    sprs.resize(thread, NULL);
    for(int i = 0; i < thread; ++i) sprs[i] = new spl_res_t();
}

void sbc_biom_t::idx2file(){
    // load whitelist
    util::makeListFromFileByLine(wlist, warray);
    // build index
    bcs2mmidx();
    bcs2ididx();
    // write index
    uint32_t t = 0;
    kstring_t* m = (kstring_t*)calloc(1, sizeof(kstring_t));
    kstring_t* s = (kstring_t*)calloc(1, sizeof(kstring_t));
#ifdef GET_SBC_DEBUG
    ksprintf(s, "#idxtype,idxkey,idxval,mismatch,count,oriseq,mutseq,chkdist\n");
#else
    ksprintf(s, "#idxtype,idxkey,idxval,mismatch,count,oriseq,mutseq\n");
#endif
    // mm
    for(size_t i = 0; i < wmindex.size(); ++i){
        if(get_barcode_mismatch(wmindex[i]) != MISSIONBIO_BARCODE_INDEX_DEFAULT){
            m->l = 0;
            t = i;
            for(int k = 0; k < MISSIONBIO_BARCODE_LENGTH; ++k){
                kputc("ACGT"[t & 0x3], m);
                t >>= 2;
            }
            util::reverse(m->s, m->l);
#ifdef GET_SBC_DEBUG
            ksprintf(s, "0,%ld,%d,%d,%d,%s,%s,%d\n",
#else
            ksprintf(s, "0,%ld,%d,%d,%d,%s,%s\n",
#endif
                        i,
                        wmindex[i],
                        get_barcode_mismatch(wmindex[i]),
                        get_barcode_bccount(wmindex[i]),
                        warray[get_bardoce_bcindex(wmindex[i])].c_str(),
#ifdef GET_SBC_DEBUG
                        m->s,
                        edit_distance_dp(m->s, m->l, warray[get_bardoce_bcindex(wmindex[i])].c_str(), warray[get_bardoce_bcindex(wmindex[i])].size())
#else
                        m->s
#endif
                        );

        }
    }
    // ins
    for(size_t i = 0; i < wiindex.size(); ++i){
        if(get_barcode_mismatch(wiindex[i]) != MISSIONBIO_BARCODE_INDEX_DEFAULT){
            m->l = 0;
            t = i;
            for(int k = 0; k < MISSIONBIO_BARCODE_LENGTH+1; ++k){
                kputc("ACGT"[t & 0x3], m);
                t >>= 2;
            }
            util::reverse(m->s, m->l);
#ifdef GET_SBC_DEBUG
            ksprintf(s, "1,%ld,%d,%d,%d,%s,%s,%d\n",
#else
            ksprintf(s, "1,%ld,%d,%d,%d,%s,%s\n",
#endif
                        i,
                        wiindex[i],
                        get_barcode_mismatch(wiindex[i]),
                        get_barcode_bccount(wiindex[i]),
                        warray[get_bardoce_bcindex(wiindex[i])].c_str(),
#ifdef GET_SBC_DEBUG
                        m->s,
                        edit_distance_dp(m->s, m->l, warray[get_bardoce_bcindex(wiindex[i])].c_str(), warray[get_bardoce_bcindex(wiindex[i])].size())
#else
                        m->s
#endif
                    );

        }
    }
    // del
    for(size_t i = 0; i < wdindex.size(); ++i){
        if(get_barcode_mismatch(wdindex[i]) != MISSIONBIO_BARCODE_INDEX_DEFAULT){
            m->l = 0;
            t = i;
            for(int k = 0; k < MISSIONBIO_BARCODE_LENGTH-1; ++k){
                kputc("ACGT"[t & 0x3], m);
                t >>= 2;
            }
            util::reverse(m->s, m->l);
#ifdef GET_SBC_DEBUG
            ksprintf(s, "2,%ld,%d,%d,%d,%s,%s,%d\n",
#else
            ksprintf(s, "2,%ld,%d,%d,%d,%s,%s\n",
#endif
                        i,
                        wdindex[i],
                        get_barcode_mismatch(wdindex[i]),
                        get_barcode_bccount(wdindex[i]),
                        warray[get_bardoce_bcindex(wdindex[i])].c_str(),
#ifdef GET_SBC_DEBUG
                        m->s,
                        edit_distance_dp(m->s, m->l, warray[get_bardoce_bcindex(wdindex[i])].c_str(), warray[get_bardoce_bcindex(wdindex[i])].size())
#else
                        m->s
#endif
                    );

        }
    }
    FILE* fout = fopen(owidx, "w");
    fwrite(s->s, sizeof(char), s->l, fout);
    fclose(fout);
    if(s->s) free(s->s);
    free(s);
    if(m->s) free(m->s);
    free(m);
}

void sbc_biom_t::file2idx(){
    wmindex.resize(1UL << (2*MISSIONBIO_BARCODE_LENGTH), MISSIONBIO_BARCODE_INDEX_DEFAULT);
    wiindex.resize(1UL << (2*(MISSIONBIO_BARCODE_LENGTH+1)), MISSIONBIO_BARCODE_INDEX_DEFAULT);
    wdindex.resize(1UL << (2*(MISSIONBIO_BARCODE_LENGTH-1)), MISSIONBIO_BARCODE_INDEX_DEFAULT);
    std::unordered_map<std::string, uint32_t> mbc;
    util::LineReader lr(iwidx);
    std::string line;
    std::vector<std::string> vstr;
    lr.getline(line);
    uint32_t idxtype, idxkey, idxval, mismatch;
    while(lr.getline(line)){
        util::split(line, vstr, ",");
        idxtype = atol(vstr[0].c_str());
        idxkey = atol(vstr[1].c_str());
        idxval = atol(vstr[2].c_str());
        mismatch = atol(vstr[3].c_str());
        if(mismatch == 0){ // original seq found
            mbc[vstr[5]] = get_bardoce_bcindex(idxval);
        }
        if(idxtype == 0){
            wmindex[idxkey] = idxval;
        }else if(idxtype == 1){
            wiindex[idxkey] = idxval;
        }else if(idxtype == 2){
            wdindex[idxkey] = idxval;
        }
        
    }
    warray.resize(mbc.size());
    for(auto& p: mbc) warray[p.second] = p.first;
}

void sbc_biom_t::bcs2mmidx(){
    wmindex.resize(1UL << (2*MISSIONBIO_BARCODE_LENGTH), MISSIONBIO_BARCODE_INDEX_DEFAULT);
    uint32_t fmask = (1UL << (2*MISSIONBIO_BARCODE_LENGTH))-1;
    uint32_t tval = 0, tidx = 0;
    // original seq index
    for(size_t i = 0; i < warray.size(); ++i){
        tval = 0;
        for(size_t j = 0; j < warray[i].size(); ++j){
            tval <<= 2;
            tval |= nuc_to_2bit[(int)warray[i][j]];
        }
        set_barcode_index(tidx, i, 1, 0);
        wmindex[tval] = tidx;
    }
    uint32_t orib = 0;
    uint32_t nval = 0;
    std::vector<uint32_t> smasks;
    smasks.resize(MISSIONBIO_BARCODE_LENGTH, 0);
    for(size_t i = 0; i < smasks.size(); ++i) smasks[i] = fmask & (~(3UL << 2 * i));
    // 1 mismatch seq
    for(size_t i = 0; i < warray.size(); ++i){
        tval = 0;
        for(size_t j = 0; j < warray[i].size(); ++j){
            tval <<= 2;
            tval |= nuc_to_2bit[(int)warray[i][j]];
        }
        for(size_t k = 0; k < warray[i].size(); ++k){
            orib = nuc_to_2bit[(int)warray[i][k]];
            size_t j = warray[i].size()-1-k;
            for(uint32_t b = 0; b < 4; ++b){
                if(b != orib){ // change the bit
                    nval = tval & smasks[j];
                    nval |= (b << 2 * j);
                    tidx = 0;
                    int cnt = get_barcode_bccount(wmindex[nval]);
                    size_t sid = get_bardoce_bcindex(wmindex[nval]);
                    int xmm = get_barcode_mismatch(wmindex[nval]);
                    if(xmm == MISSIONBIO_BARCODE_INDEX_DEFAULT){ // new
                        set_barcode_index(tidx, i, 1, 1);
                        wmindex[nval] = tidx;
                    }else{
                        if(xmm == 1 && sid != i){ // existing
                            set_barcode_index(tidx, i, (1+cnt), 1);
                            wmindex[nval] = tidx;
                        }
                    }
                }
            }
        }
    }
    // 2 mismatch seq
    uint32_t oric = 0;
    for(size_t i = 0; i < warray.size(); ++i){
        tval = 0;
        for(size_t j = 0; j < warray[i].size(); ++j){
            tval <<= 2;
            tval |= nuc_to_2bit[(int)warray[i][j]];
        }
        for(size_t j = 0; j < warray[i].size(); ++j){
            for(size_t k = j+1; k < warray[i].size(); ++k){
                size_t xj = warray[i].size()-1-j;
                size_t xk = warray[i].size()-1-k;
                orib = nuc_to_2bit[(int)warray[i][j]];
                oric = nuc_to_2bit[(int)warray[i][k]];
                for(uint8_t b = 0; b < 4; ++b){
                    for(uint8_t c = 0; c < 4; ++c){
                        if(b != orib && c != oric){
                            nval = tval & smasks[xj] & smasks[xk];
                            nval |= (b << 2 * xj);
                            nval |= (c << 2 * xk);
                            tidx = 0;
                            int cnt = get_barcode_bccount(wmindex[nval]);
                            size_t sid = get_bardoce_bcindex(wmindex[nval]);
                            int xmm = get_barcode_mismatch(wmindex[nval]);
                            if(xmm == MISSIONBIO_BARCODE_INDEX_DEFAULT){ // new
                                set_barcode_index(tidx, i, 1, 2);
                                wmindex[nval] = tidx;
                            }else{
                                if(xmm == 2 && sid != i){
                                    set_barcode_index(tidx, i, (1+cnt), 2);
                                    wmindex[nval] = tidx;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void sbc_biom_t::bcs2ididx(){
    wiindex.resize(1UL << (2*(MISSIONBIO_BARCODE_LENGTH+1)), MISSIONBIO_BARCODE_INDEX_DEFAULT);
    wdindex.resize(1UL << (2*(MISSIONBIO_BARCODE_LENGTH-1)), MISSIONBIO_BARCODE_INDEX_DEFAULT);
    uint32_t fmask = (1UL << (2*MISSIONBIO_BARCODE_LENGTH))-1;
    uint32_t tval = 0, tidx = 0, nval = 0, lval = 0, rval = 0;
    std::vector<uint32_t> lmasks, rmasks;
    lmasks.resize(MISSIONBIO_BARCODE_LENGTH, 0);
    rmasks.resize(MISSIONBIO_BARCODE_LENGTH, 0);
    for(int i = 0; i < MISSIONBIO_BARCODE_LENGTH; ++i){
        lmasks[i] = (fmask << (2 * i)) & fmask;
        rmasks[i] = (lmasks[i] ^ fmask) & fmask;
    }
    for(size_t i = 0; i < warray.size(); ++i){
        tval = 0;
        for(size_t j = 0; j < warray[i].size(); ++j){
            tval <<= 2;
            tval |= nuc_to_2bit[(int)warray[i][j]];
        }
        // 1bp ins
        for(size_t j = 0; j <= warray[i].size(); ++j){
            if(j < warray[i].size()){
                lval = tval & lmasks[j];
                rval = tval & rmasks[j];
            }else{
                lval = 0;
                rval = tval;
            }
            for(uint8_t b = 0; b < 4; ++b){
                nval = (lval << 2) | (b << 2 *j);
                nval |= rval;
                int cnt = get_barcode_bccount(wiindex[nval]);
                size_t sid = get_bardoce_bcindex(wiindex[nval]);
                int xmm = get_barcode_mismatch(wiindex[nval]);
                if(xmm == MISSIONBIO_BARCODE_INDEX_DEFAULT){ // new
                    tidx = 0;
                    set_barcode_index(tidx, i, 1, 1);
                    wiindex[nval] = tidx;
                }else{ // already exist
                    if(sid != i){// from different barcode, accumulate count
                        tidx = 0;
                        set_barcode_index(tidx, i, (cnt+1), 1);
                        wiindex[nval] = tidx;
                    }
                }
            }
        }
        // 1bp del
        for(size_t j = 0; j < warray[i].size(); ++j){
            lval = tval & lmasks[j];
            rval = tval & rmasks[j];
            nval = ((lval >> 2) & lmasks[j]) | rval;
            int cnt = get_barcode_bccount(wdindex[nval]);
            size_t sid = get_bardoce_bcindex(wdindex[nval]);
            int xmm = get_barcode_mismatch(wdindex[nval]);
            if(xmm == MISSIONBIO_BARCODE_INDEX_DEFAULT){ // new
                tidx = 0;
                set_barcode_index(tidx, i, 1, 1);
                wdindex[nval] = tidx;
            }else{ // already exist
                if(sid != i){// from different barcode, accumulate count
                    tidx = 0;
                    set_barcode_index(tidx, i, (cnt+1), 1);
                    wdindex[nval] = tidx;
                }
            }
        }
    }
}

void sbc_biom_t::split_all(){
    init();
    ReadGenerator* rg = new ReadGenerator(inr1, inr2); // read generator
    rg->init(maxp, maxr);
    ReadPack** rp =(ReadPack**)calloc(2, sizeof(ReadPack*));
    rg->start();
    std::vector<std::pair<int32_t, int32_t>> vpidx;
    std::vector<std::future<void>> rets;
    while(rg->getPack(rp)){
        int n = util::divideVecIdx(rp[0]->n, thread, vpidx);
        rets.clear();
        for(int i = 0; i < n; ++i){
            rets.push_back(tpl->enqueue(&sbc_biom_t::split_rng, this, rp, vpidx[i].first, vpidx[i].second, i));
        }
        for(auto& e: rets) e.get();
        rp[0]->r = false;
        rp[0]->w = true;
        if(inr2){
            rp[1]->r = false;
            rp[1]->w = true;
        }
    }
    delete rg; rg = NULL;
    free(rp); rp = NULL;
    // close raw writer
    if(wr1){ delete wr1; wr1 = NULL; }
    if(wr2){ delete wr2; wr2 = NULL; }
    // report
    for(size_t i = 1; i < sprs.size(); ++i){
        sprs[0]->dropr += sprs[i]->dropr;
        sprs[0]->gotr += sprs[i]->gotr;
        sprs[0]->totr += sprs[i]->totr;
        for(int ci = 0; ci < BARCODE_DROP_STATUS_ARR_LEN; ++ci){
            sprs[0]->r2drop[ci] += sprs[i]->r2drop[ci];
        }
        for(int ci = 0; ci < BARCODE_MATCH_STATUS_ARR_LEN; ++ci){
            sprs[0]->b1mm[ci] += sprs[i]->b1mm[ci];
            sprs[0]->b2mm[ci] += sprs[i]->b2mm[ci];
            sprs[0]->bbmm[ci] += sprs[i]->bbmm[ci];
        }
        for(int oi = 0; oi <= MISSIONBIO_BARCODE_CONST_MAX_OFFSET; ++oi){
            sprs[0]->offc[oi] += sprs[i]->offc[oi];
        }
        for(int ci = 0; ci < MISSIONBIO_BARCODE_COUNT; ++ci){
            for(int cj = 0; cj < MISSIONBIO_BARCODE_COUNT; ++cj){
                sprs[0]->rcmap[ci][cj] += sprs[i]->rcmap[ci][cj];
            }
        }
    }
    FILE* ofp = fopen(outtsv.c_str(), "w");
    kstring_t ks = {0, 0, 0};
    for(int ci = 0; ci < MISSIONBIO_BARCODE_COUNT; ++ci){
        for(int cj = 0; cj < MISSIONBIO_BARCODE_COUNT; ++cj){
            if(sprs[0]->rcmap[ci][cj]){
                ksprintf(&ks, "%d\t%s%s\t%lld\n",
                              bcidx[ci][cj],
                              warray[ci].c_str(), 
                              warray[cj].c_str(), 
                              sprs[0]->rcmap[ci][cj]);
            }
        }
    }
    if(ks.s){
        fwrite(ks.s, sizeof(char), ks.l, ofp);
        free(ks.s);
    }
    fclose(ofp);
    sprs[0]->res2jsn(jks);
}

void sbc_biom_t::raw2clean_all(){
    // initialize
    for(int i = 0; i < thread; ++i){
        sprs[i]->dropr = sprs[i]->gotr = sprs[i]->totr = 0;
        if(i > 0){
            for(int j = 0; j < MISSIONBIO_BARCODE_COUNT; ++j){
                memset(sprs[i]->rcmap, 0, MISSIONBIO_BARCODE_COUNT*sizeof(int64_t));
            }
        }
    }
    wr1 = new Writer(outc1);
    if(inr2) wr2 = new Writer(outc2);
    ReadGenerator* rg = new ReadGenerator(outr1.c_str(), outr2.c_str()); // read generator
    rg->init(maxp, maxr);
    ReadPack** rp =(ReadPack**)calloc(2, sizeof(ReadPack*));
    rg->start();
    std::vector<std::pair<int32_t, int32_t>> vpidx;
    std::vector<std::future<void>> rets;
    while(rg->getPack(rp)){
        int n = util::divideVecIdx(rp[0]->n, thread, vpidx);
        rets.clear();
        for(int i = 0; i < n; ++i){
            rets.push_back(tpl->enqueue(&sbc_biom_t::raw2clean_rng, this, rp, vpidx[i].first, vpidx[i].second, i));
        }
        for(auto& e: rets) e.get();
        rp[0]->r = false;
        rp[0]->w = true;
        if(inr2){
            rp[1]->r = false;
            rp[1]->w = true;
        }
    }
    delete rg; rg = NULL;
    free(rp); rp = NULL;
    // close raw writer
    if(wr1){ delete wr1; wr1 = NULL; }
    if(wr2){ delete wr2; wr2 = NULL; }
    // report
    for(size_t i = 1; i < sprs.size(); ++i){
        sprs[0]->dropr += sprs[i]->dropr;
        sprs[0]->gotr += sprs[i]->gotr;
        sprs[0]->totr += sprs[i]->totr;
    }
    ksprintf(jks, ",\n");
    ksprintf(jks, "  \"BarcodeSupportFilteringStat\": {\n");
    ksprintf(jks, "      \"ReadsWithHighAbundanceCellBarcode\": %lld,\n", sprs[0]->gotr); 
    ksprintf(jks, "      \"ReadsWithLowAbundanceCellBarcode\": %lld,\n", sprs[0]->dropr); 
    ksprintf(jks, "      \"HighAbundanceCellBarcodeRate\": %lf\n  }\n", sprs[0]->gotr ? (double)sprs[0]->gotr/(double)sprs[0]->totr : .0);
    ksprintf(jks, "}\n");
    FILE* fw = fopen(outjsn.c_str(), "w");
    fwrite(jks->s, sizeof(char), jks->l, fw);
    fclose(fw);
}

void sbc_biom_t::split_rng(ReadPack** rp, int beg, int end, int idx){
    kstring_t* s1 = br1[idx];
    kstring_t* s2 = br2[idx];
    kstring_t* ds1 = dbr1[idx];
    kstring_t* ds2 = dbr2[idx];
    spl_res_t* spr = sprs[idx];
    s1->l = s2->l = ds1->l = ds2->l = 0;
    krec1_t* r1 = NULL;
    krec1_t* r2 = NULL;
    char* s1p = NULL;
    char* s2p = NULL;
    bool add = true;
    int maxbc1l = 0, maxbc2l = 0, offn = 0, zzn = 0;
    int32_t bc1hid = -1, bc2hid = -1, droprx = -1, dropry = -1;
    int32_t bc1mm = 0, bc2mm = 0, bc1ins = 0, bc2ins = 0, bc1del = 0, bc2del = 0;
    for(int i = beg; i < end; ++i){
        bc1mm = bc2mm = bc1ins = bc2ins = bc1del = bc2del = 0;
        bc1hid = bc2hid = droprx = dropry = -1;
        ++spr->totr;
        add = false;
        r1 = rp[0]->reads[i];
        if(inr2) r2 = rp[1]->reads[i];
        if(r1->seq.l < (size_t)minl){
            ++spr->dropr;
            ++spr->r2drop[BARCODE_DROP_BY_READ_TOO_SHORT];
            if(outdrop){
                ksprintf(ds1, "@%s %s\n%s\n+\n%s\n", 
                              r1->name.s, 
                              MISSIONBIO_BARCODE_DROP_REASON[BARCODE_DROP_BY_READ_TOO_SHORT], 
                              r1->seq.s, 
                              r1->qual.s);
                if(r2) ksprintf(ds2, "@%s\n%s\n+\n%s\n", r2->name.s, r2->seq.s, r2->qual.s);
            }
            continue;
        }
        s1p = util::strnstr(r1->seq.s+sbeg1, cs1s, cs1l, srng1); 
        if(s1p){// found const sequence1
            offn = 0;
            zzn = MIN(s1p-r1->seq.s-bc1l, MISSIONBIO_BARCODE_CONST_MAX_OFFSET);
            for(int xxn = 0; xxn < zzn; ++xxn){
                if(*(s1p-1-xxn) == MISSIONBIO_BARCODE_CONST_MAX_OFFSTR[xxn]){
                    offn = xxn+1;
                }else{
                    break;
                }
            }
            s1p -= offn;
            s2p = util::strnstr(r1->seq.s+sbeg2, cs2s, cs2l, srng2);
            if(s2p){// found const sequence2
                maxbc1l = s1p-r1->seq.s; // barcode1 max length
                maxbc2l = s2p-s1p-cs1l-offn; // barcode2 length
                if(maxbc2l>bc2l+1 || maxbc2l+1<bc2l ||
                   maxbc1l+1<bc1l || maxbc1l>bc1l+1){ // too long or too short
                    ++spr->r2drop[BARCODE_DROP_BY_INDEL_TOO_MANY];
                    if(outdrop){
                        ksprintf(ds1, "@%s %s\n%s\n+\n%s\n", 
                                      r1->name.s, 
                                      MISSIONBIO_BARCODE_DROP_REASON[BARCODE_DROP_BY_INDEL_TOO_MANY], 
                                      r1->seq.s, 
                                      r1->qual.s);
                        if(r2) ksprintf(ds2, "@%s\n%s\n+\n%s\n", r2->name.s, r2->seq.s, r2->qual.s);
                    }
                }else{
                    uint32_t bc2int = 0UL;
                    for(int bi = 0; bi < maxbc2l; ++bi){
                        bc2int <<= 2;
                        bc2int |= nuc_to_2bit[(int)s1p[cs1l+offn+bi]];
                    }
                    uint32_t hitr = 0;
                    if(maxbc2l == bc2l){
                        hitr = wmindex[bc2int];
                        if(get_barcode_mismatch(hitr) == MISSIONBIO_BARCODE_INDEX_DEFAULT){
                            if(dropry < 0) dropry = BARCODE_DROP_BY_NINEBP_MM_TOO_MANY;
                        }else{
                            if(get_barcode_bccount(hitr) > 1){
                                if(dropry < 0) dropry = BARCODE_DROP_BY_NINEBP_HIT_TOO_MANY;
                            }else{
                                bc2hid = get_bardoce_bcindex(hitr);
                                bc2mm = get_barcode_mismatch(hitr);
                            }
                        }
                    }else{
                        if(maxbc2l == bc2l-1) hitr = wdindex[bc2int];
                        else hitr = wiindex[bc2int];
                        if(get_barcode_mismatch(hitr) == MISSIONBIO_BARCODE_INDEX_DEFAULT){
                            if(dropry < 0) dropry = BARCODE_DROP_BY_INDEL_TOO_MANY;
                        }else{
                            if(get_barcode_bccount(hitr) > 1){
                                if(dropry < 0) dropry = BARCODE_DROP_BY_INDEL_HIT_TOO_MANY;
                            }else{
                                bc2hid = get_bardoce_bcindex(hitr);
                                bc2mm = get_barcode_mismatch(hitr);
                                if(maxbc2l == bc2l-1) bc2del = bc2mm;
                                else bc2ins = bc2mm;
                                bc2mm = 0;
                            }
                        }
                    }
                    if(bc2hid >= 0){ // barcode2 found, then sliding window to find barcode1
                        uint32_t s2mval = 0UL, s2ival = 0UL, s2dval = 0UL;
                        for(int bi = 0; bi < maxbc1l; ++bi){
                            if(bi < bc1l-1){
                                s2dval <<= 2;
                                s2dval |= nuc_to_2bit[(int)r1->seq.s[bi]];
                            }
                            if(bi < bc1l){
                                s2mval <<= 2;
                                s2mval |= nuc_to_2bit[(int)r1->seq.s[bi]];
                            }
                            if(bi < bc1l+1){
                                s2ival <<= 2;
                                s2ival |= nuc_to_2bit[(int)r1->seq.s[bi]];
                            }
                        }
                        // try mm first
                        droprx = -1;
                        if(maxbc1l >= bc1l){
                            hitr = wmindex[s2mval];
                            if(get_barcode_mismatch(hitr) != MISSIONBIO_BARCODE_INDEX_DEFAULT){
                               if(get_barcode_bccount(hitr) == 1){
                                   bc1hid = get_bardoce_bcindex(hitr);
                                   bc1mm = get_barcode_mismatch(hitr);
                               }else{
                                   droprx = BARCODE_DROP_BY_NINEBP_HIT_TOO_MANY;
                               }
                            }else{
                                droprx = BARCODE_DROP_BY_NINEBP_MM_TOO_MANY;
                            }
                        }
                        // try ins then
                        if(bc1hid < 0 && maxbc1l > bc1l){
                            hitr = wiindex[s2ival];
                            if(get_barcode_mismatch(hitr) != MISSIONBIO_BARCODE_INDEX_DEFAULT){
                               if(get_barcode_bccount(hitr) == 1){
                                   bc1hid = get_bardoce_bcindex(hitr);
                                   bc1ins = get_barcode_mismatch(hitr);
                               }else{
                                   if(droprx < 0) droprx = BARCODE_DROP_BY_INDEL_HIT_TOO_MANY;
                               }
                            }else{
                                if(droprx < 0) droprx = BARCODE_DROP_BY_INDEL_TOO_MANY;
                            }
                        }
                        // try del then
                        if(bc1hid < 0){
                            hitr = wdindex[s2dval];
                            if(get_barcode_mismatch(hitr) != MISSIONBIO_BARCODE_INDEX_DEFAULT){
                                if(get_barcode_bccount(hitr) == 1){
                                    bc1hid = get_bardoce_bcindex(hitr);
                                    bc1del = get_barcode_mismatch(hitr);
                                }else{
                                    if(droprx < 0) droprx = BARCODE_DROP_BY_INDEL_HIT_TOO_MANY;
                                }
                            }else{
                                if(dropr < 0) droprx = BARCODE_DROP_BY_INDEL_TOO_MANY;
                            }
                        }
                        if(bc1hid < 0){
                            ++spr->r2drop[droprx];
                            if(outdrop){
                                ksprintf(ds1, "@%s %s\n%s\n+\n%s\n", 
                                              r1->name.s, 
                                              MISSIONBIO_BARCODE_DROP_REASON[droprx], 
                                              r1->seq.s, 
                                              r1->qual.s);
                                if(r2) ksprintf(ds2, "@%s\n%s\n+\n%s\n", r2->name.s, r2->seq.s, r2->qual.s);
                            }
                        }else{
                            // stat
                            if(bc1mm){
                                ++spr->b1mm[bc1mm];
                            }else if(bc1del){
                                ++spr->b1mm[BARCODE_DEL1_CNTARR_INDEX];
                            }else if(bc2ins){
                                ++spr->b1mm[BARCODE_INS1_CNTARR_INDEX];
                            }else{
                                ++spr->b1mm[BARCODE_MISMATCH0_CNTARR_INDEX];
                            }
                            if(bc2mm){
                                ++spr->b2mm[bc2mm];
                            }else if(bc2del){
                                ++spr->b2mm[BARCODE_DEL1_CNTARR_INDEX];
                            }else if(bc2ins){
                                ++spr->b2mm[BARCODE_INS1_CNTARR_INDEX];
                            }else{
                                ++spr->b2mm[BARCODE_MISMATCH0_CNTARR_INDEX];
                            }
                            ++spr->bbmm[bc1mm+bc2mm+bc1ins+bc2ins+bc1del+bc2del];
                            ++spr->offc[offn];
                            // out read1/2 to buff
                            ksprintf(s1, "@%s ", r1->name.s);
                            ksprintf(s1, "%d_%d", bc1hid, bc2hid);
                            ksprintf(s1, "\n%s\n+\n%s\n", s2p+cs2l, r1->qual.s+(s2p-r1->seq.s)+cs2l);
                            if(r2) ksprintf(s2, "@%s\n%s\n+\n%s\n", r2->name.s, r2->seq.s, r2->qual.s);
                            add = true;
                            spr->rcmap[bc1hid][bc2hid] +=1;
                        }
                    }else{
                        ++spr->r2drop[dropry];
                        if(outdrop){
                            ksprintf(ds1, "@%s %s\n%s\n+\n%s\n", 
                                          r1->name.s, 
                                          MISSIONBIO_BARCODE_DROP_REASON[dropry], 
                                          r1->seq.s, 
                                          r1->qual.s);
                            if(r2) ksprintf(ds2, "@%s\n%s\n+\n%s\n", r2->name.s, r2->seq.s, r2->qual.s);
                        }
                    }
                }
            }else{
                ++spr->r2drop[BARCODE_DROP_BY_CONSTSEQ2_NOT_FOUND];
                if(outdrop){
                    ksprintf(ds1, "@%s %s\n%s\n+\n%s\n", 
                                  r1->name.s, 
                                  MISSIONBIO_BARCODE_DROP_REASON[BARCODE_DROP_BY_CONSTSEQ2_NOT_FOUND], 
                                  r1->seq.s, 
                                  r1->qual.s);
                    if(r2) ksprintf(ds2, "@%s\n%s\n+\n%s\n", r2->name.s, r2->seq.s, r2->qual.s);
                }
            }
        }else{
            ++spr->r2drop[BARCODE_DROP_BY_CONSTSEQ1_NOT_FOUND];
            if(outdrop){
                ksprintf(ds1, "@%s %s\n%s\n+\n%s\n", 
                              r1->name.s, 
                              MISSIONBIO_BARCODE_DROP_REASON[BARCODE_DROP_BY_CONSTSEQ1_NOT_FOUND],
                              r1->seq.s, 
                              r1->qual.s);
                if(r2) ksprintf(ds2, "@%s\n%s\n+\n%s\n", r2->name.s, r2->seq.s, r2->qual.s);
            }
        }
        if(add) ++spr->gotr;
        else ++spr->dropr;
    }
    std::lock_guard<std::mutex> lg(lock);
    if(s1->l) wr1->write(s1->s, s1->l);
    if(s2->l) wr2->write(s2->s, s2->l);
    if(ds1->l) dw1->write(ds1->s, ds1->l);
    if(ds2->l) dw2->write(ds2->s, ds2->l);
}

void sbc_biom_t::raw2clean_rng(ReadPack** rp, int beg, int end, int idx){
    kstring_t* s1 = br1[idx];
    kstring_t* s2 = br2[idx];
    s1->l = s2->l = 0;
    spl_res_t* spr = sprs[idx];
    krec1_t* r1 = NULL;
    krec1_t* r2 = NULL;
    char* endptr = NULL;
    long int bc1hid = 0, bc2hid = 0;
    for(int i = beg; i < end; ++i){
        ++spr->totr;
        r1 = rp[0]->reads[i];
        bc1hid = strtol(r1->comment.s, &endptr, 10);
        bc2hid = strtol(endptr+1, &endptr, 10);
        if(sprs[0]->rcmap[bc1hid][bc2hid] > minbc){
            ++spr->gotr;
            ksprintf(s1, "@%s %d", r1->name.s, bcidx[bc1hid][bc2hid]);
            ksprintf(s1, "\n%s\n+\n%s\n", r1->seq.s, r1->qual.s);
            if(inr2){
                r2 = rp[1]->reads[i];
                ksprintf(s2, "@%s %d", r2->name.s, bcidx[bc1hid][bc2hid]);
                ksprintf(s2, "\n%s\n+\n%s\n", r2->seq.s, r2->qual.s);
            }
        }else{
            ++spr->dropr;
        }
    }
    std::lock_guard<std::mutex> lg(lock);
    if(s1->l) wr1->write(s1->s, s1->l);
    if(s2->l) wr2->write(s2->s, s2->l);
}

void sbc_biom_usage(sbc_biom_t* opt, char* arg0){
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: %s [options]\n\n", arg0);
    fprintf(stderr, "Options: -i FILE input FASTQ read1\n");
    fprintf(stderr, "         -I FILE input FASTQ read2\n");
    fprintf(stderr, "         -w FILE cellbarcode whitelist index file\n");
    fprintf(stderr, "         -W FILE cellbarcode whitelist sequence file\n");
    fprintf(stderr, "         -o FILE output directory [%s]\n", opt->outdir.c_str());
    fprintf(stderr, "         -m INT  high abundance cellbarcode threshold [%d]\n", opt->minbc);
    fprintf(stderr, "         -t INT  working thread number [%d]\n", opt->thread);
    fprintf(stderr, "         -P INT  max packs to store reads [%d]\n", opt->maxp);
    fprintf(stderr, "         -R INT  max reads in each pack [%d]\n", opt->maxr);
    fprintf(stderr, "         -d      output dropped reads if set\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "ps: cell barcode must be in read1 and in the flowing structure \n");
    fprintf(stderr, "    [CellBarcode1][Offset][ConstSeq1][CellBarcode2][ConstSeq2]\n");
    fprintf(stderr, "    CellBarcode1: %dbp\n", MISSIONBIO_BARCODE_LENGTH);
    fprintf(stderr, "    Offset      : <= %dbp\n", MISSIONBIO_BARCODE_CONST_MAX_OFFSET);
    fprintf(stderr, "    ConstSeq1   : AGTACGTACGAGTC [14bp]\n");
    fprintf(stderr, "    CellBarcode2: %dbp\n", MISSIONBIO_BARCODE_LENGTH);
    fprintf(stderr, "    ConstSeq2   : GTACTCGCAGTAGTC [15bp]\n");
    fprintf(stderr, "\n");
}

int sbc_biom_main(int argc, char** argv){
    sbc_biom_t opt;
    if(argc == 1){
        sbc_biom_usage(&opt, argv[0]);
        return 0;
    }
    int c = -1;
    while((c = getopt(argc, argv, "i:I:w:W:o:m:t:P:R:dh")) >= 0){
        switch(c){
            case 'i': opt.inr1 = optarg; break;
            case 'I': opt.inr2 = optarg; break;
            case 'w': opt.iwidx = optarg; break;
            case 'W': opt.wlist = optarg; break;
            case 'o': opt.outdir = optarg; break;
            case 'm': opt.minbc = atoi(optarg); break;
            case 't': opt.thread = atoi(optarg); break;
            case 'P': opt.maxp = atoi(optarg); break;
            case 'R': opt.maxr = atoi(optarg); break;
            case 'd': opt.outdrop = true; break;
            case 'h': sbc_biom_usage(&opt, argv[0]); return 0; break;
            default: break;
        }
    }
    if(opt.valid4spl()){
        opt.split_all();
        opt.raw2clean_all();
        return 0;
    }else{
        return 1;
    }
}

void sbc_bidx_usage(char* arg0){
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: %s [options]\n\n", arg0);
    fprintf(stderr, "Options: -i FILE cellbarcode whitelist sequence list file \n");
    fprintf(stderr, "         -o FILE cellbarcode sequence index file\n");
    fprintf(stderr, "\n");
}

int sbc_bidx_main(int argc, char** argv){
    if(argc == 1){
        sbc_bidx_usage(argv[0]);
        return 0;
    }
    sbc_biom_t opt;
    int c = -1;
    while((c = getopt(argc, argv, "i:o:h")) >= 0){
        switch(c){
            case 'o': opt.owidx = optarg; break;
            case 'i': opt.wlist = optarg; break;
            case 'h': sbc_bidx_usage(argv[0]); return 0; break;
            default: break;
        }
    }
    if(opt.valid4idx()){
        opt.idx2file();
        return 0;
    }else{
        return 1;
    }
}
