#include "hapcnt.h"

bool hap_opt_t::valid() {
    if(inbam.empty() || (!util::exists(inbam))) {
        fprintf(stderr, "valid input bam must be provided\n");
        return false;
    }
    if(ref.empty() || (!util::exists(ref))) {
        fprintf(stderr, "valid reference must be provided\n");
        return false;
    }
    if(alt.empty() || (!util::exists(alt))) {
        fprintf(stderr, "valid donor seq must be provided\n");
        return false;
    }
    if(beg < 0 || end < 0 || beg > end) {
        fprintf(stderr, "beg and end must be positive and beg < end\n");
        return false;
    }
    return true;
}

void hap_opt_t::bm2ex(bam1_t* b){
    int oi, ol, oplen_x = 0, oplen_m = 0;
    int32_t s = b->core.n_cigar+2;
    int32_t p = 0;
    uint32_t* nc = (uint32_t*)malloc(s*sizeof(uint32_t));
    uint32_t* oc = bam_get_cigar(b);
    int qpos = 0;
    int rpos = b->core.pos;
    int opt = 0;
    for(uint32_t i = 0; i < b->core.n_cigar; ++i){
        oi = bam_cigar_op(oc[i]);
        ol = bam_cigar_oplen(oc[i]);
        if(oi == BAM_CMATCH){
            for(int j = 0; j < ol; ++j){
                if(seq_nt16_int[bam_seqi(bam_get_seq(b), qpos+j)] != refints[rpos+j]){
                    if(oplen_m) nc = kswge_add_cigar(nc, &p, &s, oplen_m, '=');
                    oplen_m = 0;
                    ++oplen_x;
                }else{
                    if(oplen_x) nc = kswge_add_cigar(nc, &p, &s, oplen_x, 'X');
                    oplen_x = 0;
                    ++oplen_m;
                }
            }
            qpos += ol;
            rpos += ol;
        }else{
            if(oplen_m) nc = kswge_add_cigar(nc, &p, &s, oplen_m, '=');
            if(oplen_x) nc = kswge_add_cigar(nc, &p, &s, oplen_x, 'X');
            oplen_m = oplen_x = 0;
            nc = kswge_add_cigar(nc, &p, &s, ol, BAM_CIGAR_STR[oi]);
            opt = bam_cigar_type(oi);
            if(opt & 0x1) qpos += ol;
            if(opt & 0x2) rpos += ol;
        }
    }
    if(oplen_m) nc = kswge_add_cigar(nc, &p, &s, oplen_m, '=');
    if(oplen_x) nc = kswge_add_cigar(nc, &p, &s, oplen_x, 'X');
    int32_t new_datal = b->l_data + (p<<2)-(b->core.n_cigar<<2);
    uint8_t* ndata = (uint8_t*)malloc(new_datal*sizeof(uint8_t));
    memcpy(ndata, bam_get_qname(b), b->core.l_qname);
    memcpy(ndata+b->core.l_qname, (uint8_t*)nc, p<<2);
    memcpy(ndata+b->core.l_qname+(p<<2), bam_get_seq(b), b->l_data-(b->core.l_qname+(b->core.n_cigar<<2)));
    free(b->data);
    b->data = ndata;
    b->core.n_cigar = p;
    b->m_data = new_datal;
    b->l_data = new_datal;
    free(nc);
}

void hap_opt_t::update() {
    if(!util::exists(outdir)) util::makedir(outdir);
    outtsv = util::joinpath(outdir, outtsv);
    outhtml = util::joinpath(outdir, outhtml);
    outjsn = util::joinpath(outdir, outjsn);
    outntsv = util::joinpath(outdir, outntsv);
    outatsv = util::joinpath(outdir, outatsv);
    outstsv = util::joinpath(outdir, outstsv);
    hapb2tsv = util::joinpath(outdir, hapb2tsv);
    hapb1tsv = util::joinpath(outdir, hapb1tsv);
    outsnv = util::joinpath(outdir, outsnv);
    outaac = util::joinpath(outdir, outaac);
    outmfaac = util::joinpath(outdir, outmfaac);
    outmfnuc = util::joinpath(outdir, outmfnuc);
    outrcntnall2ref = util::joinpath(outdir, outrcntnall2ref);
    outrcntaall2ref = util::joinpath(outdir, outrcntaall2ref);
    outrcntnndef2ref = util::joinpath(outdir, outrcntnndef2ref);
    outrcntandef2ref = util::joinpath(outdir, outrcntandef2ref);
    outrcntnydef2ref = util::joinpath(outdir, outrcntnydef2ref);
    outrcntaydef2ref = util::joinpath(outdir, outrcntaydef2ref);
    outrcntnall2alt = util::joinpath(outdir, outrcntnall2alt);
    outrcntaall2alt = util::joinpath(outdir, outrcntaall2alt);
    outrcntnndef2alt = util::joinpath(outdir, outrcntnndef2alt);
    outrcntandef2alt = util::joinpath(outdir, outrcntandef2alt);
    outrcntnydef2alt = util::joinpath(outdir, outrcntnydef2alt);
    outrcntaydef2alt = util::joinpath(outdir, outrcntaydef2alt);
    --beg;
    umaxh = maxh;
    dropmask = (KSW_FHIERRLOWQ | KSW_FPRIMDIMER | KSW_FMISPRIMINN | KSW_FMISPRIMOUT | KSW_FLOWSURPT | KSW_FMANYVARS);
    if(l2s){
        rbms.resize(maxh+1, RecombBitMaskList());
        for(int i = 1; i <= maxh; ++i) {
            rbms[i].resize(i+1);
            for(int j = 1; j <= i; ++j) {
                util::combSet(i, j, rbms[i][j]);
            }
        }
    }else{
        maxh = umaxh = end - beg + 1;
        rbms.resize(maxh + 1, RecombBitMaskList());
        for(int i = 1; i <= maxh; ++i){
            rbms[i].resize(1);
            rbms[i][0].resize(1);
            rbms[i][0][0].reserve(i);
            for(int j = 0; j < i; ++j){
                rbms[i][0][0].push_back(j);
            }
        }
    }
    nthn = otn;
    if(gtn > nthn) nthn = gtn;
    subl = end-beg+1;
    sual = subl/3;
    foc2mnv();
    for(auto& e: outmsnv){
        if(e.size()) e = util::joinpath(outdir, e);
    }
}

void hap_opt_t::allocmem() {
    nuccnt = (int64_t**)malloc(5*sizeof(int64_t*));
    nucttt = (int64_t**)malloc(5*sizeof(int64_t*));
    bnuccnt = (int64_t**)malloc(37*sizeof(int64_t*));
    bnucttt = (int64_t**)malloc(37*sizeof(int64_t*));
    snuccnt = (int64_t**)malloc(5*sizeof(int64_t*));
    aaccnt = (int64_t**)malloc(ALL_AA_CNT*sizeof(int64_t*));
    for(int i = 0; i < 5; ++i) nuccnt[i] = (int64_t*)calloc(refl, sizeof(int64_t));
    for(int i = 0; i < 5; ++i) nucttt[i] = (int64_t*)calloc(refl, sizeof(int64_t));
    for(int i = 0; i < 5; ++i) snuccnt[i] = (int64_t*)calloc(5, sizeof(int64_t));
    for(int i = 0; i < ALL_AA_CNT; ++i) aaccnt[i] = (int64_t*)calloc(sual, sizeof(int64_t));
    for(int i = 0; i < 37; ++i) bnuccnt[i] = (int64_t*)calloc(37, sizeof(int64_t));
    for(int i = 0; i < 37; ++i) bnucttt[i] = (int64_t*)calloc(37, sizeof(int64_t));
    ttdepn = (int64_t*)calloc(refl, sizeof(int64_t));
    ttdepa = (int64_t*)calloc(sual, sizeof(int64_t));
    ttmutn = (int64_t*)calloc(refl, sizeof(int64_t));
    ttmuta = (int64_t*)calloc(sual, sizeof(int64_t));
    mutn2rcntall2ref = (int64_t*)calloc(refl+1, sizeof(int64_t));
    muta2rcntall2ref = (int64_t*)calloc(sual+1, sizeof(int64_t));
    mutn2rcntndef2ref = (int64_t*)calloc(refl+1, sizeof(int64_t));
    muta2rcntndef2ref = (int64_t*)calloc(sual+1, sizeof(int64_t));
    mutn2rcntydef2ref = (int64_t*)calloc(refl+1, sizeof(int64_t));
    muta2rcntydef2ref = (int64_t*)calloc(sual+1, sizeof(int64_t));
    mutn2rcntall2alt = (int64_t*)calloc(refl+1, sizeof(int64_t));
    muta2rcntall2alt = (int64_t*)calloc(sual+1, sizeof(int64_t));
    mutn2rcntndef2alt = (int64_t*)calloc(refl+1, sizeof(int64_t));
    muta2rcntndef2alt = (int64_t*)calloc(sual+1, sizeof(int64_t));
    mutn2rcntydef2alt = (int64_t*)calloc(refl+1, sizeof(int64_t));
    muta2rcntydef2alt = (int64_t*)calloc(sual+1, sizeof(int64_t));
    ks = (kstring_t*)calloc(1, sizeof(kstring_t));
    ms = (kstring_t*)calloc(1, sizeof(kstring_t));
    rs = (kstring_t*)calloc(1, sizeof(kstring_t));
}

void hap_opt_t::buildam(){
    int idx = 0;
    refais.clear();
    for(int i = beg; i+2 < end; i += 3){
        idx = codon2aminoid(far+i, ismt);
        refais.push_back(idx);
        refas.push_back(codon_1letter_names[idx][0]);
        idx = codon2aminoid(faa+i, ismt);
        altas.push_back(codon_1letter_names[idx][0]);
    }
    aas2type[refas] = AAS_TYPE_SYNONYM;
    aas2type[altas] = AAS_TYPE_ONLYALT;
    if(refas == altas){
        aas2type[refas] |= AAS_TYPE_ONLYALT;
        aas2type[altas] |= AAS_TYPE_SYNONYM;
    }
    allaltc = 0;
    for(size_t i = 0; i < altas.size(); ++i){
        if(altas[i] != refas[i]){
            ++allaltc;
            xaas.push_back(altas[i]);
            xras.push_back(refas[i]);
            xapos.push_back(i);
            xamark.push_back(altas[i]);
        }else{
            xamark.push_back('X');
        }
    }
}

void hap_opt_t::annoaas(){
    int idx = 0;
    for(auto& e: brl){
        e->aas.clear();
        for(size_t i = 0, xdx = 0; i+2 < e->seq.size(); i += 3, ++xdx){
            idx = codon2aminoid(e->seq.c_str()+i, ismt);
            aaccnt[idx][xdx] += e->cnt;
            e->aas.push_back(codon_1letter_names[idx][0]);
        }
        e->ma2r = e->ma2a = 0;
        for(size_t j = 0; j < e->aas.size(); ++j){
            if(e->aas[j] != refas[j]) ++e->ma2r;
            if(e->aas[j] != altas[j]) ++e->ma2a;
        }
        auto iter = aas2type.find(e->aas);
        if(iter != aas2type.end()){
            e->aat |= iter->second;
            if(!(e->aat & AAS_TYPE_SYNONYM)){
                e->aat |= AAS_TYPE_NONSYNO;
                e->oriaa = xras;
                e->mutaa = xaas;
                e->apos = xapos;
                xaonlym += e->cnt;
            }else{
                e->oriaa = xras;
                e->mutaa = xras;
                e->apos = xapos;
                xasyn += e->cnt;
            }
        }else{
            e->aat |= AAS_TYPE_NONSYNO;
            int xvc = 0;
            for(size_t j = 0; j < e->aas.size(); ++j){
                if(e->aas[j] == altas[j] && e->aas[j] != refas[j]) ++xvc;
                if(e->aas[j] != refas[j]){
                    e->apos.push_back(j);
                    e->mutaa.push_back(e->aas[j]);
                    e->oriaa.push_back(refas[j]);
                }
            }
            if(xvc == allaltc){
                e->aat |= AAS_TYPE_WITHALT;
                xawithm += e->cnt;
            }else{
                xaother += e->cnt;
            }
        }
    }
}

void hap_opt_t::hapcnt1b(bam1_t* b){
    int32_t cc = 1;
    int64_t vs = 0;
    if(!nget){
        cc = bam_aux2i(bam_aux_get(b, "CC"));
        vs = bam_aux2i(bam_aux_get(b, "VS"));
        if(!(vs & KSW_FREPSEQR)) return;
        if(vs & dropmask) return;
    }else{
        uint8_t* nmd = bam_aux_get(b, "NM");
        if(nmd && bam_aux2i(nmd) == 0) vs = KSW_FREFTYPE;
        else vs = KSW_FMAYVARSEQ;
        bm2ex(b);
    }
    ttr += cc;
    if(vs & KSW_FREFTYPE){
        trr += cc;
        for(int32_t p = b->core.pos; p < bam_endpos(b); ++p) nuccnt[nuc_to_3bit[(int)far[p]]][p] += cc;
        for(int32_t p = 0; p < sual; ++p) aaccnt[refais[p]][p] += cc;
        mutn2rcntall2ref[0] += cc;
        muta2rcntall2ref[0] += cc;
        mutn2rcntndef2ref[0] += cc;
        muta2rcntndef2ref[0] += cc;
        mutn2rcntall2alt[mc8nuc] += cc;
        muta2rcntall2alt[mc8aac] += cc;
        mutn2rcntndef2alt[mc8nuc] += cc;
        muta2rcntndef2alt[mc8aac] += cc;
        return;
    }
    tar += cc;
    if(vs & KSW_FMAYVARSEQ){ tmr += cc; }
    uint32_t* cigar = NULL;
    int oplen, opint, idl = 0;
    int rpos = 0, qpos = 0;
    int qbeg = 0, qend = 0;
    bpos.clear();
    rpos = b->core.pos;
    qpos = qbeg = 0;
    cigar = bam_get_cigar(b);
    char bchr;
    int withdmutn = 0;
    int n3bi = 0;
    uint8_t b2bi = 0;
    for(int i = 0; i < 5; ++i) memset(nucttt[i], 0, refl*sizeof(int64_t));
    for(int i = 0; i < 37; ++i) memset(bnucttt[i], 0, 37*sizeof(int64_t));
    for(uint32_t i = 0; i < b->core.n_cigar; ++i) {
        opint = bam_cigar_op(cigar[i]);
        oplen = bam_cigar_oplen(cigar[i]);
        switch(opint) {
            case BAM_CDEL:
                if(rpos >= beg || rpos < end) idl += oplen;
                rpos += oplen;
                break;
            case BAM_CINS:
                if(rpos < beg) qbeg += oplen;
                if(rpos < end-1) qend += oplen;
                qpos += oplen;
                if(rpos >= beg || rpos < end) idl += oplen;
                break;
            case BAM_CSOFT_CLIP:
                if(rpos < beg) qbeg += oplen;
                if(rpos < end-1) qend += oplen;
                if(rpos > beg || rpos < end) idl += oplen;
                qpos += oplen;
                break;
            case BAM_CEQUAL:
                if(rpos < beg) qbeg += MIN(beg-rpos, oplen);
                if(rpos < end-1) qend += MIN(end-1-rpos, oplen);
                for(int p = 0; p < oplen; ++p){
                    n3bi = nuc_to_3bit[(int)far[rpos+p]];
                    nuccnt[n3bi][rpos+p] += cc;
                    nucttt[n3bi][rpos+p] += cc;
                }
                rpos += oplen;
                qpos += oplen;
                break;
            case BAM_CDIFF:
                b2bi = 0;
                if(rpos >= beg && rpos < end) {
                    for(int ol = 0; ol < oplen; ++ol) {
                        bchr = seq_nt16_str[bam_seqi(bam_get_seq(b), qpos + ol)];
                        if(oplen == 2){
                            b2bi <<= 3;
                            b2bi |= nuc_to_3bit[(int)bchr];
                        }
                        if(bchr == faa[rpos+ol]) ++withdmutn;
                        bseq[rpos+ol] = bchr;
                        bpos.push_back(rpos+ol);
                        n3bi = seq_nt16_int[bam_seqi(bam_get_seq(b), qpos + ol)];
                        nuccnt[n3bi][rpos+ol] += cc;
                        nucttt[n3bi][rpos+ol] += cc;
                    }
                }else{
                    for(int ol = 0; ol < oplen; ++ol) {
                        n3bi = seq_nt16_int[bam_seqi(bam_get_seq(b), qpos + ol)];
                        nuccnt[n3bi][rpos+ol] += cc;
                        if(oplen == 2){
                            bchr = seq_nt16_str[bam_seqi(bam_get_seq(b), qpos + ol)];
                            b2bi <<= 3;
                            b2bi |= nuc_to_3bit[(int)bchr];
                        }
                        nucttt[n3bi][rpos+ol] += cc;
                    }
                }
                if(oplen == 2){
                    b2bi &= 0x3F;
                    bnuccnt[bnucref[rpos]][b2bi] += cc;
                    bnucttt[bnucref[rpos]][b2bi] += cc;
                }
                if(rpos < beg) qbeg += MIN(beg-rpos, oplen);
                if(rpos < end-1) qend += MIN(end-1-rpos, oplen);
                rpos += oplen;
                qpos += oplen;
                break;
            default:
                break;
        }
    }
    if(idl){
        tfr += cc;
        for(int i = 0; i < 5; ++i){
            for(int j = 0; j < refl; ++j){
                nuccnt[i][j] -=  nucttt[i][j];
            }
        }
        for(int i = 0; i < 37; ++i){
            for(int j = 0; j < 37; ++j){
                bnuccnt[i][j] -= bnucttt[i][j];
            }
        }
    }else if(bpos.empty()){
        for(int32_t p = 0; p < sual; ++p) aaccnt[refais[p]][p] += cc;
        mutn2rcntall2ref[0] += cc;
        muta2rcntall2ref[0] += cc;
        mutn2rcntndef2ref[0] += cc;
        muta2rcntndef2ref[0] += cc;
        mutn2rcntall2alt[mc8nuc] += cc;
        muta2rcntall2alt[mc8aac] += cc;
        mutn2rcntndef2alt[mc8nuc] += cc;
        muta2rcntndef2alt[mc8aac] += cc;
        tnr += cc;
    }else if(bpos.size() > umaxh) txr += cc;
    if(idl == 0 &&  bpos.size() && bpos.size() <= umaxh) {
        for(auto& e: rbms[bpos.size()]) {
            for(auto& f: e) {
                ks->l = 0;
                ms->l = 0;
                rs->l = 0;
                tpos.clear();
                for(auto& pi: f) {
                    ksprintf(ks, "%d%c", bpos[pi], bseq[bpos[pi]]);
                    ksprintf(ms, "%c", bseq[bpos[pi]]);
                    ksprintf(rs, "%c", far[bpos[pi]]);
                    tpos.push_back(bpos[pi]);
                }
                if(ks->l) {
                    if(tpos.size() == bpos.size()){ // full link snp
                        bam_rec_t* br = new bam_rec_t();
                        br->aat = withdmutn ? NUC_WITHALT : 0;
                        br->cnt = cc;
                        br->mn2a = bpos.size()-withdmutn;
                        br->mn2r = bpos.size();
                        br->muts = ms->s;
                        br->refs = rs->s;
                        br->pos = tpos;
                        br->seq.resize(qend-qbeg+1, '\0');
                        int ttbl = 0;
                        for(int bbl = qbeg; bbl <= qend; ++bbl){
                            br->seq[ttbl++] = seq_nt16_str[bam_seqi(bam_get_seq(b), bbl)];
                        }
                        brl.push_back(br);
                    }
                    std::string key = ks->s;
                    auto iter = hapm.find(key);
                    if(iter == hapm.end()) {
                        hap_rec_t* ht = new hap_rec_t();
                        ht->cnt = cc;
                        ht->muts = ms->s;
                        ht->refs = rs->s;
                        ht->pos = tpos;
                        hapm[key] = ht;
                    } else {
                        iter->second->cnt += cc;
                    }
                }
            }
        }
    }
}

void hap_opt_t::hapcnt8bf(){
    // compute hapcnt
    util::loginfo(stderr, "beg parse hap from bam");
    samFile* fp = sam_open(inbam.c_str(), "r");
    bam_hdr_t* hdr = sam_hdr_read(fp);
    bam1_t* b = bam_init1();
    hts_idx_t* idx = sam_index_load(fp, inbam.c_str());
    tid = sam_hdr_name2tid(hdr, contig.c_str());
    refl = hdr->target_len[tid];
    int32_t msize = refl * sizeof(char);
    bseq = (char*)malloc(msize);
    bpos.resize(refl, 0);
    tpos.resize(refl, 0);
    allocmem();
    fai = fai_load(ref.c_str());
    int fal = 0;
    far = faidx_fetch_seq(fai, contig.c_str(), 0, refl, &fal);
    fai_destroy(fai);
    refints = kswge_seq2ints(far, refl);
    fai = fai_load(alt.c_str());
    faa = faidx_fetch_seq(fai, contig.c_str(), 0, refl, &fal);
    fai_destroy(fai);
    fai = NULL;
    ref2altp();
    buildam();
    // update bnucref
    ref2bnuc();
    // hapm to hapl
    mc8nuc = 0;
    for(auto& e: aps){ if(e != 4) ++mc8nuc;};
    mc8aac = allaltc;
    hts_itr_t* itr = sam_itr_queryi(idx, tid, beg, end);
    while(sam_itr_next(fp, itr, b) >= 0){ hapcnt1b(b); }
    // relaese mem
    sam_close(fp);
    bam_hdr_destroy(hdr);
    bam_destroy1(b);
    hts_idx_destroy(idx);
    hts_itr_destroy(itr);
    util::loginfo(stderr, "end parse hap from bam");
    // update tcr
    tcr = ttr - tfr -txr;
    util::loginfo(stderr, "beg sort hap by spr");
    hapl.resize(hapm.size(), NULL);
    int64_t vid = 0;
    for(auto& e: hapm) hapl[vid++] = e.second;
    std::sort(hapl.begin(), hapl.end(), hap_rec_sorter());
    util::loginfo(stderr, "end sort hap by spr");
}

void hap_opt_t::hapcnt8bv(std::vector<bam1_t*>& bv, int blen){
    // compute hapcnt
    int32_t msize = blen * sizeof(char);
    refl = blen;
    refints = kswge_seq2ints(far, refl);
    bseq = (char*)malloc(msize);
    bpos.resize(refl, 0);
    tpos.resize(refl, 0);
    allocmem();
    buildam();
    // update bnucref
    ref2bnuc();
    // hapm to hapl
    mc8nuc = 0;
    for(auto& e: aps){ if(e != 4) ++mc8nuc;};
    mc8aac = allaltc;
    for(auto& b: bv) hapcnt1b(b);
    // update tcr
    tcr = ttr - tfr -txr;
    hapl.resize(hapm.size(), NULL);
    int64_t vid = 0;
    for(auto& e: hapm) hapl[vid++] = e.second;
    std::sort(hapl.begin(), hapl.end(), hap_rec_sorter());
}

void hap_opt_t::hap2tsv() {
    // haptype
    std::ofstream fw(outtsv);
    hap_rec_t::hap2head(fw);
    int64_t tcm = 0;
    for(auto& h: hapl) tcm += h->cnt;
    if(otn > 0){
        for(int i = 0; i < otn && i < (int)hapl.size(); ++i){
            if(hapl[i]->cnt > minr) hapl[i]->hap2rec(fw, ttr, tar, tcr, tcm);
        }
    }else{
        for(auto& hap: hapl){
            if(hap->cnt > minr) hap->hap2rec(fw, ttr, tar, tcr, tcm);
        }
    }
    fw.close();
}

void hap_opt_t::bam2tsv() {
    merge_and_sort_bam_rec8seq(brl);
    // summary nuc mut cnt of each seq
    for(auto& bm: brl){
        if(bm->aat & NUC_WITHALT){
            mutn2rcntydef2ref[bm->mn2r] += bm->cnt;
            mutn2rcntydef2alt[bm->mn2a] += bm->cnt;
        }else{
            mutn2rcntndef2ref[bm->mn2r] += bm->cnt;
            mutn2rcntndef2alt[bm->mn2a] += bm->cnt;
        }
        mutn2rcntall2ref[bm->mn2r] += bm->cnt;
        mutn2rcntall2alt[bm->mn2a] += bm->cnt;
    }
    // output all rc dist
    int64_t rydef = 0, rndef = 0;
    for(int i = 0; i <= subl; ++i){ rydef += mutn2rcntydef2ref[i]; rndef += mutn2rcntndef2ref[i];}
    mc2tsv(mutn2rcntydef2ref, subl, rydef, outrcntnydef2ref);
    mc2tsv(mutn2rcntydef2alt, subl, rydef, outrcntnydef2alt);
    mc2tsv(mutn2rcntndef2ref, subl, rndef, outrcntnndef2ref);
    mc2tsv(mutn2rcntndef2alt, subl, rndef, outrcntnndef2alt);
    mc2tsv(mutn2rcntall2ref, subl, tcr, outrcntnall2ref);
    mc2tsv(mutn2rcntall2alt, subl, tcr, outrcntnall2alt);
    // output seq cnt tsv
    int64_t tcm = 0;
    for(auto& b: brl) tcm += b->cnt;
    std::ofstream fw(outntsv);
    bam_rec_t::bam2head(fw);
    if(otn > 0){
        for(int i = 0; i < otn && i < (int)brl.size(); ++i){
            if(brl[i]->cnt > minr) brl[i]->bam2rec(fw, ttr, tar, tcr, tcm);
        }
    }else{
        for(auto& bam: brl){
            if(bam->cnt > minr) bam->bam2rec(fw, ttr, tar, tcr, tcm);
        }
    }
    fw.close();
}

void hap_opt_t::aas2tsv() {
    annoaas();
    merge_and_sort_bam_rec8aas(brl);
    // summary aa mut cnt of each seq
    int withrefaa = 0;
    int64_t tcm = 0;
    for(auto& aas: brl){
        if(aas->aat & (AAS_TYPE_ONLYALT | AAS_TYPE_WITHALT)){
            muta2rcntydef2ref[aas->ma2r] += aas->cnt;
            muta2rcntydef2alt[aas->ma2a] += aas->cnt;
        }else{
            muta2rcntndef2ref[aas->ma2r] += aas->cnt;
            muta2rcntndef2alt[aas->ma2a] += aas->cnt;
        }
        muta2rcntall2ref[aas->ma2r] += aas->cnt;
        muta2rcntall2alt[aas->ma2a] += aas->cnt;
        if(aas->aat & AAS_TYPE_SYNONYM) withrefaa = 1;
        else tcm += aas->cnt;
    }
    tat = brl.size();
    if((!withrefaa) && trr) ++tat;
    // out all ac dist
    int64_t rydef = 0, rndef = 0;
    for(int i = 0; i <= sual; ++i){ rydef += muta2rcntydef2ref[i]; rndef += muta2rcntndef2ref[i]; }
    mc2tsv(muta2rcntydef2ref, sual, rydef, outrcntaydef2ref);
    mc2tsv(muta2rcntydef2alt, sual, rydef, outrcntaydef2alt);
    mc2tsv(muta2rcntndef2ref, sual, rndef, outrcntandef2ref);
    mc2tsv(muta2rcntndef2alt, sual, rndef, outrcntandef2alt);
    mc2tsv(muta2rcntall2ref, sual, tcr, outrcntaall2ref);
    mc2tsv(muta2rcntall2alt, sual, tcr, outrcntaall2alt);
    // output aaseq tsv
    std::ofstream fw(outatsv);
    bam_rec_t::aas2head(fw);
    if(otn > 0){
        for(int i = 0; i < otn && i < (int)brl.size(); ++i){
            if(brl[i]->cnt > minr) brl[i]->aas2rec(fw, ttr, tar, tcr, tcm);
        }
    }else{
        for(auto& aas: brl){
            if(aas->cnt > minr) aas->aas2rec(fw, ttr, tar, tcr, tcm);
        }
    }
    fw.close();
}

void hap_opt_t::mc2tsv(int64_t* vec, int32_t len, int64_t sum, const std::string& outf){
    std::ofstream fw(outf);
    // row1
    for(int32_t i = 0; i < len; ++i) fw << i << "\t";
    fw << len << "\n";
    // row2
    for(int32_t i = 0; i < len; ++i){
        if(sum) fw << (double)vec[i]/(double)sum << "\t";
        else fw << .0 << "\t";
    }
    if(sum) fw << (double)vec[len]/(double)sum << "\n";
    else fw << .0 << "\n";
    // row3
    for(int32_t i = 0; i < len; ++i) fw << vec[i] << "\t";
    fw << vec[len] << "\n";
    fw.close();
}

void hap_opt_t::bias2tsv(){
    // summary 1base bias
    for(int i = 0; i < refl; ++i){
        for(uint8_t b = 0; b <= 4; ++b){
            if(refints[i] != b){
                snuccnt[refints[i]][b] += nuccnt[b][i];
            }
        }
    }
    // output 1base bias
    std::ofstream fw(hapb1tsv);
    int64_t ttbc = 0;
    fw << "Ref\\Alt";
    for(uint8_t b = 0; b <= 4; ++b) fw << "\t" << bit3_to_nuc[b];
    fw << "\n";
    for(uint8_t b = 0; b <= 4; ++b){
        fw << bit3_to_nuc[b];
        ttbc = 0;
        for(uint8_t k = 0; k <= 4; ++k){
            if(b != k) ttbc += snuccnt[b][k];
        }
        for(uint8_t k = 0; k <= 4; ++k){
            if(ttbc) fw << "\t" << (double)snuccnt[b][k]/(double)ttbc;
            else fw << "\t.0";
        }
        fw << "\n";
    }
    fw.close();
    // output 2base bias
    fw.open(hapb2tsv);
    fw << "Ref\\Alt";
    for(uint8_t b = 0; b <= 4; ++b){
        for(uint8_t k = 0; k <= 4; ++k) fw << "\t" << bit3_to_nuc[b] << bit3_to_nuc[k];
    }
    fw << "\n";
    uint8_t ridb = 0, aidb = 0;
    for(uint8_t b = 0; b <= 4; ++b){
        for(uint8_t k = 0; k <= 4; ++k){
            ttbc = 0;
            fw << bit3_to_nuc[b] << bit3_to_nuc[k];
            ridb = (b << 3) | k;
            for(uint8_t j = 0; j < 37; ++j){
                ttbc += bnuccnt[ridb][j];
            }
            for(uint8_t xb = 0; xb <= 4; ++xb){
                for(uint8_t xk = 0; xk <= 4; ++xk){
                    aidb = (xb << 3) | xk;
                    if(ttbc) fw << "\t" << (double)bnuccnt[ridb][aidb]/(double)ttbc;
                    else fw << "\t.0";
                }
            }
            fw << "\n";
        }
    }
    fw.close();
}

void hap_opt_t::hap2jsn(){
    std::ofstream fw(outjsn);
    fw << "{\n";
    fw << " \"TotalReads\":" << ttr << ",\n";
    fw << " \"RefReads\":" << trr << ",\n";
    fw << " \"NonRefReads\":" << tar << ",\n";
    fw << " \"ReadsWithDesiredNucMut\":" << tmr << ",\n";
    fw << " \"ReadsWithInDelOrClipInRng\":" << tfr << ",\n";
    fw << " \"ReadsWithTooManySNPInRng\":" << txr << ",\n";
    fw << " \"TotalHaplotype8Nuc\":";
    if(trr) fw << hapm.size() + 1 << ",\n";
    else fw << hapm.size() << ",\n";
    fw << " \"TotalHaplotype8Aac\":" << tat << ",\n";
    fw << " \"TotalHaplotypeSNV\":" << hapm.size() << ",\n";
    fw << " \"ReadsValid4HapAna\":" << tcr << ",\n";
    fw << " \"ReadsWithoutSNPInRng\":" << tnr << ",\n";
    fw << " \"ReadsWithSynonymAAMut\":" << xasyn << ",\n";
    fw << " \"ReadsWithOnlyDesiredAAMut\":" << xaonlym << ",\n";
    fw << " \"ReadsWithDesiredAAMut\":" << xawithm << ",\n";
    fw << " \"ReadsWithOtherAAMut\":" << xaother << "\n";
    fw << "}\n";
    fw.close();
}

void hap_opt_t::stat2tsv(){
    std::ofstream fw(outstsv);
    fw << "TotalReads\t";
    fw << "RefReads\t";
    fw << "NonRefReads\t";
    fw << "ReadsWithDesiredNucMut\t";
    fw << "ReadsWithInDelOrClipInRng\t";
    fw << "ReadsWithTooManySNPInRng\t";
    fw << "TotalHaplotype8Nuc\t";
    fw << "TotalHaplotype8Aac\t";
    fw << "TotalHaplotypeSNV\t";
    fw << "ReadsValid4HapAna\t";
    fw << "ReadsWithoutSNPInRng\t";
    fw << "ReadsWithSynonymAAMut\t";
    fw << "ReadsWithOnlyDesiredAAMut\t";
    fw << "ReadsWithDesiredAAMut\t";
    fw << "ReadsWithOtherAAMut\n";

    fw << ttr << "\t";
    fw << trr << "\t";
    fw << tar << "\t"; 
    fw << tmr << "\t";
    fw << tfr << "\t";
    fw << txr << "\t";
    if(trr) fw << hapm.size() + 1 << "\t";
    else fw << hapm.size() << "\t";
    fw << tat << "\t";
    fw << hapm.size() << "\t";
    fw << tcr << "\t";
    fw << tnr << "\t";
    fw << xasyn << "\t";
    fw << xaonlym << "\t";
    fw << xawithm << "\t";
    fw << xaother << "\n";

    fw.close();
}

void hap_opt_t::naa2dep(){
    if(!cntdep){
        for(int i = 0; i < refl; ++i){
            for(uint8_t b = 0; b <= 4; ++b){
                ttdepn[i] += nuccnt[b][i];
                if(refints[i] != b) ttmutn[i] += nuccnt[b][i];
            }
        }
        for(int i = 0; i < sual; ++i){
            for(int j = 0; j < ALL_AA_CNT; ++j){
                ttdepa[i] += aaccnt[j][i];
                if(refas[i] != codon_1letter_names[j][0]) ttmuta[i] += aaccnt[j][i];
            }
        }
        cntdep = 1;
    }
}

void hap_opt_t::snv2tsv(){
    if(!cntdep) naa2dep();
    std::ofstream fw(outsnv);
    fw << "Pos\tRef\tIsOn\tA\tC\tG\tT\tN\tAllDep\tMutDep\n";
    for(int i = 0; i < refl; ++i){
        fw << i+1 << "\t";
        fw << far[i] << "\t";
        fw << kswge_int2nt[aps[i]] << "\t";
        fw << nuccnt[nuc_to_3bit[(int)'A']][i] << "\t";
        fw << nuccnt[nuc_to_3bit[(int)'C']][i] << "\t";
        fw << nuccnt[nuc_to_3bit[(int)'G']][i] << "\t";
        fw << nuccnt[nuc_to_3bit[(int)'T']][i] << "\t";
        fw << nuccnt[nuc_to_3bit[(int)'N']][i] << "\t";
        fw << ttdepn[i] << "\t";
        fw << ttmutn[i] << "\n";
    }
    fw.close();
    // output nuc mut freq in range
    fw.open(outmfnuc);
    for(int i = beg; i < end-1; ++i){
        if(aps[i] == 4){
            fw << far[i] << "\t";
        }else{
            fw << far[i] << "(On)" << "\t";
        }
    }
    if(aps[end-1] == 4){
        fw << far[end-1] << "\n";
    }else{
        fw << far[end-1] << "(On)" << "\n";
    }
    for(int i = beg; i < end-1; ++i){
        if(ttmutn[i]) fw << (double)ttmutn[i]/(double)ttdepn[i] << "\t";
        else fw << .0 << "\t";
    }
    if(ttmutn[end-1]) fw << (double)ttmutn[end-1]/(double)ttdepn[end-1] << "\n";
    else fw << .0 << "\n";
    fw.close();
}

void hap_opt_t::aac2tsv(){
    if(!cntdep) naa2dep();
    std::ofstream fw(outaac);
    fw << "Pos\tRef\tIsOn\tAllDep\tMutDep";
    for(int i = 0; i < ALL_AA_CNT; ++i) fw << "\t" << codon_1letter_names[i];
    fw << "\n";
    for(int i = 0; i < sual; ++i){
        fw << i+1 << "\t";
        fw << refas[i] << "\t";
        fw << xamark[i] << "\t";
        fw << ttdepa[i] << "\t";
        fw << ttmuta[i];
        for(int j = 0; j < ALL_AA_CNT; ++j) fw << "\t" << aaccnt[j][i];
        fw << "\n";
    }
    fw.close();
    // output aac mut freq in range
    fw.open(outmfaac);
    for(int i = 0; i < sual-1; ++i){
        if(xamark[i] == 'X') fw << refas[i] << "\t";
        else fw << refas[i] <<  "(On)" << "\t";
    }
    if(xamark[sual-1] == 'X') fw << refas[sual-1] << "\n";
    else fw << refas[sual-1] << "(On)" << "\n";
    for(int i = 0; i < sual-1; ++i){
        if(ttmuta[i]) fw << (double)ttmuta[i]/(double)ttdepa[i] << "\t";
        else fw << .0 << "\t";
    }
    if(ttmuta[sual-1]) fw << (double)ttmuta[sual-1]/(double)ttdepa[sual-1] << "\n";
    else fw << .0 << "\n";
    fw.close();
}

void hap_opt_t::ref2bnuc(){
    bnucref = (uint8_t*)calloc(refl, sizeof(uint8_t));
    uint8_t txbi = nuc_to_3bit[(int)far[0]];
    for(int i = 1; i < refl; ++i){
        txbi <<= 3;
        txbi |= nuc_to_3bit[(int)far[i]];
        bnucref[i-1] = txbi & 0x3F;
    }
}

void hap_opt_t::ref2altp(){
    int match = 5, mismatch = 4, gapopen = 25, gapext = 0;
    int rlen = strlen(far);
    int qlen = strlen(faa);
    if(qlen != rlen){
        fprintf(stderr, "Error: ref[%d] and mut[%d] is not the same length!\n", rlen, qlen);
        exit(EXIT_FAILURE);
    }
    int8_t* score_mat = kswge_gen_smat(match, mismatch);
    uint8_t* rsfints = kswge_seq2ints(far, rlen);
    uint8_t* donorintsf = kswge_seq2ints(faa, qlen);
    char* donorrev = util::revComp2NewSeq((char*)faa, qlen);
    uint8_t* donorintsr = kswge_seq2ints(donorrev, qlen);
    kswr_t *retf = kswge_semi_global(qlen, donorintsf, rlen, rsfints, 5, score_mat, gapopen, gapext);
    kswr_t *retr = kswge_semi_global(qlen, donorintsr, rlen, rsfints, 5, score_mat, gapopen, gapext);
    kswr_t *bestaln = retr;
    uint8_t* bestints = donorintsr;
    if(retf->score > retr->score){
        bestaln = retf;
        bestints = donorintsf;
    }
    kswge_mark_mismatch(bestaln, rsfints, rlen, bestints, qlen);
    aps.resize(rlen+1, 4);
    char opchr;
    int oplen;
    int rpos = bestaln->tb;
    int qpos = bestaln->qb;
    for(int c = 0; c < bestaln->ncigar; ++c){
        opchr = kswge_cigar_opchr(bestaln->cigar[c]);
        oplen = kswge_cigar_oplen(bestaln->cigar[c]);
        if(opchr == 'I'){
            fprintf(stderr, "Error, ins found in alt aginist ref\n");
            exit(EXIT_FAILURE);
            qpos += oplen;
        }else if(opchr == 'M' || opchr == 'X' || opchr == '='){
            if(opchr == 'X'){
                for(int i = 1; i < oplen; ++i) aps[rpos] = bestints[qpos];
            }
            qpos += oplen;
            rpos += oplen;
        }else if(opchr == 'D'){
            fprintf(stderr, "Error, del found in alt aginist ref\n");
            exit(EXIT_FAILURE);
            rpos += oplen;
        }else if(opchr == 'S'){
            fprintf(stderr, "Error, clip found in alt aginist ref\n");
            exit(EXIT_FAILURE);
            qpos += oplen;
        }
    }
    // free up
    free(score_mat); score_mat = NULL;
    free(rsfints); rsfints = NULL;
    free(donorintsf); donorintsf = NULL;
    free(donorrev); donorrev = NULL;
    free(donorintsr); donorintsr = NULL;
    kswr_destroy(retf); retf = NULL;
    kswr_destroy(retr); retr = NULL;
}

void hap_opt_t::foc2mnv(){
    mnv.resize(16, 0);
    mmv.resize(5, 0);
    outmsnv.resize(16, "");
    outmstr.resize(16, NULL);
    if(extsnv.empty()) return;
    std::vector<std::string> vstr;
    util::split(extsnv, vstr, ",");
    std::vector<std::string> mts;
    uint8_t mval = 0, oval = 0, sval = 0;
    for(auto& p: vstr){
        if(p.empty()) continue;
        util::split(p, mts, ":");
        if(mts.size() < 2){
            fprintf(stderr, "Error format[%s] of muts, must be \"REF:ALT\"\n", p.c_str());
            exit(EXIT_FAILURE);
        }
        oval = nuc_to_3bit[(int)mts[0][0]];
        sval = nuc_to_3bit[(int)mts[1][0]];
        if(oval == 4){
            fprintf(stderr, "Error ref nucleotide provided[%c], must be A/C/G/T\n", mts[0][0]);
            exit(EXIT_FAILURE);
        }
        if(sval == 4){
            fprintf(stderr, "Error alt nucleotide provided[%c], must be A/C/G/T\n", mts[1][0]);
            exit(EXIT_FAILURE);
        }
        mval = (oval << 2) | sval;
        mmv[oval] = 1;
        mnv[mval] = 1;
        outmsnv[mval] = mts[0] + mts[1] + ".tsv";
        if(outmstr[mval] == NULL){
            outmstr[mval] = (kstring_t*)calloc(1, sizeof(kstring_t));
        }
    }
    // add target
    for(size_t i = 0; i < aps.size(); ++i){
        if(aps[i] == 4) continue;
        oval = nuc_to_3bit[(int)far[i]];
        sval = aps[i];
        mval = (oval << 2) | sval;
        mmv[oval] = 1;
        mnv[mval] = 1;
        outmsnv[mval] = std::string(1,far[i]) + std::string(1, kswge_int2nt[aps[i]])  + ".tsv";
        if(outmstr[mval] == NULL){
            outmstr[mval] = (kstring_t*)calloc(1, sizeof(kstring_t));
        }
    }
    // marker of acgtn
    outpstr.resize(4, NULL);
    for(int i = 0; i < 4; ++i) outpstr[i] = (kstring_t*)calloc(1, sizeof(kstring_t));
}

void hap_opt_t::snv2foc(){
    if(!cntdep) naa2dep();
    uint8_t oval = 0, sval = 0, mval = 0;
    for(int i = 0; i < refl; ++i){
        oval = nuc_to_3bit[(int)far[i]];
        if(!mmv[oval]) continue;
        if(aps[i] == 4){
            ksprintf(outpstr[oval], "%d\t", i+1);
        }else{
            ksprintf(outpstr[oval], "%d(On)\t", i+1);
        }
        for(uint8_t j = 0; j < 4; ++j){
            sval = j;
            mval = (oval << 2) | j;
            if(mnv[mval]){
                ksprintf(outmstr[mval], "%lf\t", ttdepn[i] ? (double)nuccnt[sval][i]/(double)ttdepn[i] : .0);
            }
        }
    }
    for(size_t i = 0; i < outpstr.size(); ++i){
        if(outpstr[i]){
            if(outpstr[i]->s){
                outpstr[i]->s[outpstr[i]->l-1] = '\n';
            }
        }
    }
    FILE* fp = NULL;
    for(size_t i = 0; i < outmstr.size(); ++i){
        if(outmstr[i]){
            oval = ((uint8_t)i >> 2) & 0x3;
            fp = fopen(outmsnv[i].c_str(), "w");
            fwrite(outpstr[oval]->s, sizeof(char), outpstr[oval]->l, fp);
            outmstr[i]->l -= 1;
            fwrite(outmstr[i]->s, sizeof(char), outmstr[i]->l, fp);
            fclose(fp);
        }
    }
}

void hapc_usage(hap_opt_t* opt, char* arg0) {
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: %s %s [options]\n\n", PACKAGE_NAME, arg0);
    fprintf(stderr, "Options: -i FILE input bam\n");
    fprintf(stderr, "         -c STR  contig name\n");
    fprintf(stderr, "         -b INT  beg pos to count\n");
    fprintf(stderr, "         -e INT  end pos to count\n");
    fprintf(stderr, "         -r FILE reference nucleotide fasta file\n");
    fprintf(stderr, "         -a FILE alternate nucleotide fasta file\n");
    fprintf(stderr, "         -7 STR  extra snv to focus on, format(REF:ALT,REF:ALT2,...)\n");
    fprintf(stderr, "         -m INT  max linked snv in one haplotype [%d]\n", opt->maxh);
    fprintf(stderr, "         -s INT  min support reads for haplotype output [%d]\n", opt->minr);
    fprintf(stderr, "         -9 FILE js/css cdn configure json file\n");
    fprintf(stderr, "         -n INT  topN top plot [%d]\n", opt->gtn);
    fprintf(stderr, "         -N INT  topN top table(negative value for all) [%d]\n", opt->otn);
    fprintf(stderr, "         -o STR  output directory [%s]\n", opt->outdir.c_str());
    fprintf(stderr, "         -l      long haplotype snv accumulate to short haplotype[memory might explodes]\n");
    fprintf(stderr, "         -M      nucleotide sequence is from mitochondrion\n");
    fprintf(stderr, "         -3      input bam was not from 'getools caledit' if set\n");
    fprintf(stderr, "\n");
}

int hapc_main(int argc, char** argv) {
    hap_opt_t opt;
    if(argc == 1) {
        hapc_usage(&opt, argv[0]);
        return 0;
    }
    int c = -1;
    while((c = getopt(argc, argv, "i:c:b:e:m:o:r:a:9:s:n:7:N:3Mlh")) >= 0) {
        switch(c) {
            case 'i': opt.inbam = optarg; break;
            case 'c': opt.contig = optarg; break;
            case 'b': opt.beg = atoi(optarg); break;
            case 'e': opt.end = atoi(optarg); break;
            case 'm': opt.maxh = atoi(optarg); break;
            case '3': opt.nget = true; break;
            case 'M': opt.ismt = true; break;
            case 'r': opt.ref = optarg; break;
            case 'a': opt.alt = optarg; break;
            case '9': opt.jscdn = optarg; break;
            case 'n': opt.gtn = atoi(optarg); break;
            case 'N': opt.otn = atoi(optarg); break;
            case '7': opt.extsnv = optarg; break;
            case 'o': opt.outdir = optarg; break;
            case 's': opt.minr = atoi(optarg); break;
            case 'l': opt.l2s = true; break;
            case 'h': hapc_usage(&opt, argv[0]); return 0; break;
            default: break;
        }
    }
    if(opt.valid()) {
        util::loginfo(stderr, "beg update args");
        opt.update();
        util::loginfo(stderr, "end update args");
        util::loginfo(stderr, "beg count haps");
        opt.hapcnt8bf();
        util::loginfo(stderr, "end count haps");
        util::loginfo(stderr, "beg output haps");
        opt.hap2tsv();
        opt.bam2tsv();
        opt.aas2tsv();
        opt.hap2jsn();
        opt.stat2tsv();
        opt.hap2html();
        opt.snv2tsv();
        opt.aac2tsv();
        opt.snv2foc();
        opt.bias2tsv();
        util::loginfo(stderr, "end output haps");
        return 0;
    } else {
        return 1;
    }
}

void out_rbms(RecombBitMaskSet& rbms){
    for(size_t i = 0; i < rbms.size(); ++i){
        fprintf(stderr, "---n:%ld---\n", i);
        for(size_t j = 0; j < rbms[i].size(); ++j){
            fprintf(stderr, "c(%ld, %ld)\n", i, j);
            for(size_t k = 0; k < rbms[i][j].size(); ++k){
                fprintf(stderr, "[");
                for(size_t l = 0; l < rbms[i][j][k].size(); ++l){
                    fprintf(stderr, "%d,", rbms[i][j][k][l]);
                }
                fprintf(stderr, "]\n");
            }
        }
    }
}

void merge_and_sort_bam_rec8seq(BamRecList& brl){
    if(brl.empty()) return;
    std::sort(brl.begin(), brl.end(), bam_rec_seq_sorter());
    size_t i = 1, j = 0;
    bam_rec_t* pbr = brl[j];
    bam_rec_t* pcr = NULL;
    for(;i < brl.size(); ++i){
        pcr = brl[i];
        if(pcr->seq == pbr->seq){
            pbr->cnt += pcr->cnt;
            delete pcr;
            pcr = NULL;
        }else{
            ++j;
            brl[j] = brl[i];
            if(j < i) brl[i] = NULL;
            pbr = brl[j];
        }
    }
    brl.resize(j+1);
    std::sort(brl.begin(), brl.end(), bam_rec_cnt_sorter());
}

void merge_and_sort_bam_rec8aas(BamRecList& brl){
    if(brl.empty()) return;
    std::sort(brl.begin(), brl.end(), bam_rec_aas_sorter());
    size_t i = 1, j = 0;
    bam_rec_t* pbr = brl[j];
    bam_rec_t* pcr = NULL;
    for(;i < brl.size(); ++i){
        pcr = brl[i];
        if(pcr->aas == pbr->aas){
            pbr->cnt += pcr->cnt;
            delete pcr;
            pcr = NULL;
        }else{
            ++j;
            brl[j] = brl[i];
            if(j < i) brl[i] = NULL;
            pbr = brl[j];
        }
    }
    brl.resize(j+1);
    std::sort(brl.begin(), brl.end(), bam_rec_cnt_sorter());
}
