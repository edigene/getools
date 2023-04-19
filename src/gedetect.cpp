#include "gedetect.h"

bam_hdr_t* GEDetector::get_bam_hdr(){
    kstring_t* ks = (kstring_t*)calloc(1, sizeof(kstring_t));
    kputs("@SQ\tSN:", ks);
    kputs(name, ks);
    ksprintf(ks, "\tLN:%d\n", rlen);
    bam_hdr_t* h = sam_hdr_parse(ks->l, ks->s);
    h->l_text = ks->l;
    h->text = ks->s;
    free(ks);
    if(hap){
        delete hap;
        hap = NULL;
    }
    return h;
}

kswr_t* GEDetector::align(char* qseq, int qlen){
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
        if(ret->nvar == 0 && ret->lsc < 10 && ret->tsc < 10) ret->smask |= KSW_FREFTYPE;
        else if(ret->ndel || ret->nins){ // cas aware realign of indel
            int indel_left_cls = 0;
            int indel_right_cls = 0;
            int indel_in_cls = 0;
            int rpos = ret->tb;
            char opchr = '\0';
            int oplen = 0;
            for(int i = 0; i < ret->ncigar; ++i){
                opchr = kswge_cigar_opchr(ret->cigar[i]);
                oplen = kswge_cigar_oplen(ret->cigar[i]);
                switch(opchr){
                    case 'I':
                        if(rpos < clsbeg) ++indel_left_cls;
                        else if(rpos > clsend) ++indel_right_cls;
                        else ++indel_in_cls;
                        break;
                    case 'D':
                        if(rpos < clsbeg) ++indel_left_cls;
                        else if(rpos > clsend) ++indel_right_cls;
                        else ++indel_in_cls;
                        if((!opt->usesem) && oplen > maxdel) ret->smask |= KSW_FPRIMDIMER;
                        rpos += oplen;
                        break;
                    case 'X': case '=':
                        rpos += oplen;
                        break;
                    default: break;
                }
            }
            if(indel_left_cls > 0){ // right align these indels
                kswge_right_align(ret, qints, refints, 0, clsbeg);
            }
            if(indel_right_cls > 0){ // left align these indels
                kswge_left_align(ret, qints, refints, clsend, ret->te);
            }
            // update indel in sgrna range
            indel_left_cls = 0;
            indel_right_cls = 0;
            indel_in_cls = 0;
            rpos = ret->tb;
            for(int i = 0; i < ret->ncigar; ++i){
                switch(kswge_cigar_opchr(ret->cigar[i])){
                    case 'I':
                        if(rpos < clsbeg) ++indel_left_cls;
                        else if(rpos > clsend) ++indel_right_cls;
                        else ++indel_in_cls;
                        break;
                    case 'D':
                        if(rpos < clsbeg) ++indel_left_cls;
                        else if(rpos > clsend) ++indel_right_cls;
                        else ++indel_in_cls;
                        rpos += kswge_cigar_oplen(ret->cigar[i]);
                        break;
                    case 'X': case '=':
                        rpos += kswge_cigar_oplen(ret->cigar[i]);
                        break;
                    default: break;
                }
            }
            if(indel_in_cls){
                int32_t ocigarlen = ret->ncigar;
                uint32_t *ocigar = (uint32_t*)malloc(ocigarlen * sizeof(uint32_t));
                memcpy(ocigar, ret->cigar, sizeof(uint32_t) * ocigarlen);
                int32_t ovar = ret->nvar;
                int32_t omm = ret->nmm;
                int32_t oins = ret->nins;
                int32_t odel = ret->ndel;
                bool rtaln = right_realn_worked(ret->smask);
                kswge_right_align(ret, qints, refints, sgrbeg, sgrend);
                int new_indel_in_cls = 0;
                int rpos = ret->tb;
                for(int i = 0; i < ret->ncigar; ++i){
                    switch(kswge_cigar_opchr(ret->cigar[i])){
                        case 'I':
                            if(rpos >= clsbeg && rpos <= clsend) ++new_indel_in_cls;
                            break;
                        case 'D':
                            if(rpos >= clsbeg && rpos <= clsend) ++new_indel_in_cls;
                            rpos += kswge_cigar_oplen(ret->cigar[i]);
                            break;
                        case 'X': case '=':
                            rpos += kswge_cigar_oplen(ret->cigar[i]);
                            break;
                        default: break;
                    }
                }
                if(new_indel_in_cls < indel_in_cls){
                    free(ret->cigar);
                    ret->cigar = ocigar;
                    ret->ncigar = ocigarlen;
                    ret->nvar = ovar;
                    ret->nmm = omm;
                    ret->nins = oins;
                    ret->ndel = odel;
                    if(!rtaln) ret->smask &= (~KSW_FRALN);
                }else{
                    free(ocigar);
                }
            }
        }
        // check donor
        if(donorseq){
            kswr_t* dret = kswge_semi_global(qlen, qints, donorlen, donorints, 5, score_mat, gapopen, gapext);
            if(dret){
                kswge_mark_mismatch(dret, donorints, donorlen, qints, qlen);
                if(dret->nvar == 0 && ret->lsc < donorbeg && ret->tsc < donorlen - donorendd){ // exact recomb
                    ret->smask |= (KSW_FRECEXACT | KSW_FRECALLHIT | KSW_FRECANYHIT);
                }else{
                    int dextramm = 0;
                    int ddreammc = 0;
                    int rpos = ret->tb;
                    int qpos = 0;
                    char opchr;
                    int oplen;
                    int i, ddmm;
                    size_t xj;
                    for(int c = 0; c < ret->ncigar; ++c){
                        opchr = kswge_cigar_opchr(ret->cigar[c]);
                        oplen = kswge_cigar_oplen(ret->cigar[c]);
                        if(rpos > donorenda) break;
                        if(opchr == 'I'){
                            if(rpos >= donorbeg){
                                if(donorins[rpos].size() != (size_t)oplen) ++dextramm;
                                else{
                                    for(i = 0; i < oplen; ++i){
                                        if(qseq[qpos+i] != donorins[rpos][i]) break;
                                    }
                                    if(i != oplen) ++dextramm;
                                    else ++ddreammc;
                                }
                            }
                            qpos += oplen;
                        }else if(opchr == 'M' || opchr == 'X' || opchr == '='){
                            if(rpos >= donorbeg){
                                if(opchr == 'X'){
                                    if(oplen == 1){
                                        if(donorsnv[rpos] == qints[qpos]) ++ddreammc; // single idea snp
                                        else ++dextramm; // single unwanted snp
                                    }else{
                                        ddmm = 0;
                                        for(i = 0; i < oplen; ++i){
                                            if(donorsnv[rpos+i] == qints[qpos+i]){
                                                ++ddreammc; // single idea snp
                                                ++ddmm;
                                            }
                                            if(donordis[rpos+i].size()){
                                                if(donordis[rpos+i].size() <= (size_t)(oplen-i)){
                                                    for(xj = 0; xj < donordis[rpos+i].size(); ++xj){
                                                       if(qseq[qpos+i+xj] != donordis[rpos+i][xj]) break;
                                                    }
                                                    if(xj != donordis[rpos+i].size()) ++dextramm;
                                                    else{
                                                        ++ddreammc;
                                                        ddmm += donordis[rpos+i].size();
                                                    }
                                                }else{
                                                    ++dextramm;
                                                }
                                            }
                                        }
                                        if(dextramm == 0 && ddmm < oplen) ++dextramm;
                                    }
                                }
                            }
                            rpos += oplen;
                            qpos += oplen;
                        }else if(opchr == 'D'){
                            if(rpos >= donorbeg){
                                if(donordel[rpos] != oplen) ++dextramm;
                                else ++ddreammc;
                            }
                            rpos += oplen;
                        }else if(opchr == 'S'){
                            qpos += oplen;
                        }
                    }
                    if(dextramm == 0 && ddreammc == donormutcnt){
                        ret->smask |= (KSW_FRECEXACT | KSW_FRECALLHIT | KSW_FRECANYHIT);
                    }else if(ddreammc == donormutcnt){
                        ret->smask |= (KSW_FRECALLHIT | KSW_FRECANYHIT);
                    }else if(ddreammc){
                        ret->smask |= KSW_FRECANYHIT;
                    }
                }
            }
            if(dret) kswr_destroy(dret);
        }
    }else{
        ret = (kswr_t*)calloc(1, sizeof(kswr_t));
        ret->smask = 0;
        ret->cigar = NULL;
    }
    free(qints);
    return ret;
}

bam1_t* GEDetector::ksw2bam(kswr_t* ret, char* qseq, int qlen, char* qname, int qnlen, int mask){
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

void GEDetector::add_var(bam1_t* b){
    vars.push_back(b);
}

void GEDetector::mark_misprimeout(){
    uint32_t* cigar = NULL;
    int opint = 0, oplen = 0, qpos = 0, rpos = 0;
    int8_t *qmm = (int8_t*)malloc((2*rlen) * sizeof(int8_t));
    int8_t *rmm = (int8_t*)malloc((2*rlen) * sizeof(int8_t));
    int qmatch  = 0, rmatch = 0, scoff = 0, bequal = 0, bdiff = 0;
    KSW_FTYPE vsv = 0;
    bool ismisp = false, seedgt30 = false, seedgt60 = false;
    for(auto& b: vars){
        if(b->core.l_qseq > rlen) continue; // invalid read length
        qpos = 0; rpos = b->core.pos;
        cigar = bam_get_cigar(b);
        memset(qmm, 0, 2*rlen);
        memset(rmm, 0, 2*rlen);
        ismisp = seedgt30 = seedgt60 = false;
        bequal = bdiff = 0;
        // stat match/mism
        for(uint32_t i = 0; i < b->core.n_cigar; ++i){
            opint = bam_cigar_op(cigar[i]);
            oplen = bam_cigar_oplen(cigar[i]);
            switch(opint){
                case BAM_CEQUAL:
                    if(b->core.flag & BAM_FREVERSE){
                        for(int j = 0; j < oplen; ++j){
                            qmm[b->core.l_qseq-1-qpos-j] = 1;
                            rmm[rpos+j] = 1;
                        }
                    }else{
                        for(int j = 0; j < oplen; ++j){
                            qmm[qpos+j] = 1;
                            rmm[rpos+j] = 1;
                        }
                    }
                    bequal += oplen;
                    if(oplen > 30) seedgt30 = true;
                    if(oplen > 60) seedgt60 = true;
                    qpos += oplen;
                    rpos += oplen;
                    break;
                case BAM_CDIFF:
                    if(b->core.flag & BAM_FREVERSE){
                        for(int j = 0; j < oplen; ++j){
                            qmm[b->core.l_qseq-1-qpos-j] = 1;
                            rmm[rpos+j] = 1;
                        }
                    }else{
                        for(int j = 0; j < oplen; ++j){
                            qmm[qpos+j] = 1;
                            rmm[rpos+j] = 1;
                        }
                    }
                    bdiff += oplen;
                    qpos += oplen;
                    rpos += oplen;
                    break;
                case BAM_CINS:
                    qpos += oplen;
                    break;
                case BAM_CDEL:
                    rpos += oplen;
                    break;
                case BAM_CSOFT_CLIP:
                    qpos += oplen;
                    break;
                default:
                    break;
            }
        }
        if(seedgt60) continue;
        if(seedgt30 && (double)bequal/double(bequal+bdiff) > 0.8) continue;
        vsv = bam_aux2i(bam_aux_get(b, "VS"));
        // now do the test
        if(vsv & KSW_FMERGED){
            qmatch = 0;
            scoff = 0;
            while(!qmm[scoff] && scoff < b->core.l_qseq) ++scoff; // skipping leading softclip
            qmatch = scoff;
            while(qmm[qmatch] && qmatch < b->core.l_qseq) ++qmatch;
            qmatch -= scoff;
            rmatch = b->core.pos;
            while(rmm[rmatch++] && rmatch < rpos);
            rmatch -= (b->core.pos+1);
            if(qmatch < mpmlen || rmatch < mpmlen) ismisp = true;
            if(!ismisp){
                qmatch = rmatch = 0;
                scoff = b->core.l_qseq-1;
                while(scoff >= 0 && !qmm[scoff]) --scoff;
                qmatch = scoff;
                while(qmatch >= 0 && qmm[qmatch]) --qmatch;
                qmatch = scoff - qmatch;
                rmatch = rpos-1;
                while(rmm[rmatch--] && rmatch >= b->core.pos);
                rmatch = rpos-2-rmatch;
                if(qmatch < mpmlen || rmatch < mpmlen) ismisp = true;
            }
        }else{
            qmatch = 0;
            scoff = 0;
            while(!qmm[scoff] && scoff < b->core.l_qseq) ++scoff;
            qmatch = scoff;
            while(qmm[qmatch] && qmatch < b->core.l_qseq) ++qmatch;
            qmatch -= scoff;
            if(qmatch < mpmlen) ismisp = true;
            else{
                if(b->core.flag & BAM_FREVERSE){
                    rmatch = rpos-1;
                    while(rmm[rmatch--] && rmatch >= b->core.pos);
                    rmatch = rpos-2-rmatch;
                }else{
                    rmatch = b->core.pos;
                    while(rmm[rmatch++] && rmatch < rpos);
                    rmatch -= b->core.pos;
                }
                if(rmatch < mpmlen) ismisp = true;
            }
        }
        if(ismisp){
            vsv |= KSW_FMISPRIMOUT;
            bam_aux_update_int(b, "VS", vsv);
        }
    }
    free(qmm); qmm = NULL;
    free(rmm); rmm = NULL;
}

void GEDetector::mark_misprimeinn(){
    int same = 0, diff = 0, opint = 0, oplen = 0, opni = 0, opnl = 0, onni = 0, onnl = 0;
    uint32_t* cigar = NULL;
    bool ismisp = false;
    KSW_FTYPE vs = 0;
    for(auto& b: vars){
        same = diff = opint = oplen = opni = opnl = onni = onnl = 0;
        cigar = bam_get_cigar(b);
        ismisp = false;
        for(uint32_t i = 0; i < b->core.n_cigar;){
            opint = bam_cigar_op(cigar[i]);
            oplen = bam_cigar_oplen(cigar[i]);
            if(opint == BAM_CINS || opint == BAM_CDEL){// got I or D
                same = diff = 0;
                while(i+1 < b->core.n_cigar){
                    opni = bam_cigar_op(cigar[i+1]);
                    opnl = bam_cigar_oplen(cigar[i+1]);
                    if(opni == BAM_CEQUAL) same += opnl;
                    else if(opni == BAM_CDIFF) diff += opnl;
                    else break;
                    ++i;
                }
                ++i;
                onnl = 0;
                if(same || diff){// got M
                    if(i < b->core.n_cigar){
                        onni = bam_cigar_op(cigar[i]);
                        if((onni == BAM_CDEL || onni == BAM_CINS) && onni != opint){// got D or I
                            onnl = bam_cigar_oplen(cigar[i]);
                        }
                    }
                }
                if(oplen && onnl){// IMD got
                    if((onnl > 9 && oplen > 9 && (same+diff < 6 || (double)same/(double)(same+diff) < 0.8)) ||
                       ((onnl > 29 || oplen > 29) && (same+diff < 6 || (double)same/(double)(same+diff) < 0.8))){
                        ismisp = true;
                        break;
                   }
                }
            }else{
                ++i;
            }
        }
        if(ismisp){
            vs = bam_aux2i(bam_aux_get(b, "VS"));
            vs |= KSW_FMISPRIMINN;
            bam_aux_update_int(b, "VS", vs);
        }
    }
}

void GEDetector::fixpe(){
    KSW_FTYPE vsv1 = 0, vsv2 = 0;
    bool upd1 = false, upd2 = false;
    // fix pe partner
    for(size_t i = 0; i < vars.size(); ){
        upd1 = upd2 = false;
        vsv1 = bam_aux2i(bam_aux_get(vars[i], "VS"));
        if(vsv1 & KSW_FMERGED) ++i;
        else{
            vsv2 = bam_aux2i(bam_aux_get(vars[i+1], "VS"));
            // fix KSW_FMISPRIMINN
            if(vsv1 & KSW_FMISPRIMINN){
                if(!(vsv2 & KSW_FMISPRIMINN)){
                    vsv2 |= KSW_FMISPRIMINN;
                    upd2 = true;
                }
            }else if(vsv2 & KSW_FMISPRIMINN){
                vsv1 |= KSW_FMISPRIMINN;
                upd1 = true;
            }
            // fix KSW_FMISPRIMOUT
            if(vsv1 & KSW_FMISPRIMOUT){
                if(!(vsv2 & KSW_FMISPRIMOUT)){
                    vsv2 |= KSW_FMISPRIMOUT;
                    upd2 = true;
                }
            }else if(vsv2 & KSW_FMISPRIMOUT){
                vsv1 |= KSW_FMISPRIMOUT;
                upd1 = true;
            }
            // fix KSW_FPRIMDIMER
            if(vsv1 & KSW_FPRIMDIMER){
                if(!(vsv2 & KSW_FPRIMDIMER)){
                    vsv2 |= KSW_FPRIMDIMER;
                    upd2 = true;
                }
            }else if(vsv2 & KSW_FPRIMDIMER){
                vsv1 |= KSW_FPRIMDIMER;
                upd1 = true;
            }
            // fix KSW_FHIERRLOWQ
            if(vsv1 & KSW_FHIERRLOWQ){
                if(!(vsv2 & KSW_FHIERRLOWQ)){
                    vsv2 |= KSW_FHIERRLOWQ;
                    upd2 = true;
                }
            }else if(vsv2 & KSW_FHIERRLOWQ){
                vsv1 |= KSW_FHIERRLOWQ;
                upd1 = true;
            }
            // fix KSW_FMANYVARS
            if(vsv1 & KSW_FMANYVARS){
                if(!(vsv2 & KSW_FMANYVARS)){
                    vsv2 |= KSW_FMANYVARS;
                    upd2 = true;
                }
            }else if(vsv2 & KSW_FMANYVARS){
                vsv1 |= KSW_FMANYVARS;
                upd1 = true;
            }
            // update if needed
            if(upd1) bam_aux_update_int(vars[i], "VS", vsv1);
            if(upd2) bam_aux_update_int(vars[i+1], "VS", vsv2);
            i += 2;
        }
    }
}

void GEDetector::cal_clspos(){
    if(clsbeg == clsend){ clspos = clsbeg; return; } // if provided as a single pos, do not detect
    int* depth = (int*)calloc(rlen+1, sizeof(int));
    uint8_t* ccd = NULL;
    int64_t ccv = 0;
    uint32_t* cigar = NULL;
    int opint, oplen, rpos;
    for(auto& b: vars){
        if(b->core.tid < 0) continue;
        rpos = b->core.pos;
        cigar = bam_get_cigar(b);
        ccd = bam_aux_get(b, "CC");
        ccv = bam_aux2i(ccd);
        for(uint32_t i = 0; i < b->core.n_cigar; ++i){
            oplen = bam_cigar_oplen(cigar[i]);
            opint = bam_cigar_op(cigar[i]);
            switch(opint){
                case BAM_CINS: case BAM_CSOFT_CLIP: 
                    break;
                case BAM_CEQUAL: case BAM_CDIFF:
                    for(int j = 0; j < oplen; ++j) depth[rpos+j] += ccv;
                    rpos += oplen;
                    break;
                case BAM_CDEL:
                    rpos += oplen;
                    break;
                default:
                    break;
            }
        }
    }
    clspos = clsbeg;
    int mindep = depth[clspos];
    for(int i = clsbeg + 1; i <= clsend; ++i){
        if(depth[i] < mindep){
            mindep = depth[i];
            clspos = i;
        }
    }
    free(depth);
}

void GEDetector::stat_var(bam1_t* b, bool toadd){
    int cc = bam_aux2i(bam_aux_get(b, "CC"));
    int ins = 0, del = 0, din = 0, snv = 0, scl = 0;
    uint32_t* cigar = bam_get_cigar(b);
    int opint, oplen;
    int rpos = b->core.pos;
    int qpos = 0;
    for(size_t i = 0; i < b->core.n_cigar; ++i){
        opint = bam_cigar_op(cigar[i]);
        oplen = bam_cigar_oplen(cigar[i]);
        switch(opint){
            case BAM_CINS:
                ++ins;
                ins_len_dist->statvar(oplen, cc, toadd);
                ins_pos_dist->statvar(rpos, cc, toadd);
                ins_lps_dist->statvar(rpos, cc * oplen, toadd);
                qpos += oplen;
                break;
            case BAM_CDEL:
                ++del;
                del_len_dist->statvar(oplen, cc, toadd);
                del_lps_dist->statvar(rpos, cc * oplen, toadd);
                for(int dp = 0; dp < oplen; ++dp) del_pos_dist->statvar(rpos+dp, cc, toadd);
                rpos += oplen;
                break;
            case BAM_CDIFF:
                if(oplen > 1){// delins
                    ++din;
                    din_len_dist->statvar(oplen, cc, toadd);
                    din_lps_dist->statvar(rpos, cc * oplen, toadd);
                    for(int dp = 0; dp < oplen; ++dp) din_pos_dist->statvar(rpos+dp, cc, toadd);
                }else{
                    ++snv;
                    snv_pos_dist->statvar(rpos, cc, toadd);
                }
                if(toadd){
                    for(int p = 0; p < oplen; ++p) nccnt[seq_nt16_int[bam_seqi(bam_get_seq(b), qpos + p)]][rpos+p] += cc;
                }else{
                    for(int p = 0; p < oplen; ++p) nccnt[seq_nt16_int[bam_seqi(bam_get_seq(b), qpos + p)]][rpos+p] -= cc;
                }
                rpos += oplen;
                qpos += oplen;
                break;
            case BAM_CEQUAL:
                if(toadd){
                    for(int p = 0; p < oplen; ++p) nccnt[seq_nt16_int[bam_seqi(bam_get_seq(b), qpos + p)]][rpos+p] += cc;
                }else{
                    for(int p = 0; p < oplen; ++p) nccnt[seq_nt16_int[bam_seqi(bam_get_seq(b), qpos + p)]][rpos+p] -= cc;
                }
                rpos += oplen;
                qpos += oplen;
                break;
            case BAM_CSOFT_CLIP:
                ++scl;
                scl_len_dist->statvar(oplen, cc, toadd);
                scl_pos_dist->statvar(rpos, cc, toadd);
                qpos += oplen;
                break;
            default:
                break;
        }
    }
    tot_cnt_dist->statvar(bam_aux2i(bam_aux_get(b, "VC")), cc, toadd);
    ins_cnt_dist->statvar(ins, cc, toadd);
    del_cnt_dist->statvar(del, cc, toadd);
    din_cnt_dist->statvar(din, cc, toadd);
    snv_cnt_dist->statvar(snv, cc, toadd);
    scl_cnt_dist->statvar(scl, cc, toadd);
    if(!ins) ins_len_dist->statvar(0, cc, toadd);
    if(!del) del_len_dist->statvar(0, cc, toadd);
    if(!din) din_len_dist->statvar(0, cc, toadd);
    if(!scl) scl_len_dist->statvar(0, cc, toadd);
}

void GEDetector::var_count(bam1_t* b1, bam1_t* b2){
    std::vector<int> varm1(rlen, 0);
    std::vector<int> varm2(rlen, 0);
    uint8_t* vsd1 = bam_aux_get(b1, "VS");
    int64_t vsv1 = bam_aux2i(vsd1);
    // stat r1
    int r1beg = b1->core.pos, r1end = b1->core.pos;
    uint32_t* cigar = bam_get_cigar(b1);
    int rpos = b1->core.pos;
    int oplen = 0;
    int opint = 0;
    for(size_t i = 0; i < b1->core.n_cigar; ++i){
        opint = bam_cigar_op(cigar[i]);
        oplen = bam_cigar_oplen(cigar[i]);
        switch(opint){
            case BAM_CDEL:
                if(rpos <= vvend && rpos + oplen - 1 >= vvbeg) vsv1 |= KSW_FDELINRNG;
                if(varmask & KSW_FDELINRNG) vsv1 |= KSW_FVARINSEQ;
                varm1[rpos] = -oplen;
                rpos += oplen;
                break;
            case BAM_CINS:
                if(rpos > vvbeg && rpos <= vvend) vsv1 |= KSW_FINSINRNG;
                if(varmask & KSW_FINSINRNG) vsv1 |= KSW_FVARINSEQ;
                varm1[rpos] = oplen;
                break;
            case BAM_CDIFF:
                if(rpos <= vvend && rpos + oplen - 1 >= vvbeg){
                    if(oplen == 1){
                        vsv1 |= KSW_FSNVINRNG;
                        if(varmask & KSW_FSNVINRNG) vsv1 |= KSW_FVARINSEQ;
                    }else{
                        vsv1 |= KSW_FDIINRNG;
                        if(varmask & KSW_FDIINRNG) vsv1 |= KSW_FVARINSEQ;
                    }
                }
                for(int p = 0; p < oplen; ++p) varm1[rpos+p] = rlen;
                rpos += oplen;
                break;
            case BAM_CEQUAL:
                rpos += oplen;
                break;
            default:
                break;
        }
    }
    r1end = rpos;
    if(b1->core.pos <= sgrbeg && r1end > sgrend) vsv1 |= KSW_FSPANSGR;
    // stat r2 and overlap diff penalty
    uint8_t* vsd2 = bam_aux_get(b2, "VS");
    int64_t vsv2 = bam_aux2i(vsd2);
    int32_t r2beg = b2->core.pos, r2end = b2->core.pos;
    cigar = bam_get_cigar(b2);
    rpos = b2->core.pos;
    for(size_t i = 0; i < b2->core.n_cigar; ++i){
        opint = bam_cigar_op(cigar[i]);
        oplen = bam_cigar_oplen(cigar[i]);
        switch(opint){
            case BAM_CDEL:
                if(rpos <= vvend && rpos + oplen - 1 >= vvbeg) vsv2 |= KSW_FDELINRNG;
                if(varmask & KSW_FDELINRNG) vsv2 |= KSW_FVARINSEQ;
                varm2[rpos] = -oplen;
                rpos += oplen;
                break;
            case BAM_CINS:
                if(rpos > vvbeg && rpos <= vvend) vsv2 |= KSW_FINSINRNG;
                if(varmask & KSW_FINSINRNG) vsv2 |= KSW_FVARINSEQ;
                varm2[rpos] = oplen;
                break;
            case BAM_CDIFF:
                if(rpos <= vvend && rpos + oplen - 1 >= vvbeg){
                    if(oplen == 1){
                        vsv2 |= KSW_FSNVINRNG;
                        if(varmask & KSW_FSNVINRNG) vsv2 |= KSW_FVARINSEQ;
                    }else{
                        vsv2 |= KSW_FDIINRNG;
                        if(varmask & KSW_FDIINRNG) vsv2 |= KSW_FVARINSEQ;
                    }
                }
                for(int p = 0; p < oplen; ++p) varm2[rpos+p] = rlen;
                rpos += oplen;
                break;
            case BAM_CEQUAL:
                rpos += oplen;
                break;
            default:
                break;
        }
    }
    r2end = rpos;
    if(b2->core.pos <= sgrbeg && r2end > sgrend) vsv2 |= KSW_FSPANSGR;
    // stat overlap
    int olen = MIN(r2end, r1end) - MAX(r1beg, r2beg);
    if(olen > 0){
        vsv1 |= KSW_FPAIROLP;
        vsv2 |= KSW_FPAIROLP;
        b1->core.flag |= BAM_FPROPER_PAIR;
        b2->core.flag |= BAM_FPROPER_PAIR;
    }
    // continuous IDX count as one biological meaningful event
    int cnt1 = 0, cnt2 = 0;
    for(size_t i = 0; i < varm1.size();){
        if(varm1[i] != 0){
            ++cnt1;
            while(i < varm1.size() && varm1[i] != 0) ++i;
        }else{
            ++i;
        }
    }
    for(size_t i = 0; i < varm2.size();){
        if(varm2[i] != 0){
            ++cnt2;
            while(i < varm2.size() && varm2[i] != 0) ++i;
        }else{
            ++i;
        }
    }
    bam_aux_update_int(b1, "VC", cnt1);
    bam_aux_update_int(b2, "VC", cnt2);
    if(cnt1 > maxvar) vsv1 |= KSW_FMANYVARS;
    if(cnt2 > maxvar) vsv2 |= KSW_FMANYVARS;
    if(cnt1 > maxvc) maxvc = cnt1;
    if(cnt2 > maxvc) maxvc = cnt2;
    /***************************FIXME TO CHOOSE AN OPTIMAL REP***************************/
    // choose representative one
    if((vsv1 & KSW_FSPANSGR) && (vsv2 & KSW_FSPANSGR)){// both span sgrna
        if(vsv1 & varmask){
            vsv1 |= KSW_FREPSEQR;
            if((vsv1 & KSW_FIDLOCMP) || greedy) vsv1 |= KSW_FMAYVARSEQ;
        }else if(vsv2 & varmask){
            vsv2 |= KSW_FREPSEQR;
            if((vsv2 & KSW_FIDLOCMP) || greedy) vsv2 |= KSW_FMAYVARSEQ;
        }else{
            if(vsv1 & KSW_FSNVINRNG){
                vsv1 |= KSW_FREPSEQR;
            }else if(vsv2 & KSW_FSNVINRNG){
                vsv2 |= KSW_FREPSEQR;
            }else{
                vsv1 |= KSW_FREPSEQR;
            }
        }
    }else if(vsv1 & KSW_FSPANSGR){ // read1 span sgrna
        vsv1 |= KSW_FREPSEQR;
        if(vsv1 & varmask){
            if((vsv1 & KSW_FIDLOCMP) || greedy) vsv1 |= KSW_FMAYVARSEQ;
        }
    }else if(vsv2 & KSW_FSPANSGR){ // read2 span sgrna
        vsv2 |= KSW_FREPSEQR;
        if(vsv2 & varmask){
            if((vsv2 & KSW_FIDLOCMP) || greedy) vsv2 |= KSW_FMAYVARSEQ;
        }
    }else{ // neither span sgrna
        vsv1 |= KSW_FREPSEQR;
    }
    /***************************FIXME TO CHOOSE AN OPTIMAL REP***************************/
    // donor rec
    if(vsv1 & recmask) vsv1 |= (KSW_FMAYVARSEQ | KSW_FVARINSEQ);
    if(vsv2 & recmask) vsv2 |= (KSW_FMAYVARSEQ | KSW_FVARINSEQ);
    bam_aux_update_int(b1, "VS", vsv1);
    bam_aux_update_int(b2, "VS", vsv2);
}

void GEDetector::mark_var(){
    uint16_t srm = (BAM_FPAIRED | BAM_FREAD1 | BAM_FREAD2);
    for(uint32_t i = 0; i < vars.size();){
        if(vars[i]->core.flag & srm){
            var_count(vars[i], vars[i+1]);
            i += 2;
        }else{
            var_count(vars[i]);
            ++i;
        }
    }
}

void GEDetector::var_count(bam1_t* b){
    if(b->core.n_cigar == 0){
        bam_aux_update_int(b, "VC", 0);
        return;
    }
    uint32_t* cigar = bam_get_cigar(b);
    std::vector<int> varm(b->core.n_cigar, 0);
    int32_t rpos = b->core.pos;
    uint8_t* vsd = bam_aux_get(b, "VS");
    int64_t vsv = bam_aux2i(vsd);
    int opint = 0;
    int oplen = 0;
    for(size_t i = 0; i < b->core.n_cigar; ++i){
        opint = bam_cigar_op(cigar[i]);
        oplen = bam_cigar_oplen(cigar[i]);
        switch(opint){
            case BAM_CDEL:
                varm[i] = 1;
                if(rpos <= vvend && rpos + oplen - 1 >= vvbeg) vsv |= KSW_FDELINRNG;
                if(varmask & KSW_FDELINRNG) vsv |= KSW_FVARINSEQ;
                rpos += oplen;
                break;
            case BAM_CINS:
                varm[i] = 1;
                if(rpos > vvbeg && rpos <= vvend) vsv |= KSW_FINSINRNG;
                if(varmask & KSW_FINSINRNG) vsv |= KSW_FVARINSEQ;
            case BAM_CDIFF:
                if(rpos <= vvend && rpos + oplen - 1 >= vvbeg){
                    if(oplen == 1){
                        vsv |= KSW_FSNVINRNG;
                        if(varmask & KSW_FSNVINRNG) vsv |= KSW_FVARINSEQ;
                    }else{
                        vsv |= KSW_FDIINRNG;
                        if(varmask & KSW_FDIINRNG) vsv |= KSW_FVARINSEQ;
                    }
                }
                rpos += oplen;
                varm[i] = 1;
                break;
            case BAM_CEQUAL:
                rpos += oplen;
            default:
                break;
        }
    }
    if(b->core.pos <= sgrbeg && rpos > sgrend) vsv |= KSW_FSPANSGR;
    vsv |= KSW_FREPSEQR;
    if((vsv & KSW_FSPANSGR) && (vsv & varmask)) vsv |= KSW_FMAYVARSEQ;
    // continuous IDX count as one biological meaningful event
    int cnt = 0;
    for(size_t i = 0; i < varm.size();){
        if(varm[i] == 1){
            ++cnt;
            while(i < varm.size() && varm[i] == 1) ++i;
        }else{
            ++i;
        }
    }
    bam_aux_update_int(b, "VC", cnt);
    if(cnt > maxvc) maxvc = cnt;
    if(cnt > maxvar) vsv |= KSW_FMANYVARS;
    // donor rec
    if(vsv & recmask) vsv |= (KSW_FMAYVARSEQ | KSW_FVARINSEQ);
    bam_aux_update_int(b, "VS", vsv);
}

void GEDetector::kcluster(){
    if(maxvc == 0) return; // all reftype, no need
    if(vars.size() < mincpnt) return; // not enough data
    int scale = 100;
    std::vector<int> memids(1 << 14, -1);
    std::vector<int> scores; scores.reserve(vars.size());
    std::vector<int> varcnt; varcnt.reserve(vars.size());
    std::vector<int> seqcnt; seqcnt.reserve(vars.size());
    std::vector<uint16_t> varids; varids.reserve(vars.size());
    uint8_t* vcd = NULL;
    uint8_t* asd = NULL;
    uint8_t* vsd = NULL;
    int vcv = 0, asv = 0, vsv = 0, asf = 0, vcf = 0;
    uint16_t vci = 0, asi = 0, idx = 0;
    for(size_t i = 0; i < vars.size(); ++i){
        vcd = bam_aux_get(vars[i], "VC");
        asd = bam_aux_get(vars[i], "AS");
        vcv = bam_aux2i(vcd);
        asv = bam_aux2i(asd);
        asf = (double)asv/(double)(match * vars[i]->core.l_qseq) * scale;
        if(asf < 0) asf = 0;
        if(asf > scale) asf = scale;
        vcf = (double)vcv/(double)maxvc * scale;
        if(vcf < 0) vcf = 0;
        if(vcf > scale) vcf = scale;
        vci = vcf; asi = asf;
        idx = (vci << 7) | asi; // do not bit shift on singed ints
        if(memids[idx] < 0){// new points
            memids[idx] = scores.size();
            scores.push_back(asf);
            varcnt.push_back(vcf);
            seqcnt.push_back(bam_aux2i(bam_aux_get(vars[i], "CC")));
        }else{
            seqcnt[memids[idx]] += bam_aux2i(bam_aux_get(vars[i], "CC"));
        }
        varids.push_back(memids[idx]);
    }
    // k-means cluster
    std::vector<std::array<int, 2>> points; points.reserve(scores.size());
    std::vector<int> badm(scores.size(), 0);
    for(size_t i = 0; i < scores.size(); ++i) points.push_back({scores[i], varcnt[i]});
    std::tuple<std::vector<std::array<int, 2>>, std::vector<uint32_t>> cls2ret = kmeans_lloyd(points, 2, maxciter);
    std::tuple<std::vector<std::array<int, 2>>, std::vector<uint32_t>> cls3ret = kmeans_lloyd(points, 3, maxciter);
    double kcoef2 = silhouette_coef(cls2ret, points);
    double kcoef3 = silhouette_coef(cls3ret, points);
    if(kcoef3 > kcoef2){ // filter available, mask filtered points
        cfworked = true;
        int ms = std::get<0>(cls3ret)[0][0];
        uint32_t mi = 0;
        for(uint32_t i = 0; i < 3; ++i){
            if(std::get<0>(cls3ret)[i][0] < ms){
                ms = std::get<0>(cls3ret)[i][0];
                mi = i;
            }
        }
        int mc = std::get<0>(cls3ret)[mi][1];
        uint32_t mj = mi;
        for(uint32_t i = 0; i < 3; ++i){
            if(std::get<0>(cls3ret)[i][1] > mc){
                mc = std::get<0>(cls3ret)[i][1];
                mj = i;
            }
        }
        if(mi != mj) return; // not perfect lowQ/highV cluster
        auto& lv = std::get<1>(cls3ret);
        for(size_t i = 0; i < lv.size(); ++i){
            if(lv[i] == mi) badm[i] = 1;
        }
        for(uint32_t i = 0; i < vars.size(); ++i){
            if(badm[varids[i]]){
                vsd = bam_aux_get(vars[i], "VS");
                vsv = bam_aux2i(vsd);
                vsv |= KSW_FHIERRLOWQ;
                bam_aux_update_int(vars[i], "VS", vsv);
            }
        }
        if(opt->tpdetail){
            // construct rscript
            kstring_t* pstr = (kstring_t*)calloc(1, sizeof(kstring_t));
            ksprintf(pstr, "score_v <- c(");
            for(uint32_t i = 0; i < points.size(); ++i){
                ksprintf(pstr, "%d,", points[i].at(0));
            }
            pstr->s[pstr->l-1] = ')';
            kputc('\n', pstr);

            ksprintf(pstr, "vacnt_v <- c(");
            for(uint32_t i = 0; i < points.size(); ++i){
                ksprintf(pstr, "%d,", points[i].at(1));
            }
            pstr->s[pstr->l-1] = ')';
            kputc('\n', pstr);

            ksprintf(pstr, "group_v <- c(");
            for(size_t i = 0; i < lv.size(); ++i){
                ksprintf(pstr, "%d,", lv[i]+1);
            }
            pstr->s[pstr->l-1] = ')';
            kputc('\n', pstr);

            ksprintf(pstr, "color_v <- c(");
            const char* optcol[3] = {"green", "purple", "blue"};
            const char* selcol[3] = {NULL, NULL, NULL};
            selcol[mi] = "red";
            for(size_t i = 0; i < lv.size(); ++i){
                if(lv[i] == mi){
                    ksprintf(pstr, "\"red\",");
                }else{
                    ksprintf(pstr, "\"%s\",", optcol[lv[i]]);
                    selcol[lv[i]] = optcol[lv[i]];
                }
            }
            pstr->s[pstr->l-1] = ')';
            kputc('\n', pstr);

            ksprintf(pstr, "pdf(file=\"%s/%s.cluster.pdf\", width=8, height=8, onefile = F)\n", opt->pltdir.c_str(), name);
            ksprintf(pstr, "plot(x=score_v, y=vacnt_v, xlim=c(0,100), ylim=c(0, 100), type=\"p\", pch=20, col=color_v, ");
            ksprintf(pstr, "xlab=\"Scaled Alignment Score\", ylab=\"Scaled Variant Count\", main=\"K-Means Cluster of Sequences\")\n");
            ksprintf(pstr, "legend(\"topleft\", col=c(\"%s\",\"%s\",\"%s\"), pch=20, legend=c(", selcol[0], selcol[1], selcol[2]);
            for(uint32_t i = 0; i < 3; ++i){
                if(i == mi){
                    ksprintf(pstr, "\"drop\",");
                }else{
                    ksprintf(pstr, "\"keep\",");
                }
            }
            pstr->s[pstr->l-1] = ')';
            ksprintf(pstr, ")\n");
            ksprintf(pstr, "points(x=c(%d, %d, %d), y=c(%d, %d, %d), pch=25, bg=c(\"%s\",\"%s\",\"%s\"), col=c(\"%s\",\"%s\",\"%s\"))\n",
                            std::get<0>(cls3ret)[0][0], std::get<0>(cls3ret)[1][0], std::get<0>(cls3ret)[2][0],
                            std::get<0>(cls3ret)[0][1], std::get<0>(cls3ret)[1][1], std::get<0>(cls3ret)[2][1],
                            selcol[0], selcol[1], selcol[2],
                            selcol[0], selcol[1], selcol[2]);
            int64_t ccumc[3] = {0, 0, 0};
            for(size_t i = 0; i < varids.size(); ++i) ccumc[lv[varids[i]]] += bam_aux2i(bam_aux_get(vars[i], "CC"));

            ksprintf(pstr, "text(x=c(%d, %d, %d), y=c(%d, %d, %d), labels=c(%lld, %lld, %lld), cex=0.5, pos=1)\n",
                            std::get<0>(cls3ret)[0][0], std::get<0>(cls3ret)[1][0], std::get<0>(cls3ret)[2][0],
                            std::get<0>(cls3ret)[0][1], std::get<0>(cls3ret)[1][1], std::get<0>(cls3ret)[2][1],
                            ccumc[0], ccumc[1], ccumc[2]);
            ksprintf(pstr, "gbinf <- dev.off()\n");

            // write rscript
            kstring_t tmpstr = {0, 0, 0};
            ksprintf(&tmpstr, "%s/%s.cluster.R", opt->pltdir.c_str(), name);
            FILE* fp = fopen(tmpstr.s, "w");
            fwrite(pstr->s, sizeof(char), pstr->l, fp); free(pstr->s); free(pstr);
            fclose(fp);

            // do plot
            kstring_t rcmd = {0, 0, 0};
            char* rspath = util::which("Rscript");
            if(rspath){
                ksprintf(&rcmd, "%s %s", rspath, tmpstr.s);
                kputc('\0', &rcmd);
                system(rcmd.s);
                free(rspath);
                free(rcmd.s);
            }
            free(tmpstr.s);
        }
    }
    // javascript output
    {
        kstring_t* s = kmcjs;
        std::string subsect = "Sequence cluster";
        std::string divName = util::replace(subsect, " ", "_");
        std::string title = "kmeans cluster sequences into 2 and 3 clusters";
        ksprintf(s, "<div class='subsection_title'><a title='click to hide/show' onclick=showOrHide('%s')>%s</a></div>\n", divName.c_str(), subsect.c_str());
        ksprintf(s, "<div id='%s'>\n", divName.c_str());
        ksprintf(s, "<div class='sub_section_tips'>kmeans cluster sequences into 2(k2) and 3(k3) clusters, cluster centers are star shape in black, points information will be shown on mouth over(ar: alignment score ratio, vc: variant count, sc: sequence count, cs: cluster sequence count).</div>\n");
        ksprintf(s, "<div class='figure' id='plot_%s'></div>\n", divName.c_str());
        ksprintf(s, "</div>\n");
        ksprintf(s, "\n<script type=\"text/javascript\">\n");
        std::vector<std::string> optcol = {"green", "purple", "blue"};
        std::vector<std::string> optpch = {"x", "o", "5", "17"}; // diamond, dot, triangle, star
        // k2
        std::string jsnstr;
        {
            jsnstr.append("var k2c = {\n");
            jsnstr.append(" x: [");
            for(uint32_t i = 0; i < points.size(); ++i) jsnstr.append(std::to_string(points[i].at(0)) + ",");
            for(int i = 0; i < 2; ++i) jsnstr.append(std::to_string(std::get<0>(cls2ret)[i][0]) + ",");
            jsnstr.append("],\n");
            jsnstr.append("  y: [");
            for(uint32_t i = 0; i < points.size(); ++i) jsnstr.append(std::to_string(points[i].at(1)) + ",");
            for(int i = 0; i < 2; ++i) jsnstr.append(std::to_string(std::get<0>(cls2ret)[i][1]) + ",");
            jsnstr.append("],\n");
            jsnstr.append("  name: 'k2',\n");
            jsnstr.append("  type: 'scatter',\n");
            jsnstr.append("  text: [");
            for(uint32_t i = 0; i < points.size(); ++i){
                jsnstr.append("'ar:" + std::to_string((double)points[i].at(0)/(double)scale) + "<br>");
                jsnstr.append("vc:" + std::to_string((double)(points[i].at(1)*maxvc)/(double)scale) + "<br>");
                jsnstr.append("sc:" + std::to_string(seqcnt[i]) + "',");
            }
            int64_t ccumc[3] = {0, 0, 0};
            auto& lv = std::get<1>(cls2ret);
            for(size_t i = 0; i < varids.size(); ++i) ccumc[lv[varids[i]]] += bam_aux2i(bam_aux_get(vars[i], "CC"));
            for(int i = 0; i < 2; ++i) jsnstr.append("'cs: " + std::to_string(ccumc[i]) + "',");
            jsnstr.append("],\n");
            jsnstr.append("  domain: {row: 0, col: 0,},\n");
            jsnstr.append("  marker: {\n");
            jsnstr.append("    color: [");
            for(size_t i = 0; i < lv.size(); ++i) jsnstr.append("'" + optcol[lv[i]] + "',");
            for(int i = 0; i < 2; ++i) jsnstr.append("'black',");
            jsnstr.append("],\n");
            jsnstr.append("    symbol: [");
            for(size_t i = 0; i < lv.size(); ++i) jsnstr.append("'" + optpch[lv[i]] + "',");
            for(int i = 0; i < 2; ++i) jsnstr.append("'" + optpch[3] + "',");
            jsnstr.append("],\n");
            jsnstr.append("  },\n");
            jsnstr.append("  mode: 'markers',\n");
            jsnstr.append("  xaxis: 'x1',\n");
            jsnstr.append("  yaxis: 'y1',\n");
            jsnstr.append("};\n");
        }
        // k3
        {
            jsnstr.append("var k3c = {\n");
            jsnstr.append(" x: [");
            for(uint32_t i = 0; i < points.size(); ++i) jsnstr.append(std::to_string(points[i].at(0)) + ",");
            for(int i = 0; i < 3; ++i) jsnstr.append(std::to_string(std::get<0>(cls3ret)[i][0]) + ",");
            jsnstr.append("],\n");
            jsnstr.append("  y: [");
            for(uint32_t i = 0; i < points.size(); ++i) jsnstr.append(std::to_string(points[i].at(1)) + ",");
            for(int i = 0; i < 3; ++i) jsnstr.append(std::to_string(std::get<0>(cls3ret)[i][1]) + ",");
            jsnstr.append("],\n");
            jsnstr.append("  text: [");
            for(uint32_t i = 0; i < points.size(); ++i){
                jsnstr.append("'ar:" + std::to_string((double)points[i].at(0)/(double)scale) + "<br>");
                jsnstr.append("vc:" + std::to_string((double)(points[i].at(1)*maxvc)/(double)scale) + "<br>");
                jsnstr.append("sc:" + std::to_string(seqcnt[i]) + "',");
            }
            int64_t ccumc[3] = {0, 0, 0};
            auto& lv = std::get<1>(cls3ret);
            for(size_t i = 0; i < varids.size(); ++i) ccumc[lv[varids[i]]] += bam_aux2i(bam_aux_get(vars[i], "CC"));
            for(int i = 0; i < 3; ++i) jsnstr.append("'cs: " + std::to_string(ccumc[i]) + "',");
            jsnstr.append("],\n");
            jsnstr.append("  name: 'k3',\n");
            jsnstr.append("  type: 'scatter',\n");
            jsnstr.append("  domain: {row: 0, col: 1,},\n");
            jsnstr.append("  marker: {\n");
            jsnstr.append("    color: [");
            for(size_t i = 0; i < lv.size(); ++i) jsnstr.append("'" + optcol[lv[i]] + "',");
            for(int i = 0; i < 3; ++i) jsnstr.append("'black',");
            jsnstr.append("],\n");
            jsnstr.append("    symbol: [");
            for(size_t i = 0; i < lv.size(); ++i) jsnstr.append("'" + optpch[lv[i]] + "',");
            for(int i = 0; i < 3; ++i) jsnstr.append("'" + optpch[3] + "',");
            jsnstr.append("],\n");
            jsnstr.append("  },\n");
            jsnstr.append("  mode: 'markers',\n");
            jsnstr.append("  xaxis: 'x2',\n");
            jsnstr.append("  yaxis: 'y2',\n");
            jsnstr.append("};\n");
        }
        jsnstr.append("var data = [k2c, k3c];\n");
        // layout
        {
            jsnstr.append("var layout={\n");
            jsnstr.append("  title: \"clustering of seuqnces, silhouette(k2:");
            jsnstr.append(std::to_string(kcoef2) + ",");
            jsnstr.append("k3:");
            jsnstr.append(std::to_string(kcoef3) + ")");
            jsnstr.append("\",\n");
            jsnstr.append("  grid: {rows: 1, columns: 2, pattern: 'independent', roworder: 'top to bottom'},\n");
            jsnstr.append("  xaxis1: {range:[-1,101], title:\"scaled alignment score\",},\n");
            jsnstr.append("  yaxis1: {range:[-1,101], title:\"scaled variant count\",},\n");
            jsnstr.append("  xaxis2: {range:[-1,101], title:\"scaled alignment score\",},\n");
            jsnstr.append("  yaxis2: {range:[-1,101], title:\"scaled variant count\",},\n");
            jsnstr.append("  showlegend: true,\n");
            jsnstr.append("  width: 1280,\n");
            jsnstr.append("  height: 600,\n");
            jsnstr.append("};\n");
        }
        // config
        {
            jsnstr.append("var config = {\n");
            jsnstr.append("  toImageButtonOptions: {\n");
            jsnstr.append("    format: 'svg',\n");
            jsnstr.append("     filename: '" + divName + "_of_" + std::string(name) + "',\n");
            jsnstr.append("     height: 600,\n");
            jsnstr.append("     width: 1280,\n");
            jsnstr.append("     scale: 1,\n");
            jsnstr.append("  }\n");
            jsnstr.append("};\n");
        }
        jsnstr.append("Plotly.newPlot('plot_" + divName + "', data, layout, config);\n");
        ksprintf(s, "%s", jsnstr.c_str());
        ksprintf(s, "</script>\n");
    }
}

void GEDetector::af_filter(){
    if(vars.empty()) return;
    int64_t minsrpt = totcnt * minaf;
    uint8_t* vsd = NULL;
    uint8_t* ccd = NULL;
    int64_t vsv = 0, ccv = 0;
    for(size_t i = 0; i < vars.size(); ++i){
        vsd = bam_aux_get(vars[i], "VS");
        vsv = bam_aux2i(vsd);
        ccd = bam_aux_get(vars[i], "CC");
        ccv = bam_aux2i(ccd);
        if(vsv & dropmask) continue;
        if(vsv & KSW_FREPSEQR){
            if(vsv & KSW_FMAYVARSEQ){
                if(ccv < minsrpt){
                    vsv |= KSW_FLOWFREQ;
                    stat_var(vars[i], false);
                    edicnt -= ccv;
                    othcnt += ccv;
                    if(vsv & KSW_FRECANYHIT) recpct -= ccv;
                    if(vsv & KSW_FRECALLHIT) reccnt -= ccv;
                    if(vsv & KSW_FRECEXACT) recent -= ccv;
                    bam_aux_update_int(vars[i], "VS", vsv);
                }
            }
        }
    }
    totcnt = refcnt + edicnt + othcnt;
    edieff = recpef = receff = reeeff = muteff = .0;
    if(totcnt > 0){
        edieff = (double)edicnt/(double)totcnt;
        recpef = (double)recpct/(double)totcnt;
        receff = (double)reccnt/(double)totcnt;
        reeeff = (double)recent/(double)totcnt;
        muteff = (double)mutcnt/(double)totcnt;
    }
}

int64_t GEDetector::get_ttr(const char* bam){
    int64_t ret = 0;
    samFile* fp = sam_open(bam, "r");
    bam_hdr_t* h = sam_hdr_read(fp);
    hts_idx_t* idx = sam_index_load(fp, bam);
    int tid = sam_hdr_name2tid(h, name);
    hts_itr_t* itr = sam_itr_queryi(idx, tid, 0, h->target_len[tid]);
    bam1_t* b = bam_init1();
    uint8_t* vsd = NULL;
    int64_t vsv = 0;
    while(sam_itr_next(fp, itr, b) >= 0){
        vsd = bam_aux_get(b, "VS");
        vsv = bam_aux2i(vsd);
        if(vsv & KSW_FREPSEQR){
            if(!(vsv & dropmask)) ret += bam_aux2i(bam_aux_get(b, "CC"));
        }
    }
    sam_close(fp);
    bam_destroy1(b);
    bam_hdr_destroy(h);
    hts_idx_destroy(idx);
    hts_itr_destroy(itr);
    return ret;
}

void GEDetector::cal_edieff(){
    uint8_t* vsd = NULL;
    uint8_t* ccd = NULL;
    int64_t vsv = 0, ccv = 0;
    for(size_t i = 0; i < vars.size(); ++i){
        vsd = bam_aux_get(vars[i], "VS");
        vsv = bam_aux2i(vsd);
        ccd = bam_aux_get(vars[i], "CC");
        ccv = bam_aux2i(ccd);
        if(ccv < mincnt){
            vsv |= KSW_FLOWSURPT;
            bam_aux_update_int(vars[i], "VS", vsv);
        }
        if(vsv & KSW_FREPSEQR){
            if(vsv & dropmask){
                dropcnt += ccv;
                if(vsv & KSW_FHIERRLOWQ){
                    droprcnt[dropschem.KSW_FHIERRLOWQ_IDX] += ccv;
                }else if(vsv & KSW_FPRIMDIMER){
                    droprcnt[dropschem.KSW_FPRIMDIMER_IDX] += ccv;
                }else if(vsv & KSW_FMISPRIMINN){
                    droprcnt[dropschem.KSW_FMISPRIMINN_IDX] += ccv;
                }else if(vsv & KSW_FMISPRIMOUT){
                    droprcnt[dropschem.KSW_FMISPRIMOUT_IDX] += ccv;
                }else if(vsv & KSW_FLOWSURPT){
                    droprcnt[dropschem.KSW_FLOWSURPT_IDX] += ccv;
                }else if(vsv & KSW_FMANYVARS){
                    droprcnt[dropschem.KSW_FMANYVARS_IDX] += ccv;
                }
                continue;
            }
            if(vsv & KSW_FMAYVARSEQ){
                edicnt += ccv;
                stat_var(vars[i], true);
            }else if(vsv & KSW_FREFTYPE){
                refcnt += ccv;
                stat_var(vars[i], true);
            }else{
                othcnt += ccv;
                stat_var(vars[i], true);
            }
            if(vsv & KSW_FVARINSEQ) mutcnt += ccv;
            if(vsv & KSW_FRECANYHIT) recpct += ccv;
            if(vsv & KSW_FRECALLHIT) reccnt += ccv;
            if(vsv & KSW_FRECEXACT) recent += ccv;
        }
    }
    totcnt = refcnt + edicnt + othcnt;
    edieff = recpef = receff = reeeff =  muteff = .0;
    if(totcnt > 0){
        edieff = (double)edicnt/(double)totcnt;
        recpef = (double)recpct/(double)totcnt;
        receff = (double)reccnt/(double)totcnt;
        reeeff = (double)recent/(double)totcnt;
        muteff = (double)mutcnt/(double)totcnt;
    }
    allcnt = refcnt + edicnt + othcnt + dropcnt;
}

void GEDetector::cal_scedff(){
    uint8_t* vsd = NULL;
    uint8_t* ccd = NULL;
    uint8_t* cbd = NULL;
    int64_t vsv = 0, ccv = 0;
    char* cbv = NULL;
    sc_ge_t* pscg = NULL;
    for(size_t i = 0; i < vars.size(); ++i){
        vsd = bam_aux_get(vars[i], "VS");
        vsv = bam_aux2i(vsd);
        ccd = bam_aux_get(vars[i], "CC");
        ccv = bam_aux2i(ccd);
        cbd = bam_aux_get(vars[i], CELL_BARCODE_ID_TAG);
        if(cbd) cbv = bam_aux2Z(cbd);
        else cbv = NULL;
        if(!cbv) continue;
        char* pcbv = strdup(cbv);
        int smax = 0, *sofs = 0;
        int snf = ksplit_core(pcbv, '|', &smax, &sofs);
        if(ccv < mincnt){
            vsv |= KSW_FLOWSURPT;
            bam_aux_update_int(vars[i], "VS", vsv);
        }
        for(int sni = 0; sni < snf; ++sni){
            char* psstr = pcbv+sofs[sni];
            char* pcstr = strchr(psstr, '-');
            int cbcnt = atoi(pcstr+1);
            *pcstr = '\0';
            auto iter = scretm.find(psstr);
            if(iter == scretm.end()){
                pscg = (sc_ge_t*)calloc(1, sizeof(sc_ge_t));
                scretm[psstr] = pscg;
            }else{
                pscg = iter->second;
            }
            if(vsv & KSW_FREPSEQR){
                if(vsv & dropmask){
                    pscg->droped += cbcnt;
                    continue;
                }
                if(vsv & KSW_FMAYVARSEQ){
                    pscg->edicnt += cbcnt;
                }else if(vsv & KSW_FREFTYPE){
                    pscg->refcnt += cbcnt;
                }else{
                    pscg->othcnt += cbcnt;
                }
                if(vsv & KSW_FRECANYHIT){ pscg->recpct += cbcnt; }
                if(vsv & KSW_FRECALLHIT){ pscg->reccnt += cbcnt; }
                if(vsv & KSW_FRECEXACT){ pscg->recent += cbcnt; }
            }
        }
        free(pcbv);
        free(sofs);
    }
}

void GEDetector::summary(){
    if(summarized) return;
    // cnt
    tot_cnt_dist->summary();
    ins_cnt_dist->summary();
    del_cnt_dist->summary();
    din_cnt_dist->summary();
    snv_cnt_dist->summary();
    scl_cnt_dist->summary();
    // len
    ins_len_dist->summary();
    del_len_dist->summary();
    din_len_dist->summary();
    scl_len_dist->summary();
    // pos
    ins_pos_dist->summary();
    del_pos_dist->summary();
    din_pos_dist->summary();
    snv_pos_dist->summary();
    scl_pos_dist->summary();
    // refrate
    for(int i = 0; i < rlen; ++i){
        int64_t ttb = 0;
        for(int j = 0; j < 5; ++j) ttb += nccnt[j][i];
        int64_t rtb = nccnt[nuc_to_3bit[(int)ref[i]]][i];
        refcov[i] = rtb;
        altcov[i] = ttb - rtb;
        if(ttb) refrate[i] = (double)rtb/(double)ttb;
        else refrate[i] = .0;
    }
    summarized = true;
}

void GEDetector::reportJSON(kstring_t* s, const char* dh, const char* dm){
    if(!summarized) summary();
    // summary
    ksprintf(s, "%s\"EditSummary\": {\n", dh);
    ksprintf(s, "%s%s\"EstimatedCleavagePos\": %d,\n", dh, dm, clspos);
    ksprintf(s, "%s%s\"QualifiedSeqCount\": %lld,\n", dh, dm, totcnt);
    ksprintf(s, "%s%s\"RefSeqCount\": %lld,\n", dh, dm, refcnt);
    ksprintf(s, "%s%s\"EditSeqCount\": %lld,\n", dh, dm, edicnt);
    ksprintf(s, "%s%s\"RecombSeqCount\": %lld,\n", dh, dm, reccnt);
    ksprintf(s, "%s%s\"RecombExactSeqCount\": %lld,\n", dh, dm, recent);
    ksprintf(s, "%s%s\"RecombAnySeqCount\": %lld,\n", dh, dm, recpct);
    ksprintf(s, "%s%s\"OtherSeqCount\": %lld,\n", dh, dm, othcnt);
    ksprintf(s, "%s%s\"DropSeqCount\": %lld,\n", dh, dm, dropcnt);
    ksprintf(s, "%s%s\"EditRate\": %lf,\n", dh, dm, edieff);
    ksprintf(s, "%s%s\"RecombRate\": %lf,\n", dh, dm, receff);
    ksprintf(s, "%s%s\"RecombExactRate\": %lf,\n", dh, dm, reeeff);
    ksprintf(s, "%s%s\"RecombAnyRate\": %lf,\n", dh, dm, recpef);
    ksprintf(s, "%s%s\"MutsRate\": %lf\n", dh, dm, muteff);
    ksprintf(s, "%s},\n", dh);
    // dropcnt
    ksprintf(s, "%s\"QCFailDetails\": {\n", dh);
    for(int i = 0; i < TOTAL_DROP_REASON; ++i){
        ksprintf(s, "%s%s\"%s\":[%lld, %f],\n", dh, dm, drop_reasons_str[i], droprcnt[i], dropcnt > 0 ? (double)droprcnt[i]/(double)dropcnt : .0);
    }
    s->s[s->l-2] = '\n';
    s->l -= 1;
    ksprintf(s, "%s},\n", dh);
    // count
    tot_cnt_dist->reportJSON(s, dh, dm, "TotalVarCount");
    ins_cnt_dist->reportJSON(s, dh, dm, "InsertionCount");
    del_cnt_dist->reportJSON(s, dh, dm, "DeletionCount");
    din_cnt_dist->reportJSON(s, dh, dm, "DelinsCount");
    snv_cnt_dist->reportJSON(s, dh, dm, "SNVCount");
    scl_cnt_dist->reportJSON(s, dh, dm, "SoftClipCount");
    // len
    ins_len_dist->reportJSON(s, dh, dm, "InsertionLength");
    del_len_dist->reportJSON(s, dh, dm, "DeletionLength");
    din_len_dist->reportJSON(s, dh, dm, "DelinsLength");
    scl_len_dist->reportJSON(s, dh, dm, "SoftClipLength");
    // pos
    ins_pos_dist->reportJSON(s, dh, dm, "InsertionPosition");
    del_pos_dist->reportJSON(s, dh, dm, "DeletionPosition");
    din_pos_dist->reportJSON(s, dh, dm, "DelinsPosition");
    snv_pos_dist->reportJSON(s, dh, dm, "SNVPosition");
    scl_pos_dist->reportJSON(s, dh, dm, "SoftClipPosition");
    int tmplv = 0;
    // lps -> ins
    ksprintf(s, "%s\"PoswiseMeanInsertionLen\": [", dh);
    for(int c = 0; c < ins_lps_dist->vcdist_aloc; ++c){
        if(ins_pos_dist->vcdist_array[c]) tmplv = (int)((double)ins_lps_dist->vcdist_array[c]/(double)ins_pos_dist->vcdist_array[c]);
        else tmplv = 0;
        ksprintf(s, "%d,", tmplv);
    }
    s->s[s->l-1] = ']';
    ksprintf(s, ",\n");
    // lps -> del
    ksprintf(s, "%s\"PoswiseMeanDeletionLen\": [", dh);
    for(int c = 0; c < del_lps_dist->vcdist_aloc; ++c){
        if(del_pos_dist->vcdist_array[c]) tmplv = (int)((double)del_lps_dist->vcdist_array[c]/(double)del_pos_dist->vcdist_array[c]);
        else tmplv = 0;
        ksprintf(s, "%d,", tmplv);
    }
    s->s[s->l-1] = ']';
    ksprintf(s, ",\n");
    // lps -> delins
    ksprintf(s, "%s\"PoswiseMeanDelinsLen\": [", dh);
    for(int c = 0; c < din_lps_dist->vcdist_aloc; ++c){
        if(din_pos_dist->vcdist_array[c]) tmplv = (int)((double)din_lps_dist->vcdist_array[c]/(double)din_pos_dist->vcdist_array[c]);
        else tmplv = 0;
        ksprintf(s, "%d,", tmplv);
    }
    s->s[s->l-1] = ']';
    ksprintf(s, ",\n");
    // ACGTN
    ksprintf(s, "%s\"ACGTNDist\": [", dh);
    for(int c = 0; c < rlen; ++c){
        ksprintf(s, "[%lld,%lld,%lld,%lld,%lld],", nccnt[0][c], nccnt[1][c], nccnt[2][c], nccnt[3][c], nccnt[4][c]);
    }
    s->s[s->l-1] = ']';
    ksprintf(s, ",\n");
    // RefRate
    ksprintf(s, "%s\"RefRate\": [", dh);
    for(int c = 0; c < rlen; ++c){
        ksprintf(s, "%lf,", refrate[c]);
    }
    s->s[s->l-1] = ']';
    ksprintf(s, ",\n");
    // REFIDX
    ksprintf(s, "%s\"RefNCIdx\": [", dh);
    for(int c = 0; c < rlen; ++c){
         ksprintf(s, "%d,", seq_nt16_int[seq_nt16_table[(int)ref[c]]]);
    }
    s->s[s->l-1] = ']';
}

void GEDetector::tsvHead(kstring_t* s){
    ksprintf(s, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s", "EstimatedCleavagePos", "QualifiedSeqCount", "RefSeqCount", "EditSeqCount", "RecombSeqCount", "RecombExactSeqCount",  "RecombAnySeqCount","OtherSeqCount", "DropSeqCount", "EditRate", "RecombRate", "RecombExactRate",  "RecombAnyRate", "MutsRate");

}

void GEDetector::tsvBody(kstring_t* s){
    if(!summarized) summary();
    ksprintf(s, "%d\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lld\t%lf\t%lf\t%lf\t%lf\t%lf", clspos, totcnt, refcnt, edicnt, reccnt, recent, recpct, othcnt, dropcnt, edieff, receff, reeeff, recpef, muteff); 
}
