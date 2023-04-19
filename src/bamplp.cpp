#include "bamplp.h"

void mplp_aux_destroy(mplp_aux_t* ma){
    if(!ma) return;
    if(ma->fp){ sam_close(ma->fp); };
    if(ma->hdr){ bam_hdr_destroy(ma->hdr); }
    if(ma->sidx){ hts_idx_destroy(ma->sidx); }
    if(ma->iter){ hts_itr_destroy(ma->iter); }
};

void mplp_pileup_init(mplp_pileup_t* mp, int _n){
    mp->n = _n;
    mp->n_plp = (int*)calloc(mp->n, sizeof(int));
    mp->plp = (const bam_pileup1_t**)calloc(mp->n, sizeof(bam_pileup1_t*));
}

void mplp_pileup_destroy(mplp_pileup_t* mp){
    if(!mp) return;
    if(mp->n_plp) free(mp->n_plp);
    if(mp->plp) free(mp->plp);
}

int mplp_func(void *data, bam1_t* b){
    mplp_aux_t* ma = (mplp_aux_t*)data;
    int ret = 0, skip = 0;
    do{
        if(ma->iter){// from file
            ret = sam_itr_next(ma->fp, ma->iter, b);
            if(ret < 0) break;
        }else{
            if(ma->idx >= ma->bams->size()){
                ret = -1;
                break;
            }else{
                if(bam_copy1(b, ma->bams->at(ma->idx)) == NULL){
                    ret = -1;
                    ++ma->idx;
                    break;
                }else{
                    ret = 0;
                    ++ma->idx;
                }
            }
        }
        if((b->core.flag & BAM_FUNMAP) || (bam_aux2i(bam_aux_get(b, "VS")) & ma->dropmask)) skip = 1;
        else skip = 0;
    }while(skip);
    return ret;
}

int pileup_constructor(void *data, const bam1_t *b, bam_pileup_cd *cd){
    bam_info_aux_t *ba = (bam_info_aux_t*)calloc(1, sizeof(bam_info_aux_t));
    ba->bi = NULL;
#ifdef BAM_CNT_ONE
    ba->cc = 1;
#else
    ba->cc = bam_aux2i(bam_aux_get(b, "CC"));
#endif
    uint8_t* bid = bam_aux_get(b, CELL_BARCODE_ID_TAG);
    if(bid){
        ba->bi = strdup(bam_aux2Z(bid));
        ba->nfs = ksplit_core(ba->bi, '|', &ba->mfs, &ba->ofs);
        ba->cbc = (int*)calloc(ba->nfs, sizeof(int));
        for(int sni = 0; sni < ba->nfs; ++sni){
            char* pscb = ba->bi+ba->ofs[sni];
            char* pccb = strchr(pscb, '-');
            ba->cbc[sni] = atoi(pccb+1);
            *pccb = '\0';
        }
    }
    cd->p = ba;
    return 0;
}

int pileup_destructor(void *data, const bam1_t *b, bam_pileup_cd *cd){
    bam_info_aux_t *ba = (bam_info_aux_t*)cd->p;
    if(ba->bi) free(ba->bi);
    if(ba->ofs) free(ba->ofs);
    if(ba->cbc) free(ba->cbc);
    free(ba);
    return 0;
}

void plp_call(int tid, int pos, const bam_pileup1_t* plp, int n_plp, call_ret_t* cr, uint16_t vtag){
    cr->init();
    cr->tid = tid;
    cr->pos = pos;
    int isdi = 0;
    for(int i = 0; i < n_plp; ++i){
        const bam_pileup1_t *p = plp + i;
        bam_info_aux_t* ba = (bam_info_aux_t*)p->cd.p;
        isdi = 0;
        if(p->indel > 0){// insertion after this pos
            cr->vm |= INS_CALLED;
            std::string ins;
            for(int j = 0; j <= p->indel; ++j){
                ins.append(1, seq_nt16_str[bam_seqi(bam_get_seq(p->b), p->qpos+j)]);
            }
            auto iter = cr->ins.find(ins);
            if(iter == cr->ins.end()){
                cr->ins[ins] = ba->cc;
            }else{
                iter->second += ba->cc;
            }
            if(ba->bi && (vtag & GEVAR_INS)){
                auto iter = cr->sc_ins.find(ins);
                if(iter == cr->sc_ins.end()){
                    cr->sc_ins[ins] = std::unordered_map<std::string, int>();
                    for(int sni = 0; sni < ba->nfs; ++sni){
                        cr->sc_ins[ins][ba->bi+ba->ofs[sni]] = ba->cbc[sni];
                    }
                }else{
                    for(int sni = 0; sni < ba->nfs; ++sni){
                        iter->second[ba->bi+ba->ofs[sni]] += ba->cbc[sni];
                    }
                }
                cr->sc_called[SC_INS_CALLED] = 1;
            }
        }else if(p->indel < 0){ // deletion beginning from this pos
            cr->vm |= DEL_CALLED;
            auto iter = cr->del.find(-p->indel);
            if(iter == cr->del.end()){
                cr->del[-p->indel] = ba->cc;
            }else{
                iter->second += ba->cc;
            }
            if(ba->bi && (vtag & GEVAR_DEL)){
                auto iter = cr->sc_del.find(-p->indel);
                if(iter == cr->sc_del.end()){
                    cr->sc_del[-p->indel] = std::unordered_map<std::string, int>();
                    for(int sni = 0; sni < ba->nfs; ++sni){
                        cr->sc_del[-p->indel][ba->bi+ba->ofs[sni]] = ba->cbc[sni];
                    }
                }else{
                    for(int sni = 0; sni < ba->nfs; ++sni){
                        iter->second[ba->bi+ba->ofs[sni]] += ba->cbc[sni];
                    }
                }
                cr->sc_called[SC_DEL_CALLED] = 1;
            }
        }else{ // check delins
            int opint = bam_cigar_op(bam_get_cigar(p->b)[p->cigar_ind]);
            int oplen = bam_cigar_oplen(bam_get_cigar(p->b)[p->cigar_ind]);
            if(opint == BAM_CDIFF && oplen > 1){ // delins, make sure only call once
                isdi = 1;
                int qpos = 0;
                for(int icig=0; icig<p->cigar_ind; icig++){
                    int opi = bam_cigar_op(bam_get_cigar(p->b)[icig]);
                    int opl = bam_cigar_oplen(bam_get_cigar(p->b)[icig]);
                    if(opi == BAM_CMATCH || opi == BAM_CEQUAL || 
                       opi == BAM_CDIFF || opi == BAM_CINS || opi == BAM_CSOFT_CLIP){
                        qpos += opl;
                    }
                }
                if(qpos == p->qpos){// first time
                    cr->vm |= DIN_CALLED;
                    std::string dis;
                    for(int j = 0; j < oplen; ++j){
                        dis.append(1, seq_nt16_str[bam_seqi(bam_get_seq(p->b), p->qpos+j)]);
                    }
                    auto iter = cr->din.find(dis);
                    if(iter == cr->din.end()){
                        cr->din[dis] = ba->cc;
                    }else{
                        iter->second += ba->cc;
                    }
                    if(ba->bi && (vtag & GEVAR_DIN)){
                        auto iter = cr->sc_din.find(dis);
                        if(iter == cr->sc_din.end()){
                            cr->sc_ins[dis] = std::unordered_map<std::string, int>();
                            for(int sni = 0; sni < ba->nfs; ++sni){
                                cr->sc_ins[dis][ba->bi+ba->ofs[sni]] = ba->cbc[sni];
                            }
                        }else{
                            for(int sni = 0; sni < ba->nfs; ++sni){
                                iter->second[ba->bi+ba->ofs[sni]] += ba->cbc[sni];
                            }
                        }
                        cr->sc_called[SC_DIN_CALLED] = 1;
                    }
                }
            }
        }
        if(!(p->is_del || p->is_refskip) && bam_get_qual(p->b)[p->qpos] > 0) cr->depth += ba->cc;
        if(p->is_del || p->is_refskip || isdi) continue; // non snv
        if(bam_get_qual(p->b)[p->qpos] > 0){
            cr->vm |= SNP_CALLED;
            int bid = seq_nt16_int[bam_seqi(bam_get_seq(p->b), p->qpos)];
            cr->snp[bid] += ba->cc;
            if(ba->bi && (vtag & GEVAR_SNV)){
                if(cr->sc_snp.empty()) cr->sc_snp.resize(5, std::unordered_map<std::string, int>());
                for(int sni = 0; sni < ba->nfs; ++sni){
                    auto xsiter = cr->sc_snp[bid].find(ba->bi+ba->ofs[sni]);
                    if(xsiter == cr->sc_snp[bid].end()){
                        cr->sc_snp[bid][ba->bi+ba->ofs[sni]] = ba->cbc[sni];
                    }else{
                        xsiter->second += ba->cbc[sni];
                    }
                }
                cr->sc_called[SC_SNP_CALLED] = 1;
            }
        }
    }
}

int call2bcf1(call_ret_t* cr, bcf_hdr_t* hdr, std::vector<bcf1_t*>& recs, char* ref){
    if(!(cr->vm & VAR_CALLED)) return -1; // no calling
    if(cr->vm & SNP_CALLED){
        bcf1_t* rec = bcf_init1();
        rec->rid = cr->tid; // CHROM
        rec->pos = cr->pos; // POS
        rec->qual = .0; // QUAL
        cr->s.l = 0;
        char rch = toupper(ref[cr->pos]);
        for(int i = 0; i < 4; ++i){
            if("ACGTN"[i] == rch){
                cr->ridx = i;
                break;
            }
        }
        kputc(rch, &cr->s); // REF
        int val = 0;
        int nalt = 0;
        for(int i = 0; i < 4; ++i){ // ALT
            if(cr->snp[i] && i != cr->ridx){
                cr->ad[nalt] = cr->snp[i];
                cr->af[nalt] = (float)(cr->snp[i])/(float)(cr->totcnt);
                kputc(',', &cr->s);
                kputc("ACGT"[i], &cr->s);
                ++nalt;
            }
        }
        bcf_update_alleles_str(hdr, rec, cr->s.s);
        if(nalt) val = SNP_CALLED;
        // INFO
        bcf_update_info_int32(hdr, rec, "VT", &val, 1);
        // FORMAT
        rec->n_sample = 1;
        bcf_update_format_float(hdr, rec, "AF", cr->af, rec->n_sample*nalt); // AF
        bcf_update_format_int32(hdr, rec, "AD", cr->ad, rec->n_sample*nalt); // AD
        bcf_update_format_int32(hdr, rec, "DP", &cr->depth, rec->n_sample); // DP
        if(cr->sc_called[SC_SNP_CALLED]){
            cr->ts->l = 0;
            for(int i = 0; i < 4; ++i){
                for(auto& p: cr->sc_snp[i]){
                    ksprintf(cr->ts, "%s%s-%d|", SC_ID_PREFIX, p.first.c_str(), p.second);
                }
                kputc(',', cr->ts);
            }
            const char* scs = cr->ts->s;
            bcf_update_format_string(hdr, rec, "SC", &scs, rec->n_sample); // SC
        }
        recs.push_back(rec);
    }
    if(cr->vm & DEL_CALLED){
        for(auto& e: cr->del){
            int64_t dl = e.first;
            int64_t dc = e.second;
            bcf1_t* rec = bcf_init1();
            rec->rid = cr->tid; // CHROM
            rec->pos = cr->pos; // POS
            rec->qual = .0; // QUAL
            cr->s.l = 0;
            // REF
            kputc(ref[cr->pos], &cr->s);
            for(int j = 1; j <= dl; ++j) kputc(ref[cr->pos+j], &cr->s);
            // ALT
            kputc(',', &cr->s); kputc(ref[cr->pos], &cr->s);
            bcf_update_alleles_str(hdr, rec, cr->s.s);
            // INFO
            int val = DEL_CALLED; 
            bcf_update_info_int32(hdr, rec, "VT", &val, 1);
            // FORMAT
            rec->n_sample = 1;
            cr->ad[0] = dc;
            cr->af[0] = (double)dc/(double)cr->totcnt;
            bcf_update_format_float(hdr, rec, "AF", cr->af, rec->n_sample); // AF
            bcf_update_format_int32(hdr, rec, "AD", cr->ad, rec->n_sample); // AD
            bcf_update_format_int32(hdr, rec, "DP", &cr->depth, rec->n_sample); // DP
            if(cr->sc_called[SC_DEL_CALLED]){
                cr->ts->l = 0;
                for(auto& p: cr->sc_del[dl]){
                    ksprintf(cr->ts, "%s%s-%d|", SC_ID_PREFIX, p.first.c_str(), p.second);
                }
                const char* scs = cr->ts->s;
                bcf_update_format_string(hdr, rec, "SC", &scs, rec->n_sample); // SC
            }
            recs.push_back(rec);
        }
    }
    if(cr->vm & INS_CALLED){
        for(auto& e: cr->ins){
            bcf1_t* rec = bcf_init1();
            rec->rid = cr->tid; // CHROM
            rec->pos = cr->pos; // POS
            rec->qual = .0; // QUAL
            cr->s.l = 0;
            // REF
            kputc(ref[cr->pos], &cr->s);
            // ALT
            kputc(',', &cr->s); 
            kputs(e.first.c_str(), &cr->s);
            bcf_update_alleles_str(hdr, rec, cr->s.s);
            // INFO
            int val = INS_CALLED; 
            bcf_update_info_int32(hdr, rec, "VT", &val, 1);
            // FORMAT
            rec->n_sample = 1;
            cr->ad[0] = e.second;
            cr->af[0] = (double)e.second/(double)cr->totcnt;
            bcf_update_format_float(hdr, rec, "AF", cr->af, rec->n_sample); // AF
            bcf_update_format_int32(hdr, rec, "AD", cr->ad, rec->n_sample); // AD
            bcf_update_format_int32(hdr, rec, "DP", &cr->depth, rec->n_sample); // DP
            if(cr->sc_called[SC_INS_CALLED]){
                cr->ts->l = 0;
                for(auto& p: cr->sc_ins[e.first]){
                    ksprintf(cr->ts, "%s%s-%d|", SC_ID_PREFIX, p.first.c_str(), p.second);
                }
                const char* scs = cr->ts->s;
                bcf_update_format_string(hdr, rec, "SC", &scs, rec->n_sample); // SC
            }
            recs.push_back(rec);
        }
    }
    if(cr->vm & DIN_CALLED){
        for(auto& e: cr->din){
            bcf1_t* rec = bcf_init1();
            rec->rid = cr->tid; // CHROM
            rec->pos = cr->pos; // POS
            rec->qual = .0; // QUAL
            cr->s.l = 0;
            // REF
            for(size_t j = 0; j < e.first.size(); ++j) kputc(ref[cr->pos+j], &cr->s);
            // ALT
            kputc(',', &cr->s);
            kputs(e.first.c_str(), &cr->s);
            bcf_update_alleles_str(hdr, rec, cr->s.s);
            // INFO
            int val = DIN_CALLED; 
            bcf_update_info_int32(hdr, rec, "VT", &val, 1);
            // FORMAT
            rec->n_sample = 1;
            cr->ad[0] = e.second;
            cr->af[0] = (double)e.second/(double)cr->totcnt;
            bcf_update_format_float(hdr, rec, "AF", cr->af, rec->n_sample); // AF
            bcf_update_format_int32(hdr, rec, "AD", cr->ad, rec->n_sample); // AD
            bcf_update_format_int32(hdr, rec, "DP", &cr->depth, rec->n_sample); // DP
            if(cr->sc_called[SC_DIN_CALLED]){
                cr->ts->l = 0;
                for(auto& p: cr->sc_din[e.first]){
                    ksprintf(cr->ts, "%s%s-%d|", SC_ID_PREFIX, p.first.c_str(), p.second);
                }
                const char* scs = cr->ts->s;
                bcf_update_format_string(hdr, rec, "SC", &scs, rec->n_sample); // SC
            }
            recs.push_back(rec);
        }
    }
    return 0;
}

int call2bcf2(call_ret_t* cr1, call_ret_t* cr2, bcf_hdr_t* hdr, std::vector<bcf1_t*>& recs, char* ref){
    if(!(cr1->vm & VAR_CALLED) && !(cr2->vm & VAR_CALLED)) return -1; // no calling
    if((cr1->vm & SNP_CALLED) || (cr2->vm & SNP_CALLED)){
        bcf1_t* rec = bcf_init1();
        rec->rid = cr1->tid; // CHROM
        rec->pos = cr1->pos; // POS
        rec->qual = .0; // QUAL
        cr1->s.l = 0;
        char rch = toupper(ref[cr1->pos]);
        for(int i = 0; i < 4; ++i){
            if("ACGTN"[i] == rch){
                cr1->ridx = i;
                cr2->ridx = i;
                break;
            }
        }
        kputc(rch, &cr1->s); // REF
        int nalt = 0;
        int val = 0;
        int8_t acgtn[5] = {0, 0, 0, 0, 0};
        for(int i = 0; i < 4; ++i){ // ALT
            if(i != cr1->ridx && (cr1->snp[i] || cr2->snp[i])){
                acgtn[i] = 1;
                kputc(',', &cr1->s);
                kputc("ACGT"[i], &cr1->s);
                ++nalt;
            }
        }
        if(nalt) val = SNP_CALLED;
        int ncct = 0;
        for(int i = 0; i < 4; ++i){
            if(acgtn[i]){
                cr1->ad[ncct] = cr1->snp[i];
                cr1->ad[ncct+nalt] = cr2->snp[i];
                cr1->af[ncct] = (double)cr1->ad[ncct]/(double)cr1->totcnt;
                cr1->af[ncct+nalt] = (double)cr1->ad[ncct+nalt]/(double)cr2->totcnt;
                ++ncct;
            }
        }
        bcf_update_alleles_str(hdr, rec, cr1->s.s);
        // INFO
        bcf_update_info_int32(hdr, rec, "VT", &val, 1);
        // FORMAT
        rec->n_sample = 2;
        cr1->dp[0] = cr1->depth; cr1->dp[1] = cr2->depth;
        cr1->tr[0] = cr1->totcnt; cr1->tr[1] = cr2->totcnt;
        bcf_update_format_float(hdr, rec, "AF", cr1->af, rec->n_sample*nalt); // AF
        bcf_update_format_int32(hdr, rec, "AD", cr1->ad, rec->n_sample*nalt); // AD
        bcf_update_format_int32(hdr, rec, "DP", cr1->dp, rec->n_sample); // DP
        if(cr1->sc_called[SC_SNP_CALLED] && cr2->sc_called[SC_SNP_CALLED]){
            cr1->ts->l = 0; cr2->ts->l = 0;
            for(int i = 0; i < 4; ++i){
                for(auto& p: cr1->sc_snp[i]){
                    ksprintf(cr1->ts, "%s%s-%d|", SC_ID_PREFIX, p.first.c_str(), p.second);
                }
                kputc(',', cr1->ts);
            }
            for(int i = 0; i < 4; ++i){
                for(auto& p: cr2->sc_snp[i]){
                    ksprintf(cr2->ts, "%s%s-%d|", SC_ID_PREFIX, p.first.c_str(), p.second);
                }
                kputc(',', cr2->ts);
            }
            const char* scs[2] = {(const char*)cr1->ts->s, (const char*)cr2->ts->s};
            bcf_update_format_string(hdr, rec, "SC", scs, rec->n_sample); // SC
        }
        recs.push_back(rec);
    }
    if((cr1->vm & DEL_CALLED) || (cr2->vm & DEL_CALLED)){
        std::unordered_set<int64_t> alld;
        for(auto& e: cr1->del) alld.insert(e.first);
        for(auto& e: cr2->del) alld.insert(e.first);
        for(auto& e: alld){
            int64_t dl = e;
            bcf1_t* rec = bcf_init1();
            rec->rid = cr1->tid; // CHROM
            rec->pos = cr1->pos; // POS
            rec->qual = .0; // QUAL
            cr1->s.l = 0;
            // REF
            kputc(ref[cr1->pos], &cr1->s);
            for(int j = 1; j <= dl; ++j) kputc(ref[cr1->pos+j], &cr1->s);
            // ALT
            kputc(',', &cr1->s); kputc(ref[cr1->pos], &cr1->s);
            bcf_update_alleles_str(hdr, rec, cr1->s.s);
            // INFO
            int val = DEL_CALLED; 
            bcf_update_info_int32(hdr, rec, "VT", &val, 1);
            // FORMAT
            auto iter1 = cr1->del.find(dl);
            if(iter1 == cr1->del.end()){
                cr1->ad[0] = 0;
                cr1->af[0] = .0;
            }else{
                cr1->ad[0] = iter1->second;
                cr1->af[0] = (double)cr1->ad[0]/(double)cr1->totcnt;
            }
            auto iter2 = cr2->del.find(dl);
            if(iter2 == cr2->del.end()){
                cr1->ad[1] = 0;
                cr1->af[1] = .0;
            }else{
                cr1->ad[1] = iter2->second;
                cr1->af[1] = (double)cr1->ad[1]/(double)cr2->totcnt;
            }
            rec->n_sample = 2;
            cr1->dp[0] = cr1->depth; cr1->dp[1] = cr2->depth;
            cr1->tr[0] = cr1->totcnt; cr1->tr[1] = cr2->totcnt;
            bcf_update_format_float(hdr, rec, "AF", cr1->af, rec->n_sample); // AF
            bcf_update_format_int32(hdr, rec, "AD", cr1->ad, rec->n_sample); // AD
            bcf_update_format_int32(hdr, rec, "DP", cr1->dp, rec->n_sample); // DP
            if(cr1->sc_called[SC_DEL_CALLED] && cr2->sc_called[SC_DEL_CALLED]){
                cr1->ts->l = 0; cr2->ts->l = 0;
                for(auto& p: cr1->sc_del[dl]){
                    ksprintf(cr1->ts, "%s%s-%d|", SC_ID_PREFIX, p.first.c_str(), p.second);
                }
                for(auto& p: cr2->sc_del[dl]){
                    ksprintf(cr2->ts, "%s%s-%d|", SC_ID_PREFIX, p.first.c_str(), p.second);
                }
                const char* scs[2] = {(const char*)cr1->ts->s, (const char*)cr2->ts->s};
                bcf_update_format_string(hdr, rec, "SC", scs, rec->n_sample); // SC
            }
            recs.push_back(rec);
        }
    }
    if((cr1->vm & INS_CALLED) || (cr2->vm & INS_CALLED)){
        std::unordered_set<std::string> allins;
        for(auto& e: cr1->ins) allins.insert(e.first);
        for(auto& e: cr2->ins) allins.insert(e.first);
        for(auto& e: allins){
            bcf1_t* rec = bcf_init1();
            rec->rid = cr1->tid; // CHROM
            rec->pos = cr1->pos; // POS
            rec->qual = .0; // QUAL
            cr1->s.l = 0;
            // REF
            kputc(ref[cr1->pos], &cr1->s);
            // ALT
            kputc(',', &cr1->s); kputc(ref[cr1->pos], &cr1->s);
            kputs(e.c_str(), &cr1->s);
            bcf_update_alleles_str(hdr, rec, cr1->s.s);
            // INFO
            int val = INS_CALLED; 
            bcf_update_info_int32(hdr, rec, "VT", &val, 1);
            // FORMAT
            rec->n_sample = 2;
            auto iter1 = cr1->ins.find(e);
            if(iter1 == cr1->ins.end()){
                cr1->ad[0] = 0;
                cr1->af[0] = .0;
            }else{
                cr1->ad[0] = iter1->second;
                cr1->af[0] = (double)cr1->ad[0]/(double)cr1->totcnt;
            }
            auto iter2 = cr2->ins.find(e);
            if(iter2 == cr2->ins.end()){
                cr1->ad[1] = 0;
                cr1->af[1] = .0;
            }else{
                cr1->ad[1] = iter2->second;
                cr1->af[1] = (double)cr1->ad[1]/(double)cr1->totcnt;
            }
            cr1->tr[0] = cr1->totcnt; cr1->tr[1] = cr2->totcnt;
            bcf_update_format_float(hdr, rec, "AF", cr1->af, rec->n_sample); // AF
            bcf_update_format_int32(hdr, rec, "AD", cr1->ad, rec->n_sample); // AD
            bcf_update_format_int32(hdr, rec, "DP", cr1->dp, rec->n_sample); // DP
            if(cr1->sc_called[SC_INS_CALLED] && cr2->sc_called[SC_INS_CALLED]){
                cr1->ts->l = 0; cr2->ts->l = 0;
                for(auto& p: cr1->sc_ins[e]){
                    ksprintf(cr1->ts, "%s%s-%d|", SC_ID_PREFIX, p.first.c_str(), p.second);
                }
                for(auto& p: cr2->sc_ins[e]){
                    ksprintf(cr2->ts, "%s%s-%d|", SC_ID_PREFIX, p.first.c_str(), p.second);
                }
                const char* scs[2] = {(const char*)cr1->ts->s, (const char*)cr2->ts->s};
                bcf_update_format_string(hdr, rec, "SC", scs, rec->n_sample); // SC
            }
            recs.push_back(rec);
        }
    }
    if((cr1->vm & DIN_CALLED) || (cr2->vm & DIN_CALLED)){
        std::unordered_set<std::string> alldin;
        for(auto& e: cr1->din) alldin.insert(e.first);
        for(auto& e: cr2->din) alldin.insert(e.first);
        for(auto& e: alldin){
            bcf1_t* rec = bcf_init1();
            rec->rid = cr1->tid; // CHROM
            rec->pos = cr1->pos; // POS
            rec->qual = .0; // QUAL
            cr1->s.l = 0;
            // REF
            for(size_t j = 0; j < e.size(); ++j) kputc(ref[cr1->pos+j], &cr1->s);
            // ALT
            kputc(',', &cr1->s);
            kputs(e.c_str(), &cr1->s);
            bcf_update_alleles_str(hdr, rec, cr1->s.s);
            // INFO
            int val = DIN_CALLED; 
            bcf_update_info_int32(hdr, rec, "VT", &val, 1);
            // FORMAT
            rec->n_sample = 2;
            auto iter1 = cr1->din.find(e);
            if(iter1 == cr1->din.end()){
                cr1->ad[0] = 0;
                cr1->af[0] = .0;
            }else{
                cr1->ad[0] = iter1->second;
                cr1->af[0] = (double)cr1->ad[0]/(double)cr1->totcnt;
            }
            auto iter2 = cr2->din.find(e);
            if(iter2 == cr2->din.end()){
                cr1->ad[1] = 0;
                cr1->af[1] = .0;
            }else{
                cr1->ad[1] = iter2->second;
                cr1->af[1] = (double)cr1->ad[0]/(double)cr2->totcnt;
            }
            cr1->tr[0] = cr1->totcnt; cr1->tr[1] = cr2->totcnt;
            bcf_update_format_float(hdr, rec, "AF", cr1->af, rec->n_sample); // AF
            bcf_update_format_int32(hdr, rec, "AD", cr1->ad, rec->n_sample); // AD
            bcf_update_format_int32(hdr, rec, "DP", cr1->dp, rec->n_sample); // DP
            if(cr1->sc_called[SC_DIN_CALLED] && cr2->sc_called[SC_DIN_CALLED]){
                cr1->ts->l = 0; cr2->ts->l = 0;
                for(auto& p: cr1->sc_din[e]){
                    ksprintf(cr1->ts, "%s%s-%d|", SC_ID_PREFIX, p.first.c_str(), p.second);
                }
                for(auto& p: cr2->sc_din[e]){
                    ksprintf(cr2->ts, "%s%s-%d|", SC_ID_PREFIX, p.first.c_str(), p.second);
                }
                const char* scs[2] = {(const char*)cr1->ts->s, (const char*)cr2->ts->s};
                bcf_update_format_string(hdr, rec, "SC", scs, rec->n_sample); // SC
            }
            recs.push_back(rec);
        }
    }
    return 0;
}
