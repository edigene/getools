#include "gedetect.h"

void GEDetector::cbam2vars(){
    samFile* fp = sam_open(cbam.c_str(), "r");
    bam_hdr_t* h = sam_hdr_read(fp);
    bam1_t* b = bam_init1();
    hts_idx_t* idx = sam_index_load(fp, cbam.c_str());
    int32_t tid = sam_hdr_name2tid(h, name);
    hts_itr_t* itr = sam_itr_queryi(idx, tid, 0, rlen);
    cins.resize(rlen+1, std::unordered_map<std::string, int32_t>());
    cdis.resize(rlen+1, std::unordered_map<std::string, int32_t>());
    cdel.resize(rlen+1, std::vector<int32_t>());
    for(int32_t i = 0; i < rlen+1; ++i) cdel[i].resize(rlen+1, 0);
    csnv.resize(rlen+1, std::vector<int32_t>());
    for(int32_t i = 0; i < rlen+1; ++i) csnv[i].resize(5, 0);
    int rpos = 0, qpos = 0, opint = 0, oplen = 0;
    uint32_t* cigar = NULL;
    uint8_t* vsd = NULL;
    uint8_t* seq = NULL;
    uint8_t* ccd = NULL;
    int64_t vsv = 0;
    int64_t ccv = 0;
    ctrltt = 0;
    std::string vseq; vseq.reserve(rlen);
    while(sam_itr_next(fp, itr, b) >= 0){
        vsd = bam_aux_get(b, "VS");
        ccd = bam_aux_get(b, "CC");
        vsv = bam_aux2i(vsd);
        ccv = bam_aux2i(ccd);
        if(vsv & dropmask) continue;
        if(vsv & KSW_FREPSEQR) ctrltt += ccv;
        if((vsv & KSW_FREPSEQR) && (vsv & KSW_FMAYVARSEQ)){// rep and var
            rpos = b->core.pos;
            qpos = 0;
            cigar = bam_get_cigar(b);
            seq = bam_get_seq(b);
            for(uint32_t i = 0; i < b->core.n_cigar; ++i){
                opint = bam_cigar_op(cigar[i]);
                oplen = bam_cigar_oplen(cigar[i]);
                switch(opint){
                    case BAM_CINS:
                        if(rpos > vvbeg && rpos <= vvend){
                            vseq.clear();
                            for(int p = 0; p < oplen; ++p) vseq.push_back(seq_nt16_str[bam_seqi(seq, qpos+p)]);
                            auto iter = cins[rpos].find(vseq);
                            if(iter == cins[rpos].end()) cins[rpos][vseq] = ccv;
                            else iter->second += ccv;
                        }
                        qpos += oplen;
                        break;
                    case BAM_CDEL:
                        if(rpos <= vvend && rpos + oplen - 1 >= vvbeg) cdel[rpos][oplen] += ccv;
                        rpos += oplen;
                        break;
                    case BAM_CDIFF:
                        if(rpos <= vvend && rpos + oplen - 1 >= vvbeg){
                            if(oplen == 1){
                                csnv[rpos][seq_nt16_int[bam_seqi(seq, qpos)]] += ccv;
                            }else{
                                vseq.clear();
                                for(int p = 0; p < oplen; ++p) vseq.push_back(seq_nt16_str[bam_seqi(seq, qpos+p)]);
                                auto iter = cdis[rpos].find(vseq);
                                if(iter == cdis[rpos].end()) cdis[rpos][vseq] = ccv;
                                else iter->second += ccv;
                            }
                        }
                        rpos += oplen;
                        qpos += oplen;
                        break;
                    case BAM_CEQUAL:
                        rpos += oplen;
                        qpos += oplen;
                        break;
                    case BAM_CSOFT_CLIP:
                        qpos += oplen;
                        break;
                    default:
                        break;
                }
            }
        }
    }
    sam_close(fp);
    bam_hdr_destroy(h);
    bam_destroy1(b);
    hts_idx_destroy(idx);
    hts_itr_destroy(itr);
}

void GEDetector::markctrlref(){
    int rpos = 0, qpos = 0, opint = 0, oplen = 0;
    uint32_t* cigar = NULL;
    uint8_t* vsd = NULL;
    uint8_t* seq = NULL;
    int64_t vsv = 0;
    std::string vseq; vseq.reserve(rlen);
    KSW_FTYPE rev = ~(KSW_FMAYVARSEQ | KSW_FINSINRNG | KSW_FDELINRNG | KSW_FSNVINRNG | KSW_FDIINRNG);
    for(auto& b: vars){
        vsd = bam_aux_get(b, "VS");
        vsv = bam_aux2i(vsd);
        if(vsv & dropmask) continue;
        if((vsv & KSW_FREPSEQR) && (vsv & KSW_FMAYVARSEQ)){// rep and var
            vsv &= rev;
            rpos = b->core.pos;
            qpos = 0;
            cigar = bam_get_cigar(b);
            seq = bam_get_seq(b);
            for(uint32_t i = 0; i < b->core.n_cigar; ++i){
                opint = bam_cigar_op(cigar[i]);
                oplen = bam_cigar_oplen(cigar[i]);
                switch(opint){
                    case BAM_CINS:
                        if(rpos > vvbeg && rpos <= vvend){
                            vseq.clear();
                            for(int p = 0; p < oplen; ++p) vseq.push_back(seq_nt16_str[bam_seqi(seq, qpos+p)]);
                            auto iter = cins[rpos].find(vseq);
                            if(iter == cins[rpos].end() || iter->second < minctrlc) vsv |= KSW_FINSINRNG;
                        }
                        qpos += oplen;
                        break;
                    case BAM_CDEL:
                        if(rpos <= vvend && rpos + oplen - 1 >= vvbeg){
                            if(cdel[rpos][oplen] < minctrlc) vsv |= KSW_FDELINRNG;
                        }
                        rpos += oplen;
                        break;
                    case BAM_CDIFF:
                        if(rpos <= vvend && rpos + oplen - 1 >= vvbeg){
                            if(oplen == 1){
                                if(csnv[rpos][seq_nt16_int[bam_seqi(seq, qpos)]] < minctrlc) vsv |= KSW_FSNVINRNG;
                            }else{
                                vseq.clear();
                                for(int p = 0; p < oplen; ++p) vseq.push_back(seq_nt16_str[bam_seqi(seq, qpos+p)]);
                                auto iter = cdis[rpos].find(vseq);
                                if(iter == cdis[rpos].end() || iter->second < minctrlc) vsv |= KSW_FDIINRNG;
                            }
                        }
                        rpos += oplen;
                        qpos += oplen;
                        break;
                    case BAM_CEQUAL:
                        rpos += oplen;
                        qpos += oplen;
                        break;
                    case BAM_CSOFT_CLIP:
                        qpos += oplen;
                        break;
                    default:
                        break;
                }
            }
            if((vsv & KSW_FSPANSGR) && (vsv & varmask)) vsv |= KSW_FMAYVARSEQ;
            bam_aux_update_int(b, "VS", vsv);
        }
    }
}
