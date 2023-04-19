#include "geplot.h"

void gep_opt_t::plot(){
    // get bams
    samFile* fp = sam_open(inbam, "r");
    hts_idx_t* idx = sam_index_load(fp, inbam);
    if(!idx){
        fprintf(stderr, "input bam must be sorted by coordinats\n");
        exit(EXIT_FAILURE);
    }
    bamhdr = sam_hdr_read(fp);
    int tid = sam_hdr_name2tid(bamhdr, name);
    hts_itr_t* itr= sam_itr_queryi(idx, tid, 0, sam_hdr_tid2len(bamhdr, tid));
    bam1_t* b = bam_init1();
    std::vector<bam1_t*> vars;
    while(sam_itr_next(fp, itr, b) >= 0){
        vars.push_back(b);
        b = bam_init1();
    }
    // plot
    plot(vars);
    // release src
    sam_close(fp);
    hts_idx_destroy(idx);
    sam_itr_destroy(itr);
    for(auto& e: vars){
        bam_destroy1(e);
    }
    bam_destroy1(b);
}

void gep_opt_t::plot(const std::vector<bam1_t*>& vbs){
    std::vector<var_rec_t*> vars;
    uint8_t* dcn = NULL;
    int64_t vsv = 0;
    int64_t ttr = 0;
    int64_t ttv = 0;
    for(auto& b: vbs){
        vsv = bam_aux2i(bam_aux_get(b, "VS"));
        if((vsv & KSW_FREPSEQR) && (!(vsv & dropm))){
            var_rec_t* var = (var_rec_t*)calloc(1, sizeof(var_rec_t));
            var->b = b;
            var->cc = bam_aux2i(bam_aux_get(b, "CC"));
            var->vs = vsv;
            dcn = bam_aux_get(b, CELL_BARCODE_COUNT_TAG);
            if(dcn) var->cn = bam_aux2i(dcn);
            vars.push_back(var);
            ttr += var->cc;
            macnt += var->cc;
            ++ttv;
            if(vsv & KSW_FMAYVARSEQ) edcnt += var->cc;
        }
    }
    allrs = ttr;
    std::sort(vars.begin(), vars.end(), var_sort_t());
    if(vars.empty()){ empty = true; return; }
    // find REF to put in the first place
    if(!(vars[0]->vs & KSW_FREFTYPE)){
        for(size_t i = 1; i < vars.size(); ++i){
            if(vars[i]->vs & KSW_FREFTYPE){
                std::swap(vars[0], vars[i]);
                break;
            }
        }
    }
    // output variant details
    kstring_t* ks = (kstring_t*)calloc(1, sizeof(kstring_t)); // variants tsv buffer
    ksprintf(ks, "Freq\tCount\tWidth\tStart\tEnd\tType\tREF\tALT\tTmark\n");
    std::vector<int> rb;
    std::vector<int> re;
    std::vector<std::string> rs;
    std::vector<std::string> as;
    kstring_t* ts = (kstring_t*)calloc(1, sizeof(kstring_t)); // temp buffer
    ts->l = 0; ts->m = rlen; ts->s = (char*)malloc(rlen*sizeof(char));
    uint32_t* cigar = NULL;
    int rpos = 0, qpos = 0, oplen = 0, opint = 0;
    kstring_t* ss = (kstring_t*)calloc(1, sizeof(kstring_t)); // topn buffer
    bool ltntop = false;
    int topc = 0;
    for(int64_t i = 0; i < ttv; ++i){
        vars[i]->af = (double)vars[i]->cc/(double)ttr;
        ltntop = false; if((topc < topn + 1) && (vars[i]->vs & topm)){ ltntop = true; ++topc; }
        ksprintf(ks, "%lf\t", vars[i]->af); // Freq
        ksprintf(ks, "%d\t", vars[i]->cc); // Count
        rb.clear(), re.clear(), rs.clear(), as.clear();
        cigar = bam_get_cigar(vars[i]->b);
        rpos = vars[i]->b->core.pos;
        qpos = 0;
        if(ltntop) ss->l = 0;
        for(uint32_t c = 0; c < vars[i]->b->core.n_cigar; ++c){
            opint = bam_cigar_op(cigar[c]);
            oplen = bam_cigar_oplen(cigar[c]);
            switch(opint){
                case BAM_CINS:
                    if(rpos > vbeg && rpos <= vend){
                        rb.push_back(rpos);
                        re.push_back(rpos);
                        rs.push_back("-");
                        ts->l = 0;
                        for(int e = 0; e < oplen; ++e){
                            kputc(seq_nt16_str[bam_seqi(bam_get_seq(vars[i]->b), qpos + e)], ts);
                        }
                        kputc('\0', ts);
                        as.push_back(ts->s);
                    }
                    if(ltntop && rpos > refbeg && rpos <= refend){
                        vars[i]->ipos.push_back(rpos - refbeg - 1);
                        vars[i]->ins += oplen;
                    }
                    qpos += oplen;
                    break;
                case BAM_CDEL:
                    if(ltntop && (rpos + oplen - 1 >= refbeg && rpos <= refend)){
                        for(int d = MAX(refbeg, rpos); d <= MIN(rpos+oplen-1, refend); ++d) kputc('-', ss);
                        vars[i]->del += oplen;
                    }
                    if(rpos + oplen -1 >= vbeg && rpos <= vend){
                        rb.push_back(rpos);
                        re.push_back(rpos+oplen-1);
                        ts->l = 0;
                        kputsn(amplicon + rpos, oplen, ts);
                        kputc('\0', ts);
                        rs.push_back(ts->s);
                        as.push_back("-");
                    }
                    rpos += oplen;
                    break;
                case BAM_CEQUAL:
                    if(ltntop && (rpos + oplen - 1 >= refbeg && rpos <= refend)){
                        for(int d = MAX(refbeg, rpos)-rpos; d <= MIN(rpos+oplen-1, refend)-rpos; ++d){
                            kputc(seq_nt16_str[bam_seqi(bam_get_seq(vars[i]->b), qpos + d)], ss);
                        }
                    }
                    rpos += oplen;
                    qpos += oplen;
                    break;
                case BAM_CDIFF:
                    if(ltntop && (rpos + oplen - 1 >= refbeg && rpos <= refend)){
                        for(int d = MAX(refbeg, rpos)-rpos; d <= MIN(rpos+oplen-1, refend)-rpos; ++d){
                            kputc(seq_nt16_str[bam_seqi(bam_get_seq(vars[i]->b), qpos + d)], ss);
                        }
                    }
                    if(rpos + oplen - 1 >= refbeg && rpos <= refend){
                        rb.push_back(rpos);
                        re.push_back(rpos+oplen-1);
                        ts->l = 0;
                        kputsn(amplicon + rpos, oplen, ts);
                        kputc('\0', ts);
                        rs.push_back(ts->s);
                        ts->l = 0;
                        for(int d = 0; d < oplen; ++d){
                            char base = seq_nt16_str[bam_seqi(bam_get_seq(vars[i]->b), qpos + d)];
                            kputc(base, ts);
                        }
                        kputc('\0', ts);
                        as.push_back(ts->s);
                    }
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
        if(abs(vars[i]->ins + vars[i]->del) % 3) fscnt += vars[i]->cc;
        if(ltntop){
            vars[i]->seq = strndup(ss->s, ss->l);
            vars[i]->lseq = ss->l;
        }
        // output var
        if(rb.empty()){// REF in this range
            ksprintf(ks, "0\t0\t0\tREF\t\t\t0\n");
        }else{// choose only one var, delins first
            std::vector<std::string> tvrs;
            std::vector<std::string> tvas;
            std::vector<std::string> tvty;
            std::vector<int> tvrb;
            std::vector<int> tvre;
            uint32_t vbeg = 0;
            uint32_t vend = 0;
            uint32_t vpos = 1;
            re.push_back(rlen + rlen); // magic
            rb.push_back(rlen + rlen); // magic
            for(vpos = 1; vpos < rb.size(); ++vpos){
                if(rb[vpos] != re[vpos-1] + 1){// output previous
                    if(vbeg != vend){// indel
                        std::string trs;
                        std::string tas;
                        for(uint32_t tvi = vbeg; tvi <= vend; ++tvi){
                            if(rs[tvi] != "-") trs.append(rs[tvi]);
                            if(as[tvi] != "-") tas.append(as[tvi]);
                        }
                        tvrs.push_back(trs);
                        tvas.push_back(tas);
                        tvrb.push_back(rb[vbeg]);
                        tvre.push_back(re[vend]);
                        tvty.push_back("DELINS");
                    }else{
                        tvrs.push_back(rs[vbeg]);
                        tvas.push_back(as[vbeg]);
                        tvrb.push_back(rb[vbeg]);
                        tvre.push_back(re[vbeg]);
                        if(rs[vbeg] == rs[vend]){
                            if(rs[vbeg] == "-") tvty.push_back("INS");
                            else if(as[vbeg] == "-") tvty.push_back("DEL");
                            else if(rs[vbeg].size() == as[vbeg].size()){
                                if(rs[vbeg].size() == 1) tvty.push_back("SNV");
                                else tvty.push_back("DELINS");
                            }
                        }else{
                            if(as[vbeg] == "-") tvty.push_back("DEL");
                            else tvty.push_back("DELINS");
                        }
                    }
                }
            }
            int32_t delins = -1;
            int32_t vardel = -1;
            int32_t varins = -1;
            int32_t varsnv = -1;
            int outvar = -1, tvmsk = 0;
            for(uint32_t j = 0; j < tvty.size(); ++j){
                if(tvty[j] == "DELINS"){
                    delins = j;
                }else if(tvty[j] == "INS"){
                    varins = j;
                }else if(tvty[j] == "DEL"){
                    vardel = j;
                }else if(tvty[j] == "SNV"){
                    varsnv = j;
                }
            }
            if(delins >= 0){
                outvar = delins;
                tvmsk = GEVAR_DIN;
            }else if(varins >= 0){
                outvar = varins;
                tvmsk = GEVAR_INS;
            }else if(vardel >= 0){
                outvar = vardel;
                tvmsk = GEVAR_DEL;
            }else if(varsnv >= 0){
                outvar = varsnv;
                tvmsk = GEVAR_SNV;
            }
            ksprintf(ks, "%d\t%d\t%d\t%s\t%s\t%s\t%d\n",
                         tvre[outvar] - tvrb[outvar] + 1,
                         tvrb[outvar] - sgrbeg, 
                         tvre[outvar] - sgrbeg, 
                         tvty[outvar].c_str(), 
                         tvrs[outvar].c_str(),
                         tvas[outvar].c_str(), 
                         tvmsk);
        }
    }
    if(outtsv){
        gzFile fp = gzopen(outtsv, "wb");
        gzwrite(fp, ks->s, sizeof(char)*ks->l);
        gzclose(fp);
    }
    free(ks->s); free(ks);
    free(ts->s); free(ts);
    // plot
    plot(vars);
}

void gep_opt_t::plot(const std::vector<var_rec_t*>& vrs){
    kstring_t* seqs = (kstring_t*)calloc(1, sizeof(kstring_t)); // sequence bases
    ksprintf(seqs, "seqs <- c(\n");
    kstring_t* seqc = (kstring_t*)calloc(1, sizeof(kstring_t)); // sequence count
    ksprintf(seqc, "seqc <- c(");
    kstring_t* seqf = (kstring_t*)calloc(1, sizeof(kstring_t)); // sequence freq
    ksprintf(seqf, "seqf <- c(");
    kstring_t* seqw = (kstring_t*)calloc(1, sizeof(kstring_t)); // sequence frameshift count
    ksprintf(seqw, "seqw <- c(");
    kstring_t* seqi = (kstring_t*)calloc(1, sizeof(kstring_t)); // ins pos list of each seq
    ksprintf(seqi, "seqi <- list(");
    kstring_t* ampt = (kstring_t*)calloc(1, sizeof(kstring_t)); // first sequence translation
    ksprintf(ampt, "ampt <- c(\n");
    kstring_t* fseq = (kstring_t*)calloc(1, sizeof(kstring_t)); // first sequence seq
    if(genhtvar){
        hvars.reserve(topn+1);
        html_var_t* rhv = (html_var_t*)calloc(1, sizeof(html_var_t));
        rhv->seq = (kstring_t*)calloc(1, sizeof(kstring_t));
        ksprintf(rhv->seq, "%s", amplicon);
        hvars.push_back(rhv);
    }
    int n = 0;
    // output topn BAM
    samFile *bofp = NULL;
    if(outbam){
        bofp = sam_open(outbam, "wb");
        assert(sam_hdr_write(bofp, bamhdr) >= 0);
    }
    for(size_t i = 0; i < vrs.size(); ++i){
        if(n == topn + 1) break;
        if(vrs[i]->vs & topm){
            ++n;
            // seqs
            if(bofp) assert(sam_write1(bofp, bamhdr, vrs[i]->b) >= 0);
            kputc('"', seqs);
            int extb = 0;
            for(int c = 0; c < vrs[i]->b->core.pos - refbeg; ++c){ kputc('N', seqs);  ++extb; }// fill pre with N
            kputsn(vrs[i]->seq, vrs[i]->lseq, seqs);
            for(int c = 0; c <  refgot - extb - vrs[i]->lseq; ++c) kputc('N', seqs); // fill suf with N
            kputs("\",\n", seqs);
            if(genhtvar){// html seqs
                html_var_t* hv = (html_var_t*)calloc(1, sizeof(html_var_t));
                hv->seq = (kstring_t*)calloc(1, sizeof(kstring_t));
                hv->cigar = (kstring_t*)calloc(1, sizeof(kstring_t));
                bam2hv(vrs[i]->b, hv);
                hvars.push_back(hv);
            }
            if(n == 1){
                kputsn(seqs->s + 12, refgot, fseq);
                // ampt
                for(int i = 0; i < 3; ++i){
                    char* aa = nc2aa(fseq->s + i, fseq->l - i, ismt);
                    kputc('"', ampt);
                    kputs(aa, ampt);
                    kputs("\",\n", ampt);
                    free(aa);
                }
                for(int i = 0; i < 3; ++i){
                    char* aa = ncr2aa(fseq->s, fseq->l - i, ismt);
                    kputc('"', ampt);
                    kputs(aa, ampt);
                    kputs("\",\n", ampt);
                    free(aa);
                }
                ampt->l -= 1;
                ampt->s[ampt->l-1] = ')';
            }
            // seqc
            ksprintf(seqc, "%d,", vrs[i]->cc);
            // seqf
            ksprintf(seqf, "%lf,", vrs[i]->af);
            // seqw
            ksprintf(seqw, "%d,", vrs[i]->ins - vrs[i]->del);
            // seqi
            kputs("c(", seqi);
            for(auto& e: vrs[i]->ipos){
                ksprintf(seqi, "%d,", e);
            }
            if(seqi->s[seqi->l-1] == ',') seqi->s[seqi->l-1] = ')';
            else kputc(')', seqi);
            kputs(",", seqi);
        }
    }
    if(bofp) sam_close(bofp);
    seqs->l -= 1;
    seqs->s[seqs->l-1] = ')';
    seqc->s[seqc->l-1] = ')';
    seqf->s[seqf->l-1] = ')';
    seqw->s[seqw->l-1] = ')';
    seqi->s[seqi->l-1] = ')';
    // commands
    kstring_t* cmds = (kstring_t*)calloc(1, sizeof(kstring_t));
    ksprintf(cmds, "spn <- \"%s\"\n", name);
    ksprintf(cmds, "sgrbeg <- %d\n", sgrbeg - refbeg + 1);
    ksprintf(cmds, "sgrend <- %d\n", sgrend - refbeg + 1);
    ksprintf(cmds, "fspct <- %d\n", (int)((double)fscnt * 100.0 / (double)macnt));
    ksprintf(cmds, "edpct <- %d\n", (int)((double)edcnt * 100.0 / (double)macnt));
    ksprintf(cmds, "mapct <- %d\n", (int)((double)macnt * 100.0 / (double)macnt));
    ksprintf(cmds, "trsn <- c(\"1st,5'->3'\", \"2nd,5'->3'\",\n");
    ksprintf(cmds, "\"3rd,5'->3'\", \"1st,3'<-5'\", \"2nd,3'<-5'\", \"3rd,3'<-5'\")\n");
    ksprintf(cmds, "amplen <- nchar(seqs[1])\n");
    ksprintf(cmds, "trsto <- c(0, 1, 2, 0, 1, 2)\n");
    ksprintf(cmds, "trsto[4] <- amplen %%%% 3\n");
    ksprintf(cmds, "trsto[5] <- (amplen-1) %%%% 3\n");
    ksprintf(cmds, "trsto[6] <- (amplen-2) %%%% 3\n");
    ksprintf(cmds, "seqn <- c(\"amplicon\")\n");
    ksprintf(cmds, "for(i in 1:(length(seqs)-1)) seqn <- append(seqn, i)\n");
    ksprintf(cmds, "tots <- length(seqs)\n");
    ksprintf(cmds, "aa_line <- 6\n");
    ksprintf(cmds, "nc_line <- tots\n");
    ksprintf(cmds, "frame_width <- 2\n");
    ksprintf(cmds, "frame_blkpre <- 10\n");
    ksprintf(cmds, "frame_blksuf <- 3\n");
    ksprintf(cmds, "pref_width <- frame_width + frame_blkpre + frame_blksuf\n");
    ksprintf(cmds, "freq_blkpre <- 5\n");
    ksprintf(cmds, "freq_width <- 5\n");
    ksprintf(cmds, "count_width <- 5\n");
    ksprintf(cmds, "F_width <- 5\n");
    ksprintf(cmds, "fdist_width <- freq_width + count_width + F_width\n");
    ksprintf(cmds, "pctrgb <- col2rgb(\"#7FCDCD\")\n");
    ksprintf(cmds, "colnc <- c(\"#EA4335\", \"#4285F4\", \"#FBBC05\", \"#34A853\",\n");
    ksprintf(cmds, "           \"#EA4335\", \"#4285F4\", \"#FBBC05\", \"#34A853\", \"#FFFFFF\")\n");
    ksprintf(cmds, "names(colnc) <- c(\"A\", \"C\", \"G\", \"T\", \"a\", \"c\", \"g\", \"t\", \"-\")\n");
    ksprintf(cmds, "colaa <- c(\"#E15D44\", \"#98B4D4\", \"#D65076\", \"#BC243C\", \"#009B77\",\n");
    ksprintf(cmds, "           \"#D65076\", \"#BC243C\", \"#E15D44\", \"#D65076\", \"#009B77\",\n");
    ksprintf(cmds, "           \"#009B77\", \"#98B4D4\", \"#009B77\", \"#009B77\", \"#009B77\",\n");
    ksprintf(cmds, "           \"#E15D44\", \"#E15D44\", \"#009B77\", \"#009B77\", \"#009B77\",\n");
    ksprintf(cmds, "           rep(\"#DFCFBE\", 10), \"#FFFFFF\")\n");
    ksprintf(cmds, "names(colaa) <- c(\"A\", \"R\", \"N\", \"D\", \"C\",\n");
    ksprintf(cmds, "                  \"Q\", \"E\", \"G\", \"H\", \"I\",\n");
    ksprintf(cmds, "                  \"L\", \"K\", \"M\", \"F\", \"P\",\n");
    ksprintf(cmds, "                  \"S\", \"T\", \"W\", \"Y\", \"V\",\n");
    ksprintf(cmds, "                  \"U\", \"O\", \"B\", \"J\", \"Z\",\n");
    ksprintf(cmds, "                  \"X\", \"*\", \"-\", \"+\", \".\")\n");
    ksprintf(cmds, "gwidth <- 8\n");
    ksprintf(cmds, "gheight <- 6\n");
    ksprintf(cmds, "if(length(seqs) <= 15){\n");
    ksprintf(cmds, "  gwidth <- 6\n");
    ksprintf(cmds, "  gheight <- 4\n");
    ksprintf(cmds, "}\n");
    ksprintf(cmds, "if(length(seqs) > 50){\n");
    ksprintf(cmds, "  gwidth <- 12\n");
    ksprintf(cmds, "  gheight <- 9\n");
    ksprintf(cmds, "}\n");
    ksprintf(cmds, "gytx <- gheight/gwidth\n");
    ksprintf(cmds, "pdf(file=\"%s\", width = gwidth, height = gheight, onefile = F)\n", outpdf);
    ksprintf(cmds, "par(cex=1, mar=c(0, 0, 0, 0))\n");
    ksprintf(cmds, "plot.new()\n");
    ksprintf(cmds, "strwh <- c(strwidth(\"M\", cex=0.5), strheight(\"M\", cex=1))\n");
    ksprintf(cmds, "wunit <- strwh[1]\n");
    ksprintf(cmds, "hunit <- strwh[2]\n");
    ksprintf(cmds, "ofsty <- 0\n");
    ksprintf(cmds, "ofstx <- 0\n");
    ksprintf(cmds, "miny<- 0\n");
    ksprintf(cmds, "maxy <- (tots + 7) * (hunit + ofsty)\n");
    ksprintf(cmds, "minx <- 0\n");
    ksprintf(cmds, "maxx <- (frame_width + frame_blkpre + frame_blksuf + amplen + freq_blkpre + fdist_width) * (wunit + ofstx)\n");
    ksprintf(cmds, "plot(x=c(minx, maxx), y=c(miny-hunit, maxy), axes = F, type = \"n\", xlab = \"\", ylab = \"\", frame.plot = F)\n");
    ksprintf(cmds, "# first round, plot nucleotide sequence\n");
    ksprintf(cmds, "xbeg <- 0\n");
    ksprintf(cmds, "ybeg <- -hunit\n");
    ksprintf(cmds, "yend <- hunit\n");
    ksprintf(cmds, "ymid <- 0.5 * (yend - ybeg)\n");
    ksprintf(cmds, "xstep <- wunit + ofstx\n");
    ksprintf(cmds, "ystep <- hunit + ofsty\n");
    ksprintf(cmds, "xmid <- 0.5 * wunit\n");
    ksprintf(cmds, "snxpos <- (frame_blkpre + 0.5 * frame_width) * wunit\n");
    ksprintf(cmds, "for(i in length(seqs):1){\n");
    ksprintf(cmds, "  ybeg <- ybeg + ystep\n");
    ksprintf(cmds, "  yend <- ybeg + ystep\n");
    ksprintf(cmds, "  ymid <- ybeg + 0.5 * ystep\n");
    ksprintf(cmds, "  text(x=snxpos, y=ymid, labels=seqn[i], cex = 0.5, pos = 2)\n");
    ksprintf(cmds, "  xbeg <- pref_width * wunit\n");
    ksprintf(cmds, "  for(n in 1:nchar(seqs[i])){\n");
    ksprintf(cmds, "    rect(xbeg + (n-1) * xstep, ybeg, xbeg + n * xstep, yend, col=colnc[substr(seqs[i], n, n)], border=colnc[substr(seqs[i], n, n)])\n");
    ksprintf(cmds, "    text(xbeg + (n-0.65) * xstep, ymid, substr(seqs[i], n, n), cex = 0.3)\n");
    ksprintf(cmds, "  }\n");
    ksprintf(cmds, "  for(n in seqi[[i]]){\n");
    ksprintf(cmds, "    points(xbeg + n * xstep, ybeg, type=\"p\", pch=17, cex=0.5, lwd=0.5)\n");
    ksprintf(cmds, "  }\n");
    ksprintf(cmds, "}\n");
    ksprintf(cmds, "# second round(a), plot underline of nucleotide sequence\n");
    ksprintf(cmds, "xbeg <- 0\n");
    ksprintf(cmds, "ybeg <- -hunit\n");
    ksprintf(cmds, "yend <- hunit\n");
    ksprintf(cmds, "ymid <- 0.5 * (yend - ybeg)\n");
    ksprintf(cmds, "xstep <- wunit + ofstx\n");
    ksprintf(cmds, "ystep <- hunit + ofsty\n");
    ksprintf(cmds, "xmid <- 0.5 * wunit\n");
    ksprintf(cmds, "for(i in length(seqs):1){\n");
    ksprintf(cmds, "  ybeg <- ybeg + ystep\n");
    ksprintf(cmds, "  yend <- ybeg + ystep\n");
    ksprintf(cmds, "  ymid <- ybeg + 0.5 * ystep\n");
    ksprintf(cmds, "  xbeg <- pref_width * wunit\n");
    ksprintf(cmds, "  if(i != length(seqs)) segments(xbeg, ybeg, xbeg + nchar(seqs[i]) * xstep, ybeg, col=\"white\", lty=2, lwd=0.2)\n");
    ksprintf(cmds, "  else{# plot nc relative pos\n");
    ksprintf(cmds, "    relpc <- c()\n");
    ksprintf(cmds, "    for(j in 1:sgrbeg){\n");
    ksprintf(cmds, "      if((sgrbeg - j)%%%%5 == 0){\n");
    ksprintf(cmds, "        relpc <- append(relpc, j)\n");
    ksprintf(cmds, "      }\n");
    ksprintf(cmds, "    }\n");
    ksprintf(cmds, "    for(j in (sgrbeg+1):nchar(seqs[1])){\n");
    ksprintf(cmds, "      if((j-sgrbeg)%%%%5 == 0){\n");
    ksprintf(cmds, "        relpc <- append(relpc, j)\n");
    ksprintf(cmds, "      }\n");
    ksprintf(cmds, "    }\n");
    ksprintf(cmds, "    for(j in relpc){\n");
    ksprintf(cmds, "      mpx <- xbeg + (j -1) * xstep\n");
    ksprintf(cmds, "      mpyd <- ybeg - ystep\n");
    ksprintf(cmds, "      mpyu <- ybeg - 0.5 * ystep\n");
    ksprintf(cmds, "      text(mpx+0.5*xstep, mpyu, \"|\", cex = 0.5, adj=0.5)\n");
    ksprintf(cmds, "      text(mpx+0.5*xstep, mpyd, j - sgrbeg, cex = 0.3, adj=0.5)\n");
    ksprintf(cmds, "    }\n");
    ksprintf(cmds, "    pybm <- ybeg - 2 * ystep\n");
    ksprintf(cmds, "    pxbt <- xbeg + length(relpc) * 2.5 * xstep\n");
    ksprintf(cmds, "    text(pxbt, pybm, \"Relative Nucleotide Position\", cex=0.5, adj=0.5)\n");
    ksprintf(cmds, "  }\n");
    ksprintf(cmds, "}\n");
    ksprintf(cmds, "# secround(b), plot sgrna marker\n");
    ksprintf(cmds, "sgybeg <- ybeg\n");
    ksprintf(cmds, "sgyend <- ybeg + ystep\n");
    ksprintf(cmds, "sgxbeg <- pref_width * wunit + (sgrbeg - 1) * xstep\n");
    ksprintf(cmds, "sgxend <- sgxbeg + (sgrend - sgrbeg + 1) * xstep\n");
    ksprintf(cmds, "rect(sgxbeg, sgybeg, sgxend, sgyend, col=NA, border=\"black\")\n");
    ksprintf(cmds, "# third round, plot bottom left legend\n");
    ksprintf(cmds, "ypos1 <- 0.5 * yend\n");
    ksprintf(cmds, "text(0, ypos1, spn, cex = 0.5, pos = 2, srt=90)\n");
    ksprintf(cmds, "# forth round, plot aa translation\n");
    ksprintf(cmds, "ypos1 <- ybeg\n");
    ksprintf(cmds, "ybeg <- ybeg + ystep\n");
    ksprintf(cmds, "pwunit <- 3 * wunit\n");
    ksprintf(cmds, "pxstep <- 3 * wunit\n");
    ksprintf(cmds, "for(i in length(trsn):1){\n");
    ksprintf(cmds, "  ybeg <- ybeg + ystep\n");
    ksprintf(cmds, "  yend <- ybeg + ystep\n");
    ksprintf(cmds, "  ymid <- ybeg + 0.5 * ystep\n");
    ksprintf(cmds, "  text(x=snxpos, y=ymid, labels=trsn[i], cex = 0.5, pos = 2)\n");
    ksprintf(cmds, "  xbeg <- (pref_width + trsto[i]) * wunit\n");
    ksprintf(cmds, "  for(n in 1:nchar(ampt[i])){\n");
    ksprintf(cmds, "    rect(xbeg + (n - 1)* pxstep, ybeg, xbeg + n * pxstep, yend, col=colaa[substr(ampt[i], n, n)], border=\"white\")\n");
    ksprintf(cmds, "    text(xbeg + (n - 0.5) * pxstep, ymid, substr(ampt[i], n, n), cex = 0.3, adj=0.5)\n");
    ksprintf(cmds, "  }\n");
    ksprintf(cmds, "}\n");
    ksprintf(cmds, "# fifth round, plot top left legend\n");
    ksprintf(cmds, "yltbeg <- ypos1 + ystep\n");
    ksprintf(cmds, "yltend <- yend\n");
    ksprintf(cmds, "text(0, yltbeg + 0.5*(yltend-yltbeg) + 2 * hunit, \"Frame\", cex = 0.5, pos = 2, srt=90)\n");
    ksprintf(cmds, "# sixth round, plot freq table\n");
    ksprintf(cmds, "xbeg <- (pref_width + amplen + freq_blkpre) * wunit\n");
    ksprintf(cmds, "ybeg <- -hunit\n");
    ksprintf(cmds, "yend <- hunit\n");
    ksprintf(cmds, "ymid <- 0.5 * (yend - ybeg)\n");
    ksprintf(cmds, "xstep <- wunit + ofstx\n");
    ksprintf(cmds, "ystep <- hunit + ofsty\n");
    ksprintf(cmds, "xmid <- 0.5 * wunit\n");
    ksprintf(cmds, "freq_width <- freq_width * wunit\n");
    ksprintf(cmds, "count_width <- count_width * wunit\n");
    ksprintf(cmds, "F_width <- F_width * wunit\n");
    ksprintf(cmds, "pbeg <- (pref_width + amplen + freq_blkpre) * wunit\n");
    ksprintf(cmds, "sxbeg <- pbeg - 0.5 * freq_width\n");
    ksprintf(cmds, "sxend <- pbeg + fdist_width * xstep\n");
    ksprintf(cmds, "for(i in length(seqs):1){\n");
    ksprintf(cmds, "  xbeg <- pbeg\n");
    ksprintf(cmds, "  ybeg <- ybeg + ystep\n");
    ksprintf(cmds, "  yend <- ybeg + ystep\n");
    ksprintf(cmds, "  ymid <- ybeg + 0.5 * ystep\n");
    ksprintf(cmds, "  # Freq\n");
    ksprintf(cmds, "  xend <- xbeg + freq_width\n");
    ksprintf(cmds, "  reccol <- rgb(pctrgb[1], pctrgb[2], pctrgb[3], alpha = 205*seqf[i], maxColorValue=205)\n");
    ksprintf(cmds, "  rect(xbeg, ybeg, xend, yend, col=reccol, border=reccol)\n");
    ksprintf(cmds, "  text(xbeg + 0.5 * freq_width, ymid, seqf[i], cex = 0.3, adj=0.5)\n");
    ksprintf(cmds, "  # Count\n");
    ksprintf(cmds, "  xbeg <- xend\n");
    ksprintf(cmds, "  xend <- xbeg + count_width\n");
    ksprintf(cmds, "  rect(xbeg, ybeg, xend, yend, col=\"white\", border=\"white\")\n");
    ksprintf(cmds, "  text(xbeg + 0.5 * count_width, ymid, seqc[i], cex = 0.3, adj=0.5)\n");
    ksprintf(cmds, "  # Width\n");
    ksprintf(cmds, "  xbeg <- xend\n");
    ksprintf(cmds, "  xend <- xbeg + F_width\n");
    ksprintf(cmds, "  reccol <- \"#FF0000\"\n");
    ksprintf(cmds, "  if(seqw[i] > 0) reccol <- \"#FF0000\"\n");
    ksprintf(cmds, "  else if(seqw[i] < 0) reccol <- \"#79C753\"\n");
    ksprintf(cmds, "  else reccol = \"white\"\n");
    ksprintf(cmds, "  rect(xbeg, ybeg, xend, yend, col=reccol, border=reccol)\n");
    ksprintf(cmds, "  text(xbeg + 0.5 * F_width, ymid, seqw[i], cex = 0.3, adj=0.5)\n");
    ksprintf(cmds, "  # upper line\n");
    ksprintf(cmds, "  if(i != length(seqs)) segments(pbeg, ybeg, xend, ybeg, col=\"black\", lty=2, lwd=0.2)\n");
    ksprintf(cmds, "}\n");
    ksprintf(cmds, "ybeg <- yend\n");
    ksprintf(cmds, "yend <- ybeg + ystep\n");
    ksprintf(cmds, "ymid <- ybeg + 0.5 * ystep\n");
    ksprintf(cmds, "if(T){ # title\n");
    ksprintf(cmds, "  # Freq\n");
    ksprintf(cmds, "  xbeg <- pbeg\n");
    ksprintf(cmds, "  xend <- xbeg + freq_width\n");
    ksprintf(cmds, "  ymid <- ybeg + 0.5 * ystep\n");
    ksprintf(cmds, "  rect(xbeg, ybeg, xend, yend, col=\"white\", border=\"white\")\n");
    ksprintf(cmds, "  text(xbeg + 0.5 * freq_width, ymid, \"Freq\", cex = 0.3, adj=0.5)\n");
    ksprintf(cmds, "  # Count\n");
    ksprintf(cmds, "  xbeg <- xend\n");
    ksprintf(cmds, "  xend <- xbeg + count_width\n");
    ksprintf(cmds, "  rect(xbeg, ybeg, xend, yend, col=\"white\", border=\"white\")\n");
    ksprintf(cmds, "  text(xbeg + 0.5 * count_width, ymid, \"Count\", cex = 0.3, adj=0.5)\n");
    ksprintf(cmds, "  # Width\n");
    ksprintf(cmds, "  xbeg <- xend\n");
    ksprintf(cmds, "  xend <- xbeg + F_width\n");
    ksprintf(cmds, "  rect(xbeg, ybeg, xend, yend, col=\"white\", border=\"white\")\n");
    ksprintf(cmds, "  text(xbeg + 0.5 * F_width, ymid, \"F\", cex = 0.3, adj=0.5)\n");
    ksprintf(cmds, "}\n");
    ksprintf(cmds, "segments(pbeg, ybeg, xend, ybeg, col=\"black\", lty=2, lwd=0.2)\n");
    ksprintf(cmds, "if(T){ #verticle line\n");
    ksprintf(cmds, "  vybeg <- 0\n");
    ksprintf(cmds, "  vyend <- yend\n");
    ksprintf(cmds, "  xbeg <- pbeg\n");
    ksprintf(cmds, "  xend <- xbeg + freq_width\n");
    ksprintf(cmds, "  vxpos <- xend\n");
    ksprintf(cmds, "  segments(vxpos, vybeg, vxpos, vyend, col=\"black\", lty=2, lwd=0.2)\n");
    ksprintf(cmds, "  xbeg <- xend\n");
    ksprintf(cmds, "  xend <- xbeg + count_width\n");
    ksprintf(cmds, "  vxpos <- xend\n");
    ksprintf(cmds, "  segments(vxpos, vybeg, vxpos, vyend, col=\"black\", lty=2, lwd=0.2)\n");
    ksprintf(cmds, "}\n");
    ksprintf(cmds, "# seventh round, upper right corner\n");
    ksprintf(cmds, "ybeg <- yend + ystep\n");
    ksprintf(cmds, "if(T){\n");
    ksprintf(cmds, "  # F\n");
    ksprintf(cmds, "  xbeg <- pbeg\n");
    ksprintf(cmds, "  xend <- xbeg + freq_width\n");
    ksprintf(cmds, "  yend <- ybeg + ystep\n");
    ksprintf(cmds, "  ymid <- ybeg + 0.5 * ystep\n");
    ksprintf(cmds, "  rect(xbeg, ybeg, xend, yend, col=\"white\", border=\"white\")\n");
    ksprintf(cmds, "  text(xbeg + 0.5 * freq_width, ymid, \"F\", cex = 0.3, adj=0.5)\n");
    ksprintf(cmds, "  recbeg <- xbeg + freq_width\n");
    ksprintf(cmds, "  recend <- recbeg + (sxend - recbeg) * fspct / 100\n");
    ksprintf(cmds, "  rect(recbeg, ybeg, recend, yend, col=\"#79C753\", border=\"#79C753\")\n");
    ksprintf(cmds, "  nmadj <- 1\n");
    ksprintf(cmds, "  if(fspct < 10) nmadj <- 0\n");
    ksprintf(cmds, "  text(recend, ymid, fspct, cex = 0.3, adj=nmadj)\n");
    ksprintf(cmds, "  # Edited\n");
    ksprintf(cmds, "  ybeg <- yend\n");
    ksprintf(cmds, "  yend <- ybeg + ystep\n");
    ksprintf(cmds, "  ymid <- ybeg + 0.5 * ystep\n");
    ksprintf(cmds, "  rect(xbeg, ybeg, xend, yend, col=\"white\", border=\"white\")\n");
    ksprintf(cmds, "  text(xbeg + 0.5 * freq_width, ymid, \"Edited\", cex = 0.3, adj=0.5)\n");
    ksprintf(cmds, "  recbeg <- xbeg + freq_width\n");
    ksprintf(cmds, "  recend <- recbeg + (sxend - recbeg) * edpct / 100\n");
    ksprintf(cmds, "  rect(recbeg, ybeg, recend, yend, col=\"#7FCDCD\", border=\"#7FCDCD\")\n");
    ksprintf(cmds, "  nmadj <- 1\n");
    ksprintf(cmds, "  if(edpct < 10) nmadj <- 0\n");
    ksprintf(cmds, "  text(recend, ymid, edpct, cex = 0.3, adj=nmadj)\n");
    ksprintf(cmds, "  # Match\n");
    ksprintf(cmds, "  ybeg <- yend\n");
    ksprintf(cmds, "  yend <- ybeg + ystep\n");
    ksprintf(cmds, "  ymid <- ybeg + 0.5 * ystep\n");
    ksprintf(cmds, "  rect(xbeg, ybeg, xend, yend, col=\"white\", border=\"white\")\n");
    ksprintf(cmds, "  text(xbeg + 0.5 * freq_width, ymid, \"Match\", cex = 0.3, adj=0.5)\n");
    ksprintf(cmds, "  recbeg <- xbeg + freq_width\n");
    ksprintf(cmds, "  recend <- recbeg + (sxend - recbeg) * mapct / 100\n");
    ksprintf(cmds, "  rect(recbeg, ybeg, recend, yend, col=\"#898E8C\", border=\"#898E8C\")\n");
    ksprintf(cmds, "  nmadj <- 1\n");
    ksprintf(cmds, "  if(edpct < 10) nmadj <- 0\n");
    ksprintf(cmds, "  text(recend, ymid, mapct, cex = 0.3, adj=nmadj)\n");
    ksprintf(cmds, "  # ticks\n");
    ksprintf(cmds, "  ybeg <- yend\n");
    ksprintf(cmds, "  yend <- ybeg + 0.5 * ystep\n");
    ksprintf(cmds, "  ymid <- ybeg + 0.5 * ystep\n");
    ksprintf(cmds, "  reclen <- recend - recbeg\n");
    ksprintf(cmds, "  ticklen <- reclen / 4\n");
    ksprintf(cmds, "  for(i in 1:5){\n");
    ksprintf(cmds, "    sgx <- recbeg + ticklen * (i-1)\n");
    ksprintf(cmds, "    segments(sgx, ybeg, sgx, yend, col=\"black\")\n");
    ksprintf(cmds, "  }\n");
    ksprintf(cmds, "  ybeg <- ymid\n");
    ksprintf(cmds, "  yend <- ybeg + 0.5 * ystep\n");
    ksprintf(cmds, "  ymid <- ybeg + 0.5 * ystep\n");
    ksprintf(cmds, "  for(i in 1:5){\n");
    ksprintf(cmds, "    stx <- recbeg + ticklen * (i-1)\n");
    ksprintf(cmds, "    text(stx, ymid,(i-1)*25, cex=0.25, adj=0.5)\n");
    ksprintf(cmds, "  }\n");
    ksprintf(cmds, "  ybeg <- yend\n");
    ksprintf(cmds, "  yend <- ybeg + ystep\n");
    ksprintf(cmds, "  ymid <- ybeg + 0.5 * ystep\n");
    ksprintf(cmds, "  pxbeg <- recbeg + ticklen * 2\n");
    ksprintf(cmds, "  text(pxbeg, ymid, \"[%%]\", cex=0.3, adj=0.5)\n");
    ksprintf(cmds, "}\n");
    ksprintf(cmds, "gbinf <- dev.off()\n");
    // output to Rscript and free res
    if(outrsc && outpdf){
        FILE* fp = fopen(outrsc, "w");
        fwrite(seqs->s, sizeof(char), seqs->l, fp);
        fputc('\n', fp); free(seqs->s); free(seqs);
        fwrite(seqc->s, sizeof(char), seqc->l, fp);
        fputc('\n', fp); free(seqc->s); free(seqc);
        fwrite(seqf->s, sizeof(char), seqf->l, fp);
        fputc('\n', fp); free(seqf->s); free(seqf);
        fwrite(seqw->s, sizeof(char), seqw->l, fp);
        fputc('\n', fp); free(seqw->s); free(seqw);
        fwrite(seqi->s, sizeof(char), seqi->l, fp);
        fputc('\n', fp); free(seqi->s); free(seqi);
        fwrite(ampt->s, sizeof(char), ampt->l, fp);
        fputc('\n', fp); free(ampt->s); free(ampt);
        fwrite(cmds->s, sizeof(char), cmds->l, fp);
        // output remaining scripts
        fclose(fp);
        // tru ty run
        cmds->l = 0;
        if(n){
            char* rspath = util::which("Rscript");
            if(rspath){
                ksprintf(cmds, "%s %s", rspath, outrsc);
                kputc('\0', cmds);
                system(cmds->s);
                free(cmds->s); free(cmds);
                free(rspath);
            }
        }else{
            remove(outrsc);
        }
    }
}

void gep_opt_t::parse_cfg(){
    util::LineReader* lr = new util::LineReader(config);
    std::vector<std::string> vstr;
    std::string line;
    while(lr->getline(line)){
        util::split(line, vstr, "\t");
        if(strcmp(vstr[0].c_str(), name) == 0){
            amplicon = strdup(vstr[3].c_str());
            rlen = vstr[3].size();
            sgrbeg = atoi(vstr[5].c_str());
            sgrend = atoi(vstr[6].c_str());
            break;
        }
    }
    if(amplicon == NULL){
        fprintf(stderr, "Error configure file, does no contain the desired amplicon: %s\n", name);
        exit(EXIT_FAILURE);
    }
    delete lr;
}

bool gep_opt_t::valid_opt(){
    bool valid = true;
    if(!inbam){
        fprintf(stderr, "input bam must be provided and sorted by coordinates");
        valid = false;
    }
    if(!config){
        fprintf(stderr, "configure file must be provided");
        valid = false;
    }
    if(!name){
        fprintf(stderr, "amplicon name must be provided");
        valid = false;
    }
    if(!outdir){
        fprintf(stderr, "output directory must be provided");
        valid = false;
    }
    return valid;
}

void gep_opt_t::update_iof(){
    // io
    kstring_t ks = {0, 0, 0};
    ksprintf(&ks, "%s/%s.top%d.var.pdf", outdir, name, topn);
    outpdf = strndup(ks.s, ks.l);
    ks.l = 0;
    ksprintf(&ks, "%s/%s.tsv.gz", outdir, name);
    outtsv = strndup(ks.s, ks.l);
    ks.l = 0;
    ksprintf(&ks, "%s/%s.top%d.var.bam", outdir, name, topn);
    outbam = strndup(ks.s, ks.l);
    ks.l = 0;
    ksprintf(&ks, "%s/%s.top%d.var.R", outdir, name, topn);
    outrsc = strndup(ks.s, ks.l);
    free(ks.s);
}

void gep_opt_t::update_amp(){
    // coord
    refbeg = MAX(0, sgrbeg - flklen);
    refend = MIN(rlen-1, sgrend + flklen);
    refgot = refend - refbeg + 1;
    if(clsbpos < 0) clsbeg = sgrbeg;
    else clsbeg = MAX(sgrbeg, sgrend - clsbpos + 1);
    if(clsepos < 0) clsend = sgrend;
    else clsend = MAX(clsbeg, sgrend - clsepos + 1);
    vbeg = MAX(0, clsbeg - cutlen);
    vend = MIN(clsend + cutlen, refend);
    hrbeg = refbeg;
    hrend = refend;
    if(refgot != 100){
        int extrpl = (100-refgot)/2;
        hrbeg = MAX(0, refbeg - extrpl);
        hrend = MIN(rlen-1, refend + extrpl);
    }
}

void gep_opt_t::bam2hv(bam1_t* b, html_var_t* v){
    for(int c = 0; c < b->core.pos; ++c) kputc('N', v->seq);
    uint32_t* cigar = bam_get_cigar(b);
    int fs = 0;
    int qpos = 0, rpos = b->core.pos, oplen = 0, opint = 0;
    for(uint32_t i = 0; i < b->core.n_cigar; ++i){
        opint = bam_cigar_op(cigar[i]);
        oplen = bam_cigar_oplen(cigar[i]);
        ksprintf(v->cigar, "%d%c", oplen, (char)bam_cigar_opchr(cigar[i]));
        switch(opint){
            case BAM_CSOFT_CLIP:
                {
                    qpos += oplen;
                    break;
                }
            case BAM_CINS:
                {
                    fs += oplen;
                    assert(rpos > 0);
                    v->ipos.push_back(rpos-1);
                    v->seq->s[v->seq->l-1] = 'I';
                    kstring_t* is = (kstring_t*)calloc(1, sizeof(kstring_t));
                    for(int p = -1; p < oplen; ++p){
                        kputc(seq_nt16_str[bam_seqi(bam_get_seq(b), qpos+p)], is);
                    }
                    v->iseq.push_back(is);
                    is = NULL;
                    qpos += oplen;
                    break;
                }
            case BAM_CDEL:
                {
                    fs -= oplen;
                    for(int p = 0; p < oplen; ++p){
                        kputc('-', v->seq);
                    }
                    rpos += oplen;
                    break;
                }
            case BAM_CEQUAL: case BAM_CDIFF:
                {
                    for(int p = 0; p < oplen; ++p){
                        kputc(seq_nt16_str[bam_seqi(bam_get_seq(b), qpos+p)], v->seq);
                    }
                    rpos += oplen;
                    qpos += oplen;
                    break;
                }
            default:
                break;
        }
    }
    for(int c = 0; c < rlen - rpos; ++c) kputc('N', v->seq);
    v->cc = bam_aux2i(bam_aux_get(b, "CC"));
    v->af = (double)v->cc/(double)allrs;
    v->fscnt = fs;
    uint8_t* dcn = bam_aux_get(b, CELL_BARCODE_COUNT_TAG);
    if(dcn) v->cn = bam_aux2i(dcn);
    else v->cn = 0;
}

void geplot_usage(gep_opt_t* go, char* arg0){
    fprintf(stderr, "\nUsage: %s %s [options]\n\n", PACKAGE_NAME, arg0);
    fprintf(stderr, "Options: -i  FILE  input bam generated by getools\n");
    fprintf(stderr, "         -o  DIR   output directory\n");
    fprintf(stderr, "         -c  FILE  configure file of getools\n");
    fprintf(stderr, "         -r  STR   contig name to plot\n");
    fprintf(stderr, "         -b  INT   cutsite buffer length around sgRNA [%d]\n", go->cutlen);
    fprintf(stderr, "         -s  INT   cleavage beg position(relative to sgRNA_PAM end, negative for sgRNA_PAM beg) [%d]\n", go->clsbpos);
    fprintf(stderr, "         -e  INT   cleavage end position(relative to sgRNA_PAN end, negative for sgRNA_PAM end) [%d]\n", go->clsepos);
    fprintf(stderr, "         -t  INT   number of top variants to plot [%d]\n", go->topn);
    fprintf(stderr, "         -f  INT   flank length around sgRNA to plot [%d]\n", go->flklen);
    fprintf(stderr, "         -m        amplicon is from mitochondrion\n\n");
}

int geplot_main(int argc, char** argv){
    gep_opt_t *go = new gep_opt_t();
    if(argc == 1){
        geplot_usage(go, argv[0]);
        return 0;
    }
    int c = 0;
    while((c = getopt(argc, argv, "i:o:c:r:b:s:e:t:f:mh")) >= 0){
        switch(c){
            case 'i': go->inbam = optarg; break;
            case 'o': go->outdir = optarg; break;
            case 'c': go->config = optarg; break;
            case 'r': go->name = optarg; break;
            case 'b': go->cutlen = atoi(optarg); break;
            case 's': go->clsbpos = atoi(optarg); break;
            case 'e': go->clsepos = atoi(optarg); break;
            case 't': go->topn = atoi(optarg); break;
            case 'f': go->flklen = atoi(optarg); break;
            case 'm': go->ismt = true; break;
            case 'h': geplot_usage(go, argv[0]); return 0; break;
            default: break;
        }
    }
    if(!go->valid_opt()){
        return 1;
    }
    go->parse_cfg();
    go->update_iof();
    go->plot();
    return 0;
}
