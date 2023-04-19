#include "pscnt.h"

pscnt_opt_t::pscnt_opt_t(){
}

pscnt_opt_t::~pscnt_opt_t(){
    if(nuccnt){
        for(int i = 0; i < 5; ++i){
            free(nuccnt[i]);
        }
        free(nuccnt);
        nuccnt = NULL;
    }
    if(fai){
        fai_destroy(fai);
        fai = NULL;
    }
}

bool pscnt_opt_t::valid(){
    if(inbam.empty() || !util::exists(inbam)){
        fprintf(stderr, "input bam must be provided and exists\n");
        return false;
    }
    if(inref.empty()){
        fprintf(stderr, "input reference must be provided and exists\n");
    }
    if(len <= 0){
        fprintf(stderr, "pattern length must be positive\n");
        return false;
    }
    if(pos < 0){
        fprintf(stderr, "position must be positive\n");
        return false;
    }
    return true;
}

void pscnt_opt_t::init(){
    if(!util::exists(outdir)) util::makedir(outdir);
    outf = util::joinpath(outdir, outf);
    outs = util::joinpath(outdir, outs);
    outx = util::joinpath(outdir, outx);
    fai = fai_load(inref.c_str());
}

void pscnt_opt_t::fill_pam(){
    if(inpaf.size() && util::exists(inpaf)){
        std::ifstream fr(inpaf);
        std::string line;
        std::vector<std::string> vstr;
        while(std::getline(fr, line)){
            util::split(line, vstr, "\t");
            if(vstr.size() > 1){
                util::str2upper(vstr[1]);
                pam[vstr[1]] = vstr[0];
            }
        }
        fr.close();
    }
}

void pscnt_opt_t::ana(){
#ifdef PSCNT_DEBUG
        fprintf(stderr, "qname\tpos\t(qbeg,qend,rpos)\tpat\tcigar\tseq\n");
#endif
    fill_pam();
    samFile* fp = sam_open(inbam.c_str(), "r");
    bam_hdr_t* hdr = sam_hdr_read(fp);
    bam1_t* b = bam_init1();
    pcnt.resize(hdr->n_targets);
    ttvr.resize(hdr->n_targets, 0);
    ttpr.resize(hdr->n_targets, 0);
    KSW_FTYPE dropmask = (KSW_FHIERRLOWQ | KSW_FPRIMDIMER | KSW_FMISPRIMINN | KSW_FMISPRIMOUT | KSW_FLOWSURPT | KSW_FMANYVARS);
    uint16_t skipr = BAM_FUNMAP | BAM_FSUPPLEMENTARY | BAM_FSECONDARY | BAM_FQCFAIL;
    if(skr == 1) skipr |= BAM_FREAD1;
    if(skr == 2) skipr |= BAM_FREAD2;
    int32_t cc = 1;
    int64_t vs = 0;
    uint32_t* cigar = NULL;
    int rpos = 0, mm = 0, idl = 0, oplen = 0, opint = 0, ml = 0;
    int bpos = pos-1;
    int epos = bpos+len;
    int qbeg = 0, qend = 0;
    char* rseq = NULL;
    uint8_t* rints = NULL;
    int ltid = -1, flen = 0;
    std::string buf(len, '\0');
    while(sam_read1(fp, hdr, b) >= 0){
        if(b->core.flag & skipr) continue;
        if(!ngt){
            vs = bam_aux2i(bam_aux_get(b, "VS"));
            if(vs & dropmask) continue;
            if(!(vs & KSW_FREPSEQR)) continue;
            cc = bam_aux2i(bam_aux_get(b, "CC"));
        }
        if(ltid != b->core.tid){
            if(rseq){ free(rseq); rseq = NULL; }
            if(rints){ free(rints); rints = NULL; }
            rseq = faidx_fetch_seq(fai, hdr->target_name[b->core.tid], 0, hdr->target_len[b->core.tid], &flen);
            if(!rseq){
                fprintf(stderr, "Error load refseq [%s]\n", hdr->target_name[b->core.tid]);
                exit(EXIT_FAILURE);
            }else{
                rints = kswge_seq2ints(rseq, flen);
            }
        }
        ttvr[b->core.tid] += cc;
        rpos = b->core.pos;
        mm = ml = idl = 0;
        qbeg = qend = 0;
        cigar = bam_get_cigar(b);
        for(uint32_t i = 0; i < b->core.n_cigar; ++i){
            if(b->core.pos >= epos){
                qbeg = 1; qend = 0;
                break;
            }
            opint = bam_cigar_op(cigar[i]);
            oplen = bam_cigar_oplen(cigar[i]);
            switch(opint){
                case BAM_CDEL:
                    if(rpos <= bpos && rpos + oplen >= epos){
                        qbeg = 1; qend = 0;
                    }
                    if(rpos <= epos-1 && rpos + oplen >= epos){
                        --qend;
                    }
                    idl += oplen;
                    rpos += oplen;
                    break;
                case BAM_CINS:
                    if(rpos <= bpos){ qbeg += oplen; }
                    if(rpos < epos){ qend += oplen; }
                    idl += oplen;
                    break;
                case BAM_CSOFT_CLIP:
                    if(qbeg == 0) qbeg += oplen;
                    if(qend == 0) qend += oplen;
                    break;
                case BAM_CEQUAL:
                    if(rpos <= bpos) qbeg += MIN(bpos-rpos, oplen);
                    if(rpos < epos) qend += MIN(epos-1-rpos, oplen);
                    ml += oplen;
                    rpos += oplen;
                    break;
                case BAM_CDIFF:
                    if(rpos <= bpos) qbeg += MIN(bpos-rpos, oplen);
                    if(rpos < epos) qend += MIN(epos-1-rpos, oplen);
                    for(int prp = 0; prp < oplen; ++prp){
                        if(rints[rpos+prp] != 4) ++mm;
                    }
                    rpos += oplen;
                    break;
                default:
                    break;
            }
        }
        if(qbeg <= qend && 
           rpos >= epos && 
           ml > 0 && 
           ((!nidl) || idl == 0) && 
           ((!npidl) || (qend-qbeg+1) == len) &&
           (double)(idl+mm)/(double)(idl+mm+ml) < 1 - ppid){
            ttpr[b->core.tid] += cc;
            buf.clear();
            for(int i = qbeg; i <= qend; ++i){
                buf.push_back(seq_nt16_str[bam_seqi(bam_get_seq(b), i)]);
            }
#ifdef PSCNT_DEBUG
            if((int)buf.size() != len){
                fprintf(stderr, "%s\t%lld\t(%d,%d;%d)\t%s\t", 
                        bam_get_qname(b), 
                        b->core.pos,
                        qbeg, qend,
                        rpos,
                        buf.c_str());
                for(size_t i = 0; i < b->core.n_cigar; ++i){
                    fprintf(stderr, "%d%c", bam_cigar_oplen(cigar[i]), bam_cigar_opchr(cigar[i]));
                }
                fprintf(stderr, "\t");
                for(int i = 0; i < b->core.l_qseq; ++i){
                    fprintf(stderr, "%c", seq_nt16_str[bam_seqi(bam_get_seq(b), i)]);
                }
                fprintf(stderr, "\n");
            }
#endif
            auto iter = pcnt[b->core.tid].find(buf);
            if(iter == pcnt[b->core.tid].end()){
                pcnt[b->core.tid].insert({buf, cc});
            }else{
                iter->second += cc;
            }
        }else{
#ifdef PSCNT_DEBUG
            fprintf(stderr, "%s\t%lld\t(%d,%d;%d)\t-\t", 
                    bam_get_qname(b), 
                    b->core.pos,
                    qbeg, qend,
                    rpos);
            for(size_t i = 0; i < b->core.n_cigar; ++i){
                fprintf(stderr, "%d%c", bam_cigar_oplen(cigar[i]), bam_cigar_opchr(cigar[i]));
            }
            fprintf(stderr, "\t");
            for(int i = 0; i < b->core.l_qseq; ++i){
                fprintf(stderr, "%c", seq_nt16_str[bam_seqi(bam_get_seq(b), i)]);
            }
            fprintf(stderr, "\n");
#endif
        }
    }
    std::vector<int64_t> allop(hdr->n_targets, 0);
    std::vector<int64_t> allap(hdr->n_targets, 0);
    // output
    std::ofstream fw(outf);
    fw << "Contig\tPattern\tCount\tFreq\tAnno\n";
    for(size_t i = 0; i < pcnt.size(); ++i){
        for(auto& e: pcnt[i]){
            fw << hdr->target_name[i] << "\t" << e.first << "\t" << e.second << "\t" << (double)e.second/(double)ttpr[i] << "\t";
            auto iter = pam.find(e.first);
            if(iter != pam.end()){
                fw << iter->second << "\n";
                allap[i] += e.second;
            }else{
                fw << "-\n";
                allop[i] += e.second;
            }
        }
    }
    fw.close();
    // nuccnt
    fw.open(outx);
    fw << "Contig\tPos";
    for(int i = 0; i < len; ++i) fw << "\t" << i+1;
    fw << "\n";
    int64_t ttrp = 0;
    nuccnt = (int64_t**)malloc(len*sizeof(int64_t*));
    for(int i = 0; i < len; ++i) nuccnt[i] = (int64_t*)calloc(5, sizeof(int64_t));
    for(size_t i = 0; i < pcnt.size(); ++i){
        for(int pp = 0; pp < len; ++pp) memset(nuccnt[pp], 0, 5 * sizeof(int64_t));
        ttrp = 0;
        for(auto& e: pcnt[i]){
            for(int bi = 0; bi < len; ++bi){
                nuccnt[bi][nuc_to_3bit[(int)e.first[bi]]] += e.second;
            }
            ttrp += e.second;
        }
        for(int b = 0; b < 5; ++b){
            fw << hdr->target_name[i] << "\t" << kswge_int2nt[b];
            for(int pp = 0; pp < len; ++pp) fw << "\t" << (double)nuccnt[pp][b]/(double)ttrp;
            fw << "\n";
        }
    }
    fw.close();
    fw.open(outs);
    fw << "Contig\tTotalReads\tWithAnnoPat\tWithOtherPat\tNoPatFound\n";
    for(size_t i = 0; i < ttvr.size(); ++i){
        fw << hdr->target_name[i] << "\t" << ttvr[i] << "\t" << allap[i] << "\t" << allop[i] << "\t";
        fw << ttvr[i]-allap[i]-allop[i] << "\n";
    }
    fw.close();
    // release
    sam_close(fp);
    bam_hdr_destroy(hdr);
    bam_destroy1(b);
    if(rseq){ free(rseq); rseq = NULL; }
    if(rints){ free(rints); rints = NULL; }
}

void pscnt_usage(pscnt_opt_t* opt, char* arg0){
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: %s %s [options]\n\n", PACKAGE_NAME, arg0);
    fprintf(stderr, "Options: -i FILE  input bam(default generated by getools caledit)\n");
    fprintf(stderr, "         -f FILE  reference file used to generate input bam\n");
    fprintf(stderr, "         -a FILE  pattern annotation file\n");
    fprintf(stderr, "         -o FILE  output directory [%s]\n", opt->outdir.c_str());
    fprintf(stderr, "         -p INT   amplican begin position to look for pattern (1base included) [%d]\n", opt->pos);
    fprintf(stderr, "         -l INT   pattern length [%d]\n", opt->len);
    fprintf(stderr, "         -s FLOAT identity percentage outside of pattern [%f]\n", opt->ppid);
    fprintf(stderr, "         -r INT   skip read(1:read1,2:read2) [%d]\n", opt->skr);
    fprintf(stderr, "         -x       skip reads with indel if set\n");
    fprintf(stderr, "         -X       skip pattern with indel if set\n");
    fprintf(stderr, "         -n       input bam was not from 'getools caledit' if set\n");
    fprintf(stderr, "\n");
}

int pscnt_main(int argc, char** argv){
    pscnt_opt_t opt;
    if(argc == 1){
        pscnt_usage(&opt, argv[0]);
        return 0;
    }
    int c = 0;
    while((c = getopt(argc, argv, "i:f:a:o:p:l:s:r:nxXh")) >= 0){
        switch(c){
            case 'i': opt.inbam = optarg; break;
            case 'f': opt.inref = optarg; break;
            case 'a': opt.inpaf = optarg; break;
            case 'o': opt.outdir = optarg; break;
            case 'p': opt.pos = atoi(optarg); break;
            case 'l': opt.len = atoi(optarg); break;
            case 's': opt.ppid = atof(optarg); break;
            case 'r': opt.skr = atoi(optarg); break;
            case 'x': opt.nidl = true; break;
            case 'X': opt.npidl = true; break;
            case 'n': opt.ngt = true; break;
            case 'h': pscnt_usage(&opt, argv[0]); return 0; break;
            default: break;
        }
    }
    if(opt.valid()){
        opt.init();
        opt.ana();
        return 0;
    }else{
        return 1;
    }
}
