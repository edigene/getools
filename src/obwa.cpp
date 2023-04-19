#include "obwa.h"

void obwa_usage(char* arg0){
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: %s [options]\n\n", arg0);
    fprintf(stderr, "Options: -r FILE reference fasta file\n");
    fprintf(stderr, "         -i FILE reference index file\n");
    fprintf(stderr, "         -q STR  query sequence 1\n");
    fprintf(stderr, "         -Q STR  query sequence 2\n");
    fprintf(stderr, "         -m INT  match match score [1]\n");
    fprintf(stderr, "         -M INT  mismatch score [4]\n");
    fprintf(stderr, "         -o INT  gap open score [6]\n");
    fprintf(stderr, "         -e INT  gap extend score [6]\n");
    fprintf(stderr, "         -c INT  5' and 3' clip score [5]\n");
    fprintf(stderr, "         -s INT  seed length [19]\n");
    fprintf(stderr, "         -t INT  min output score [30]\n");
    fprintf(stderr, "         -f      only build forward strand index [false]\n");
    fprintf(stderr, "\n");
}

int obwa_main(int argc, char** argv){
    if(argc == 1){
        obwa_usage(argv[0]);
        return 0;
    }
    char* rfile = NULL;
    char* index = NULL;
    char* qseq1 = NULL;
    char* qseq2 = NULL;
    int match = 1;
    int mism = 4;
    int gapo = 6;
    int gape = 6;
    int clip = 5;
    int c = -1;
    int seedl = 19;
    int minscr = 30;
    while((c = getopt(argc, argv, "r:i:q:Q:m:M:o:s:e:c:t:h")) >= 0){
        switch(c){
            case 'r': rfile = optarg; break;
            case 'i': index = optarg; break;
            case 'q': qseq1 = optarg; break;
            case 'Q': qseq2 = optarg; break;
            case 'm': match = atoi(optarg); break;
            case 'M': mism = atoi(optarg); break;
            case 'o': gapo = atoi(optarg); break;
            case 'e': gape = atoi(optarg); break;
            case 'c': clip = atoi(optarg); break;
            case 's': seedl = atoi(optarg); break;
            case 't': minscr = atoi(optarg); break;
            case 'h': obwa_usage(argv[0]); return 0; break;
            default: break;
        }
    }
       if(rfile == NULL && index == NULL){
        fprintf(stderr, "index file or ref file must be provided\n");
        return 1;
    }
    if(qseq1 == NULL){
        fprintf(stderr, "first query must be provided\n");
        return 1;
    }
    // get index
    obwa_t* obwa = new obwa_t();
    obwa->setMismatchPenalty(mism);
    obwa->setGapOpenPenalty(gapo);
    obwa->setGapExtendPenalty(gape);
    obwa->set3PrimeClipPenalty(clip);
    obwa->set5PrimeClipPenalty(clip);
    obwa->setMinOutScore(minscr);
    obwa->setMatchScore(match); // after all penalty set
    obwa->fillScoreMatrix(); // after all score/penalty set
    obwa->setMinSeedLength(seedl);
    if(rfile){
        std::vector<krec1_t*> recs;
        gzFile ifp = gzopen(rfile, "r");
        krec1_t* rec = krec1_init();
        kseq1_t* seq = kseq1_init(ifp);
        while(kseq1_read(seq, rec) >= 0){
            recs.push_back(rec);
            rec = krec1_init();
        }
        krec1_destroy(rec);
        kseq1_destroy(seq);
        gzclose(ifp);
        if(recs.empty()){
            fprintf(stderr, "No reference seuqnces in %s, aborted\n", rfile);
            return 1;
        }
        obwa->buildIndex(recs);
    }else{
        obwa->loadIndex(index);
    }
    if(!obwa->bwaidx){
        fprintf(stderr, "Index construction failed\n");
        return 1;
    }else{
        fprintf(stderr, "Index construction finished\n");
    }
    // align work
    std::vector<mem_aln_t> rets;
    std::vector<char*> seqs;
    seqs.push_back(qseq1);
    if(qseq2) seqs.push_back(qseq2);
    for(auto& qseq: seqs){
        int qlen = strlen(qseq);
        if(qlen < obwa->getMinSeedLength()) obwa->setMinSeedLength(qlen);
#ifdef ALN_SEQ_TO_MEM
        obwa->alignSeq2mem(qseq, qlen, rets);
#else
        obwa->alignSeq2set(qseq, qlen, rets, 4);
#endif
        bntseq_t* bs = obwa->bwaidx->bns;
        fprintf(stderr, "qseq alignment results\n");
        int i = 0;
        for(auto& b: rets){
            fprintf(stderr, "#    : %d\n", i++);
            fprintf(stderr, "rid  : %d\n", b.rid);
            fprintf(stderr, "rname: %s\n", bs->anns[b.rid].name);
            fprintf(stderr, "rlen : %d\n", bs->anns[b.rid].len);
            fprintf(stderr, "rpos : %lld\n", b.pos);
            fprintf(stderr, "NM   : %d\n", b.NM);
            fprintf(stderr, "cigar: ");
            for(int i = 0; i < b.n_cigar; ++i){
                fprintf(stderr, "%d%c", b.cigar[i] >> 4, "MIDSH"[b.cigar[i] & 0xf]);
            }
            fprintf(stderr, "\n");
            fprintf(stderr, "score: %d\n", b.score);
            fprintf(stderr, "mapQ : %d\n", b.mapq);
            fprintf(stderr, "flag : %d\n", b.flag);
            fprintf(stderr, "isrev: %d\n", b.is_rev);
            fprintf(stderr, "XA   : %s\n", b.XA);
            free(b.cigar);
            free(b.XA);
        }
    }
    // cleanup
    delete obwa;
    return 0;
}

#ifdef OBWA_MAIN_FUN
int main(int argc, char** argv){
    return obwa_main(argc, argv);
}
#endif
