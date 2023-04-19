#include "options.h"

Options::Options(){
    soft = new Software();
    outdir = "./out/";
    hapdir = "hapcnt";
    thread = 4;
    f4bout = GEOUT_SIMP_BAM;
    maxpack = 3;
    maxreads = 100000;
    compression = 1;
    maxmm = 2;
    maxbo = 0;
    maxoff = 2;
    maxpl = 0;
    maxml = 0;
    minlen = 1;
    droplib = "DROP";
    tpl = NULL;
    tpdetail = false;
    useHash = false;
    droppre = false;
    outspl = false;
    sematch = false;
    usesem = false;
    adjfr = false;
    ostrand = false;
    recstat = false;
    isq64 = false;
    fq2bam = false;
    fq2cal = false;
    fq2spl = false;
    memlow = false;
    dmatch = false;
    memone = 0;
    s8aln = false;
    dosccal = false;
    preidx = NULL;
    ampidx = NULL;
    bcfh = NULL;
    bamh = NULL;
    skipfr1 = 0;
    skipfr2 = 0;
}

Options::~Options(){
    cleanup();
}

void Options::cleanup(){
    if(soft){ delete soft; soft = NULL; };
    if(tpl){ delete tpl; tpl = NULL; }
    for(auto& e: amplicons){
        if(e){ delete e; e = NULL; }
    }
    for(auto& e: prl){
        if(e){ delete e; e = NULL; }
    }
    if(preidx) delete preidx;
    if(ampidx) delete ampidx;
    if(bcfh){ bcf_hdr_destroy(bcfh); bcfh = NULL; }
    if(bamh){ bam_hdr_destroy(bamh); bamh = NULL; }
}

void Options::update(int argc, char** argv){
    outdir = util::abspath(outdir);
    if(!util::exists(outdir)){
        util::makedir(outdir);
    }
    if(!fq2spl){
        if(tpdetail){
            pltdir = outdir + "/detail/";
            if(!util::exists(pltdir)){
                util::makedir(pltdir);
            }
        }
        tmpdir = outdir + "/.tmp";
        if(memlow){
            util::makedir(tmpdir);
            memone = MIN(thread, memone);
        }
    }
    htmldir = outdir + "/" + hmo.subdir;
    if(fq2cal || fq2spl){
        if(!util::exists(htmldir)) util::makedir(htmldir);
    }
    jsr = util::abspath(outdir + "/split.json");
    tsr = util::abspath(outdir + "/" + GETOOLS_SPLITTSV);
    tsm = util::abspath(outdir + "/" + GETOOLS_SPLITCMBTSV);
    hmr = util::abspath(outdir + "/report.html");
    jss = util::abspath(outdir + "/sublib.json");
    tss = util::abspath(outdir + "/" + GETOOLS_SUBLIBTSV);
    txx_se = util::abspath(outdir + "/fbdiff_se.tsv");
    txx_pe = util::abspath(outdir + "/fbdiff_pe.tsv");
    soft->update(argc, argv);
    hmo.tnfigh = 20 * edo.topn + 100;
    if(edo.hapcnt){
        hapdir = util::joinpath(outdir, hapdir);
        util::makedir(hapdir);
    }
}

bool Options::parseFq2bamCfg(){
    // parse prefix cfg
    std::ifstream fr(prf);
    std::string line;
    std::vector<std::string> vstr;
    std::map<std::string, int> smps;
    int count = 0, idx = 0, bfo = -1, bro = -1;
    while(std::getline(fr, line)){
        util::split(line, vstr, "\t");
        if(vstr.size() < 4){
            fprintf(stderr, "ERROR: input configure file format wrong, it must be at least 4-column TSV format\n");
            fprintf(stderr, "the 5 required fields are: name,fwd_primer,rev_primer,amplicon\n");
            fr.close();
            return false;
        }
        util::nuc2upper(vstr[1]); // fwd primer
        util::nuc2upper(vstr[2]); // rev primer
        util::nuc2upper(vstr[3]);
        // check fwd, rev
        if(!checkps(vstr[1], vstr[2], vstr[3], bfo, bro)){
            if(checkps(vstr[2], vstr[1], vstr[3], bfo, bro)){
                std::swap(vstr[1], vstr[2]);
            }else{
                if(!s8aln){
                    fprintf(stderr, "ERROR: error primer sequences of amplicon: %s\n", vstr[0].c_str());
                    fr.close();
                    return false;
                }
            }
        }
        vblens.insert(vstr[1].size());
        vblens.insert(vstr[2].size());
        if(maxpl < vstr[1].size()) maxpl = vstr[1].size();
        if(maxpl < vstr[2].size()) maxpl = vstr[2].size();
        auto iter = smps.find(vstr[0]);
        if(iter == smps.end()){
            idx = count;
            smps[vstr[0]] = count++;
            samples.push_back(vstr[0]);
            amplicon_t* amp = new amplicon_t();
            amp->aname = vstr[0];
            amp->aseq = vstr[3];
            amplicons.push_back(amp);
        }else{
            idx = iter->second;
        }
        prefix_t* pr = new prefix_t(vstr[1], vstr[2], '+', idx);
        pr->fbo= bfo;
        pr->rbo = bro; 
        pr->seidx = false;
        prl.push_back(pr); prs.push_back(pr);
        prefix_t* rp = new prefix_t(vstr[2], vstr[1], '-', idx);
        rp->fbo = bro;
        rp->rbo = bfo;
        rp->seidx = false;
        prl.push_back(rp);
    }
    fr.close();
    std::sort(prl.begin(), prl.end(), prefix_sort_t());
    // set drop lib
    dropidx = samples.size();
    samples.push_back(droplib);
    if(samples.size() > MAX_SAMPLE_IN_ONE_LIB){
        fprintf(stderr, "Too many sample/contigs to be split(%ld), max is: %d\n", samples.size(), MAX_SAMPLE_IN_ONE_LIB);
        return false;
    }
    return true;
}

bool Options::parsePreCfg(){
    // parse prefix cfg
    std::ifstream fr(prf);
    std::string line;
    std::vector<std::string> vstr;
    std::map<std::string, int> smps;
    int count = 0, idx = 0, bfo = -1, bro = -1;
    while(std::getline(fr, line)){
        util::split(line, vstr, "\t");
        if(vstr.size() < 5){
            fprintf(stderr, "ERROR: input configure file format wrong, it must be at least 5-column TSV format\n");
            fprintf(stderr, "the 5 required fields are: name,fwd_primer,rev_primer,amplicon,sgRNA_PAM\n");
            fprintf(stderr, "the 6th optional filed is: donor\n");
            fr.close();
            return false;
        }
        util::nuc2upper(vstr[1]); // fwd primer
        util::nuc2upper(vstr[2]); // rev primer
        util::nuc2upper(vstr[3]);
        util::nuc2upper(vstr[4]);
        // check fwd, rev
        if(!checkps(vstr[1], vstr[2], vstr[3], bfo, bro)){
            if(checkps(vstr[2], vstr[1], vstr[3], bfo, bro)){
                std::swap(vstr[1], vstr[2]);
            }else{
                fprintf(stderr, "ERROR: error primer sequences of amplicon: %s\n", vstr[0].c_str());
                fr.close();
                return false;
            }
        }
        vblens.insert(vstr[1].size());
        vblens.insert(vstr[2].size());
        if(maxpl < vstr[1].size()) maxpl = vstr[1].size();
        if(maxpl < vstr[2].size()) maxpl = vstr[2].size();
        auto iter = smps.find(vstr[0]);
        if(iter == smps.end()){
            idx = count;
            smps[vstr[0]] = count++;
            samples.push_back(vstr[0]);
            amplicon_t* amp = new amplicon_t();
            amp->aname = vstr[0];
            amp->aseq = vstr[3];
            amp->sgrseq = vstr[4];
            if(vstr.size() >= 6){ // donor provided
                amp->donor = vstr[5];
                util::nuc2upper(amp->donor);
                recstat = true;
            }
            // adjust sgRNA sequence same as amplicon
            size_t sgbeg = amp->aseq.find(amp->sgrseq);
            if(sgbeg == std::string::npos){// fwd strand not found
                util::revComp2OriSeq(amp->sgrseq); // reverse compliment sgrna
                sgbeg = amp->aseq.find(amp->sgrseq); // another try
            }
            if(sgbeg == std::string::npos){// fwd/rev both not found, something fatal 
                fprintf(stderr, "sgRNA_PAM of %s not found in its amplicon sequence forward or reverse strand, program aborted.\n", amp->aname.c_str());
                delete amp;
                fr.close();
                return false;
            }
            amp->sgrbeg = sgbeg;
            amp->sgrend = sgbeg + amp->sgrseq.size() - 1;
            amp->maxdel = amp->aseq.size();
            amp->fplen = vstr[1].size()-bfo;
            amp->rvlen = vstr[2].size()-bro;
            if(vstr[1].size() + vstr[2].size() + edo.minextram < amp->aseq.size()){
                amp->maxdel = amp->aseq.size() - (vstr[1].size() + vstr[2].size() + edo.minextram);
            }
            // compute clsbeg and clsend
            if(edo.clsbpos < 0){
                amp->clsbeg = amp->sgrbeg;
            }else{
                amp->clsbeg = MAX(amp->sgrbeg, amp->sgrend - edo.clsbpos + 1);
            }
            if(edo.clsepos < 0){
                amp->clsend = amp->sgrend;
            }else{
                amp->clsend = MAX(amp->sgrend - edo.clsepos + 1, amp->clsbeg);
            }
            amp->maxrpl = util::longestRepeatLen(amp->aseq.c_str(), amp->aseq.size());
            // test donor strand info and adjust it to be the same as amplicon
            if(!donor2amp(amp)){
                fprintf(stderr, "whole donor template of %s construction failed.\n", amp->aname.c_str());
                delete amp;
                fr.close();
                return false;
            }
            amplicons.push_back(amp);
        }else{
            idx = iter->second;
        }
        prefix_t* pr = new prefix_t(vstr[1], vstr[2], '+', idx);
        pr->fbo= bfo;
        pr->rbo = bro; 
        pr->seidx = false;
        prl.push_back(pr); prs.push_back(pr);
        prefix_t* rp = new prefix_t(vstr[2], vstr[1], '-', idx);
        rp->fbo = bro;
        rp->rbo = bfo;
        rp->seidx = false;
        prl.push_back(rp);
    }
    fr.close();
    std::sort(prl.begin(), prl.end(), prefix_sort_t());
    // set drop lib
    dropidx = samples.size();
    samples.push_back(droplib);
    if(samples.size() > MAX_SAMPLE_IN_ONE_LIB){
        fprintf(stderr, "Too many sample/contigs to be split(%ld), max is: %d\n", samples.size(), MAX_SAMPLE_IN_ONE_LIB);
        return false;
    }
    return true;
}

void Options::setdm(){
    dmatch = false;
    if(maxmm == 0) return;
    std::set<std::string> dms;
    for(auto& p: prs){
        if(dmatch) break;
        if(p->fseq.size()){
            if(dms.find(p->fseq) != dms.end()){
                dmatch = true;
                util::loginfo(stderr, "duplicated index exists, slow match will be used");
                break;
            }else{
                dms.insert(p->fseq);
            }
        }
        if(p->rseq.size()){
            if(dms.find(p->rseq) != dms.end()){
                dmatch = true;
                util::loginfo(stderr, "duplicated index exists, slow match will be used");
                break;
            }else{
                dms.insert(p->rseq);
            }
        }
    }
}

bool Options::parseSplCfg(){
    // parse prefix cfg
    std::ifstream fr(prf);
    std::string line;
    std::vector<std::string> vstr;
    std::map<std::string, int> smps;
    int count = 0, idx = 0;
    bool seidx = true;
    while(std::getline(fr, line)){
        util::split(line, vstr, "\t");
        if(vstr.size() < 2){
            fprintf(stderr, "input configure file format wrong, it must be at least 2-column TSV format\n");
            fprintf(stderr, "the 2 fields are: name,primer(or name,amplicon if -g set)\n");
            fprintf(stderr, "in the case your library is constructed with double index, 3-column TSV is needed\n");
            fprintf(stderr, "the 3 fields are: name,fwd_primer,rev_primer(or name,fwd_primer,rev_primer,amplicon if -g set)\n");
            fr.close();
            return false;
        }
        if(s8aln){
            auto iter = smps.find(vstr[0]);
            if(iter == smps.end()){
                idx = count;
                smps[vstr[0]] = count++;
                samples.push_back(vstr[0]);
                amplicon_t* amp = new amplicon_t();
                amp->aname = vstr[0];
                if(vstr.size() == 2) amp->aseq = vstr[1];
                else if(vstr.size() >= 4) amp->aseq = vstr[3];
                amplicons.push_back(amp);
            }else{
                idx = iter->second;
            }
        }else{
            seidx = vstr.size() == 2;
            util::nuc2upper(vstr[1]);
            if(seidx){
                vstr.push_back(util::revComp2NewSeq(vstr[1]));
            }else{
                util::nuc2upper(vstr[2]);
            }
            vblens.insert(vstr[1].size());
            if(!seidx) vblens.insert(vstr[2].size());
            if(maxpl < vstr[1].size()) maxpl = vstr[1].size();
            if((!seidx) && maxpl < vstr[2].size()) maxpl = vstr[2].size();
            auto iter = smps.find(vstr[0]);
            if(iter == smps.end()){
                idx = count;
                smps[vstr[0]] = count++;
                samples.push_back(vstr[0]);
            }else{
                idx = iter->second;
            }
            if(seidx){
                prefix_t* pr = new prefix_t(vstr[1], "", '+', idx);
                pr->seidx = seidx;
                prl.push_back(pr); prs.push_back(pr);
            }else{
                prefix_t* pr = new prefix_t(vstr[1], vstr[2], '+', idx);
                pr->seidx = seidx;
                prl.push_back(pr); prs.push_back(pr);
                prefix_t* rp = new prefix_t(vstr[2], vstr[1], '-', idx);
                rp->seidx = seidx;
                prl.push_back(rp);
            }
        }
        std::sort(prl.begin(), prl.end(), prefix_sort_t());
    }
    fr.close();
    // set drop lib
    dropidx = samples.size();
    samples.push_back(droplib);
    if(samples.size() > MAX_SAMPLE_IN_ONE_LIB){
        fprintf(stderr, "Too many sample/contigs to be split(%ld), max is: %d\n", samples.size(), MAX_SAMPLE_IN_ONE_LIB);
        return false;
    }
    return true;
}

void Options::genHash(){
    if(maxpl <= 32) useHash = true;
    prh.resize(maxpl+1);
    psh.resize(maxpl+1);
    for(auto& e: prl){
        auto iter = prh[e->fseq.size()].find(e->fhash);
        if(iter == prh[e->fseq.size()].end()) prh[e->fseq.size()][e->fhash] = {e};
        else iter->second.push_back(e);
        auto xter = psh[e->fseq.size()].find(e->fseq);
        if(xter == psh[e->fseq.size()].end()) psh[e->fseq.size()][e->fseq] = {e};
        else xter->second.push_back(e);
    }
}

void Options::init(int argc, char** argv){
    update(argc, argv);
    bool pg = false;
    if(fq2spl) pg = parseSplCfg();
    else if(fq2bam) pg = parseFq2bamCfg();
    else if(fq2cal) pg = parsePreCfg();
    if(!pg){
        cleanup();
        exit(EXIT_FAILURE);
    }
    if(samples.size() > INT16_MAX){
        fprintf(stderr, "Max amplicon allowed is %d\n", INT16_MAX);
        cleanup();
        exit(EXIT_FAILURE);
    }
    if(s8aln){
        ampidx  = new amp_idx_t();
        ampidx->init(amplicons, sematch);
    }else{
        if(maxmm){
            if(maxpl > MIN_SEED_LEN){
                setdm();
                preidx = new pre_idx_t();
                preidx->init(prs, maxmm, maxoff, sematch, droppre);
            }else{
                util::loginfo(stderr, "primers too short[%d], only use exact match", maxpl);
                maxmm = 0;
            }
        }
    }
    maxml = maxpl + maxoff;
    minlen = maxml;
    genHash();

    tpl = new ThreadPool(thread);
    if(fq2bam || fq2cal){ // bam hdr
        // bam header
        kstring_t* ks = (kstring_t*)calloc(1, sizeof(kstring_t));
        for(size_t i = 0; i < amplicons.size(); ++i){
            ksprintf(ks, "@SQ\tSN:%s\tLN:%lu\n", amplicons[i]->aname.c_str(), amplicons[i]->aseq.size());
        }
        bamh = sam_hdr_parse(ks->l, ks->s);
        bamh->l_text = ks->l;
        bamh->text = ks->s;
        free(ks);
    }
    if(fq2cal){
        // bcf header
        bcfh = bcf_hdr_init("w");
        kstring_t str = {0, 0, 0};
        // reference
        for(int i = 0; i < sam_hdr_nref(bamh); ++i){
            str.l = 0;
            int geb = MAX(0, amplicons[i]->clsbeg - edo.cutbuflen);
            int ged = MIN(amplicons[i]->clsend + edo.cutbuflen, (int)amplicons[i]->aseq.size() - 1);
            ksprintf(&str, "##contig=<ID=%s,length=%" PRId64 ",geb=%d, gee=%d>",
                     sam_hdr_tid2name(bamh, i), (int64_t)sam_hdr_tid2len(bamh, i),
                     geb, ged);
            bcf_hdr_append(bcfh, str.s);
        }
        free(str.s);
        bcf_hdr_append(bcfh, "##INFO=<ID=VT,Number=1,Type=Integer,Description=\"Variant type(1:SNP,2:INS,4:DEL,8:DIN)\">");
        bcf_hdr_append(bcfh, "##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">");
        bcf_hdr_append(bcfh, "##FORMAT=<ID=AD,Number=A,Type=Integer,Description=\"Allelic depths\">");
        bcf_hdr_append(bcfh, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Raw depth(all valid reads included)\">");
        if(dosccal) bcf_hdr_append(bcfh, "##FORMAT=<ID=SC,Number=1,Type=String,Description=\"Single cells whith have shis allele\">");
        bcf_hdr_add_sample(bcfh, "case");
        if(ctrl.size()) bcf_hdr_add_sample(bcfh, "ctrl");
        if(bcfh->dirty){
            if(bcf_hdr_sync(bcfh) < 0){
                fprintf(stderr, "Error in sync bcf header\n");
                cleanup();
                exit(EXIT_FAILURE);
            }
        }
    }
    hmo.init(hrjsn);
}

bool Options::valid(){
    if(prf.empty()){
        fprintf(stderr, "configure file must be provided\n");
        return false;
    }
    if(!util::exists(prf)){
        fprintf(stderr, "configure file '%s' does not exist\n", prf.c_str());
        return false;
    }
    if(in1.empty()){
        fprintf(stderr, "input read1 must be provided\n");
        return false;
    }
    if(!util::exists(in1)){
        fprintf(stderr, "input fastq read1 '%s'  does not existn\n", in1.c_str());
        return false;
    }
    if(in2.size() && !util::exists(in2)){
        fprintf(stderr, "input fastq read2 '%s'  does not existn\n", in2.c_str());
        return false;
    }
    if(thread < 2){
        fprintf(stderr, "threads(%d) must be more than 2\n", thread);
        return false;
    }
    if(maxpack < 2 || maxpack > 10){
        fprintf(stderr, "max packs(%d)must be in range[2, 10]\n", maxpack);
        return false;
    }
    if(compression < 0 || compression > 9){
        fprintf(stderr, "compression level(%d) must be in range[0, 9]\n", compression);
        return false;
    }
    if(maxreads < 10000 || maxreads > 1000000){
        fprintf(stderr, "max reads in each pack(%d) must be in range[10000, 100000]\n", maxreads);
        return false;
    }
    if(maxoff < 0 || maxoff > MAX_OFF_ALLOWED){
        fprintf(stderr, "ERROR: max offset(%d) before looking for match must be in range[0, 2]\n", maxoff);
        return false;
    }
    if(maxmm < 0 || maxmm > 5){
        fprintf(stderr, "ERROR: max mismatch(%d) allowed for good match must be in range[0, 5]\n", maxmm);
        return false;
    }
    if(memone < 0 || memone > thread){
        fprintf(stderr, "ERROR: max number of samples processed one timei(%d) must be positive and less than working threads(%d)\n", memone, thread);
        return false;
    }
    // check qual options
    if(lowq.maxLowQualFrac < 0 || lowq.maxLowQualFrac > 1){
        fprintf(stderr, "ERROR: max lowquality base fraction(%f) allowed muste be in range[0, 1]\n", lowq.maxLowQualFrac);
        return false;
    }
    // alignment options
    if(aln.match < 0 || aln.mismatch < 0 || aln.gapopen < 0 || aln.gapext < 0){
        fprintf(stderr, "ERROR: all alignment score should be their absolute positive value\n");
        return false;
    }
    // edit event
    if(edo.minaf < 0 || edo.minaf > 1){
        fprintf(stderr, "ERROR: min af(%f) allowed for valid edit events must be in range[0, 1]\n", edo.minaf);
        return false;
    }
    if(edo.cutbuflen < 0){
        fprintf(stderr, "ERROR: extra range(%d) around sgRNA to compute variants should be positive\n", edo.cutbuflen);
        return false;
    }
    if(edo.flklen < 0){
        fprintf(stderr, "ERROR: flank length(%d) to plot edit event must be positive\n", edo.flklen);
        return false;
    }
    if(edo.minextram < 0){
        fprintf(stderr, "ERROR: min extrai(%d) match needed for large indel must be positive\n", edo.minextram);
        return false;
    }

    if(edo.topn < 1){
        fprintf(stderr, "ERROR: topn(%d) events to plot must be positive\n", edo.topn);
        return false;
    }

    if(edo.vartypem < 0){
        fprintf(stderr, "ERROR: variant type masks must be positive\n");
        return false;
    }
    if(edo.rectypem < 0 || edo.rectypem > 3){
        fprintf(stderr, "ERROR: valid edited recombination type must be in range[0, 3]\n");
        return false;
    }
    if(maxbo < 0){
        fprintf(stderr, "ERROR: invalid max barcode length, must be positive\n");
        return false;
    }
    return true;
}

bool Options::donor2amp(amplicon_t* amp){
    if(amp->donor.empty()) return true;
    int8_t* score_mat = kswge_gen_smat(aln.match, aln.mismatch);
    uint8_t* refints = kswge_seq2ints(amp->aseq.c_str(), amp->aseq.size());
    uint8_t* donorintsf = kswge_seq2ints(amp->donor.c_str(), amp->donor.size());
    char* donorrev = util::revComp2NewSeq((char*)amp->donor.c_str(), amp->donor.size());
    uint8_t* donorintsr = kswge_seq2ints(donorrev, amp->donor.size());
    kswr_t *retf = kswge_semi_global(amp->donor.size(), donorintsf, amp->aseq.size(), refints, 5, score_mat, aln.gapopen, aln.gapext);
    kswr_t *retr = kswge_semi_global(amp->donor.size(), donorintsr, amp->aseq.size(), refints, 5, score_mat, aln.gapopen, aln.gapext);
    kswr_t *bestaln = retr;
    char* bestseq = donorrev;
    uint8_t* bestints = donorintsr;
    if(retf->score > retr->score){
        bestaln = retf;
        bestseq = (char*)amp->donor.c_str();
        bestints = donorintsf;
    }
    kswge_mark_mismatch(bestaln, refints, amp->aseq.size(), bestints, amp->donor.size());
    amp->dins.resize(amp->aseq.size()+1);
    amp->ddis.resize(amp->aseq.size()+1);
    amp->ddel.resize(amp->aseq.size()+1, 0);
    amp->dsnv.resize(amp->aseq.size()+1, -1);
    amp->dmutcnt = 0;
    kstring_t wds = {0, 0, 0};
    char opchr;
    int oplen;
    int idl = 0;
    int rpos = bestaln->tb;
    int qpos = bestaln->qb;
    for(int c = 0; c < bestaln->tb; ++c) kputc(amp->aseq[c], &wds); // leading amplicon
    for(int c = 0; c < bestaln->ncigar; ++c){
        opchr = kswge_cigar_opchr(bestaln->cigar[c]);
        oplen = kswge_cigar_oplen(bestaln->cigar[c]);
        if(opchr == 'I'){
            for(int i = 0; i < oplen; ++i){
                kputc(bestseq[qpos+i], &wds);
                amp->dins[rpos].append(1, bestseq[qpos+i]);
            }
            ++amp->dmutcnt;
            qpos += oplen;
            idl += oplen;
        }else if(opchr == 'M' || opchr == 'X' || opchr == '='){
            for(int i = 0; i < oplen; ++i) kputc(bestseq[qpos+i], &wds);
            if(opchr == 'X'){
                if(oplen == 1) amp->dsnv[rpos] = bestints[qpos];
                else{
                    for(int i = 0; i < oplen; ++i) amp->ddis[rpos].append(1, bestseq[qpos+i]);
                }
                ++amp->dmutcnt;
            }
            qpos += oplen;
            rpos += oplen;
        }else if(opchr == 'D'){
            amp->ddel[rpos] = oplen;
            ++amp->dmutcnt;
            rpos += oplen;
            idl += oplen;
        }else if(opchr == 'S'){
            qpos += oplen;
        }
    }
    int extral = amp->aseq.size() - rpos;
    for(size_t c = rpos; c < amp->aseq.size(); ++c) kputc(amp->aseq[c], &wds); // tailing amplicon
    amp->donor = wds.s;
    amp->dbeg = bestaln->tb;
    amp->denda = bestaln->te;
    amp->dendd = wds.l - extral - 1;
    // free up
    free(score_mat); score_mat = NULL;
    free(refints); refints = NULL;
    free(donorintsf); donorintsf = NULL;
    free(donorrev); donorrev = NULL;
    free(donorintsr); donorintsr = NULL;
    kswr_destroy(retf); retf = NULL;
    kswr_destroy(retr); retr = NULL;
    if(wds.s){
        free(wds.s);
        if(edo.hapcnt && idl){
            fprintf(stderr, "Error, no indel allowed for amplicon[%s] hapcnt analysis\n", amp->aname.c_str());
            return false;
        }else{
            return true;
        }
    }else{
        return false;
    }
}

bool Options::checkps(const std::string& pfwd, const std::string& prev, const std::string& amps, int& bfo, int &bro){
    bfo = bro = -1;
    for(int i = 0; i <= maxbo; ++i){
        if(amps.find(pfwd.substr(i)) == 0){
            bfo = i;
            break;
        }
    }
    if(bfo < 0) return false;
    else{
        std::string rprev = util::revComp2NewSeq(prev);
        for(int i = 0; i <= maxbo; ++i){
            if(amps.rfind(rprev.substr(0, rprev.size()-i)) == amps.size()-rprev.size()+i){
                bro = i;
                break;
            }
        }
        return bro >= 0;
    }
}
