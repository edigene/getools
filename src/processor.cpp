#include "processor.h"

Processor::Processor(Options* opt){
    mOptions = opt;
    mBamHeader = mOptions->bamh;
    mBcfHeader = mOptions->bcfh;
    mIsPE = mOptions->in2.size();
    mReadGenerator = NULL;
    if(mIsPE) mReadGenerator = new ReadGenerator(mOptions->in1.c_str(), mOptions->in2.c_str());
    else mReadGenerator = new ReadGenerator(mOptions->in1.c_str());
    mReadGenerator->init(mOptions->maxpack, mOptions->maxreads);
    mSpliters = (Spliter**)calloc(mOptions->thread, sizeof(Spliter*));
    for(int i = 0; i < mOptions->thread; ++i){
        mSpliters[i] = new Spliter(mOptions);
        mSpliters[i]->init();
    }
    if(mOptions->outspl){
        mWriter1 = (Writer**)calloc(mOptions->samples.size(), sizeof(Writer*));
        for(size_t i = 0; i < mOptions->samples.size(); ++i) mWriter1[i] = NULL;
        mWriter2 = NULL;
        if(mIsPE){
            mWriter2 = (Writer**)calloc(mOptions->samples.size(), sizeof(Writer*));
            for(size_t i = 0; i < mOptions->samples.size(); ++i) mWriter2[i] = NULL;
        }
    }
    mPeMergers = NULL;
    mDeRepers = NULL;
    mGEDetectors = NULL;
    mPAligners = NULL;
    mScAna = NULL;
    if(mOptions->fq2bam){
        mPAligners = (PAligner**)calloc(mOptions->amplicons.size(), sizeof(PAligner*));
        for(size_t i = 0; i < mOptions->amplicons.size(); ++i){
            mPAligners[i] = new PAligner(mOptions->amplicons[i], mOptions);
            mPAligners[i]->rtid = i;
        }
        mOutRef = mOptions->outdir + "/ref.fa";
        mOutTotBam = mOptions->outdir + "/split.bam";
    }
    if(mOptions->fq2cal){
        mPeMergers = (PairMerger**)calloc(mOptions->amplicons.size(), sizeof(PairMerger*));
        for(size_t i = 0; i < mOptions->amplicons.size(); ++i){
            mPeMergers[i] = new PairMerger(mOptions);
            if(mPeMergers[i]->opt_fastq_minovlen < mOptions->amplicons[i]->maxrpl){
                mPeMergers[i]->opt_fastq_minovlen = mOptions->amplicons[i]->maxrpl;
            }
            if(mPeMergers[i]->opt_olpm_minolplen < mOptions->amplicons[i]->maxrpl){
                mPeMergers[i]->opt_olpm_minolplen = mOptions->amplicons[i]->maxrpl;
            }
        }
        mDeRepers = (DeReper**)calloc(mOptions->amplicons.size(), sizeof(DeReper*));
        for(size_t i = 0; i < mOptions->amplicons.size(); ++i){
            mDeRepers[i] = new DeReper(opt);
        }
        mGEDetectors = (GEDetector**)calloc(mOptions->amplicons.size(), sizeof(GEDetector*));
        for(size_t i = 0; i < mOptions->amplicons.size(); ++i){
            mGEDetectors[i] = new GEDetector(mOptions->amplicons[i], mOptions);
            mGEDetectors[i]->rtid = i;
        }
        mTopnHtml.resize(mOptions->samples.size());
        for(size_t i = 0; i < mOptions->samples.size(); ++i){
            mTopnHtml[i] = (kstring_t*)calloc(1, sizeof(kstring_t));
        }
        mOutRef = mOptions->outdir + "/ref.fa";
        mOutTotBam = mOptions->outdir + "/var.bam";
        mOutFulBam = mOptions->outdir + "/ful.bam";
        mOutTotBcf = mOptions->outdir + "/var.bcf";
        mOutTotSCJSON = mOptions->outdir + "/sccal.json";
        if(mOptions->dosccal){
            mScAna = new scana_opt_t();
            mScAna->inbcf = mOutTotBcf;
            mScAna->outdir = mOptions->outdir;
            mScAna->hrjsn = mOptions->hrjsn;
            mScAna->topn = mOptions->edo.topn;
        }
    }
    mOutHTMLs.resize(mOptions->samples.size());
    for(size_t i = 0; i < mOptions->samples.size(); ++i){
        mOutHTMLs[i] = mOptions->htmldir + "/" + mOptions->samples[i] + ".html";
    }
    mQcStats1 = mQcStats2 = NULL;
    mSplitResults = NULL;
}

Processor::~Processor(){
    if(mReadGenerator) delete mReadGenerator;
    if(mSpliters){
        for(int i = 0; i < mOptions->thread; ++i){
            delete  mSpliters[i];
        }
        free(mSpliters); mSpliters = NULL;
    }
    if(mOptions->outspl){
        if(mWriter1){
            for(size_t i = 0; i < mOptions->samples.size(); ++i){
                if(mWriter1[i]) delete mWriter1[i];
            }
            free(mWriter1); mWriter1 = NULL;
        }
        if(mWriter2){
            for(size_t i = 0; i < mOptions->samples.size(); ++i){
                if(mWriter2[i]) delete mWriter2[i];
            }
            free(mWriter2); mWriter2 = NULL;
        }
    }
    if(mPeMergers){
        for(size_t i = 0; i < mOptions->samples.size()-1; ++i){
            if(mPeMergers[i]) delete mPeMergers[i];
        }
        free(mPeMergers); mPeMergers = NULL;
    }
    if(mDeRepers){
        for(size_t i = 0; i < mOptions->samples.size()-1; ++i){
            delete mDeRepers[i];
        }
    }
    if(mGEDetectors){
        for(size_t i = 0; i < mOptions->samples.size()-1; ++i){
            delete mGEDetectors[i];
        }
    }
    if(mPAligners){
        for(size_t i = 0; i < mOptions->samples.size()-1; ++i){
            delete mPAligners[i];
        }
    }
    for(auto& e: mTopnHtml){
        if(e){
            if(e->s){ free(e->s); e->s = NULL; }
            free(e); e = NULL;
        }
    }
    if(mQcStats1){ delete mQcStats1; mQcStats1 = NULL; }
    if(mQcStats2){ delete mQcStats2; mQcStats2 = NULL; }
    if(mSplitResults){ delete mSplitResults; mSplitResults = NULL; }
    if(mScAna){ delete mScAna; mScAna = NULL; }
}

void Processor::init(){
    if(mOptions->outspl){
        if(mOptions->fq2cal || mOptions->fq2spl){
            for(size_t i = 0; i < mOptions->samples.size(); ++i){
                std::string outFileName1 = mOptions->outdir + "/" + mOptions->samples[i] + ".R1.fq";
                std::string outFileName2 = mOptions->outdir + "/" + mOptions->samples[i] + ".R2.fq";
                if(mOptions->outgz){
                    outFileName1.append(".gz");
                    outFileName2.append(".gz");
                }
                mWriter1[i] = new Writer(outFileName1, mOptions->compression);
                if(mIsPE) mWriter2[i] = new Writer(outFileName2, mOptions->compression);
            }
        }else if(mOptions->fq2bam){
            std::string outFileName1 = mOptions->outdir + "/" + mOptions->samples[mOptions->dropidx] + ".R1.fq";
            std::string outFileName2 = mOptions->outdir + "/" + mOptions->samples[mOptions->dropidx] + ".R2.fq";
            mWriter1[mOptions->dropidx] = new Writer(outFileName1, mOptions->compression);
            if(mIsPE) mWriter2[mOptions->dropidx] = new Writer(outFileName2, mOptions->compression);
        }
    }
}

void Processor::process4fq2bam(){
    // start read generator
    util::loginfo("beg split and align reads");
    mReadGenerator->start();
    // process
    std::vector<std::pair<int32_t, int32_t>> vpidx;
    std::vector<std::future<void>> rets;
    ReadPack** rps = (ReadPack**)malloc(2 * sizeof(ReadPack*));
    rps[0] = rps[1] = NULL;
    while(mReadGenerator->getPack(rps)){
        int n = util::divideVecIdx(rps[0]->n, mOptions->thread, vpidx);
        // split
        rets.clear();
        for(int i = 0; i < n; ++i){
            rets.push_back(mOptions->tpl->enqueue(&Processor::splitRange, this, mSpliters[i], rps[0], rps[1], vpidx[i].first, vpidx[i].second));
        }
        for(auto& e: rets) e.get();
        rets.clear();
        // align
        for(size_t s = 0; s < mOptions->samples.size()-1; ++s){
            rets.push_back(mOptions->tpl->enqueue(&Processor::alignFq2BAM, this, rps[0], rps[1], s));
        }
        for(auto& e: rets) e.get();
        // write DROP
        if(mOptions->outspl){
            for(int i = 0; i < n; ++i){
                mWriter1[mOptions->dropidx]->write(mSpliters[i]->mResult->mSplitRead1[mOptions->dropidx]->s,
                                                   mSpliters[i]->mResult->mSplitRead1[mOptions->dropidx]->l);
                if(mIsPE){
                    mWriter2[mOptions->dropidx]->write(mSpliters[i]->mResult->mSplitRead2[mOptions->dropidx]->s,
                                                     mSpliters[i]->mResult->mSplitRead2[mOptions->dropidx]->l);
                }
            }
        }
        // reset seq/qual
        if(mOptions->skipfr1 || mOptions->skipfr2){
            for(int i = 0; i < rps[0]->n; ++i){
                if(rps[0]->reads[i]->dwhy != MATCH_FAIL_QC){
                    rps[0]->reads[i]->seq.s -= mOptions->skipfr1;
                    rps[0]->reads[i]->seq.l += mOptions->skipfr1;
                    rps[0]->reads[i]->qual.s -= mOptions->skipfr1;
                    rps[0]->reads[i]->qual.l += mOptions->skipfr1;
                    if(mIsPE){
                        rps[1]->reads[i]->seq.s -= mOptions->skipfr2;
                        rps[1]->reads[i]->seq.l += mOptions->skipfr2;
                        rps[1]->reads[i]->qual.s -= mOptions->skipfr2;
                        rps[1]->reads[i]->qual.l += mOptions->skipfr2;
                    }
                }
            }
        }
        // release pack
        rps[0]->r = false;
        rps[0]->w = true;
        if(mIsPE){
            rps[1]->r = false;
            rps[1]->w = true;
        }
    }
    util::loginfo(stderr, "end split and align reads", mReadGenerator->incnt);
    rets.clear();
    // reconstruct threadpool to maxify performance if needed
    if(mOptions->thread < (int)(mOptions->samples.size()) - 1){
        delete mOptions->tpl; 
        mOptions->tpl = new ThreadPool(mOptions->thread + 2);
    }
    // sort BAM
    rets.clear();
    util::loginfo("beg sort BAM");
    for(size_t s = 0; s < mOptions->samples.size()-1; ++s){
        rets.push_back(mOptions->tpl->enqueue(&Processor::sortBAM4Spl, this, s));
    }
    for(auto& e: rets) e.get();
    util::loginfo("end sort BAM");
    // writre BAM
    util::loginfo("beg write BAM");
    writeBAM4Spl();
    util::loginfo("end write BAM");
    // write ref
    util::loginfo("beg write REF fasta");
    writeRef4Spl();
    util::loginfo("end write REF fasta");
    util::loginfo("beg generate reports");
    // generat report
    std::vector<SplitResult*> srs;
    std::vector<QcStat*> qcs1;
    std::vector<QcStat*> qcs2;
    for(int i = 0; i < mOptions->thread; ++i){
        srs.push_back(mSpliters[i]->mResult);
        qcs1.push_back(mSpliters[i]->mQcR1);
    }
    if(mOptions->in2.size()){
        for(int i = 0; i < mOptions->thread; ++i) qcs2.push_back(mSpliters[i]->mQcR2);
    }
    // report summary/QC 
    mSplitResults = SplitResult::mergeResult(srs);
    mQcStats1 = QcStat::merge(qcs1);
    mQcStats2 = QcStat::merge(qcs2);
    Reporter* rptr = new Reporter(mOptions);
    rptr->report(mSplitResults, mQcStats1, mQcStats2, NULL);
    // report qc/stat of each split result
    reportTSV();
    reportFq2BamDiffTSV();
    reportJSON();
    util::loginfo("end generate reports");
}

void Processor::process4caledit1(){
    // start read generator
    util::loginfo("beg split, merge, derep and align reads");
    mReadGenerator->start();
    // process
    std::vector<std::pair<int32_t, int32_t>> vpidx;
    std::vector<std::future<void>> rets;
    ReadPack** rps = (ReadPack**)malloc(2 * sizeof(ReadPack*));
    rps[0] = rps[1] = NULL;
    while(mReadGenerator->getPack(rps)){
        int n = util::divideVecIdx(rps[0]->n, mOptions->thread, vpidx);
        // split
        rets.clear();
        for(int i = 0; i < n; ++i){
            rets.push_back(mOptions->tpl->enqueue(&Processor::splitRange, this, mSpliters[i], rps[0], rps[1], vpidx[i].first, vpidx[i].second));
        }
        for(auto& e: rets) e.get();
        rets.clear();
        // var call
        for(size_t s = 0; s < mOptions->samples.size()-1; ++s){
            if(mIsPE){
                if(mOptions->usesem) rets.push_back(mOptions->tpl->enqueue(&Processor::callVarMS, this, rps[0], rps[1], s));
                else rets.push_back(mOptions->tpl->enqueue(&Processor::callVarPE, this, rps[0], rps[1], s));
            }else rets.push_back(mOptions->tpl->enqueue(&Processor::callVarSE, this, rps[0], s));
        }
        for(auto& e: rets) e.get();
        // reset seq/qual
        if(mOptions->skipfr1 || mOptions->skipfr2){
            for(int i = 0; i < rps[0]->n; ++i){
                if(rps[0]->reads[i]->dwhy != MATCH_FAIL_QC){
                    rps[0]->reads[i]->seq.s -= mOptions->skipfr1;
                    rps[0]->reads[i]->seq.l += mOptions->skipfr1;
                    rps[0]->reads[i]->qual.s -= mOptions->skipfr1;
                    rps[0]->reads[i]->qual.l += mOptions->skipfr1;
                    if(mIsPE){
                        rps[1]->reads[i]->seq.s -= mOptions->skipfr2;
                        rps[1]->reads[i]->seq.l += mOptions->skipfr2;
                        rps[1]->reads[i]->qual.s -= mOptions->skipfr2;
                        rps[1]->reads[i]->qual.l += mOptions->skipfr2;
                    }
                }
            }
        }
        // release pack
        rps[0]->r = false;
        rps[0]->w = true;
        if(mIsPE){
            rps[1]->r = false;
            rps[1]->w = true;
        }
        // write fastq
        if(mOptions->outspl){
            for(int i = 0; i < n; ++i){
                rets.clear();
                for(size_t s = 0; s < mOptions->samples.size(); ++s){
                    rets.push_back(mOptions->tpl->enqueue(&Processor::writePack1, this, 
                                                          mSpliters[i]->mResult, s, mWriter1[s]));
                    if(mIsPE){
                        rets.push_back(mOptions->tpl->enqueue(&Processor::writePack2, this, 
                                                              mSpliters[i]->mResult, s, mWriter2[s]));
                    }
                }
                for(auto& e: rets) e.get();
            }
        }
    }
    util::loginfo(stderr, "end split, merge, derep and align reads", mReadGenerator->incnt);
    rets.clear();
    // reconstruct threadpool to maxify performance if needed
    if(mOptions->thread < (int)(mOptions->samples.size()) - 1){
        delete mOptions->tpl; 
        mOptions->tpl = new ThreadPool(mOptions->thread + 2);
    }
    util::loginfo("beg record duplicate count");
    // update merged count
    for(size_t s = 0; s < mOptions->samples.size()-1; ++s){
        rets.push_back(mOptions->tpl->enqueue(&Processor::bamRepCnt, this, s));
    }
    for(auto& e: rets) e.get();
    util::loginfo("end record duplicate count");
    // call edit events
    rets.clear();
    util::loginfo("beg compute edit effifiency");
    for(size_t s = 0; s < mOptions->samples.size()-1; ++s){
        rets.push_back(mOptions->tpl->enqueue(&Processor::calEdit, this, s));
    }
    for(auto& e: rets) e.get();
    util::loginfo("end compute edit effifiency");
    // plot ge
    rets.clear();
    util::loginfo("beg plot efficiency");
    for(size_t s = 0; s < mOptions->samples.size()-1; ++s){
        rets.push_back(mOptions->tpl->enqueue(&Processor::geplot, this, s));
    }
    for(auto& e: rets) e.get();
    util::loginfo("end plot efficiency");
    // sort BAM
    rets.clear();
    if(mOptions->f4bout & (GEOUT_FULL_BAM | GEOUT_SIMP_BAM)){
        util::loginfo("beg sort BAM");
        for(size_t s = 0; s < mOptions->samples.size()-1; ++s){
            rets.push_back(mOptions->tpl->enqueue(&Processor::sortBAM4Cal, this, s));
        }
        for(auto& e: rets) e.get();
        util::loginfo("end sort BAM");
        // writre BAM
        util::loginfo("beg write BAM");
        if(mOptions->f4bout & GEOUT_SIMP_BAM) writeVarBAM4Cal();
        if(mOptions->f4bout & GEOUT_FULL_BAM) writeFulBAM4Cal();
        util::loginfo("end write BAM");
    }
    // write BCF
    util::loginfo("beg write BCF");
    writeBCF();
    util::loginfo("end write BCF");
    // single cell analysis if needed
    if(mOptions->dosccal){
        util::loginfo("beg single cell analysis");
        mScAna->ana4all();
        util::loginfo("end single cell analysis");
    }
    // write ref
    util::loginfo("beg write REF fasta");
    writeRef4Cal();
    util::loginfo("end write REF fasta");
    util::loginfo("beg generate reports");
    // generat report
    std::vector<SplitResult*> srs;
    std::vector<QcStat*> qcs1;
    std::vector<QcStat*> qcs2;
    for(int i = 0; i < mOptions->thread; ++i){
        srs.push_back(mSpliters[i]->mResult);
        qcs1.push_back(mSpliters[i]->mQcR1);
    }
    if(mOptions->in2.size()){
        for(int i = 0; i < mOptions->thread; ++i) qcs2.push_back(mSpliters[i]->mQcR2);
    }
    // report summary/QC 
    mSplitResults = SplitResult::mergeResult(srs);
    mQcStats1 = QcStat::merge(qcs1);
    mQcStats2 = QcStat::merge(qcs2);
    Reporter* rptr = new Reporter(mOptions);
    rptr->report(mSplitResults, mQcStats1, mQcStats2, NULL);
    // report qc/stat of each split result
    reportTSV();
    reportJSON();
    // report html
    rets.clear();
    for(size_t i = 0; i < mOptions->samples.size() - 1; ++i){
        rets.push_back(mOptions->tpl->enqueue(&Processor::reportOneHTML, this, i));
    }
    for(auto& e: rets) e.get();
    reportAllHTML();
    // report sc report if needed
    if(mOptions->dosccal) reportSCGE();
    util::loginfo("end generate reports");
}

void Processor::process4split(){
    // start read generator
    util::loginfo("beg split fastq");
    mReadGenerator->start();
    // process
    std::vector<std::pair<int32_t, int32_t>> vpidx;
    std::vector<std::future<void>> rets;
    ReadPack** rps = (ReadPack**)malloc(2 * sizeof(ReadPack*));
    rps[0] = rps[1] = NULL;
    while(mReadGenerator->getPack(rps)){
        int n = util::divideVecIdx(rps[0]->n, mOptions->thread, vpidx);
        // split
        rets.clear();
        for(int i = 0; i < n; ++i){
            rets.push_back(mOptions->tpl->enqueue(&Processor::splitRange, this, mSpliters[i], rps[0], rps[1], vpidx[i].first, vpidx[i].second));
        }
        for(auto& e: rets) e.get();
        // reset seq/qual
        if(mOptions->skipfr1 || mOptions->skipfr2){
            for(int i = 0; i < rps[0]->n; ++i){
                if(rps[0]->reads[i]->dwhy != MATCH_FAIL_QC){
                    rps[0]->reads[i]->seq.s -= mOptions->skipfr1;
                    rps[0]->reads[i]->seq.l += mOptions->skipfr1;
                    rps[0]->reads[i]->qual.s -= mOptions->skipfr1;
                    rps[0]->reads[i]->qual.l += mOptions->skipfr1;
                    if(mIsPE){
                        rps[1]->reads[i]->seq.s -= mOptions->skipfr2;
                        rps[1]->reads[i]->seq.l += mOptions->skipfr2;
                        rps[1]->reads[i]->qual.s -= mOptions->skipfr2;
                        rps[1]->reads[i]->qual.l += mOptions->skipfr2;
                    }
                }
            }
        }
        // release pack
        rps[0]->r = false;
        rps[0]->w = true;
        if(mIsPE){
            rps[1]->r = false;
            rps[1]->w = true;
        }
        // write fastq
        if(mOptions->outspl){
            for(int i = 0; i < n; ++i){
                rets.clear();
                for(size_t s = 0; s < mOptions->samples.size(); ++s){
                    rets.push_back(mOptions->tpl->enqueue(&Processor::writePack1, this, 
                                                          mSpliters[i]->mResult, s, mWriter1[s]));
                    if(mIsPE){
                        rets.push_back(mOptions->tpl->enqueue(&Processor::writePack2, this, 
                                                              mSpliters[i]->mResult, s, mWriter2[s]));
                    }
                }
                for(auto& e: rets) e.get();
            }
        }
    }
    util::loginfo("end split fastq");
    util::loginfo("beg generate reports");
    // generat report
    std::vector<SplitResult*> srs;
    std::vector<QcStat*> qcs1;
    std::vector<QcStat*> qcs2;
    for(int i = 0; i < mOptions->thread; ++i){
        srs.push_back(mSpliters[i]->mResult);
        qcs1.push_back(mSpliters[i]->mQcR1);
    }
    if(mOptions->in2.size()){
        for(int i = 0; i < mOptions->thread; ++i) qcs2.push_back(mSpliters[i]->mQcR2);
    }
    // report summary/QC 
    mSplitResults = SplitResult::mergeResult(srs);
    mQcStats1 = QcStat::merge(qcs1);
    mQcStats2 = QcStat::merge(qcs2);
    Reporter* rptr = new Reporter(mOptions);
    rptr->report(mSplitResults, mQcStats1, mQcStats2, NULL);
    // report qc/stat of each split result
    reportTSV();
    reportJSON();
    // report html
    rets.clear();
    for(size_t i = 0; i < mOptions->samples.size() - 1; ++i){
        rets.push_back(mOptions->tpl->enqueue(&Processor::reportQCOneHTML, this, i));
    }
    for(auto& e: rets) e.get();
    reportQCAllHTML();
    // report sc report if needed
    if(mOptions->dosccal) reportSCGE();
    util::loginfo("end generate reports");
}

void Processor::process4caledit2(){
    // initialize writer
    mOptions->adjfr = false;
    mOptions->ostrand = true;
    mOptions->outspl = true;
    mOptions->outgz = true;
    std::vector<std::string> outPathsRead1;
    std::vector<std::string> outPathsRead2;
    mWriter1 = (Writer**)calloc(mOptions->samples.size(), sizeof(Writer*));
    for(size_t i = 0; i < mOptions->samples.size(); ++i) mWriter1[i] = NULL;
    mWriter2 = NULL;
    if(mIsPE){
        mWriter2 = (Writer**)calloc(mOptions->samples.size(), sizeof(Writer*));
        for(size_t i = 0; i < mOptions->samples.size(); ++i) mWriter2[i] = NULL;
    }
    mOutBams.reserve(mOptions->amplicons.size());
    mOutBcfs.reserve(mOptions->amplicons.size());
    mOutSCJSONs.reserve(mOptions->amplicons.size());
    for(size_t i = 0; i < mOptions->amplicons.size(); ++i){
        mOutBams.push_back(mOptions->tmpdir + "/" + mOptions->samples[i] + ".bam");
        mOutBcfs.push_back(mOptions->tmpdir + "/" + mOptions->samples[i] + ".bcf");
        mOutSCJSONs.push_back(mOptions->tmpdir + "/" + mOptions->samples[i] + ".sc.json");
    }
    for(size_t i = 0; i < mOptions->samples.size(); ++i){
        std::string outFileName1 = mOptions->tmpdir + "/" + mOptions->samples[i] + ".R1.fq";
        std::string outFileName2 = mOptions->tmpdir + "/" + mOptions->samples[i] + ".R2.fq";
        if(mOptions->outgz){
            outFileName1.append(".gz");
            outFileName2.append(".gz");
        }
        mWriter1[i] = new Writer(outFileName1, mOptions->compression);
        outPathsRead1.push_back(outFileName1);
        if(mIsPE && (!mOptions->usesem)){
            mWriter2[i] = new Writer(outFileName2, mOptions->compression);
            outPathsRead2.push_back(outFileName2);
        }
    }
    // initialize read generator
    util::loginfo("beg split libraries");
    mReadGenerator->start();
    // split
    std::vector<std::pair<int32_t, int32_t>> vpidx;
    std::vector<std::future<void>> rets;
    ReadPack** rps = (ReadPack**)malloc(2 * sizeof(ReadPack*));
    rps[0] = rps[1] = NULL;
    while(mReadGenerator->getPack(rps)){
        int n = util::divideVecIdx(rps[0]->n, mOptions->thread, vpidx);
        rets.clear();
        for(int i = 0; i < n; ++i){
            rets.push_back(mOptions->tpl->enqueue(&Processor::splitRange, this, 
                                                  mSpliters[i], rps[0], rps[1], vpidx[i].first, vpidx[i].second));
        }
        for(auto& e: rets) e.get();
        // reset seq/qual
        if(mOptions->skipfr1 || mOptions->skipfr2){
            for(int i = 0; i < rps[0]->n; ++i){
                if(rps[0]->reads[i]->dwhy != MATCH_FAIL_QC){
                    rps[0]->reads[i]->seq.s -= mOptions->skipfr1;
                    rps[0]->reads[i]->seq.l += mOptions->skipfr1;
                    rps[0]->reads[i]->qual.s -= mOptions->skipfr1;
                    rps[0]->reads[i]->qual.l += mOptions->skipfr1;
                    if(mIsPE){
                        rps[1]->reads[i]->seq.s -= mOptions->skipfr2;
                        rps[1]->reads[i]->seq.l += mOptions->skipfr2;
                        rps[1]->reads[i]->qual.s -= mOptions->skipfr2;
                        rps[1]->reads[i]->qual.l += mOptions->skipfr2;
                    }
                }
            }
        }
        // release pack
        rps[0]->r = false;
        rps[0]->w = true;
        if(mIsPE){
            rps[1]->r = false;
            rps[1]->w = true;
        }
        // write fastq
        for(int i = 0; i < n; ++i){
            rets.clear();
            for(size_t s = 0; s < mOptions->samples.size(); ++s){
                rets.push_back(mOptions->tpl->enqueue(&Processor::writePack1, this, 
                                                      mSpliters[i]->mResult, s, mWriter1[s]));
                if(mIsPE && (!mOptions->usesem)){
                    rets.push_back(mOptions->tpl->enqueue(&Processor::writePack2, this, 
                                                          mSpliters[i]->mResult, s, mWriter2[s]));
                }
            }
            for(auto& e: rets) e.get();
        }
    }
    for(size_t i = 0; i < mOptions->samples.size(); ++i){
        delete mWriter1[i]; mWriter1[i] = NULL;
        if(mIsPE && (!mOptions->usesem)){ delete mWriter2[i]; mWriter2[i] = NULL; }
    }
    free(mWriter1); mWriter1 = NULL;
    if(mIsPE){ free(mWriter2); mWriter2 = NULL; }
    delete mReadGenerator; mReadGenerator = NULL;
    util::loginfo("end split libraries");
    util::loginfo("beg merge, derep, align and call");
    // begin do work on each samples
    // initialize readers
    mSoloReaders = (ReadGenerator**)calloc(mOptions->samples.size()-1, sizeof(ReadGenerator*));
    for(size_t i = 0; i < mOptions->samples.size()-1; ++i){
        if(mIsPE && (!mOptions->usesem)) mSoloReaders[i] = new ReadGenerator(outPathsRead1[i].c_str(), outPathsRead2[i].c_str());
        else mSoloReaders[i] = new ReadGenerator(outPathsRead1[i].c_str());
    }
    int nworker = util::divideVecIdx(mOptions->samples.size()-1, mOptions->memone, vpidx);
    rets.clear();
    for(int i = 0; i < nworker; ++i){
        rets.push_back(mOptions->tpl->enqueue(&Processor::process4caledit3, this, vpidx[i].first, vpidx[i].second));
    }
    for(auto& e: rets) e.get();
    if(mOptions->f4bout & (GEOUT_FULL_BAM | GEOUT_SIMP_BAM)){
        if(mOptions->f4bout & GEOUT_SIMP_BAM) mergeVarBAM();
        if(mOptions->f4bout & GEOUT_FULL_BAM) mergeFulBAM();
    }
    mergeBCF();
    if(mOptions->dosccal) mergeSCJSON();
    util::loginfo("end merge, derep, align and call");
    // single cell analysis if needed
    if(mOptions->dosccal){
        util::loginfo("beg single cell analysis");
        mScAna->ana4all();
        util::loginfo("end single cell analysis");
    }
    // write ref
    util::loginfo("beg write REF fasta");
    writeRef4Cal();
    util::loginfo("end write REF fasta");
    util::loginfo("beg generate reports");
    // generat report
    std::vector<SplitResult*> srs;
    std::vector<QcStat*> qcs1;
    std::vector<QcStat*> qcs2;
    for(int i = 0; i < mOptions->thread; ++i){
        srs.push_back(mSpliters[i]->mResult);
        qcs1.push_back(mSpliters[i]->mQcR1);
    }
    if(mOptions->in2.size()){
        for(int i = 0; i < mOptions->thread; ++i) qcs2.push_back(mSpliters[i]->mQcR2);
    }
    // report summary/QC 
    mSplitResults = SplitResult::mergeResult(srs);
    mQcStats1 = QcStat::merge(qcs1);
    mQcStats2 = QcStat::merge(qcs2);
    Reporter* rptr = new Reporter(mOptions);
    rptr->report(mSplitResults, mQcStats1, mQcStats2, NULL);
    // report qc/stat of each split result
    reportTSV();
    reportJSON();
    // report html
    rets.clear();
    for(size_t i = 0; i < mOptions->samples.size() - 1; ++i){
        rets.push_back(mOptions->tpl->enqueue(&Processor::reportOneHTML, this, i));
    }
    for(auto& e: rets) e.get();
    reportAllHTML();
    util::loginfo("end generate reports");
    // cleanup
#ifndef KEEP_TMP_RESULT
    util::rmrf(mOptions->tmpdir.c_str());
#endif
}

void Processor::process4caledit3(int beg, int end){
    ReadPack** rps = (ReadPack**)malloc(2 * sizeof(ReadPack*));
    for(int i = beg ; i < end; ++i){
        mSoloReaders[i]->init(3, 10000);
        mSoloReaders[i]->start();
        rps[0] = rps[1] = NULL;
        while(mSoloReaders[i]->getPack(rps)){
            if(mIsPE && (!mOptions->usesem)){
                for(int j = 0; j < rps[0]->n; ++j){
                    rps[0]->reads[j]->sample = i;
                    rps[1]->reads[j]->sample = i;
                    rps[0]->reads[j]->off = rps[1]->reads[j]->off = 0;
                }
                if(mOptions->usesem) callVarMS(rps[0], rps[1], i);
                else callVarPE(rps[0], rps[1], i);
            }else{
                for(int j = 0; j < rps[0]->n; ++j){
                    rps[0]->reads[j]->sample = i; 
                    rps[0]->reads[j]->off = 0;
                    rps[0]->reads[j]->match = true;
                }
                callVarSE(rps[0], i);
            }
            rps[0]->r = false;
            rps[0]->w = true;
            if(mIsPE && (!mOptions->usesem)){
                rps[1]->r = false;
                rps[1]->w = true;
            }
        }
        delete mSoloReaders[i];  mSoloReaders[i] = NULL;
        bamRepCnt(i);
        calEdit(i);
        geplot(i);
        if(mOptions->f4bout & (GEOUT_FULL_BAM | GEOUT_SIMP_BAM)) writeOneBAM(i);
        writeOneBCF(i);
        if(mOptions->dosccal) writeOneSCJSON(i);
        for(auto& b: mGEDetectors[i]->vars){ bam_destroy1(b); b = NULL; }
        mGEDetectors[i]->vars.clear();
        for(auto& b: mGEDetectors[i]->bcfs){ bcf_destroy1(b); b = NULL; }
        mGEDetectors[i]->bcfs.clear();
        mLocker.lock();
        updateOneBCFHeader(i);
        mLocker.unlock();
        for(auto& b: mGEDetectors[i]->scretm){ if(b.second) { free(b.second); b.second = NULL; } } 
        mGEDetectors[i]->scretm.clear();
    }
    free(rps);
}

void Processor::splitRange(Spliter* s, ReadPack* p1, ReadPack* p2, int b, int e){
    s->init();
    for(int i = b; i < e; ++i){
        if(p1 && p2){
            s->splitRead(p1->reads[i], p2->reads[i]);
        }else{
            s->splitRead(p1->reads[i]);
        }
    }
}

void Processor::writePack1(SplitResult* r, size_t i, Writer* w){
    w->write(r->mSplitRead1[i]->s, r->mSplitRead1[i]->l);
}

void Processor::writePack2(SplitResult* r, size_t i, Writer* w){
    w->write(r->mSplitRead2[i]->s, r->mSplitRead2[i]->l);
}

void Processor::callVarPE(ReadPack* rp1, ReadPack* rp2, int i){
    kstring_t* pmses = (kstring_t*)calloc(1, sizeof(kstring_t));
    kstring_t* pmqss = (kstring_t*)calloc(1, sizeof(kstring_t));
    for(int r = 0; r < rp1->n; ++r){
        if(rp1->reads[r]->sample != i) continue;
        // merge
        if(rp1->reads[r]->strand) mPeMergers[i]->merge(rp1->reads[r], rp2->reads[r]);
        else mPeMergers[i]->merge(rp2->reads[r], rp1->reads[r]);
        if(mPeMergers[i]->ip->merged){// merge succeed
            // derep
            bool dup = mDeRepers[i]->derep_rec(mPeMergers[i]->ip->merged_sequence, 
                                               mPeMergers[i]->ip->merged_length,
                                               mPeMergers[i]->ip->merged_quality,
                                               false,
                                               rp1->reads[r]->comment.s,
                                               rp1->reads[r]->comment.l,
                                               true);
            if(dup) continue;
            kswr_t* ret = mGEDetectors[i]->align(mPeMergers[i]->ip->merged_sequence, mPeMergers[i]->ip->merged_length);
            set_merge_status(ret->smask, mPeMergers[i]->ip->reason);
            ret->smask |= KSW_FMERGED;
            bam1_t* aln = mGEDetectors[i]->ksw2bam(ret, mPeMergers[i]->ip->merged_sequence, mPeMergers[i]->ip->merged_length,
                                                   rp1->reads[r]->name.s, rp1->reads[r]->name.l, 0);
            bam_aux_update_int(aln, "VS", ret->smask);
            bam_aux_update_int(aln, "CC", mDeRepers[i]->clusters-1);
            mGEDetectors[i]->add_var(aln);
            kswr_destroy(ret);
        }else{// merge failed, align seperately
            // derep
            pmses->l = pmqss->l = 0;
            if(rp1->reads[r]->strand){
                kputsn(rp1->reads[r]->seq.s+rp1->reads[r]->off, rp1->reads[r]->seq.l-rp1->reads[r]->off, pmses);
                kputsn(rp2->reads[r]->seq.s+rp2->reads[r]->off, rp2->reads[r]->seq.l-rp2->reads[r]->off, pmses);
                kputsn(rp1->reads[r]->qual.s+rp1->reads[r]->off, rp1->reads[r]->qual.l-rp1->reads[r]->off, pmqss);
                kputsn(rp2->reads[r]->qual.s+rp2->reads[r]->off, rp2->reads[r]->qual.l-rp2->reads[r]->off, pmqss);
            }else{
                kputsn(rp2->reads[r]->seq.s+rp2->reads[r]->off, rp2->reads[r]->seq.l-rp2->reads[r]->off, pmses);
                kputsn(rp1->reads[r]->seq.s+rp1->reads[r]->off, rp1->reads[r]->seq.l-rp1->reads[r]->off, pmses);
                kputsn(rp2->reads[r]->qual.s+rp2->reads[r]->off, rp2->reads[r]->qual.l-rp2->reads[r]->off, pmqss);
                kputsn(rp1->reads[r]->qual.s+rp1->reads[r]->off, rp1->reads[r]->qual.l-rp1->reads[r]->off, pmqss);
            }
            bool dup = mDeRepers[i]->derep_rec(pmses->s, pmses->l, pmqss->s, false, rp1->reads[r]->comment.s, rp1->reads[r]->comment.l, true);
            if(dup) continue;
            kswr_t* swr1 = NULL;
            kswr_t* swr2 = NULL;
            if(rp1->reads[r]->strand){ // FR
                swr1 = mGEDetectors[i]->align(rp1->reads[r]->seq.s+rp1->reads[r]->off, rp1->reads[r]->seq.l-rp1->reads[r]->off);
                util::revComp2OriSeq(rp2->reads[r]->seq.s, rp2->reads[r]->seq.l);
                swr2 = mGEDetectors[i]->align(rp2->reads[r]->seq.s, rp2->reads[r]->seq.l-rp2->reads[r]->off);
            }else{ // RF
                swr2 = mGEDetectors[i]->align(rp2->reads[r]->seq.s+rp2->reads[r]->off, rp2->reads[r]->seq.l-rp2->reads[r]->off);
                util::revComp2OriSeq(rp1->reads[r]->seq.s, rp1->reads[r]->seq.l);
                swr1 = mGEDetectors[i]->align(rp1->reads[r]->seq.s, rp1->reads[r]->seq.l-rp1->reads[r]->off);
            }
            // smask-merge status bits
            set_merge_status(swr1->smask, mPeMergers[i]->ip->reason);
            set_merge_status(swr2->smask, mPeMergers[i]->ip->reason);
            // write bam records
            int mask1 = 0, mask2 = 0;
            mask1 |= (BAM_FREAD1 | BAM_FPAIRED);
            mask2 |= (BAM_FREAD2 | BAM_FPAIRED);
            if(!swr1->cigar){ mask1 |= BAM_FUNMAP; mask2 |= BAM_FMUNMAP; }
            if(!swr2->cigar){ mask2 |= BAM_FUNMAP; mask1 |= BAM_FMUNMAP; }
            if(rp1->reads[r]->strand){ mask1 |= BAM_FMREVERSE; mask2 |= BAM_FREVERSE; }
            else{ mask1 |= BAM_FREVERSE; mask2 |= BAM_FMREVERSE; }
            bam1_t* br1 = NULL; 
            bam1_t* br2 = NULL;
            if(rp1->reads[r]->strand){
                br1 = mGEDetectors[i]->ksw2bam(swr1, 
                                               rp1->reads[r]->seq.s+rp1->reads[r]->off, 
                                               rp1->reads[r]->seq.l-rp1->reads[r]->off,
                                               rp1->reads[r]->name.s, rp1->reads[r]->name.l, mask1);
                br2 = mGEDetectors[i]->ksw2bam(swr2, 
                                               rp2->reads[r]->seq.s,
                                               rp2->reads[r]->seq.l-rp2->reads[r]->off,
                                               rp2->reads[r]->name.s, rp2->reads[r]->name.l, mask2);
            }else{
                br1 = mGEDetectors[i]->ksw2bam(swr1, 
                                               rp1->reads[r]->seq.s, 
                                               rp1->reads[r]->seq.l-rp1->reads[r]->off,
                                               rp1->reads[r]->name.s, rp1->reads[r]->name.l, mask1);
                br2 = mGEDetectors[i]->ksw2bam(swr2, 
                                               rp2->reads[r]->seq.s+rp2->reads[r]->off,
                                               rp2->reads[r]->seq.l-rp2->reads[r]->off,
                                               rp2->reads[r]->name.s, rp2->reads[r]->name.l, mask2);
            }
            br1->core.mtid = br2->core.tid;
            br1->core.mpos = br2->core.pos;
            br2->core.mtid = br1->core.tid;
            br2->core.mpos = br1->core.pos;
            if(swr1->cigar && swr2->cigar){
                int isize = 1 + MAX(swr1->te, swr2->te) - MIN(swr1->tb, swr2->tb);
                if(br1->core.pos < br2->core.pos){
                    br1->core.isize = isize;
                    br2->core.isize = -isize;
                }else{
                    br1->core.isize = -isize;
                    br2->core.isize = isize;
                }
            }
            // mark compability of swr1 and swr2, overlap and same indel
            if((rp1->reads[r]->strand && kswge_indel_compatible(swr1, swr2, rp1->reads[r]->seq.s+rp1->reads[r]->off, rp2->reads[r]->seq.s) == 1) ||
               ((!rp1->reads[r]->strand) && kswge_indel_compatible(swr1, swr2, rp1->reads[r]->seq.s, rp2->reads[r]->seq.s+rp2->reads[r]->off) == 1)){
                swr1->smask |= KSW_FIDLOCMP;
                swr2->smask |= KSW_FIDLOCMP;
            }
            bam_aux_update_int(br1, "VS", swr1->smask);
            bam_aux_update_int(br2, "VS", swr2->smask);
            bam_aux_update_int(br1, "CC", mDeRepers[i]->clusters-1);
            bam_aux_update_int(br2, "CC", mDeRepers[i]->clusters-1);
            mGEDetectors[i]->add_var(br1);
            mGEDetectors[i]->add_var(br2);
            kswr_destroy(swr1);
            kswr_destroy(swr2);
        }
    }
    if(pmses->s) free(pmses->s);
    free(pmses); pmses = NULL;
    if(pmqss->s) free(pmqss->s);
    free(pmqss); pmqss = NULL;
}

void Processor::alignFq2BAM(ReadPack* rp1, ReadPack* rp2, int i){
    for(int r = 0; r < rp1->n; ++r){
        if(rp1->reads[r]->sample != i) continue;
        // align
        kswr_t* swr1 = NULL;
        kswr_t* swr2 = NULL;
        if(rp1->reads[r]->strand){ // FR
            swr1 = mPAligners[i]->align(rp1->reads[r]->seq.s+rp1->reads[r]->off, rp1->reads[r]->seq.l-rp1->reads[r]->off);
            util::revComp2OriSeq(rp2->reads[r]->seq.s, rp2->reads[r]->seq.l);
            swr2 = mPAligners[i]->align(rp2->reads[r]->seq.s, rp2->reads[r]->seq.l-rp2->reads[r]->off);
        }else{ // RF
            swr2 = mPAligners[i]->align(rp2->reads[r]->seq.s+rp2->reads[r]->off, rp2->reads[r]->seq.l-rp2->reads[r]->off);
            util::revComp2OriSeq(rp1->reads[r]->seq.s, rp1->reads[r]->seq.l);
            swr1 = mPAligners[i]->align(rp1->reads[r]->seq.s, rp1->reads[r]->seq.l-rp1->reads[r]->off);
        }
        // write bam records
        int mask1 = 0, mask2 = 0;
        mask1 |= (BAM_FREAD1 | BAM_FPAIRED);
        mask2 |= (BAM_FREAD2 | BAM_FPAIRED);
        if(swr1->cigar){ 
            ++mPAligners[i]->mapped;
        }else{
            mask1 |= BAM_FUNMAP; 
            mask2 |= BAM_FMUNMAP; 
            ++mPAligners[i]->unmapped;
        }
        if(swr2->cigar){
            ++mPAligners[i]->mapped;
        }else{
            mask2 |= BAM_FUNMAP; 
            mask1 |= BAM_FMUNMAP;
            ++mPAligners[i]->unmapped;
        }
        if(swr1->cigar && swr2->cigar){
            mPAligners[i]->stat_var(swr1, swr2);
        }else if(swr1->cigar){
            mPAligners[i]->stat_var(swr1, NULL);
        }else if(swr2->cigar){
            mPAligners[i]->stat_var(NULL, swr2);
        }
        if(rp1->reads[r]->strand){ mask1 |= BAM_FMREVERSE; mask2 |= BAM_FREVERSE; }
        else{ mask1 |= BAM_FREVERSE; mask2 |= BAM_FMREVERSE; }
        bam1_t* br1 = NULL; 
        bam1_t* br2 = NULL;
        if(rp1->reads[r]->strand){
            br1 = mPAligners[i]->ksw2bam(swr1, 
                                           rp1->reads[r]->seq.s+rp1->reads[r]->off, 
                                           rp1->reads[r]->seq.l-rp1->reads[r]->off,
                                           rp1->reads[r]->name.s, rp1->reads[r]->name.l, mask1);
            br2 = mPAligners[i]->ksw2bam(swr2, 
                                           rp2->reads[r]->seq.s,
                                           rp2->reads[r]->seq.l-rp2->reads[r]->off,
                                           rp2->reads[r]->name.s, rp2->reads[r]->name.l, mask2);
        }else{
            br1 = mPAligners[i]->ksw2bam(swr1, 
                                           rp1->reads[r]->seq.s, 
                                           rp1->reads[r]->seq.l-rp1->reads[r]->off,
                                           rp1->reads[r]->name.s, rp1->reads[r]->name.l, mask1);
            br2 = mPAligners[i]->ksw2bam(swr2, 
                                           rp2->reads[r]->seq.s+rp2->reads[r]->off,
                                           rp2->reads[r]->seq.l-rp2->reads[r]->off,
                                           rp2->reads[r]->name.s, rp2->reads[r]->name.l, mask2);
        }
        br1->core.mtid = br2->core.tid;
        br1->core.mpos = br2->core.pos;
        br2->core.mtid = br1->core.tid;
        br2->core.mpos = br1->core.pos;
        if(swr1->cigar && swr2->cigar){
            int isize = 1 + MAX(swr1->te, swr2->te) - MIN(swr1->tb, swr2->tb);
            if(br1->core.pos < br2->core.pos){
                br1->core.isize = isize;
                br2->core.isize = -isize;
            }else{
                br1->core.isize = -isize;
                br2->core.isize = isize;
            }
        }
        mPAligners[i]->add_var(br1);
        mPAligners[i]->add_var(br2);
        kswr_destroy(swr1);
        kswr_destroy(swr2);
    }
}

void Processor::callVarSE(ReadPack* rp1, int i){
    for(int r = 0; r < rp1->n; ++r){
        if(rp1->reads[r]->sample != i) continue;
        // derep
        bool dup = false;
        if(!rp1->reads[r]->strand){
            util::revComp2OriSeq(rp1->reads[r]->seq.s, rp1->reads[r]->seq.l);
            dup = mDeRepers[i]->derep_rec(rp1->reads[r]->seq.s, 
                                          rp1->reads[r]->seq.l-rp1->reads[r]->off, 
                                          rp1->reads[r]->qual.s,
                                          true,
                                          rp1->reads[r]->comment.s,
                                          rp1->reads[r]->comment.l,
                                          true);
        }else{
            dup = mDeRepers[i]->derep_rec(rp1->reads[r]->seq.s+rp1->reads[r]->off, 
                                          rp1->reads[r]->seq.l-rp1->reads[r]->off, 
                                          rp1->reads[r]->qual.s+rp1->reads[r]->off,
                                          false,
                                          rp1->reads[r]->comment.s,
                                          rp1->reads[r]->comment.l,
                                          true);
        }
        if(dup) continue;
        // align
        kswr_t* swf = NULL; bam1_t* aln = NULL;
        if(!rp1->reads[r]->strand){
            swf = mGEDetectors[i]->align(rp1->reads[r]->seq.s, rp1->reads[r]->seq.l-rp1->reads[r]->off);
            aln = mGEDetectors[i]->ksw2bam(swf, rp1->reads[r]->seq.s, rp1->reads[r]->seq.l-rp1->reads[r]->off,
                                                rp1->reads[r]->name.s, rp1->reads[r]->name.l, 0);
        }else{
            swf = mGEDetectors[i]->align(rp1->reads[r]->seq.s+rp1->reads[r]->off, rp1->reads[r]->seq.l-rp1->reads[r]->off);
            aln = mGEDetectors[i]->ksw2bam(swf, rp1->reads[r]->seq.s+rp1->reads[r]->off, rp1->reads[r]->seq.l-rp1->reads[r]->off,
                                                rp1->reads[r]->name.s, rp1->reads[r]->name.l, 0);
        }
        if(!rp1->reads[r]->strand) aln->core.flag |= BAM_FREVERSE;
        bam_aux_update_int(aln, "CC", mDeRepers[i]->clusters-1);
        bam_aux_update_int(aln, "VS", 0);
        mGEDetectors[i]->add_var(aln);
        kswr_destroy(swf);
    }
}

void Processor::callVarMS(ReadPack* rp1, ReadPack* rp2, int i){
    krec1_t* kr = NULL;
    for(int r = 0; r < rp1->n; ++r){
        if(rp1->reads[r]->sample != i) continue;
        // pick matched
        if(rp1->reads[r]->match) kr = rp1->reads[r];
        else kr = rp2->reads[r];
        // derep
        bool dup = false;
        if(!kr->strand){
            util::revComp2OriSeq(kr->seq.s, kr->seq.l);
            dup = mDeRepers[i]->derep_rec(kr->seq.s, 
                                          kr->seq.l-kr->off, 
                                          kr->qual.s,
                                          true,
                                          kr->comment.s,
                                          kr->comment.l,
                                          true);
        }else{
            dup = mDeRepers[i]->derep_rec(kr->seq.s+kr->off, 
                                          kr->seq.l-kr->off, 
                                          kr->qual.s+kr->off,
                                          false,
                                          kr->comment.s,
                                          kr->comment.l,
                                          true);
        }
        if(dup) continue;
        // align
        kswr_t* swf = NULL; bam1_t* aln = NULL;
        if(!kr->strand){
            swf = mGEDetectors[i]->align(kr->seq.s, kr->seq.l-kr->off);
            aln = mGEDetectors[i]->ksw2bam(swf, kr->seq.s, kr->seq.l-kr->off,
                                                kr->name.s, kr->name.l, 0);
        }else{
            swf = mGEDetectors[i]->align(kr->seq.s+kr->off, kr->seq.l-kr->off);
            aln = mGEDetectors[i]->ksw2bam(swf, kr->seq.s+kr->off, kr->seq.l-kr->off,
                                                kr->name.s, kr->name.l, 0);
        }
        if(!kr->strand) aln->core.flag |= BAM_FREVERSE;
        bam_aux_update_int(aln, "CC", mDeRepers[i]->clusters-1);
        bam_aux_update_int(aln, "VS", 0);
        mGEDetectors[i]->add_var(aln);
        kswr_destroy(swf);
    }
}

void Processor::sortBAM4Cal(size_t i){
    if(!mGEDetectors[i]->vsorted){
        std::sort(mGEDetectors[i]->vars.begin(), mGEDetectors[i]->vars.end(), BamSortByCoord());
        mGEDetectors[i]->vsorted = 1;
    }
}

void Processor::sortBAM4Spl(size_t i){
    if(!mPAligners[i]->vsorted){
        std::sort(mPAligners[i]->bams.begin(), mPAligners[i]->bams.end(), BamSortByCoord());
        mPAligners[i]->vsorted = 1;
    }
}

void Processor::writeOneBAM(int i){
    std::string out = mOutBams[i];
    samFile* ofp = sam_open(out.c_str(), "wb");
    assert(sam_hdr_write(ofp, mBamHeader) >= 0);
    for(auto& b: mGEDetectors[i]->vars){
        if(!(b->core.flag & BAM_FUNMAP)) b->core.tid = i;
        if(!(b->core.flag & BAM_FMUNMAP)) b->core.mtid = i;
        if(b->core.tid >= 0) assert(sam_write1(ofp, mBamHeader, b) >= 0);
    }
    for(auto& b: mGEDetectors[i]->vars){
        if(b->core.flag & BAM_FUNMAP){
            assert(sam_write1(ofp, mBamHeader, b) >= 0);
        }
    }
    sam_close(ofp);
}

void Processor::writeOneSCJSON(int i){
    kstring_t* ks = (kstring_t*)calloc(1, sizeof(kstring_t));
    // json begin
    const char* d1 = "    ";
    const char* d2 = "        ";
    const char* d3 = "            ";
    bool outone = false;
    GEDetector* gt = mGEDetectors[i];
    if(!gt->summarized) gt->summary();
    ksprintf(ks, "%s\"%s\": {\n", d1, gt->name);
    outone = false;
    sc_ge_t* psc = NULL;
    for(auto& p: gt->scretm){
        psc = p.second;
        if(!psc) continue;
        ksprintf(ks, "%s\"%s%s\": {\n", d2, SC_ID_PREFIX, p.first.c_str());
        ksprintf(ks, "%s\"DERO\": [%d,%d,%d,%d]\n",
                d3,
                psc->droped,
                psc->edicnt,
                psc->refcnt,
                psc->othcnt);
        ksprintf(ks, "%s},\n", d2);
        outone = true;
    }
    if(outone){
        ks->s[ks->l-2] = '\n';
        ks->s[ks->l-1] = '\0';
        --ks->l;
    }
    ksprintf(ks, "%s},\n", d1);
    // write2file
    FILE* ofp = fopen(mOutSCJSONs[i].c_str(), "w");
    fwrite(ks->s, sizeof(char), ks->l, ofp);
    fclose(ofp);
    free(ks->s);
    free(ks);
}

void Processor::mergeSCJSON(std::string outSCJsnPath){
    std::string out = mOutTotSCJSON;
    if(outSCJsnPath.size()) out = outSCJsnPath;
    // json
    kstring_t* ks = (kstring_t*)calloc(1, sizeof(kstring_t));
    // json begin
    ksprintf(ks, "{\n");
    std::string line;
    for(size_t i = 0; i < mOptions->samples.size()-1; ++i){
        std::ifstream fr(mOutSCJSONs[i]);
        while(std::getline(fr, line)){
            ksprintf(ks, "%s\n", line.c_str());
        }
        fr.close();
    }
    ks->s[ks->l-2] = '\n'; ks->s[ks->l-1] = '\0'; --ks->l;
    // json end
    ksprintf(ks, "}\n");
    // write2file
    FILE* ofp = fopen(out.c_str(), "w");
    fwrite(ks->s, sizeof(char), ks->l, ofp);
    fclose(ofp);
    free(ks->s);
    free(ks);
}

void Processor::writeVarBAM4Cal(std::string outBamPath){
    // write BAM, mapped
    std::string out = mOutTotBam;
    if(outBamPath.size())  out = outBamPath;
    samFile* ofp = sam_open(out.c_str(), "wb");
    assert(sam_hdr_write(ofp, mBamHeader) >= 0);
    for(size_t i = 0; i < mOptions->amplicons.size(); ++i){
        for(auto& b: mGEDetectors[i]->vars){
            if(!(b->core.flag & BAM_FUNMAP)){
                b->core.tid = i;
            }
            if(!(b->core.flag & BAM_FMUNMAP)){
                b->core.mtid = i;
            }
            if(b->core.tid >= 0){
                assert(sam_write1(ofp, mBamHeader, b) >= 0);
            }
        }
    }
    // write BAM, unmapped
    for(size_t i = 0; i < mOptions->amplicons.size(); ++i){
        for(auto& b: mGEDetectors[i]->vars){
            if(b->core.flag & BAM_FUNMAP){
                assert(sam_write1(ofp, mBamHeader, b) >= 0);
            }
        }
    }
    sam_close(ofp);
    assert(sam_index_build(out.c_str(), 0) == 0);
}

void Processor::writeFulBAM4Cal(std::string outBamPath){
    // write BAM, mapped
    std::string out = mOutFulBam;
    if(outBamPath.size())  out = outBamPath;
    samFile* ofp = sam_open(out.c_str(), "wb");
    assert(sam_hdr_write(ofp, mBamHeader) >= 0);
    uint8_t* cdata = NULL;
    int cval = 0;
    for(size_t i = 0; i < mOptions->amplicons.size(); ++i){
        for(auto& b: mGEDetectors[i]->vars){
            if(!(b->core.flag & BAM_FUNMAP)){
                b->core.tid = i;
            }
            if(!(b->core.flag & BAM_FMUNMAP)){
                b->core.mtid = i;
            }
            if(b->core.tid >= 0){
                // do CC time
                cdata = bam_aux_get(b, "CC");
                if(cdata){
                    cval = bam_aux2i(cdata);
                    bam_aux_update_int(b, "CC", 1);
                }else{
                    cval = 1;
                }
                for(int ci = 0; ci < cval; ++ci){
                    assert(sam_write1(ofp, mBamHeader, b) >= 0);
                }
                if(cdata){
                    bam_aux_update_int(b, "CC", cval);
                }
            }
        }
    }
    // write BAM, unmapped
    for(size_t i = 0; i < mOptions->amplicons.size(); ++i){
        for(auto& b: mGEDetectors[i]->vars){
            if(b->core.flag & BAM_FUNMAP){
                // do CC time
                cdata = bam_aux_get(b, "CC");
                if(cdata){
                    cval = bam_aux2i(cdata);
                    bam_aux_update_int(b, "CC", 1);
                }else{
                    cval = 1;
                }
                for(int ci = 0; ci < cval; ++ci){
                    assert(sam_write1(ofp, mBamHeader, b) >= 0);
                }
                if(cdata){
                    bam_aux_update_int(b, "CC", cval);
                }
            }
        }
    }
    sam_close(ofp);
    assert(sam_index_build(out.c_str(), 0) == 0);
}

void Processor::writeBAM4Spl(std::string outBamPath){
    // write BAM, mapped
    std::string out = mOutTotBam;
    if(outBamPath.size())  out = outBamPath;
    samFile* ofp = sam_open(out.c_str(), "wb");
    assert(sam_hdr_write(ofp, mBamHeader) >= 0);
    for(size_t i = 0; i < mOptions->amplicons.size(); ++i){
        for(auto& b: mPAligners[i]->bams){
            if(!(b->core.flag & BAM_FUNMAP)){
                b->core.tid = i;
            }
            if(!(b->core.flag & BAM_FMUNMAP)){
                b->core.mtid = i;
            }
            if(b->core.tid >= 0){
                assert(sam_write1(ofp, mBamHeader, b) >= 0);
            }
        }
    }
    // write BAM, unmapped
    for(size_t i = 0; i < mOptions->amplicons.size(); ++i){
        for(auto& b: mPAligners[i]->bams){
            if(b->core.flag & BAM_FUNMAP){
                assert(sam_write1(ofp, mBamHeader, b) >= 0);
            }
        }
    }
    sam_close(ofp);
    assert(sam_index_build(out.c_str(), 0) == 0);
}

void Processor::mergeVarBAM(std::string outBamPath){
    std::string umb = mOptions->tmpdir + "/_tmp_bam_.bam";
    std::string out = mOutTotBam;
    if(outBamPath.size()) out = outBamPath;
    samFile* ofp = sam_open(out.c_str(), "wb");
    assert(sam_hdr_write(ofp, mBamHeader) >= 0);
    samFile* ufp = sam_open(umb.c_str(), "wb");
    assert(sam_hdr_write(ufp, mBamHeader) >= 0);
    samFile* ifp = NULL;
    bam1_t* b = bam_init1();
    bam_hdr_t* h = NULL;
    for(size_t i = 0; i < mOptions->amplicons.size(); ++i){
        ifp = sam_open(mOutBams[i].c_str(), "r");
        h = sam_hdr_read(ifp);
        while(sam_read1(ifp, h, b) >= 0){
            if(b->core.tid < 0 || b->core.flag & BAM_FUNMAP){ // unmapped
                assert(sam_write1(ufp, mBamHeader, b) >= 0);
            }else{
                assert(sam_write1(ofp, mBamHeader, b) >= 0);
            }
        }
        sam_close(ifp);
    }
    sam_close(ufp);
    bam_hdr_destroy(h);
    ufp = sam_open(umb.c_str(), "r");
    h = sam_hdr_read(ufp);
    while(sam_read1(ufp, h, b) >= 0){
        assert(sam_write1(ofp, mBamHeader, b) >= 0);
    }
    sam_close(ufp);
    sam_close(ofp);
    assert(sam_index_build(out.c_str(), 0) == 0);
    bam_destroy1(b);
    bam_hdr_destroy(h);
}

void Processor::mergeFulBAM(std::string outBamPath){
    std::string umb = mOptions->tmpdir + "/_tmp_bam_.bam";
    std::string out = mOutFulBam;
    if(outBamPath.size()) out = outBamPath;
    samFile* ofp = sam_open(out.c_str(), "wb");
    assert(sam_hdr_write(ofp, mBamHeader) >= 0);
    samFile* ufp = sam_open(umb.c_str(), "wb");
    assert(sam_hdr_write(ufp, mBamHeader) >= 0);
    samFile* ifp = NULL;
    bam1_t* b = bam_init1();
    bam_hdr_t* h = NULL;
    uint8_t* cdata = NULL;
    int cval = 0;
    for(size_t i = 0; i < mOptions->amplicons.size(); ++i){
        ifp = sam_open(mOutBams[i].c_str(), "r");
        h = sam_hdr_read(ifp);
        while(sam_read1(ifp, h, b) >= 0){
            cdata = bam_aux_get(b, "CC");
            if(cdata){
                cval = bam_aux2i(cdata);
                bam_aux_update_int(b, "CC", 1);
            }else{
                cval = 1;
            }
            if(b->core.tid < 0 || b->core.flag & BAM_FUNMAP){ // unmapped
                for(int ci = 0; ci < cval; ++ci){
                    assert(sam_write1(ufp, mBamHeader, b) >= 0);
                }
            }else{
                for(int ci = 0; ci < cval; ++ci){
                    assert(sam_write1(ofp, mBamHeader, b) >= 0);
                }
            }
            if(cdata){
                    bam_aux_update_int(b, "CC", cval);
            }
        }
        sam_close(ifp);
    }
    sam_close(ufp);
    bam_hdr_destroy(h);
    ufp = sam_open(umb.c_str(), "r");
    h = sam_hdr_read(ufp);
    while(sam_read1(ufp, h, b) >= 0){
        assert(sam_write1(ofp, mBamHeader, b) >= 0);
    }
    sam_close(ufp);
    sam_close(ofp);
    assert(sam_index_build(out.c_str(), 0) == 0);
    bam_destroy1(b);
    bam_hdr_destroy(h);
}

void Processor::writeRef4Cal(){
    kstring_t ks = {0, 0, 0};
    for(size_t i = 0; i < mOptions->samples.size()-1; ++i){
        ksprintf(&ks, ">%s\n", mGEDetectors[i]->name);
        for(int c = 0; c < mGEDetectors[i]->sgrbeg; ++c) mGEDetectors[i]->ref[c] = tolower(mGEDetectors[i]->ref[c]);
        for(int c = mGEDetectors[i]->sgrend + 1; c < mGEDetectors[i]->rlen; ++c) mGEDetectors[i]->ref[c] = tolower(mGEDetectors[i]->ref[c]);
        kputsn(mGEDetectors[i]->ref, mGEDetectors[i]->rlen, &ks);
        kputc('\n', &ks);
    }
    FILE* fout = fopen(mOutRef.c_str(), "w");
    fwrite(ks.s, sizeof(char), ks.l, fout);
    fclose(fout);
    assert(fai_build(mOutRef.c_str()) == 0);
    free(ks.s);
}

void Processor::writeRef4Spl(){
    kstring_t ks = {0, 0, 0};
    for(size_t i = 0; i < mOptions->samples.size()-1; ++i){
        ksprintf(&ks, ">%s\n", mPAligners[i]->name);
        kputsn(mPAligners[i]->ref, mPAligners[i]->rlen, &ks);
        kputc('\n', &ks);
    }
    FILE* fout = fopen(mOutRef.c_str(), "w");
    fwrite(ks.s, sizeof(char), ks.l, fout);
    fclose(fout);
    assert(fai_build(mOutRef.c_str()) == 0);
    free(ks.s);
}

void Processor::bamRepCnt(int i){
    // get CC and qual
    int64_t* clscnt = (int64_t*)calloc(mDeRepers[i]->clusters+1, sizeof(int64_t));
    bucket_t** bs = (bucket_t**)calloc(mDeRepers[i]->clusters+1, sizeof(bucket_t*));
    for(uint64_t b  = 0; b < mDeRepers[i]->hashtablesize; ++b){
        if(mDeRepers[i]->hashtable[b].size){
            mDeRepers[i]->hashtable[b].mean_qual();
            clscnt[mDeRepers[i]->hashtable[b].id] = mDeRepers[i]->hashtable[b].size;
            bs[mDeRepers[i]->hashtable[b].id] = &(mDeRepers[i]->hashtable[b]);
        }
        if(mDeRepers[i]->hashtable[b].seq){
            free(mDeRepers[i]->hashtable[b].seq);
            mDeRepers[i]->hashtable[b].seq = NULL;
        }
        if(mDeRepers[i]->hashtable[b].ext){
            free(mDeRepers[i]->hashtable[b].ext);
            mDeRepers[i]->hashtable[b].ext = NULL;
        }
    }
    // update CC, qual, cell barcode
    uint8_t* data = NULL;
    uint8_t* qpts = NULL;
    kstring_t* cs = (kstring_t*)calloc(1, sizeof(kstring_t));
    int64_t id = 0;
    int scn = 0;
    for(auto& b: mGEDetectors[i]->vars){
        data = bam_aux_get(b, "CC");
        if(data){
            id = bam_aux2i(data);
            if(bs[id]){// has qual
                qpts = bam_get_qual(b);
                if(bs[id]->len == b->core.l_qseq || (b->core.flag & BAM_FMREVERSE)){ // whole or first part
                    for(int j = 0; j < b->core.l_qseq; ++j) qpts[j] = MAX(bs[id]->qual[j] - 33, 0);
                }else{// second part
                    for(int j = 0; j < b->core.l_qseq; ++j) qpts[j] = MAX(bs[id]->qual[bs[id]->len-b->core.l_qseq+j] - 33, 0);
                }
            }
            bam_aux_update_int(b, "CC", clscnt[id]);
            if(mOptions->dosccal){
                scn = 0;
                cs->l = 0;
                khiter_t k;
                for(k = kh_begin(bs[id]->emap); k != kh_end(bs[id]->emap); ++k){
                    if(kh_exist(bs[id]->emap, k)){
                        ksprintf(cs, "%s-%d|", (char*)kh_key(bs[id]->emap, k), kh_val(bs[id]->emap, k));
                        ++scn;
                    }
                }
                if(cs->l && cs->s){
                    bam_aux_update_str(b, CELL_BARCODE_ID_TAG, cs->l-1, cs->s);
                    bam_aux_update_int(b, CELL_BARCODE_COUNT_TAG, scn);
                }
            }
        }
    }
    if(cs->s) free(cs->s);
    free(cs);
    // free qual
    for(uint64_t b  = 0; b < mDeRepers[i]->hashtablesize; ++b){
        if(mDeRepers[i]->hashtable[b].qual){
            free(mDeRepers[i]->hashtable[b].qual);
            mDeRepers[i]->hashtable[b].qual = NULL;
        }
    }
    // free emap
    for(uint64_t b  = 0; b < mDeRepers[i]->hashtablesize; ++b){
        if(mDeRepers[i]->hashtable[b].emap){
            khiter_t k;
            for(k = kh_begin(mDeRepers[i]->hashtable[b].emap); k != kh_end(mDeRepers[i]->hashtable[b].emap); ++k){
                if(kh_exist(mDeRepers[i]->hashtable[b].emap, k)){
                    free((char*)kh_key(mDeRepers[i]->hashtable[b].emap, k));
                }
            }
            kh_destroy(es, mDeRepers[i]->hashtable[b].emap);
            mDeRepers[i]->hashtable[b].emap = NULL;
        }
    }
    // release res
    free(clscnt);
    free(bs);
}

void Processor::geplot(int i){
    gep_opt_t* geo = new gep_opt_t();
    geo->amplicon = mGEDetectors[i]->ref;
    geo->name = mGEDetectors[i]->name;
    geo->outdir = mOptions->pltdir.c_str();
    geo->rlen = mGEDetectors[i]->rlen;
    geo->topn = mOptions->edo.topn;
    geo->flklen = mOptions->edo.flklen;
    geo->sgrbeg = mGEDetectors[i]->sgrbeg;
    geo->sgrend = mGEDetectors[i]->sgrend;
    geo->refbeg = MAX(0, mGEDetectors[i]->sgrbeg - geo->flklen);
    geo->refend = MIN(geo->rlen-1, geo->sgrend + geo->flklen);
    geo->refgot = geo->refend - geo->refbeg + 1;
    geo->vbeg = geo->sgrbeg + mOptions->edo.cutbuflen;
    geo->vend = geo->sgrend + mOptions->edo.cutbuflen;
    geo->ismt = mOptions->edo.ismt;
    geo->bamhdr = sam_hdr_dup(mBamHeader);
    geo->genhtvar = true;
    geo->hfigh = mOptions->hmo.tnfigh;
    if(mOptions->tpdetail) geo->update_iof();
    geo->update_amp();
    geo->plot(mGEDetectors[i]->vars);
    geo->topn2html(mTopnHtml[i]);
    if(geo->empty) util::loginfo(stderr, mOptions->logmtx, "amplicon %s has no valid variants to plot", geo->name);
    delete geo;
}

void Processor::calEdit(int i){
    util::loginfo(stderr, mOptions->logmtx, "beg compute edit efficience of %s", mOptions->samples[i].c_str());
    mGEDetectors[i]->cal_clspos();
    mGEDetectors[i]->mark_var();
    mGEDetectors[i]->mark_misprime();
    if(mOptions->edo.cluster) mGEDetectors[i]->kcluster();
    if(mOptions->ctrl.size() && mGEDetectors[i]->vars.size()){
        mGEDetectors[i]->cbam2vars();
        mGEDetectors[i]->markctrlref();
    }
    if(mIsPE && (!mOptions->usesem)) mGEDetectors[i]->fixpe();
    mGEDetectors[i]->cal_edieff();
    if(mOptions->edo.hapcnt) mGEDetectors[i]->hapcnt();
    if(mOptions->dosccal) mGEDetectors[i]->cal_scedff();
    if(mOptions->edo.minaf > 0) mGEDetectors[i]->af_filter();
    mGEDetectors[i]->bam2bcf();
    util::loginfo(stderr, mOptions->logmtx, "end compute edit efficience of %s", mOptions->samples[i].c_str());
}

void Processor::reportJSON(){
    kstring_t* ks = (kstring_t*)calloc(1, sizeof(kstring_t));
    const char* d1 = "    ";
    const char* d2 = "        ";
    QcStat* qs1 = NULL;
    QcStat* qs2 = NULL;
    // json beg
    ksprintf(ks, "{\n");
    for(size_t i = 0; i < mOptions->samples.size() - 1; ++i){
        // beg sample
        ksprintf(ks, "\"%s\": {\n", mOptions->samples[i].c_str());
        qs1 = mSplitResults->mQCRead1[i];
        qs2 = mSplitResults->mQCRead2[i];
        // beg BasicQC
        ksprintf(ks, "%s\"BasicQC\": {\n", d1);
        qs1->reportJSON(ks, d2, d1, 1);
        if(qs2){
            ksprintf(ks, ",\n");
            qs2->reportJSON(ks, d2, d1, 2);
        }
        ksprintf(ks, "\n%s},\n", d1);
        // end BasicQC
        if(mPeMergers){ // MergeStat
            ksprintf(ks, "%s\"MergeResult\": {\n", d1);
            mPeMergers[i]->reportJSON(ks, d2, d1);
            ksprintf(ks, "\n%s},\n", d1);
        }
        if(mDeRepers){ // DedupStat
            ksprintf(ks, "%s\"DerepResult\": {\n", d1);
            mDeRepers[i]->reportJSON(ks, d2, d1);
            ksprintf(ks, "\n%s},\n", d1);
        }
        if(mGEDetectors){ // GeneEdit
            ksprintf(ks, "%s\"GeneEdit\": {\n", d1);
            mGEDetectors[i]->reportJSON(ks, d2, d1);
            ksprintf(ks, "\n%s},\n", d1);
        }
        if(mPAligners){ // Paligner
            ksprintf(ks, "%s\"PAlign\": {\n", d1);
            mPAligners[i]->reportJSON(ks, d2, d1);
            ksprintf(ks, "\n%s},\n", d1);
        }
        ks->l -= 1;
        ks->s[ks->l-1] = '\n';
        // end sample
        ksprintf(ks, "},\n");
    }
    ks->l -= 2;
    // json end
    ksprintf(ks, "\n}\n");
    // write to file
    FILE* fp = fopen(mOptions->jss.c_str(), "w");
    fwrite(ks->s, sizeof(char), ks->l, fp);
    fclose(fp);
    free(ks->s);
    free(ks);
}

void Processor::reportFq2BamDiffTSV(){
    kstring_t* ks = (kstring_t*)calloc(1, sizeof(kstring_t));
    // se
    int mxl = 0;
    for(size_t i = 0; i < mOptions->samples.size() - 1; ++i){
        mPAligners[i]->summary();
        if(mxl < mPAligners[i]->maxdstl_se) mxl = mPAligners[i]->maxdstl_se;
    }
    PAligner::tsv4DiffxHead(ks, mxl);
    for(size_t i = 0; i < mOptions->samples.size() - 1; ++i){
        mPAligners[i]->tsv4DiffxBodySE(ks);
    }
    // write to file
    FILE* fp = fopen(mOptions->txx_se.c_str(), "w");
    fwrite(ks->s, sizeof(char), ks->l, fp);
    fclose(fp);
    // pe
    ks->l = 0;
    mxl = 0;
    for(size_t i = 0; i < mOptions->samples.size() - 1; ++i){
        mPAligners[i]->summary();
        if(mxl < mPAligners[i]->maxdstl_pe) mxl = mPAligners[i]->maxdstl_pe;
    }
    PAligner::tsv4DiffxHead(ks, mxl);
    for(size_t i = 0; i < mOptions->samples.size() - 1; ++i){
        mPAligners[i]->tsv4DiffxBodyPE(ks);
    }
    // write to file
    fp = fopen(mOptions->txx_pe.c_str(), "w");
    fwrite(ks->s, sizeof(char), ks->l, fp);
    fclose(fp);
    // free
    free(ks->s);
    free(ks);
}


void Processor::reportTSV(){
    kstring_t* ks = (kstring_t*)calloc(1, sizeof(kstring_t));
    QcStat* qs1 = NULL;
    QcStat* qs2 = NULL;
    // head
    kputs("Amplicon\t", ks);
    QcStat::tsvHead(ks); kputc('\t', ks);
    if(mPeMergers){ PairMerger::tsvHead(ks); kputc('\t', ks); }
    if(mDeRepers){ DeReper::tsvHead(ks); kputc('\t', ks); }
    if(mGEDetectors){ GEDetector::tsvHead(ks); kputc('\t', ks); }
    if(mPAligners){ PAligner::tsvHead(ks); kputc('\t', ks); }
    ks->s[ks->l-1] = '\n';
    // body
    for(size_t i = 0; i < mOptions->samples.size() - 1; ++i){
        ksprintf(ks, "%s\t", mOptions->samples[i].c_str());
        qs1 = mSplitResults->mQCRead1[i];
        qs2 = mSplitResults->mQCRead2[i];
        qs1->tsvBody(ks, qs2); kputc('\t', ks);
        if(mPeMergers){ mPeMergers[i]->tsvBody(ks); kputc('\t', ks); }
        if(mDeRepers){ mDeRepers[i]->tsvBody(ks); kputc('\t', ks); }
        if(mGEDetectors){ mGEDetectors[i]->tsvBody(ks); kputc('\t', ks); }
        if(mPAligners){ mPAligners[i]->tsvBody(ks); kputc('\t', ks); }
        ks->s[ks->l-1] = '\n';
    }
    // write to file
    FILE* fp = fopen(mOptions->tss.c_str(), "w");
    fwrite(ks->s, sizeof(char), ks->l, fp);
    fclose(fp);
    free(ks->s);
    free(ks);
}

void Processor::reportSCGE(){
    // json
    kstring_t* ks = (kstring_t*)calloc(1, sizeof(kstring_t));
    // json begin
    const char* d1 = "    ";
    const char* d2 = "        ";
    const char* d3 = "            ";
    ksprintf(ks, "{\n");
    GEDetector* gt = NULL;
    bool outone = false;
    sc_ge_t* psc = NULL;
    for(size_t i = 0; i < mOptions->samples.size()-1; ++i){
        gt = mGEDetectors[i];
        if(!gt->summarized) gt->summary();
        ksprintf(ks, "%s\"%s\": {\n", d1, gt->name);
        outone = false;
        for(auto& p: gt->scretm){
            psc = p.second;
            ksprintf(ks, "%s\"%s%s\": {\n", d2, SC_ID_PREFIX, p.first.c_str());
            ksprintf(ks, "%s\"DERO\": [%d,%d,%d,%d]\n",
                    d3, 
                    psc->droped, 
                    psc->edicnt, 
                    psc->refcnt, 
                    psc->othcnt);
            ksprintf(ks, "%s},\n", d2);
            outone = true;
        }
        if(outone){
            ks->s[ks->l-2] = '\n';
            ks->s[ks->l-1] = '\0';
            --ks->l;
        }
        ksprintf(ks, "%s},\n", d1);
    }
    ks->s[ks->l-2] = '\n'; ks->s[ks->l-1] = '\0'; --ks->l;
    ksprintf(ks, "}\n");
    // json end
    FILE* ofp = fopen(mOutTotSCJSON.c_str(), "w");
    fwrite(ks->s, sizeof(char), ks->l, ofp);
    fclose(ofp);
    free(ks->s);
    free(ks);
}

