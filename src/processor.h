#ifndef READ_PROCESSOR_H_
#define READ_PROCESSOR_H_

#include "ksw4get.h"
#include "derep.h"
#include "scana.h"
#include "bamplp.h"
#include "geplot.h"
#include "writer.h"
#include "report.h"
#include "options.h"
#include "readgen.h"
#include "spliter.h"
#include "gedetect.h"
#include "paligner.h"
#include "bamsorter.h"
#include "mergepairs.h"
#include "htslib/vcf.h"
#include "htslib/faidx.h"

// class to do split work
class Processor{
    public:
        Options* mOptions; // options
        ReadGenerator* mReadGenerator; // read generator
        ReadGenerator** mSoloReaders; // read generators
        Spliter** mSpliters; // spliters to do parallel work
        Writer** mWriter1; // read1 writers
        Writer** mWriter2; // read2 writers
        PairMerger** mPeMergers; // pair end mergers
        DeReper** mDeRepers; // derep workers
        GEDetector** mGEDetectors; // gene edit detector
        PAligner** mPAligners; // pairwise sequence aligner
        SplitResult* mSplitResults; // merged split results
        QcStat* mQcStats1; // merged QC stats r1
        QcStat* mQcStats2; // merged QC stats r2
        bam_hdr_t* mBamHeader; // bam header
        bcf_hdr_t* mBcfHeader; // bcf header
        bool mIsPE; // pair end processor
        scana_opt_t* mScAna; // single cell analysis
        std::mutex mLocker; // locker for one by one process
        std::vector<std::string> mOutBams; // outbams
        std::vector<std::string> mOutBcfs; // outbcfs
        std::vector<std::string> mOutSCJSONs; // outscjsns
        std::vector<kstring_t*> mTopnHtml; // topn html str
        std::vector<std::string> mOutHTMLs; // outhtml paths
        std::string mOutTotBam, mOutFulBam, mOutTotBcf, mOutTotSCJSON, mOutRef;

       // constructor 
       Processor(Options* opt);

       // destructor
       ~Processor();

       // init
       void init();

       void process4caledit1(); // caledit, big mem
       void process4caledit2(); // caledit, low mem
       void process4caledit3(int beg, int end); // process batch of samples in range
       void process4split(); // split fastq
       void process4fq2bam(); // split fastq to BAM

       // split a range
       void splitRange(Spliter* s, ReadPack* p1, ReadPack* p2, int b, int e);

       // call variant of one sample
       void callVarSE(ReadPack* rp1, int i);
       void callVarPE(ReadPack* rp1, ReadPack* rp2, int i);
       void callVarMS(ReadPack* rp1, ReadPack* rp2, int i);

       // align one sample
       void alignFq2BAM(ReadPack* rp1, ReadPack* rp2, int i);

       // write a range of fastq
       void writePack1(SplitResult* r, size_t i, Writer* w);
       void writePack2(SplitResult* r, size_t i, Writer* w);

       // sort a sample BAM
       void sortBAM4Cal(size_t i);
       void sortBAM4Spl(size_t i);

       // write a whole BAM
       void writeVarBAM4Cal(std::string outBamPath = "");
       void writeFulBAM4Cal(std::string outBamPath = "");
       void writeOneBAM(int i);
       void mergeVarBAM(std::string outBamPath = "");
       void mergeFulBAM(std::string outBamPath = "");

       // write fq2bam BAM
       void writeBAM4Spl(std::string outBamPath = "");

       // write a whole BCF
       void updateOneBCFHeader(int i);
       void writeBCF(std::string outBcfPath = "");
       void writeOneBCF(int i);
       void mergeBCF(std::string outBcfPath = "");

       // write a whole Single Cell
       void writeSCJSON(std::string outSCJsnPath = "");
       void writeOneSCJSON(int i);
       void mergeSCJSON(std::string outSCJsnPath = "");

       // output reference fasta with index
       void writeRef4Cal();
       void writeRef4Spl();

       // update derep bam count
       void bamRepCnt(int i);

       // identify relative real edit events
       void calEdit(int i);

       // graph
       void geplot(int i);

       // report
       void reportSCGE();
       void reportJSON();
       void reportTSV();
       void reportFq2BamDiffTSV();
       void reportAllHTML();
       void reportQCAllHTML();
       void reportOneHTML(int i);
       void reportQCOneHTML(int i);
       void printHeader(kstring_t* s);
       void printFooter(kstring_t* s);
};
#endif
