#ifndef SPLITER_H
#define SPLITER_H

#include <string>
#include <fstream>
#include <vector>
#include "util.h"
#include "krec.h"
#include "qcstat.h"
#include "options.h"
#include "splitresult.h"

// class to process splitting of a read
class Spliter{
    public:
        Options* mOpt; // pointer to options
        pre_idx_t* mAlnSpliter; // preidx alignment based spliter
        amp_idx_t* mAmpSpliter; // amplicon alignment based spliter
        SplitResult* mResult; // pointer to split result
        QcStat* mQcR1; // read1 QC
        QcStat* mQcR2; // read2 QC

    public:
        // constructor
        Spliter(Options* opt);

        // destructor
        ~Spliter();

        // init
        void init();

        
        /** exact match a prefix against a read
         * @param r pointer to read struct
         * @param p pattern to be matched
         */
        bool exactMatch(krec1_t* r, const std::string& p);
        
        /** exact check read/prefix match
         * @param r pointer to read struct
         * @return prefix matched
         */
        prefix_t* exactMatch(krec1_t* r);

        /** exact check read/prefix match
         * @param r1 pointer to read1
         * @param r2 pointer to read2
         */
        prefix_t* exactMatch(krec1_t* r1, krec1_t* r2);

        /** split a read and get its corresponding sample file
         * @param r pointer to a Read struct
         * @param m mismatch of r against matched prefix
         */
        void splitRead(krec1_t* r);

        /** split a read pair and get their corresponding sample file
         * @param r1 pointer to a Read object
         * @param r2 pointer to a Read object
         */
        void splitRead(krec1_t* r1, krec1_t* r2);
};
        
#endif
