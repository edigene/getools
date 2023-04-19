#ifndef GE_FQ2BAM_H
#define GE_FQ2BAM_H

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <libgen.h>
#include "common.h"
#include "options.h"
#include "processor.h"

void gef2b_usage(Options* opt, char* arg0){
    fprintf(stderr, "\nUsage: %s %s [options]\n", PACKAGE_NAME, arg0);
    fprintf(stderr, "\nRead split match options:\n");
    fprintf(stderr, "         -f INT   max offset allowed before match begin [%d]\n", opt->maxoff);
    fprintf(stderr, "         -m INT   max mismatch allowed in each match [%d]\n", opt->maxmm);
    fprintf(stderr, "         -1 INT   read1 skip length before split match [%ld]\n", opt->skipfr1);
    fprintf(stderr, "         -2 INT   read2 skip length before split match [%ld]\n", opt->skipfr2);
    fprintf(stderr, "         -3 INT   max barcode length to skip during primer check [%d]\n", opt->maxbo);
    fprintf(stderr, "         -S       only one end match is okay\n");
    fprintf(stderr, "         -g       use amplicon alignment instead of primer\n");
    fprintf(stderr, "\nQuality filter options:\n");
    fprintf(stderr, "         -z INT   max N bases allowed in non-primer region [%d]\n", opt->lowq.maxN);
    fprintf(stderr, "         -q INT   min base quality as poor quality base [%d]\n", opt->lowq.minBaseQual);
    fprintf(stderr, "         -Q FLOAT max rate of poor quality bases allowed in one read [%f]\n", opt->lowq.maxLowQualFrac);
    fprintf(stderr, "         -P INT   max number of poor quality bases allowed in one read [%d]\n", opt->lowq.maxLowQualBase);
    fprintf(stderr, "\nAlignment score options:\n");
    fprintf(stderr, "         -M INT   match score[%d]\n", opt->aln.match);
    fprintf(stderr, "         -X INT   mismatch penalty [%d]\n", opt->aln.mismatch);
    fprintf(stderr, "         -O INT   gap open penalty [%d]\n", opt->aln.gapopen);
    fprintf(stderr, "         -E INT   gap extend penalty [%d]\n", opt->aln.gapext);
    fprintf(stderr, "\nInput output options:\n");
    fprintf(stderr, "         -i FILE  read1 input file\n");
    fprintf(stderr, "         -I FILE  read2 input file\n");
    fprintf(stderr, "         -c FILE  configure file(4 column TSV:[ampliconName,fwdPrimer,revPrimer,ampliconSeq])\n");
    fprintf(stderr, "         -o DIR   output directory [%s]\n", opt->outdir.c_str());
    fprintf(stderr, "         -t INT   threads used to do computation [%d]\n", opt->thread);
    fprintf(stderr, "         -p INT   max packs in buffer [%d]\n", opt->maxpack);
    fprintf(stderr, "         -r INT   max reads in each buffer pack [%d]\n", opt->maxreads);
    fprintf(stderr, "         -l INT   compression level if output gz file [%d]\n", opt->compression);
    fprintf(stderr, "         -d STR   drop library name [%s]\n", opt->droplib.c_str());
    fprintf(stderr, "         -L       output dropped library\n");
    fprintf(stderr, "         -Z       output dropped library in gzip format\n");
    fprintf(stderr, "\n");
    delete opt;
}

int gef2b_main(int argc, char** argv){
    // parse args
    Options* opt = new Options();
    if(argc == 1){
        gef2b_usage(opt, argv[0]);
        return 0;
    }
    // parse args
    const char* args = "f:m:z:q:Q:P:M:X:O:E:i:I:c:o:t:p:r:1:2:3:l:d:SgLZh";
    int c = -1;
    while((c = getopt(argc, argv, args)) >= 0){
        switch(c){
            // Primer match options
            case 'f': opt->maxoff = atoi(optarg); break;
            case 'm': opt->maxmm = atoi(optarg); break;
            case 'S': opt->sematch = true; break;
            case 'g': opt->s8aln = true; break;
            case '1': opt->skipfr1 = atoi(optarg); break;
            case '2': opt->skipfr2 = atoi(optarg); break;
            case '3': opt->maxbo = atoi(optarg); break;
            // Quality filter options
            case 'z': opt->lowq.maxN = atoi(optarg); break;
            case 'q': opt->lowq.minBaseQual = atoi(optarg); break;
            case 'Q': opt->lowq.maxLowQualFrac = atoi(optarg); break;
            case 'P': opt->lowq.maxLowQualBase = atoi(optarg); break;
            // Alignment score options
            case 'M': opt->aln.match = atoi(optarg); break;
            case 'X': opt->aln.mismatch = atoi(optarg); break;
            case 'O': opt->aln.gapopen = atoi(optarg); break;
            case 'E': opt->aln.gapext = atoi(optarg); break;
            // Input output options
            case 'i': opt->in1 = strdup(optarg); break;
            case 'I': opt->in2 = strdup(optarg); break;
            case 'c': opt->prf = strdup(optarg); break;
            case 'o': opt->outdir = strdup(optarg); break;
            case 't': opt->thread = atoi(optarg); break;
            case 'p': opt->maxpack = atoi(optarg); break;
            case 'r': opt->maxreads = atoi(optarg); break;
            case 'l': opt->compression = atoi(optarg); break;
            case 'd': opt->droplib = strdup(optarg); break;
            case 'Z': opt->outgz = true; break;
            case 'L': opt->outspl = true; break;
            case 'h': gef2b_usage(opt, argv[0]); return 0; break;
            default: break;
        }
    }
    opt->fq2bam = true;
    // check args and run
    if(opt->valid()){
        util::loginfo("beg parse args");
        if(opt->memone > 0) opt->memlow = true;
        opt->init(argc, argv);
        util::loginfo("end parse args");
        util::loginfo("beg split fastq to BAM");
        Processor *p = new Processor(opt);
        p->init();
        p->process4fq2bam();
        util::loginfo("end split fastq to BAM");
        delete p;
        delete opt;
        return 0;
    }else{
        delete opt;
        return 1;
    }
}

#endif
