#ifndef GE_SPL_H
#define GE_SPL_H

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <libgen.h>
#include "options.h"
#include "processor.h"

void gspl_usage(Options* opt, char* arg0){
    fprintf(stderr, "\nUsage: %s %s [options]\n\n", PACKAGE_NAME, arg0);
    fprintf(stderr, "Options: -i FILE  read1 input file\n");
    fprintf(stderr, "         -I FILE  read2 input file\n");
    fprintf(stderr, "         -c FILE  configure file(3 column TSV:[ampliconName,fwdPrimer,revPrimer])\n");
    fprintf(stderr, "         -o DIR   output directory [%s]\n", opt->outdir.c_str());
    fprintf(stderr, "         -z INT   max N bases allowed in non-primer region [%d]\n", opt->lowq.maxN);
    fprintf(stderr, "         -q INT   min base quality as poor quality base [%d]\n", opt->lowq.minBaseQual);
    fprintf(stderr, "         -Q FLOAT max rate of poor quality bases allowed in one read [%f]\n", opt->lowq.maxLowQualFrac);
    fprintf(stderr, "         -P INT   max number of poor quality bases allowed in one read [%d]\n", opt->lowq.maxLowQualBase);
    fprintf(stderr, "         -1 INT   read1 skip length before split match [%ld]\n", opt->skipfr1);
    fprintf(stderr, "         -2 INT   read2 skip length before split match [%ld]\n", opt->skipfr2);
    fprintf(stderr, "         -d STR   drop library name [%s]\n", opt->droplib.c_str());
    fprintf(stderr, "         -t INT   threads used to split and write [%d]\n", opt->thread);
    fprintf(stderr, "         -p INT   max packs in buffer [%d]\n", opt->maxpack);
    fprintf(stderr, "         -r INT   max reads in each buffer pack [%d]\n", opt->maxreads);
    fprintf(stderr, "         -l INT   compression level if output gz file [%d]\n", opt->compression);
    fprintf(stderr, "         -f INT   max offset allowed before match begin [%d]\n", opt->maxoff);
    fprintf(stderr, "         -m INT   max mismatch allowed in each match [%d]\n", opt->maxmm);
    fprintf(stderr, "         -9 FILE  js/css cdn configure json file\n");
    fprintf(stderr, "         -S       only one end match is okay\n");
    fprintf(stderr, "         -g       use amplicon alignment instead of primer\n");
    fprintf(stderr, "         -Z       output gz file\n");
    fprintf(stderr, "         -D       drop matched prefix seq in output file\n");
    fprintf(stderr, "         -A       adjust output r1/2 as forward/reverse reads\n\n");
    delete opt;
}

int gspl_main(int argc, char** argv){
    // parse args
    Options* opt = new Options();
    if(argc == 1){
        gspl_usage(opt, argv[0]);
        return 0;
    }
    // parse args
    const char* args = "f:m:i:I:c:o:z:q:Q:P:1:2:d:t:p:r:l:9:gZDASh";
    int c = -1;
    while((c = getopt(argc, argv, args)) >= 0){
        switch(c){
            case 'f': opt->maxoff = atoi(optarg); break;
            case 'm': opt->maxmm = atoi(optarg); break;
            case 'i': opt->in1 = strdup(optarg); break;
            case 'I': opt->in2 = strdup(optarg); break;
            case 'c': opt->prf = strdup(optarg); break;
            case 'z': opt->lowq.maxN = atoi(optarg); break;
            case 'q': opt->lowq.minBaseQual = atoi(optarg); break;
            case 'Q': opt->lowq.maxLowQualFrac = atoi(optarg); break;
            case 'P': opt->lowq.maxLowQualBase = atoi(optarg); break;
            case '1': opt->skipfr1 = atoi(optarg); break;
            case '2': opt->skipfr2 = atoi(optarg); break;
            case 'o': opt->outdir = strdup(optarg); break;
            case 'd': opt->droplib = strdup(optarg); break;
            case 't': opt->thread = atoi(optarg); break;
            case 'p': opt->maxpack = atoi(optarg); break;
            case 'r': opt->maxreads = atoi(optarg); break;
            case 'l': opt->compression = atoi(optarg); break;
            case '9': opt->hrjsn = optarg; break;
            case 'S': opt->sematch = true; break;
            case 'g': opt->s8aln = true; break;
            case 'Z': opt->outgz = true; break;
            case 'D': opt->droppre = true; break;
            case 'A': opt->adjfr = true; break;
            case 'h': gspl_usage(opt, argv[0]); return 0; break;
            default: break;
        }
    }
    opt->fq2spl = true;
    opt->outspl = true;
    // check args and run
    if(opt->valid()){
        util::loginfo("beg parse args");
        opt->init(argc, argv);
        util::loginfo("end parse args");
        util::loginfo("beg split library");
        Processor *p = new Processor(opt);
        p->init();
        p->process4split();
        util::loginfo("end split library");
        delete p;
        delete opt;
        return 0;
    }else{
        delete opt;
        return 1;
    }
}

#endif
