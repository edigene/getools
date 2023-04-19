#ifndef GE_CAL_H
#define GE_CAL_H

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <libgen.h>
#include "common.h"
#include "options.h"
#include "processor.h"

void gec_usage(Options* opt, char* arg0){
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
    fprintf(stderr, "\nEdit efficiency computation options:\n");
    fprintf(stderr, "         -a FLOAT min frequency of valid edit events [%f]\n", opt->edo.minaf);
    fprintf(stderr, "         -n INT   min support sequences of valid edit events [%d]\n", opt->edo.minseqc);
    fprintf(stderr, "         -k INT   min extra match(rather than primers) needed for valid large deletion [%d]\n", opt->edo.minextram);
    fprintf(stderr, "         -x INT   max variant count allowed for valid edit event [%d]\n", opt->edo.maxvc);
    fprintf(stderr, "         -s INT   cleavage beg position(1 based, relative to sgRNA_PAM end, negative for sgRNA_PAM beg) [%d]\n", opt->edo.clsbpos);
    fprintf(stderr, "         -e INT   cleavage end position(1 based, relative to sgRNA_PAN end, negative for sgRNA_PAM end) [%d]\n", opt->edo.clsepos);
    fprintf(stderr, "         -b INT   max window size around cleavage position range[-s, -e] valid edit events occured [%d]\n", opt->edo.cutbuflen);
    fprintf(stderr, "         -v INT   bit masks of variant type of valid edit events ([%d]:SNV,[%d]:INSERTION,[%d]:DELETION,[%d]:DELINS) [%d]\n",
                                       GEVAR_SNV, GEVAR_INS, GEVAR_DEL, GEVAR_DIN, opt->edo.vartypem);
    fprintf(stderr, "         -7 STR   extra snv to focus on, format(REF:ALT,REF:ALT2,...)\n");
    fprintf(stderr, "         -w INT   recombination type of valid edit events ([%d]:Exact,[%d]:AllAHit,[%d]:AnyHit) [%d]\n", 
                                       REC_EXACT_IS_EDIT, REC_ALLHIT_IS_EDIT, REC_ANYHIT_IS_EDIT, opt->edo.rectypem);
    fprintf(stderr, "         -j INT   min match length from ends of amplicon needed [%d]\n", opt->edo.mpmatch);
    fprintf(stderr, "         -F INT   flank length around sgRNA to plot edit event graph [%d]\n", opt->edo.flklen);
    fprintf(stderr, "         -N INT   top N most frequent events to plot in edit event graph [%d]\n", opt->edo.topn);
    fprintf(stderr, "         -K       do single cell gene edit analysis(make sure cellbarcode was put as comment in FASTQ)\n");
    fprintf(stderr, "         -G       treat incompatible indel in overlap read pairs as real edit event\n");
    fprintf(stderr, "         -C       filter out events with low alignment score and many variants\n");
    fprintf(stderr, "         -Y       input library is from mitochontron\n");
    fprintf(stderr, "         -R       only use single matched read to do analysys(if both match, use read1)\n");
    fprintf(stderr, "\nInput output options:\n");
    fprintf(stderr, "         -i FILE  read1 input file\n");
    fprintf(stderr, "         -I FILE  read2 input file\n");
    fprintf(stderr, "         -T FILE  control sample BAM generated by getools\n");
    fprintf(stderr, "         -U INT   control variant support sequences threshold to filter [%d]\n", opt->edo.ctfc);
    fprintf(stderr, "         -c FILE  configure file(5 column TSV:[ampliconName,fwdPrimer,revPrimer,ampliconSeq,sgRNA_PAM])\n");
    fprintf(stderr, "         -o DIR   output directory [%s]\n", opt->outdir.c_str());
    fprintf(stderr, "         -t INT   threads used to do computation [%d]\n", opt->thread);
    fprintf(stderr, "         -p INT   max packs in buffer [%d]\n", opt->maxpack);
    fprintf(stderr, "         -r INT   max reads in each buffer pack [%d]\n", opt->maxreads);
    fprintf(stderr, "         -u INT   max number of samples processed one time, 0 for no limit[%d]\n", opt->memone);
    fprintf(stderr, "         -l INT   compression level if output gz file [%d]\n", opt->compression);
    fprintf(stderr, "         -d STR   drop library name [%s]\n", opt->droplib.c_str());
    fprintf(stderr, "         -5 INT   bit mask to output bam([%d]:full bam, [%d]:simple bam) [%d]\n", GEOUT_FULL_BAM, GEOUT_SIMP_BAM, opt->f4bout);
    fprintf(stderr, "         -9 FILE  js/css cdn configure json file\n");
    fprintf(stderr, "         -8       hap snv count in sgrpam region\n");
    fprintf(stderr, "         -L       output split library\n");
    fprintf(stderr, "         -Z       output split library in gzip format\n");
    fprintf(stderr, "         -A       adjust output r1/2 as forward/reverse reads\n\n");
    delete opt;
}

int gec_main(int argc, char** argv){
    // parse args
    Options* opt = new Options();
    if(argc == 1){
        gec_usage(opt, argv[0]);
        return 0;
    }
    // parse args
    const char* args = "f:m:1:2:z:q:Q:P:M:X:O:E:a:k:x:b:s:e:v:w:n:j:F:N:i:I:U:T:c:o:d:5:t:p:r:l:u:9:3:7:g8RZLAGCYSKh";
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
            // Edit efficiency computation options
            case 'a': opt->edo.minaf = atof(optarg); break;
            case 'n': opt->edo.minseqc = atoi(optarg); break;
            case 'k': opt->edo.minextram = atoi(optarg); break;
            case 'x': opt->edo.maxvc = atoi(optarg); break;
            case 'b': opt->edo.cutbuflen = atoi(optarg); break;
            case 's': opt->edo.clsbpos = atoi(optarg); break;
            case 'e': opt->edo.clsepos = atoi(optarg); break;
            case 'v': opt->edo.vartypem = atoi(optarg); break;
            case 'w': opt->edo.rectypem = atoi(optarg); break;
            case 'j': opt->edo.mpmatch = atoi(optarg); break;
            case 'F': opt->edo.flklen = atoi(optarg); break;
            case 'N': opt->edo.topn = atoi(optarg); break;
            case 'G': opt->edo.greedy = true; break;
            case 'C': opt->edo.cluster = true; break;
            case 'Y': opt->edo.ismt = true; break;
            case 'K': opt->dosccal = true; break;
            case 'R': opt->usesem = true; break;
            case '7': opt->extsnv = optarg; break;
            // Input output options
            case 'i': opt->in1 = strdup(optarg); break;
            case 'I': opt->in2 = strdup(optarg); break;
            case 'T': opt->ctrl = strdup(optarg); break;
            case 'U': opt->edo.ctfc = atoi(optarg); break;
            case 'c': opt->prf = strdup(optarg); break;
            case 'o': opt->outdir = strdup(optarg); break;
            case 'd': opt->droplib = strdup(optarg); break;
            case '5': opt->f4bout = atoi(optarg); break;
            case 't': opt->thread = atoi(optarg); break;
            case 'p': opt->maxpack = atoi(optarg); break;
            case 'r': opt->maxreads = atoi(optarg); break;
            case 'l': opt->compression = atoi(optarg); break;
            case 'u': opt->memone = atoi(optarg); break;
            case '9': opt->hrjsn = optarg; break;
            case 'Z': opt->outgz = true; break;
            case 'L': opt->outspl = true; break;
            case 'A': opt->adjfr = true; break;
            case '8': opt->edo.hapcnt = true; break;
            case 'h': gec_usage(opt, argv[0]); return 0;  break;
            default: break;
        }
    }
    opt->fq2cal = true;
    // check args and run
    if(opt->valid()){
        util::loginfo("beg parse args");
        if(opt->memone > 0) opt->memlow = true;
        opt->init(argc, argv);
        util::loginfo("end parse args");
        util::loginfo("beg computation");
        Processor *p = new Processor(opt);
        p->init();
        if(opt->memlow) p->process4caledit2();
        else p->process4caledit1();
        util::loginfo("end computation");
        delete p;
        delete opt;
        return 0;
    }else{
        delete opt;
        return 1;
    }
}

#endif
