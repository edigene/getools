#include "b2g.h"
#include "gtn.h"
#include "ksw4get.h"
#include "sbc.h"
#include "scana.h"
#include "gs2as.h"
#include "pscnt.h"
#include "pscmp.h"
#include "stats.h"
#include "gecal.h"
#include "gespl.h"
#include "flags.h"
#include "derep.h"
#include "gef2b.h"
#include "getrpt.h"
#include "hapcnt.h"
#include "extrds.h"
#include "geplot.h"
#include "sgraln.h"
#include "guessbc.h"
#include "mergepairs.h"

inline int usage(){
    fprintf(stderr, "\n");
    fprintf(stderr, "Program: %s (amplicon based gene edit evaluation tools)\n", PACKAGE_NAME);
    fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
    fprintf(stderr, "Contact: %s <%s>\n\n", PACKAGE_MAINTER_NAME, PACKAGE_MAINTER_EMAIL);
    fprintf(stderr, "Usage: getools <command> [options]\n\n");
    fprintf(stderr, "Command: caledit        compute edit efficiency\n");
    fprintf(stderr, "         split          split mixed amplicon library to FASTQ\n");
    fprintf(stderr, "         fq2bam         split mixed amplicon library to BAM\n");
    fprintf(stderr, "         geplot         gene edit event plot\n");
    fprintf(stderr, "         view           view BAM from caledit\n");
    fprintf(stderr, "         stats          stat variant count and length\n");
    fprintf(stderr, "         flags          explain VS tag in BAM from caledit\n");
    fprintf(stderr, "         mergepe        merge a pair of sequence/reads\n");
    fprintf(stderr, "         mergefq        merge a pair of FASTQ files\n");
    fprintf(stderr, "         derep          dedup FASTQ file\n");
    fprintf(stderr, "         ksw            SW alignment of sequences\n");
    fprintf(stderr, "         gtncmp         consistence of topN edit events\n");
    fprintf(stderr, "         scana          single cell analysis from BCF\n");
    fprintf(stderr, "         bcf2fa         getools topN edit events to constructed fasta\n");
    fprintf(stderr, "         gs2as          convert results of bcf2fa to simulation table\n");
    fprintf(stderr, "         hapcnt         haplotype snp calling from getools BAM\n");
    fprintf(stderr, "         guessb         guess barcode from FASTQ\n");
    fprintf(stderr, "         sgraln         guess sgRNA_PAM by alignment\n");
    fprintf(stderr, "         mbspl          split mission bio single cell FASTQ\n");
    fprintf(stderr, "         mbidx          build mission bio cellbarcode whitelist index\n");
    fprintf(stderr, "         pscnt          pattern string count\n");
    fprintf(stderr, "         pscmp          pattern results compare\n");
    fprintf(stderr, "         gerpt          generate reports from batch of results\n");
    fprintf(stderr, "\n");
    return 0;
}

int main(int argc, char** argv){
    int ret = 1;
    if(argc < 2) return usage();
    if(strcmp(argv[1], "caledit") == 0) ret = gec_main(argc-1, argv+1);
    else if(strcmp(argv[1], "ksw") == 0) ret = kswge_main(argc-1, argv+1);
    else if(strcmp(argv[1], "flags") == 0) ret = flags_main(argc-1, argv+1);
    else if(strcmp(argv[1], "view") == 0) ret = extrds_main(argc-1, argv+1);
    else if(strcmp(argv[1], "stats") == 0) ret = stats_main(argc-1, argv+1);
    else if(strcmp(argv[1], "geplot") == 0) ret = geplot_main(argc-1, argv+1);
    else if(strcmp(argv[1], "mergepe") == 0) ret = mergepe_main(argc-1, argv+1);
    else if(strcmp(argv[1], "mergefq") == 0) ret = mergefq_main(argc-1, argv+1);
    else if(strcmp(argv[1], "derep") == 0) ret = derep_main(argc-1, argv+1);
    else if(strcmp(argv[1], "split") == 0) ret = gspl_main(argc-1, argv+1);
    else if(strcmp(argv[1], "fq2bam") == 0) ret = gef2b_main(argc-1, argv+1);
    else if(strcmp(argv[1], "fq2bam") == 0) ret = gef2b_main(argc-1, argv+1);
    else if(strcmp(argv[1], "gtncmp") == 0) ret = gtn_main(argc-1, argv+1);
    else if(strcmp(argv[1], "scana") == 0) ret = scana_main(argc-1, argv+1);
    else if(strcmp(argv[1], "bcf2fa") == 0) ret = b2g_main(argc-1, argv+1);
    else if(strcmp(argv[1], "gs2as") == 0) ret = gs2as_main(argc-1, argv+1);
    else if(strcmp(argv[1], "hapcnt") == 0) ret = hapc_main(argc-1, argv+1);
    else if(strcmp(argv[1], "guessb") == 0) ret = gbc_main(argc-1, argv+1);
    else if(strcmp(argv[1], "sgraln") == 0) ret = sgraln_main(argc-1, argv+1);
    else if(strcmp(argv[1], "mbspl") == 0) ret = sbc_biom_main(argc-1, argv+1);
    else if(strcmp(argv[1], "mbidx") == 0) ret = sbc_bidx_main(argc-1, argv+1);
    else if(strcmp(argv[1], "pscnt") == 0) ret = pscnt_main(argc-1, argv+1);
    else if(strcmp(argv[1], "pscmp") == 0) ret = pscmp_main(argc-1, argv+1);
    else if(strcmp(argv[1], "gerpt") == 0) ret = rpt_main(argc-1, argv+1);
    else{
        fprintf(stderr, "unrecognized command '%s'\n", argv[1]);
        return 1;
    }
    return ret;
}
