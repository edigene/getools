#include "flags.h"

void flags_usage(char* arg0){
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: %s %s <flag>\n\n", PACKAGE_NAME, arg0);
}

int flags_main(int argc, char** argv){
    int c = 0;
    while((c = getopt(argc, argv, "h")) >= 0){
        switch(c){
            case 'h':
                flags_usage(argv[0]);
                return 0;
            default:
                break;
        }
    }
    if(argc == 1){
        fprintf(stderr, "\n");
        fprintf(stderr, "about: convert getools VS flag into textual representation\n");
        fprintf(stderr, "usage: %s %s FLAGS...\n\n", PACKAGE_NAME, argv[0]);
        fprintf(stderr, "each FLAGS is a decimal int representing a combination of the following numeric flag values\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "%-10s%8s        %s\n", "Hex", "Dec", "Meaning");
        for(int i = 0; i < 86; ++i)fprintf(stderr, "-");
        fprintf(stderr, "\n");
        for(auto& e: KSW_FLAG_MAP){
            fprintf(stderr, "0x%-8x%8u        %s\n", e.first, e.first, e.second.c_str());
        }
        fprintf(stderr, "\n");
    }else{
        KSW_FTYPE i = atol(argv[1]);
        kstring_t ks = {0, 0, 0};
        ksprintf(&ks, "\n");
        ksprintf(&ks, "0x%-10x: ", i);
        if(i & KSW_FLALN) ksprintf(&ks, "%s,", "LALN");
        if(i & KSW_FRALN) ksprintf(&ks, "%s,", "RALN");
        if(i & KSW_FIDLOCMP) ksprintf(&ks, "%s,", "IDLOCMP");
        if(i & KSW_FMAYVARSEQ) ksprintf(&ks, "%s,", "MAYVARSEQ");
        if(i & KSW_FREFTYPE) ksprintf(&ks, "%s,", "REFTYPE");
        if(i & KSW_FMERGED) ksprintf(&ks, "%s,", "MERGED");
        if(i & KSW_FHIERRLOWQ) ksprintf(&ks, "%s,", "HIERRLOWQ");
        if(i & KSW_FLOWFREQ) ksprintf(&ks, "%s,", "LOWFREQ");
        if(i & KSW_FSNVINRNG) ksprintf(&ks, "%s,", "SNVINRNG");
        if(i & KSW_FINSINRNG) ksprintf(&ks, "%s,", "INSINRNG");
        if(i & KSW_FDELINRNG) ksprintf(&ks, "%s,", "DELINRNG");
        if(i & KSW_FDIINRNG) ksprintf(&ks, "%s,", "DIINRNG");
        if(i & KSW_FREPSEQR) ksprintf(&ks, "%s,", "REPSEQR");
        if(i & KSW_FSPANSGR) ksprintf(&ks, "%s,", "SPANSGR");
        if(i & KSW_FPAIROLP) ksprintf(&ks, "%s,", "PAIROLP");
        if(i & KSW_FPRIMDIMER) ksprintf(&ks, "%s,", "PRIMDIMER");
        if(i & KSW_FMISPRIMINN) ksprintf(&ks, "%s,", "MISPRIMINN");
        if(i & KSW_FMISPRIMOUT) ksprintf(&ks, "%s,", "MISPRIMOUT"); 
        if(i & KSW_FLOWSURPT) ksprintf(&ks, "%s,", "LOWSURPT");
        if(i & KSW_FMANYVARS) ksprintf(&ks, "%s,", "MANYVARS");
        if(i & KSW_FVARINSEQ) ksprintf(&ks, "%s,", "VARINSEQ");
        if(i & KSW_FRECANYHIT) ksprintf(&ks, "%s,", "ANYRECVAR");
        if(i & KSW_FRECALLHIT) ksprintf(&ks, "%s,", "ALLRECVAR");
        if(i & KSW_FRECEXACT) ksprintf(&ks, "%s,", "RECEXACT");
        ks.s[ks.l-1] = '\n';
        if(get_merge_status(i) >= MERGE_REASON_CNT){
            ksprintf(&ks, "merge status: input merge status mask %d is invalid, must be in [0, %d]\n", get_merge_status(i), MERGE_REASON_CNT);
        }else{
            ksprintf(&ks, "merge status: %s\n", reason_str[get_merge_status(i)]);
        }
        fprintf(stderr, "%s\n", ks.s);
        free(ks.s);
    }
    return 0;
}
