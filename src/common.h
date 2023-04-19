#ifndef COMMON_GETOOLS_H
#define COMMON_GETOOLS_H

#ifndef PACKAGE_NAME
#define PACKAGE_NAME "getools"
#endif

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.1.0"
#endif

#ifndef PACKAGE_MAINTER_NAME
#define PACKAGE_MAINTER_NAME "Linjun Wu"
#endif

#ifndef PACKAGE_MAINTER_EMAIL
#define PACKAGE_MAINTER_EMAIL "wljlinksly@gmail.com"
#endif

#ifndef MAX_SAMPLE_IN_ONE_LIB
#define MAX_SAMPLE_IN_ONE_LIB 65535
#endif

#ifndef MAX_OFF_ALLOWED
#define MAX_OFF_ALLOWED 255
#endif


#ifndef MATCH_FAIL_RCNT
#define MATCH_FAIL_RCNT 8
#endif

#ifndef MATCH_FAIL_UNKNOWN
#define MATCH_FAIL_UNKNOWN 0
#endif

#ifndef MATCH_FAIL_ALN
#define MATCH_FAIL_ALN 1
#endif

#ifndef MATCH_FAIL_SCORE
#define MATCH_FAIL_SCORE 2
#endif

#ifndef MATCH_FAIL_MAXOFF
#define MATCH_FAIL_MAXOFF 3
#endif

#ifndef MATCH_FAIL_MAXMM
#define MATCH_FAIL_MAXMM 4
#endif

#ifndef MATCH_FAIL_MULTIPLE
#define MATCH_FAIL_MULTIPLE 5
#endif

#ifndef MATCH_FAIL_PEFAIL
#define MATCH_FAIL_PEFAIL 6
#endif

#ifndef MATCH_FAIL_QC
#define MATCH_FAIL_QC 7
#endif

static const char* MATCH_FAIL_REASONS[MATCH_FAIL_RCNT] = {
    "Unknown", "UnMapped", "LowMapScore", "MapOffSetTooLarge", "MismatchTooMany", "MultipleMatches", "PeMatchFail", "QCFail"
};

#ifndef GETOOLS_SUBLIBTSV
#define GETOOLS_SUBLIBTSV "sublib.tsv"
#endif

#ifndef GETOOLS_SPLITCMBTSV
#define GETOOLS_SPLITCMBTSV "splcmb.tsv"
#endif

#ifndef GETOOLS_SPLITTSV
#define GETOOLS_SPLITTSV "split.tsv"
#endif

#ifndef HAPCNT_OUTTTSV
#define HAPCNT_OUTTTSV "hapcnt.tsv"
#endif

#ifndef HAPCNT_OUTAAC
#define HAPCNT_OUTAAC "aaccnt.tsv"
#endif

#ifndef HAPCNT_OUTSNV
#define HAPCNT_OUTSNV "snvcnt.tsv"
#endif

#ifndef HAPCNT_OUTSTSV
#define HAPCNT_OUTSTSV "hapstat.tsv"
#endif

#ifndef HAPBIAS1_OUTTSV
#define HAPBIAS1_OUTTSV "hapbias1.tsv"
#endif

#ifndef HAPBIAS2_OUTTSV
#define HAPBIAS2_OUTTSV "hapbias2.tsv"
#endif

#ifndef HAPCNT_OUTMFAAC
#define HAPCNT_OUTMFAAC "mfaac.tsv"
#endif

#ifndef HAPCNT_OUTMFNUC
#define HAPCNT_OUTMFNUC "mfnuc.tsv"
#endif

#ifndef HAPCNT_OUTRCNTNALL2REF
#define HAPCNT_OUTRCNTNALL2REF "rcntnall2ref.tsv"
#endif


#ifndef HAPCNT_OUTRCNTNALL2ALT
#define HAPCNT_OUTRCNTNALL2ALT "rcntnall2alt.tsv"
#endif

#ifndef HAPCNT_OUTRCNTAALL2REF
#define HAPCNT_OUTRCNTAALL2REF "rcntaall2ref.tsv"
#endif


#ifndef HAPCNT_OUTRCNTAALL2ALT
#define HAPCNT_OUTRCNTAALL2ALT "rcntaall2alt.tsv"
#endif

#ifndef HAPCNT_OUTRCNTNYDEF2REF
#define HAPCNT_OUTRCNTNYDEF2REF "rcntnydef2ref.tsv"
#endif


#ifndef HAPCNT_OUTRCNTNYDEF2ALT
#define HAPCNT_OUTRCNTNYDEF2ALT "rcntnydef2alt.tsv"
#endif

#ifndef HAPCNT_OUTRCNTAYDEF2REF
#define HAPCNT_OUTRCNTAYDEF2REF "rcntaydef2ref.tsv"
#endif


#ifndef HAPCNT_OUTRCNTAYDEF2ALT
#define HAPCNT_OUTRCNTAYDEF2ALT "rcntaydef2alt.tsv"
#endif


#ifndef HAPCNT_OUTRCNTNNDEF2REF
#define HAPCNT_OUTRCNTNNDEF2REF "rcntnndef2ref.tsv"
#endif


#ifndef HAPCNT_OUTRCNTNNDEF2ALT
#define HAPCNT_OUTRCNTNNDEF2ALT "rcntnndef2alt.tsv"
#endif

#ifndef HAPCNT_OUTRCNTANDEF2REF
#define HAPCNT_OUTRCNTANDEF2REF "rcntandef2ref.tsv"
#endif


#ifndef HAPCNT_OUTRCNTANDEF2ALT
#define HAPCNT_OUTRCNTANDEF2ALT "rcntandef2alt.tsv"
#endif

#ifndef PSCNT_PSC_TSV
#define PSCNT_PSC_TSV "pscnt.tsv"
#endif

#ifndef PSCNT_RDC_TSV
#define PSCNT_RDC_TSV "rdcnt.tsv"
#endif

#ifndef PSCNT_NXC_TSV
#define PSCNT_NXC_TSV "nxcnt.tsv"
#endif

#ifndef SCANA_AMP_ALLELE_IN_SC_TSV
#define SCANA_AMP_ALLELE_IN_SC_TSV "ampallele4sc.tsv"
#endif

#ifndef SCANA_AMP_CO_EDIT_TSV
#define SCANA_AMP_CO_EDIT_TSV "ampcoedit4sc.tsv"
#endif

#ifndef SCANA_AMP_SC_EDIT_EFF
#define SCANA_AMP_SC_EDIT_EFF "ampediteff4sc.tsv"
#endif

#ifndef SCANA_CELL_FIND_JSON
#define SCANA_CELL_FIND_JSON "cellfind4sc.json"
#endif

#ifndef SCANA_OTG_EDIT_SC
#define SCANA_OTG_EDIT_SC "targetedit4sc.tsv"
#endif

#ifndef SCANA_HTML_RESULT
#define SCANA_HTML_RESULT "scana.html"
#endif

#ifndef CELL_BARCODE_ID_TAG
#define CELL_BARCODE_ID_TAG "CB"
#endif

#ifndef CELL_BARCODE_COUNT_TAG
#define CELL_BARCODE_COUNT_TAG "CN"
#endif

#define GEVAR_REF 0x0
#define GEVAR_SNV 0x1
#define GEVAR_INS 0x2
#define GEVAR_DEL 0x4
#define GEVAR_DIN 0x8
#define GEVAR_OTH 0x9

#define REC_EXACT_IS_EDIT 1
#define REC_ALLHIT_IS_EDIT 2
#define REC_ANYHIT_IS_EDIT 3

#define GEOUT_FULL_BAM 0x1
#define GEOUT_SIMP_BAM 0x2

#ifndef SC_ID_PREFIX
#define SC_ID_PREFIX "s"
#endif

#ifndef GT_HOMO_WTWT_TYPE
#define GT_HOMO_WTWT_TYPE 0
#endif

#ifndef GT_HETERO_WTMT_TYPE
#define GT_HETERO_WTMT_TYPE 1
#endif

#ifndef GT_HOMO_MTMT_TYPE
#define GT_HOMO_MTMT_TYPE 2
#endif

#ifndef GT_HETERO_MTMT2_TYPE
#define GT_HETERO_MTMT2_TYPE 3
#endif

#ifndef GT_MIX_MTLOWFREQ_TYPE
#define GT_MIX_MTLOWFREQ_TYPE 4
#endif

#ifndef GT_MIX_MTHIGHFREQ_TYPE
#define GT_MIX_MTHIGHFREQ_TYPE 5
#endif

#ifndef GT_TYPE_CNT
#define GT_TYPE_CNT 6
#endif

static const char* GT_INT2STR_ARR[6] = { 
    "Hom Wild", "Het Wt/Mt", 
    "Hom Mt/Mt", "Het Mt/Mt2",
    "Mix Mt.LowFreq", "Mix Mt.HighFreq",
};

static const char* GT_INT2VAR_ARR[6] = { 
    "gthomw", "gthetwm", "gthomm", "gthetm", "mixmtlow", "mixmthigh",
};

#ifndef DR_METHOD_NONE
#define DR_METHOD_NONE 0
#endif

#ifndef DR_METHOD_PCA
#define DR_METHOD_PCA 1
#endif

#ifndef DR_METHOD_UMAP
#define DR_METHOD_UMAP 2
#endif

#ifndef DR_METHOD_CNT
#define DR_METHOD_CNT 2
#endif

static const char* DR_METHOD2STR_ARR[3] ={
    "None", "PCA", "UMAP",
};

static const char* DR_METHOD2AXIS_ARR[3] ={
    "None", "pc", "umap",
};

#endif
