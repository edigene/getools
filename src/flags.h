#ifndef GET_FLAGS_H
#define GET_FLAGS_H


#include "ksw4get.h"
#include "mergepairs.h"
#include "htslib/kstring.h"
#include <string>
#include <map>

const std::map<KSW_FTYPE, std::string> KSW_FLAG_MAP = {
    {KSW_FLALN, "left alignment occured"},
    {KSW_FRALN, "right alignment occured"},
    {KSW_FIDLOCMP, "indel compatible in overlap region"},
    {KSW_FMAYVARSEQ, "gene edit event candidates"},
    {KSW_FREFTYPE, "reference type sequence"},
    {KSW_FMERGED, "sequence merged from pair end reads"},
    {KSW_FHIERRLOWQ, "low alignment score and many variants"},
    {KSW_FLOWFREQ, "sequence frequency too low"},
    {KSW_FSNVINRNG, "sequence has SNV in/around sgRNA"},
    {KSW_FINSINRNG, "sequence has INS in/around sgRNA"},
    {KSW_FDELINRNG, "sequence has DEL in/around sgRNA"},
    {KSW_FDIINRNG, "sequence has DELINS in/around sgRNA"},
    {KSW_FREPSEQR, "sequence is representative"},
    {KSW_FSPANSGR, "sequence spans sgRNA"},
    {KSW_FPAIROLP, "overlaps with its partner"},
    {KSW_FPRIMDIMER, "sequence is from primer dimer"},
    {KSW_FMISPRIMINN, "sequence is from mis-priming(complex INS-Match-DEL sequence)"},
    {KSW_FMISPRIMOUT, "sequence is from mis-priming(match beyond primer too short)"},
    {KSW_FLOWSURPT, "sequence number too small"},
    {KSW_FMANYVARS, "sequence has too many variants"},
    {KSW_FVARINSEQ, "sequence has variants of desired type"},
    {KSW_FRECANYHIT, "sequence has some donor sequence variants inside"},
    {KSW_FRECALLHIT, "sequence has all donor sequence variants inside"},
    {KSW_FRECEXACT, "sequence has exact match donor sequence inside"},
};

void flags_usage(char* arg0);
int flags_main(int argc, char** argv);

#endif
