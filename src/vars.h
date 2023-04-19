#ifndef GET_VARS_H
#define GET_VARS_H

#include "util.h"
#include "htslib/sam.h"
#include <unordered_map>

typedef std::vector<std::unordered_map<std::string, int32_t>> SeqVarCntMap; // pos, vars
typedef std::vector<std::vector<int32_t>> LenVarCntVector; // pos, lengths
typedef std::vector<std::vector<int32_t>> SNPVarCntVector; // pos, muts

#endif
