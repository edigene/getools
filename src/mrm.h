#ifndef MATCH_MATCH_MAP_H
#define MATCH_MATCH_MAP_H

#include <map>

struct mm_mat_t{
    int mm = 0;
    int scr = 0;
    int off = 0;
    char strand = '+';
};

typedef std::map<int, mm_mat_t> MatResutMap;

#endif
