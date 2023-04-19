#ifndef GRP_MUT_H
#define GRP_MUT_H

#include <string>
#include <vector>
#include "htslib/kstring.h"

struct grp_mut_t{
    int nsmp;
    std::string sample;
    std::string amplicon;
    std::string group;
    std::vector<std::string> pos;
    std::vector<std::string> row;
    std::vector<double> freq;
    std::vector<int32_t> dist;

    void rec2th(kstring_t* ss, const std::string& mut, size_t maxp);
    void rec2tsvh(kstring_t* ss, const std::string& mut, size_t maxp);
    void rec2td(kstring_t* ss, const std::string& mut, size_t maxp, bool ofreq = true);
    void rec2tsvd(kstring_t* ss, const std::string& mut, size_t maxp, bool ofreq = true);

    void recm2th(kstring_t* ss, const std::string& mut, const std::string& altx, size_t maxp);
    void recm2tsvh(kstring_t* ss, const std::string& mut, const std::string& altx, size_t maxp);
    void recm2td(kstring_t* ss, const std::string& mut, size_t maxp, bool ofreq = true);
    void recm2tsvd(kstring_t* ss, const std::string& mut, size_t maxp, bool ofreq = true);
};

typedef std::vector<grp_mut_t*>  GroupMutList;

typedef std::vector<GroupMutList> GroupMutColType;

struct grp_mut_sorter{
    inline bool operator()(const grp_mut_t* gs1, const grp_mut_t* gs2) const {
        return (gs1->group < gs2->group) ||
               (gs1->group == gs2->group && gs1->sample < gs2->sample) ||
               (gs1->group == gs2->group && gs1->sample == gs2->sample && gs1->amplicon < gs2->amplicon);
    }
};

#endif
