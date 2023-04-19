#include "grpmut.h"

void grp_mut_t::rec2th(kstring_t* ss, const std::string& mut, size_t maxp) {
    ksprintf(ss, "<tr>\n");
    ksprintf(ss, "<td>grp_of_%s</td>\n", group.c_str());
    ksprintf(ss, "<td>grp_of_%s</td>\n", group.c_str());
    ksprintf(ss, "<td>%s</td>\n", group.c_str());
    ksprintf(ss, "<td>%s</td>\n", mut.c_str());
    size_t i = 0;
    for(i = 0; i < pos.size(); ++i) {
        ksprintf(ss, "<td>%s</td>\n", pos[i].c_str());
    }
    for(; i < maxp; ++i) {
        ksprintf(ss, "<td> </td>\n");
    }
    ksprintf(ss, "</tr>\n");
}

void grp_mut_t::recm2th(kstring_t* ss, const std::string& mut, const std::string& altx, size_t maxp) {
    ksprintf(ss, "<tr>\n");
    ksprintf(ss, "<td>grp_of_%s</td>\n", group.c_str());
    ksprintf(ss, "<td>grp_of_%s</td>\n", group.c_str());
    ksprintf(ss, "<td>%s</td>\n", group.c_str());
    ksprintf(ss, "<td>%s</td>\n", mut.c_str());
    ksprintf(ss, "<td>%s</td>\n", altx.c_str());
    size_t i = 0;
    for(i = 0; i < pos.size(); ++i) {
        ksprintf(ss, "<td>%s</td>\n", pos[i].c_str());
    }
    for(; i < maxp; ++i) {
        ksprintf(ss, "<td> </td>\n");
    }
    ksprintf(ss, "</tr>\n");
}

void grp_mut_t::rec2tsvh(kstring_t* ss, const std::string& mut, size_t maxp) {
    ksprintf(ss, "grp_of_%s\t", group.c_str());
    ksprintf(ss, "grp_of_%s\t", group.c_str());
    ksprintf(ss, "%s\t", group.c_str());
    ksprintf(ss, "%s\t", mut.c_str());
    size_t i = 0;
    for(i = 0; i < pos.size(); ++i) {
        ksprintf(ss, "%s\t", pos[i].c_str());
    }
    for(; i < maxp; ++i) {
        ksprintf(ss, "\t");
    }
    ss->s[ss->l-1] = '\n';
}

void grp_mut_t::recm2tsvh(kstring_t* ss, const std::string& mut, const std::string& altx, size_t maxp) {
    ksprintf(ss, "grp_of_%s\t", group.c_str());
    ksprintf(ss, "grp_of_%s\t", group.c_str());
    ksprintf(ss, "%s\t", group.c_str());
    ksprintf(ss, "%s\t", mut.c_str());
    ksprintf(ss, "%s\t", altx.c_str());
    size_t i = 0;
    for(i = 0; i < pos.size(); ++i) {
        ksprintf(ss, "%s\t", pos[i].c_str());
    }
    for(; i < maxp; ++i) {
        ksprintf(ss, "\t");
    }
    ss->s[ss->l-1] = '\n';
}

void grp_mut_t::rec2td(kstring_t* ss, const std::string& mut, size_t maxp, bool ofreq) {
    ksprintf(ss, "<tr>\n");
    ksprintf(ss, "<td>%s</td>\n", sample.c_str());
    ksprintf(ss, "<td>%s</td>\n", amplicon.c_str());
    ksprintf(ss, "<td>%s</td>\n", group.c_str());
    ksprintf(ss, "<td>%s</td>\n", mut.c_str());
    size_t i = 0;
    if(ofreq) {
        for(i = 0; i < freq.size(); ++i) {
            ksprintf(ss, "<td>%lf</td>\n", freq[i]);
        }
    } else {
        for(i = 0; i < dist.size(); ++i) {
            ksprintf(ss, "<td>%d</td>\n", dist[i]);
        }
    }
    for(; i < maxp; ++i) {
        ksprintf(ss, "<td> </td>\n");
    }
    ksprintf(ss, "</tr>\n");
}

void grp_mut_t::recm2td(kstring_t* ss, const std::string& mut, size_t maxp, bool ofreq) {
    size_t ncol = freq.size()/row.size(), i = 0, j = 0;
    for(size_t nrow = 0; nrow < row.size(); ++nrow) {
        ksprintf(ss, "<tr>\n");
        ksprintf(ss, "<td>%s</td>\n", sample.c_str());
        ksprintf(ss, "<td>%s</td>\n", amplicon.c_str());
        ksprintf(ss, "<td>%s</td>\n", group.c_str());
        ksprintf(ss, "<td>%s</td>\n", mut.c_str());
        ksprintf(ss, "<td>%s</td>\n", row[nrow].c_str());
        if(ofreq) {
            for(j = 0; j < ncol; ++j) {
                ksprintf(ss, "<td>%lf</td>\n", freq[i++]);
            }
        } else {
            for(j = 0; j < dist.size(); ++j) {
                ksprintf(ss, "<td>%d</td>\n", dist[i++]);
            }
        }
        for(; j < maxp; ++j) {
            ksprintf(ss, "<td> </td>\n");
        }
        ksprintf(ss, "</tr>\n");
    }
}

void grp_mut_t::rec2tsvd(kstring_t* ss, const std::string& mut, size_t maxp, bool ofreq) {
    ksprintf(ss, "%s\t", sample.c_str());
    ksprintf(ss, "%s\t", amplicon.c_str());
    ksprintf(ss, "%s\t", group.c_str());
    ksprintf(ss, "%s\t", mut.c_str());
    size_t i = 0;
    if(ofreq) {
        for(i = 0; i < freq.size(); ++i) {
            ksprintf(ss, "%lf\t", freq[i]);
        }
    } else {
        for(i = 0; i < dist.size(); ++i) {
            ksprintf(ss, "%d\t", dist[i]);
        }
    }
    for(; i < maxp; ++i) {
        ksprintf(ss, "\t");
    }
    ss->s[ss->l-1] = '\n';
}

void grp_mut_t::recm2tsvd(kstring_t* ss, const std::string& mut, size_t maxp, bool ofreq) {
    size_t ncol = freq.size()/row.size(), i = 0, j = 0;
    for(size_t nrow = 0; nrow < row.size(); ++nrow) {
        ksprintf(ss, "%s\t", sample.c_str());
        ksprintf(ss, "%s\t", amplicon.c_str());
        ksprintf(ss, "%s\t", group.c_str());
        ksprintf(ss, "%s\t", mut.c_str());
        ksprintf(ss, "%s\t", row[nrow].c_str());
        if(ofreq) {
            for(j = 0; j < ncol; ++j) {
                ksprintf(ss, "%lf\t", freq[i++]);
            }
        } else {
            for(j = 0; j < ncol; ++j) {
                ksprintf(ss, "%d\t", dist[i++]);
            }
        }
        for(; j < maxp; ++j) {
            ksprintf(ss, "\t");
        }
        ss->s[ss->l-1] = '\n';
    }
}
