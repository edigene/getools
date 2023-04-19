#ifndef HTML_UTIL_H
#define HTML_UTIL_H

#include "htslib/kstring.h"
#include <stdint.h>
#include <string>

namespace htmlutil{
    inline void outTableRow(kstring_t* s, const char* key, int64_t val){
        ksprintf(s, "<tr><td class='col1'>%s</td><td class='col2'>%lld</td></tr>\n", key, val);
    }

    inline void outTableRow(kstring_t* s, const char* key, const char* val){
        ksprintf(s, "<tr><td class='col1'>%s</td><td class='col2'>%s</td></tr>\n", key, val);
    }

    inline void outTableRow(kstring_t* s, std::string& key, int64_t val){
        outTableRow(s, key.c_str(), val);
    }

    inline void outTableRow(kstring_t* s, std::string& key, std::string& val){
        outTableRow(s, key.c_str(), val.c_str());
    }

    template<typename T>
    inline std::string formatNumber(T number){
        double num = static_cast<double>(number);
        std::string unit[6] = {"", "K", "M", "G", "T", "P"};
        int order = 0;
        while(num > 1000.0){
            order += 1;
            num /= 1000.0;
        }
        if(order == 0) return std::to_string(num);
        else return std::to_string(num) + " " + unit[order];
    }

    template<typename T>
    inline std::string formatNumber2(T number){
        double num = static_cast<double>(number);
        std::string unit[6] = {"", "K", "M", "G", "T", "P"};
        int order = 0;
        while(num > 1000.0){
            order += 1;
            num /= 1000.0;
        }
        if(order == 0) return std::to_string(num);
        else return std::to_string(num) + unit[order];
    }
};
#endif
