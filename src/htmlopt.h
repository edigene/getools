#ifndef HTML_OPTION_H
#define HTML_OPTION_H

#include "util.h"
#include "kson.h"
#include "comjci.h"
#include "htslib/kstring.h"
#include <unordered_map>

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

typedef std::unordered_map<std::string, std::string> Str2StrHash;

struct res_schem_js{
    // js
    std::string jszip_js = "jszip.min.js";
    std::string plotly_js = "plotly.min.js";
    std::string jquery_js = "jquery.min.js";
    std::string buttons_html5_js = "buttons.html5.min.js";
    std::string bootstrap_bundle_js = "bootstrap.bundle.min.js";
    std::string jquery_dataTables_js = "jquery.dataTables.min.js";
    std::string dataTables_buttons_js = "dataTables.buttons.min.js";
    std::string logo_bundle_js ="logojs.bundle.min.js";
    std::string druid_js = "druid.min.js";
};


struct res_schem_css{
    // css
    std::string bootstrp_css = "bootstrap.min.css";
    std::string jquery_dataTables_css = "jquery.dataTables.min.css";
    std::string buttons_dataTables_css = "buttons.dataTables.min.css";
};

// html options
struct HtmlOpt{
    // figure dimension
    int figw = 800; // default figure width
    int figh = 600; // default figure height
    int tnfigw = 800; // topn figure width
    int tnfigh = 1100; // topn figure height 20 * topn + 100
    // dir
    std::string subdir = "html"; // sub-directories of all amplicons
    // sechem
    res_schem_js schjs;
    res_schem_css schcss;
    // js address
    Str2StrHash hashjs;
    // css address
    Str2StrHash hashcss;

    // init
    void init(const std::string& cfg);
    // generate css/js lines in html
    void printJsAndCSSRes(kstring_t* s);
    // generate table export button java script section
    void printExportButtons(kstring_t* s, std::string tid, std::string fname);
    // generate general header
    void html2head(kstring_t* s, const std::string& title);
    // generate save2svg script
    void printJssave2svg(kstring_t* s);
    // generate general footer
    void html2foot(kstring_t* s);
    void html2endbody(kstring_t* s);
    void html2footonly(kstring_t* s);
};


#endif
