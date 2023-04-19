#include "getrpt.h"

rpt_opt_t::rpt_opt_t(){
}

rpt_opt_t::~rpt_opt_t(){
    if(ks){
        if(ks->s){free(ks->s); ks->s = NULL;}
        free(ks); ks = NULL;
    }
    if(ps){
        if(ps->s){free(ps->s); ps->s = NULL;}
        free(ps); ps = NULL;
    }
    if(ss){
        if(ss->s){free(ss->s); ss->s = NULL;}
        free(ss); ss = NULL;
    }
    if(ts){
        if(ts->s){free(ts->s); ts->s = NULL;}
        free(ts); ts = NULL;
    }
    for(auto& sl: amnv){
        for(auto& pl: sl){
            if(pl){
                delete pl;
                pl = NULL;
            }
        }
    }
}

bool rpt_opt_t::valid(){
    if(outdir == NULL){
        fprintf(stderr, "output directory path must be provided\n");
        return false;
    }

    if(indir == NULL){
        fprintf(stderr, "getool result directory path must be provided\n");
        return false;
    }

    if(splist == NULL && spstr.empty()){
        fprintf(stderr, "sample list path or samples must be provided\n");
        return false;
    }
    if(maxmdistl > MAX_MUTMM_DIST_LEN){
        fprintf(stderr, "too many mutation in one read set, rediculous\n");
        return false;
    }
    return true;
}

void rpt_opt_t::init(){
    ks = (kstring_t*)calloc(1, sizeof(kstring_t));
    ps = (kstring_t*)calloc(1, sizeof(kstring_t));
    ss = (kstring_t*)calloc(1, sizeof(kstring_t));
    ts = (kstring_t*)calloc(1, sizeof(kstring_t));
    if(splist) util::makeListFromFileByLine(splist, samples);
    if(spstr.size()) util::split(spstr, samples, ",", true);
    outdir = strdup(util::abspath(outdir).c_str());
    indir = strdup(util::abspath(indir).c_str());
    if(strcmp(indir, outdir) == 0) gen4web = true;
    util::makedir(outdir);
    ksprintf(ps, "%s/%s", outdir, rptd.consdir.c_str());
    util::makedir(ps->s);
    ps->l = 0;
    ksprintf(ps, "%s/%s", outdir, rptd.qcdir.c_str());
    util::makedir(ps->s);
    ps->l = 0;
    ksprintf(ps, "%s/%s", outdir, rptd.sumdir.c_str());
    util::makedir(ps->s);
    ps->l = 0;
    ksprintf(ps, "%s/%s", outdir, rptd.rptdir.c_str());
    util::makedir(ps->s);
    ps->l = 0;
    ksprintf(ps, "%s/%s", outdir, rptd.hapdir.c_str());
    util::makedir(ps->s);
    ps->l = 0;
    ksprintf(ps, "%s/%s", outdir, rptd.snvdir.c_str());
    util::makedir(ps->s);
    ps->l = 0;
    for(auto& s: samples){
        ksprintf(ps, "%s/%s/%s", outdir, rptd.rptdir.c_str(), s.c_str());
        util::makedir(ps->s);
        ps->l = 0;
    }
    if(!title) title = strdup("Gene Edit Analysis Report by getools");
    ncol = ceil(sqrt(samples.size()));
    hmo.init(jscfg);
    if(nrow <= 0) nrow = ceil(sqrt(samples.size()));
    ksprintf(ps, "%s/%s/hapcnt", indir, samples[0].c_str());
    if(util::exists(ps->s)) haprpt = true;
    ps->l = 0;
    if(extsnv.size()){
        uint8_t oval = 0, sval = 0;
        std::string mnn, xnn;
        std::set<std::string> nset;
        std::vector<std::string> vstr;
        util::split(extsnv, vstr, ",");
        std::vector<std::string> mts;
        for(auto& p: vstr){
            if(p.empty()) continue;
            util::split(p, mts, ":");
            if(mts.size() < 2){
                fprintf(stderr, "Error format[%s] of muts, must be \"REF:ALT\"\n", p.c_str());
                exit(EXIT_FAILURE);
            }
            oval = nuc_to_3bit[(int)mts[0][0]];
            sval = nuc_to_3bit[(int)mts[1][0]];
            if(oval == 4){
                fprintf(stderr, "Error ref nucleotide provided[%c], must be A/C/G/T\n", mts[0][0]);
                exit(EXIT_FAILURE);
            }
            if(sval == 4){
                fprintf(stderr, "Error alt nucleotide provided[%c], must be A/C/G/T\n", mts[1][0]);
                exit(EXIT_FAILURE);
            }
            mnn = std::string(1, std::toupper(mts[0][0])) + std::string(1, std::toupper(mts[1][0]));
            xnn = std::string(1, std::toupper(mts[0][0])) + "2" + std::string(1, std::toupper(mts[1][0]));
            if(nset.find(mnn) == nset.end()){
                nmnv.push_back(xnn);
                outmsnv.push_back(mnn + ".tsv");
            }
        }
        amnv.resize(nmnv.size(), GroupMutList());
    }
    if(haprpt){
        nfreq = {"Nuc", "Aac"};
        outmfreq = {HAPCNT_OUTMFNUC, HAPCNT_OUTMFAAC};
        nrcnt = {"NucAll2Ref", "NucYesDef2Ref", "NucNoDef2Ref", 
                 "AacAll2Ref", "AacYesDef2Ref", "AacNoDef2Ref", 
                 "NucAll2Alt", "NucYesDef2Alt", "NucNoDef2Alt", 
                 "AacAll2Alt", "AacYesDef2Alt", "AacNoDef2Alt"};
        outrcnt = {HAPCNT_OUTRCNTNALL2REF, HAPCNT_OUTRCNTNYDEF2REF, HAPCNT_OUTRCNTNNDEF2REF, 
                   HAPCNT_OUTRCNTAALL2REF, HAPCNT_OUTRCNTAYDEF2REF, HAPCNT_OUTRCNTANDEF2REF,
                   HAPCNT_OUTRCNTNALL2ALT, HAPCNT_OUTRCNTNYDEF2ALT, HAPCNT_OUTRCNTNNDEF2ALT,
                   HAPCNT_OUTRCNTAALL2ALT, HAPCNT_OUTRCNTAYDEF2ALT, HAPCNT_OUTRCNTANDEF2ALT
        };
        nbias = {"OneBase", "TwoBase"};
        outbias = {HAPBIAS1_OUTTSV, HAPBIAS2_OUTTSV};
        afreq.resize(nfreq.size());
        arnct.resize(nrcnt.size());
        abias.resize(nbias.size());
    }
    if(ampgrpf.size()){
        if(util::exists(ampgrpf)){
            util::makeMapPairFromFileByLine(ampgrpf, amp2grp);
        }else{
            fprintf(stderr, "Amplicon group file does not exist:%s\n", ampgrpf.c_str());
            exit(EXIT_FAILURE);
        }
    }
}

void rpt_opt_t::reportHTMLHeader(){
    ksprintf(ks, "<head>\n");
    ksprintf(ks, "<meta http-equiv=\"Content-Type\" content=\"text/html; charset=UTF-8\">\n");
    ksprintf(ks, "<title>%s</title>\n", title);
    ksprintf(ks, "<style type=\"text/css\">\n");
    ksprintf(ks, "BODY { font-family : ariel, monospace, sans-serif; }\n");
    ksprintf(ks, "P { font-weight: normal; font-family : ariel, monospace, sans-serif; color: black; background-color: transparent;}\n");
    ksprintf(ks, "B { font-weight: normal; color: black; background-color: transparent;}\n");
    ksprintf(ks, "A:visited { font-weight : normal; text-decoration : none; background-color : transparent; margin : 0px 0px 0px 0px; padding : 0px 0px 0px 0px; display: inline; }\n");
    ksprintf(ks, "A:link    { font-weight : normal; text-decoration : none; margin : 0px 0px 0px 0px; padding : 0px 0px 0px 0px; display: inline; }\n");
    ksprintf(ks, "A:hover   { color : #000000; font-weight : normal; text-decoration : underline; background-color : yellow; margin : 0px 0px 0px 0px; padding : 0px 0px 0px 0px; display: inline; }\n");
    ksprintf(ks, "A:active  { color : #000000; font-weight: normal; background-color : transparent; margin : 0px 0px 0px 0px; padding : 0px 0px 0px 0px; display: inline; }\n");
    ksprintf(ks, ".VERSION { font-size: small; font-family : arial, sans-serif; }\n");
    ksprintf(ks, ".NORM  { color: black;  background-color: transparent;}\n");
    ksprintf(ks, ".FIFO  { color: purple; background-color: transparent;}\n");
    ksprintf(ks, ".CHAR  { color: yellow; background-color: transparent;}\n");
    ksprintf(ks, ".DIR   { color: blue;   background-color: transparent;}\n");
    ksprintf(ks, ".BLOCK { color: yellow; background-color: transparent;}\n");
    ksprintf(ks, ".LINK  { color: aqua;   background-color: transparent;}\n");
    ksprintf(ks, ".SOCK  { color: fuchsia;background-color: transparent;}\n");
    ksprintf(ks, ".EXEC  { color: green;  background-color: transparent;}\n");
    ksprintf(ks, "</style>\n");
    ksprintf(ks, "</head>\n");
}

void rpt_opt_t::reportHTMLBody(){
    ksprintf(ks, "<body>\n");
    int i = 1;
    // Pooled Library QC
    ksprintf(ks, "<div>%d. Pooled Library QC: <a href=\"./%s/%s\">web</a> <a href=\"./%s/%s\">tsv</a></div><br>\n",
                 i, rptd.qcdir.c_str(), rptd.qcfile.c_str(), rptd.qcdir.c_str(), rptd.tsv_qcfile.c_str());
    ++i;
    if(consamp){
        // Edit Consistence
        ksprintf(ks, "<a href=\"./%s/%s.gtn.html\">%d. Edit Consistence of %s</a><br><br>\n", rptd.consdir.c_str(), consamp, i, consamp);
        ++i;
    }
    if(haprpt){
        // Hapcnt Summary
        ksprintf(ks, "<div>%d. Hap SNV Count Summary: <a href=\"./%s/%s\">web</a> <a href=\"./%s/%s\">tsv</a></div><br>\n",
                     i, rptd.hapdir.c_str(), rptd.hapfile.c_str(), rptd.hapdir.c_str(), rptd.tsv_hapfile.c_str());
        ++i;
        ksprintf(ks, "<div>%d. Whole Amplicon-wide SNV Details: <a href=\"./%s/%s\">web</a> <a href=\"./%s/%s\">tsv</a></div><br>\n",
                     i, rptd.snvdir.c_str(), rptd.snvfile.c_str(), rptd.snvdir.c_str(), rptd.tsv_snvfile.c_str());
        ++i;
        ksprintf(ks, "<div>%d. Inframe Region Aminoacid Details: <a href=\"./%s/%s\">web</a> <a href=\"./%s/%s\">tsv</a></div><br>\n",
                     i, rptd.snvdir.c_str(), rptd.aacfile.c_str(), rptd.snvdir.c_str(), rptd.tsv_aacfile.c_str());
        ++i;
        if(extsnv.size()){
            ksprintf(ks, "<div>%d. Inframe Region SNV Details: <a href=\"./%s/%s\">web</a> <a href=\"./%s/%s\">tsv</a></div><br>\n",
                         i, rptd.snvdir.c_str(), rptd.snvfocf.c_str(), rptd.snvdir.c_str(), rptd.tsv_snvfocf.c_str());
            ++i;
        }
        ksprintf(ks, "<div>%d. Inframe Region Mutation Frequency: <a href=\"./%s/%s\">web</a> <a href=\"./%s/%s\">tsv</a></div><br>\n",
                     i, rptd.snvdir.c_str(), rptd.mutfreq.c_str(), rptd.snvdir.c_str(), rptd.tsv_mutfreq.c_str());
        ++i;
        ksprintf(ks, "<div>%d. Inframe Region Mutation Count per Read Distribution: <a href=\"./%s/%s\">web</a> <a href=\"./%s/%s\">tsv</a></div><br>\n",
                     i, rptd.snvdir.c_str(), rptd.mutn8rdist.c_str(), rptd.snvdir.c_str(), rptd.tsv_mutn8rdist.c_str());
        ++i;
        ksprintf(ks, "<div>%d. Whole Amplicon-wide Mutation Bias Distribution: <a href=\"./%s/%s\">web</a> <a href=\"./%s/%s\">tsv</a></div><br>\n",
                     i, rptd.snvdir.c_str(), rptd.sbiasfile.c_str(), rptd.snvdir.c_str(), rptd.tsv_sbiasfile.c_str());
        ++i;
    }
    // Edit Report
    ksprintf(ks, "<a href=\"./\">%d. Edit Report</a><br>\n", i);
    ksprintf(ks, "<table border=1>\n");
    std::vector<std::pair<int, int>> pidx;
    util::divideVecIdx(samples.size(), nrow, pidx, false);
    for(auto& p: pidx){
        ksprintf(ks, "<tr>\n");
        for(int i = p.first; i < p.second; ++i){
            const char* s = samples[i].c_str();
            if(gen4web){
                ksprintf(ks, "<td><a href=\"./%s/report.html\">%s</a></td>\n", s, s);
            }else{
                ksprintf(ks, "<td><a href=\"./%s/%s/index.html\">%s</a></td>\n", rptd.rptdir.c_str(), s, s);
            }
        }
        ksprintf(ks, "</tr>\n");
    }
    ksprintf(ks, "</table><br>\n");
    ++i;
    // Edit Summary
    ksprintf(ks, "<a href=\"./\">%d. Edit Summary</a><br>\n", i);
    ksprintf(ks, "<table border=1>\n");
    for(auto& p: pidx){
        ksprintf(ks, "<tr>\n");
        for(int i = p.first; i < p.second; ++i){
            const char* s = samples[i].c_str();
            ksprintf(ks, "<td><a href=\"./%s/%s.html\">%s</a></td>\n", rptd.sumdir.c_str(), s, s);
        }
        ksprintf(ks, "</tr>\n");
    }
    ksprintf(ks, "</table><br>\n");
    ++i;
    // Pooled Edit QC
    ksprintf(ks, "<div>%d. Pooled Edit QC: <a href=\"./%s/%s\">web</a> <a href=\"./%s/%s\">tsv</a></div><br>\n",
                 i, rptd.qcdir.c_str(), rptd.edfile.c_str(), rptd.qcdir.c_str(), rptd.tsv_edfile.c_str());
    ksprintf(ks, "</body>\n");
}

void rpt_opt_t::reportHTML(){
    ksprintf(ks, "<!DOCTYPE HTML>\n");
    ksprintf(ks, "<html>\n");
    reportHTMLHeader();
    reportHTMLBody();
    ksprintf(ks, "</html>\n");
    ps->l = 0;
    ksprintf(ps, "%s/index.html", outdir);
    util::writestr(ks->s, ks->l, ps->s);
}

void rpt_opt_t::gen_consf(){
    // write configure file
    ss->l = 0;
    for(auto& s: samples) ksprintf(ss, "%s\t%s/%s/var.bcf\n", s.c_str(), indir, s.c_str());
    ps->l = 0;
    ksprintf(ps, "%s/.cons.psv", outdir);
    util::writestr(ss->s, ss->l, ps->s);
    // generate consitence analysis file
    gtn_opt_t opt;
    opt.inlist = ps->s;
    opt.amp = consamp;
    opt.topn = na4cons;
    opt.uniqa = uniqa;
    opt.drmethod = drmethod;
    opt.showt = showt;
    opt.calg = calg;
    opt.hw = topnW;
    opt.beg = opt.end = -1;
    opt.hrjsn = jscfg;
    ss->l = 0;
    ksprintf(ss, "%s/%s/%s", outdir, rptd.consdir.c_str(), consamp);
    opt.outpre = ss->s;
    opt.hmo.init(opt.hrjsn);
    opt.cons_ana();
    // remove temp file
    remove(ps->s);
}

void rpt_opt_t::gen_hap(){
    std::vector<std::string> amps;
    ss->l = 0;
    std::vector<std::string> xsps, xamps, fsFile;
    for(auto& smp: samples){
        ss->l = 0;
        ksprintf(ss, "%s/%s/hapcnt", indir, smp.c_str());
        amps.clear();
        util::listDir(ss->s, amps);
        for(auto& amp: amps){
            ss->l = 0;
            ksprintf(ss, "%s/%s/hapcnt/%s/%s", indir, smp.c_str(), amp.c_str(), HAPCNT_OUTSTSV);
            if(util::exists(ss->s)){
                xsps.push_back(smp);
                xamps.push_back(amp);
                fsFile.push_back(ss->s);
            }
        }
    }
    if(fsFile.empty()){
        ss->l = 0;
        fprintf(stderr, "no hapstat file\n");
        return;
    }
    ss->l = ts->l = 0;
    // html begin
    ksprintf(ss, "<!DOCTYPE html>\n");
    ksprintf(ss, "<html>\n");
    writeHTMLTableHeader(ss, "hap_summary");
    // body begin
    ksprintf(ss, "<body>\n");
    // table begin
    ksprintf(ss, "<table id='hap_summary_table'>\n");
    ksprintf(ss, "<thead>\n");
    ksprintf(ss, "<tr>\n");
    ksprintf(ss, "<th>Sample</th>\n"); ksprintf(ts, "Sample\t");
    ksprintf(ss, "<th>Amplicon</th>\n"); ksprintf(ts, "Amplicon\t");
    std::ifstream fr(fsFile[0]);
    std::string tmpline;
    std::getline(fr, tmpline);
    std::vector<std::string> vstr;
    util::split(tmpline, vstr, "\t");
    for(auto& e: vstr){ ksprintf(ts, "%s\t", e.c_str()); }
    ts->s[ts->l-1] = '\n';
    tmpline = util::replace(tmpline, "<", "&lt;");
    tmpline = util::replace(tmpline, ">", "&gt;");
    util::split(tmpline, vstr, "\t");
    for(auto& e: vstr){ ksprintf(ss, "<th>%s</th>\n", e.c_str()); }
    ksprintf(ss, "</tr>\n");
    ksprintf(ss, "</thead>\n");
    ksprintf(ss, "<tbody>\n");
    fr.close();
    for(size_t i = 0; i < fsFile.size(); ++i){
        fr.open(fsFile[i].c_str());
        std::getline(fr, tmpline);
        std::getline(fr, tmpline);
        ksprintf(ss, "<tr>\n");
        ksprintf(ss, "<td>%s</td>\n", xsps[i].c_str()); ksprintf(ts, "%s\t", xsps[i].c_str());
        ksprintf(ss, "<td>%s</td>\n", xamps[i].c_str()); ksprintf(ts, "%s\t", xamps[i].c_str());
        util::split(tmpline, vstr, "\t");
        for(auto& e: vstr){
            ksprintf(ss, "<td>%s</td>\n", e.c_str());
            ksprintf(ts, "%s\t", e.c_str());
        }
        ts->s[ts->l-1] = '\n';
        fr.close();
    }
    ksprintf(ss, "</tbody>\n");
    // table end
    ksprintf(ss, "</table>\n");
    hmo.printExportButtons(ss, "hap_summary_table", "hap.summary");
    // body end
    ksprintf(ss, "</body>\n");
    // html end
    ksprintf(ss, "</html>\n");
    ps->l = 0;
    ksprintf(ps, "%s/%s/%s", outdir, rptd.hapdir.c_str(), rptd.hapfile.c_str());
    util::writestr(ss->s, ss->l, ps->s);
    ps->l = 0;
    ksprintf(ps, "%s/%s/%s", outdir, rptd.hapdir.c_str(), rptd.tsv_hapfile.c_str());
    util::writestr(ts->s, ts->l, ps->s);
}

void rpt_opt_t::gen_snv(){
    std::vector<std::string> amps;
    ss->l = 0;
    std::vector<std::string> xsps, xamps, fsFile;
    for(auto& smp: samples){
        ss->l = 0;
        ksprintf(ss, "%s/%s/hapcnt", indir, smp.c_str());
        amps.clear();
        util::listDir(ss->s, amps);
        for(auto& amp: amps){
            ss->l = 0;
            ksprintf(ss, "%s/%s/hapcnt/%s/%s", indir, smp.c_str(), amp.c_str(), HAPCNT_OUTSNV);
            if(util::exists(ss->s)){
                xsps.push_back(smp);
                xamps.push_back(amp);
                fsFile.push_back(ss->s);
            }
        }
    }
    if(fsFile.empty()){
        ss->l = 0;
        fprintf(stderr, "no snvstat file\n");
        return;
    }
    ss->l = ts->l = 0;
    // html begin
    ksprintf(ss, "<!DOCTYPE html>\n");
    ksprintf(ss, "<html>\n");
    writeHTMLTableHeader(ss, "snv_summary");
    // body begin
    ksprintf(ss, "<body>\n");
    // table begin
    ksprintf(ss, "<table id='snv_summary_table'>\n");
    ksprintf(ss, "<thead>\n");
    ksprintf(ss, "<tr>\n");
    ksprintf(ss, "<th>Sample</th>\n"); ksprintf(ts, "Samplt\t");
    ksprintf(ss, "<th>Amplicon</th>\n"); ksprintf(ts, "Amplicon\t");
    std::ifstream fr(fsFile[0]);
    std::string tmpline;
    std::getline(fr, tmpline);
    std::vector<std::string> vstr;
    util::split(tmpline, vstr, "\t");
    for(auto& e: vstr) ksprintf(ts, "%s\t", e.c_str());
    ts->s[ts->l-1] = '\n';
    tmpline = util::replace(tmpline, "<", "&lt;");
    tmpline = util::replace(tmpline, ">", "&gt;");
    util::split(tmpline, vstr, "\t");
    for(auto& e: vstr) ksprintf(ss, "<th>%s</th>\n", e.c_str());
    ksprintf(ss, "</tr>\n");
    ksprintf(ss, "</thead>\n");
    ksprintf(ss, "<tbody>\n");
    fr.close();
    for(size_t i = 0; i < fsFile.size(); ++i){
        fr.open(fsFile[i].c_str());
        std::getline(fr, tmpline);
        while(std::getline(fr, tmpline)){
            ksprintf(ss, "<tr>\n");
            ksprintf(ss, "<td>%s</td>\n", xsps[i].c_str()); ksprintf(ts, "%s\t", xsps[i].c_str());
            ksprintf(ss, "<td>%s</td>\n", xamps[i].c_str()); ksprintf(ts, "%s\t", xamps[i].c_str());
            util::split(tmpline, vstr, "\t");
            for(auto& e: vstr){
                ksprintf(ss, "<td>%s</td>\n", e.c_str());
                ksprintf(ts, "%s\t", e.c_str());
            }
            ts->s[ts->l-1] = '\n';
        }
        fr.close();
    }
    ksprintf(ss, "</tbody>\n");
    // table end
    ksprintf(ss, "</table>\n");
    hmo.printExportButtons(ss, "snv_summary_table", "snv.summary");
    // body end
    ksprintf(ss, "</body>\n");
    // html end
    ksprintf(ss, "</html>\n");
    ps->l = 0;
    ksprintf(ps, "%s/%s/%s", outdir, rptd.snvdir.c_str(), rptd.snvfile.c_str());
    util::writestr(ss->s, ss->l, ps->s);
    ps->l = 0;
    ksprintf(ps, "%s/%s/%s", outdir, rptd.snvdir.c_str(), rptd.tsv_snvfile.c_str());
    util::writestr(ts->s, ts->l, ps->s);
}

void rpt_opt_t::gen_aac(){
    std::vector<std::string> amps;
    ss->l = 0;
    std::vector<std::string> xsps, xamps, fsFile;
    for(auto& smp: samples){
        ss->l = 0;
        ksprintf(ss, "%s/%s/hapcnt", indir, smp.c_str());
        amps.clear();
        util::listDir(ss->s, amps);
        for(auto& amp: amps){
            ss->l = 0;
            ksprintf(ss, "%s/%s/hapcnt/%s/%s", indir, smp.c_str(), amp.c_str(), HAPCNT_OUTAAC);
            if(util::exists(ss->s)){
                xsps.push_back(smp);
                xamps.push_back(amp);
                fsFile.push_back(ss->s);
            }
        }
    }
    if(fsFile.empty()){
        ss->l = 0;
        fprintf(stderr, "no aaccnt file\n");
        return;
    }
    ss->l = ts->l = 0;
    // html begin
    ksprintf(ss, "<!DOCTYPE html>\n");
    ksprintf(ss, "<html>\n");
    writeHTMLTableHeader(ss, "aac_summary");
    // body begin
    ksprintf(ss, "<body>\n");
    // table begin
    ksprintf(ss, "<table id='aac_summary_table'>\n");
    ksprintf(ss, "<thead>\n");
    ksprintf(ss, "<tr>\n");
    ksprintf(ss, "<th>Sample</th>\n"); ksprintf(ts, "Sample\t");
    ksprintf(ss, "<th>Amplicon</th>\n"); ksprintf(ts, "Amplicon\t");
    std::ifstream fr(fsFile[0]);
    std::string tmpline;
    std::getline(fr, tmpline);
    std::vector<std::string> vstr;
    util::split(tmpline, vstr, "\t");
    for(auto& e: vstr) ksprintf(ts, "%s\t", e.c_str());
    ts->s[ts->l-1] = '\n';
    tmpline = util::replace(tmpline, "<", "&lt;");
    tmpline = util::replace(tmpline, ">", "&gt;");
    util::split(tmpline, vstr, "\t");
    for(auto& e: vstr) ksprintf(ss, "<th>%s</th>\n", e.c_str());
    ksprintf(ss, "</tr>\n");
    ksprintf(ss, "</thead>\n");
    ksprintf(ss, "<tbody>\n");
    fr.close();
    for(size_t i = 0; i < fsFile.size(); ++i){
        fr.open(fsFile[i].c_str());
        std::getline(fr, tmpline);
        while(std::getline(fr, tmpline)){
            ksprintf(ss, "<tr>\n");
            ksprintf(ss, "<td>%s</td>\n", xsps[i].c_str()); ksprintf(ts, "%s\t", xsps[i].c_str());
            ksprintf(ss, "<td>%s</td>\n", xamps[i].c_str()); ksprintf(ts, "%s\t", xamps[i].c_str());
            util::split(tmpline, vstr, "\t");
            for(auto& e: vstr){
                ksprintf(ss, "<td>%s</td>\n", e.c_str());
                ksprintf(ts, "%s\t", e.c_str());
            }
            ts->s[ts->l-1] = '\n';
        }
        fr.close();
    }
    ksprintf(ss, "</tbody>\n");
    // table end
    ksprintf(ss, "</table>\n");
    hmo.printExportButtons(ss, "aac_summary_table", "aac.summary");
    // body end
    ksprintf(ss, "</body>\n");
    // html end
    ksprintf(ss, "</html>\n");
    ps->l = 0;
    ksprintf(ps, "%s/%s/%s", outdir, rptd.snvdir.c_str(), rptd.aacfile.c_str());
    util::writestr(ss->s, ss->l, ps->s);
    ps->l = 0;
    ksprintf(ps, "%s/%s/%s", outdir, rptd.snvdir.c_str(), rptd.tsv_aacfile.c_str());
    util::writestr(ts->s, ts->l, ps->s);
}

void rpt_opt_t::gen_fsnv(){
    // parsing snvs
    std::vector<std::string> amps;
    ss->l = 0;
    std::string tmpnl;
    std::vector<std::string> vstr;
    for(auto& smp: samples){
        ss->l = 0;
        ksprintf(ss, "%s/%s/hapcnt", indir, smp.c_str());
        amps.clear();
        util::listDir(ss->s, amps);
        for(auto& amp: amps){
            for(size_t i = 0; i < nmnv.size(); ++i){
                ss->l = 0;
                ksprintf(ss, "%s/%s/hapcnt/%s/%s", indir, smp.c_str(), amp.c_str(), outmsnv[i].c_str());
                if(util::exists(ss->s)){
                    grp_mut_t* pgst = new grp_mut_t();
                    pgst->sample = smp;
                    pgst->amplicon = amp;
                    auto iter = amp2grp.find(amp);
                    if(iter == amp2grp.end()){
                        fprintf(stderr, "Warning: group of %s was not provided, default will be used\n", amp.c_str());
                        pgst->group = "default";
                    }else{
                        pgst->group = iter->second;
                    }
                    util::LineReader* lr = new util::LineReader(ss->s);
                    if(!lr->getline(tmpnl)){
                        fprintf(stderr, "Error parsing position line of %s\n", ss->s);
                    }else{
                        util::split(tmpnl, pgst->pos, "\t");
                    }
                    if(!lr->getline(tmpnl)){
                        fprintf(stderr, "Error parsing frequency line of %s\n", ss->s);
                    }else{
                        util::split(tmpnl, vstr, "\t");
                        for(auto& pval: vstr){
                            pgst->freq.push_back(atof(pval.c_str()));
                        }
                    }
                    if(pgst->pos.empty() || pgst->freq.empty() || (pgst->pos.size() != pgst->freq.size())){
                        fprintf(stderr, "Error parsing file of %s\n", ss->s);
                        delete pgst;
                    }else{
                        amnv[i].push_back(pgst);
                        if(pgst->pos.size() > maxp4focsnv) maxp4focsnv = pgst->pos.size();
                    }
                }
            }
        }
    }
    for(auto& mnvl: amnv){
        std::sort(mnvl.begin(), mnvl.end(), grp_mut_sorter());
    }
    // output
    ss->l = ts->l = 0;
    // html begin
    ksprintf(ss, "<!DOCTYPE html>\n");
    ksprintf(ss, "<html>\n");
    writeHTMLTableHeader(ss, "snv_foc_summary");
    // body begin
    ksprintf(ss, "<body>\n");
    // table begin
    ksprintf(ss, "<table id='snv_foc_summary_table'>\n");
    ksprintf(ss, "<thead>\n");
    ksprintf(ss, "<tr>\n");
    ksprintf(ss, "<th>Sample</th>\n"); ksprintf(ts, "Sample\t");
    ksprintf(ss, "<th>Amplicon</th>\n"); ksprintf(ts, "Amplicon\t");
    ksprintf(ss, "<th>Group</th>\n"); ksprintf(ts, "Group\t");
    ksprintf(ss, "<th>Muts</th>\n"); ksprintf(ts, "Muts\t");
    for(size_t i = 0; i < maxp4focsnv; ++i){
        ksprintf(ss, "<th>Val%ld</th>\n", i+1);
        ksprintf(ts, "Val%ld\t", i+1);
    }
    ts->s[ts->l-1] = '\n';
    ksprintf(ss, "</tr>\n");
    ksprintf(ss, "</thead>\n");
    ksprintf(ss, "<tbody>\n");
    grp_mut_t* ptsnv = new grp_mut_t();
    for(size_t i = 0; i < amnv.size(); ++i){
        for(size_t j = 0; j < amnv[i].size(); ++j){
            if(j == 0){
                amnv[i][j]->rec2th(ss, nmnv[i], maxp4focsnv);
                amnv[i][j]->rec2td(ss, nmnv[i], maxp4focsnv);
                amnv[i][j]->rec2tsvh(ts, nmnv[i], maxp4focsnv);
                amnv[i][j]->rec2tsvd(ts, nmnv[i], maxp4focsnv);
                ptsnv->sample = "avg_of_" + amnv[i][j]->group;
                ptsnv->amplicon = "avg_of_" + amnv[i][j]->group;
                ptsnv->group = amnv[i][j]->group;
                ptsnv->freq = amnv[i][j]->freq;
                ptsnv->nsmp = 1;
            }
            if(j > 0){
                if(amnv[i][j]->group == amnv[i][j-1]->group){
                    amnv[i][j]->rec2td(ss, nmnv[i], maxp4focsnv);
                    amnv[i][j]->rec2tsvd(ts, nmnv[i], maxp4focsnv);
                    for(size_t k = 0; k < amnv[i][j]->freq.size(); ++k){
                        ptsnv->freq[k] += amnv[i][j]->freq[k];
                    }
                    ++ptsnv->nsmp;
                }else{
                    for(auto& e: ptsnv->freq) e /= ptsnv->nsmp;
                    ptsnv->rec2td(ss, nmnv[i], maxp4focsnv);
                    ptsnv->rec2tsvd(ts, nmnv[i], maxp4focsnv);
                    ptsnv->nsmp = 0;
                    amnv[i][j]->rec2th(ss, nmnv[i], maxp4focsnv);
                    amnv[i][j]->rec2td(ss, nmnv[i], maxp4focsnv);
                    amnv[i][j]->rec2tsvh(ts, nmnv[i], maxp4focsnv);
                    amnv[i][j]->rec2tsvd(ts, nmnv[i], maxp4focsnv);
                    ptsnv->sample = "avg_of_" + amnv[i][j]->group;
                    ptsnv->amplicon = "avg_of_" + amnv[i][j]->group;
                    ptsnv->group = amnv[i][j]->group;
                    ptsnv->freq = amnv[i][j]->freq;
                    ptsnv->nsmp = 1;
                }
            }
        }
        if(ptsnv->nsmp){
            ptsnv->rec2td(ss, nmnv[i], maxp4focsnv);
            ptsnv->rec2tsvd(ts, nmnv[i], maxp4focsnv);
        }
    }
    delete ptsnv; ptsnv = NULL;
    ksprintf(ss, "</tbody>\n");
    // table end
    ksprintf(ss, "</table>\n");
    hmo.printExportButtons(ss, "snv_foc_summary_table", "snv.focus.summary");
    // body end
    ksprintf(ss, "</body>\n");
    // html end
    ksprintf(ss, "</html>\n");
    ps->l = 0;
    ksprintf(ps, "%s/%s/%s", outdir, rptd.snvdir.c_str(), rptd.snvfocf.c_str());
    util::writestr(ss->s, ss->l, ps->s);
    ps->l = 0;
    ksprintf(ps, "%s/%s/%s", outdir, rptd.snvdir.c_str(), rptd.tsv_snvfocf.c_str());
    util::writestr(ts->s, ts->l, ps->s);
}

void rpt_opt_t::gen_fmutf(){
    // parsing mutfreq
    std::vector<std::string> amps;
    ss->l = 0;
    std::string tmpnl;
    std::vector<std::string> vstr;
    for(auto& smp: samples){
        ss->l = 0;
        ksprintf(ss, "%s/%s/hapcnt", indir, smp.c_str());
        amps.clear();
        util::listDir(ss->s, amps);
        for(auto& amp: amps){
            for(size_t i = 0; i < nfreq.size(); ++i){
                ss->l = 0;
                ksprintf(ss, "%s/%s/hapcnt/%s/%s", indir, smp.c_str(), amp.c_str(), outmfreq[i].c_str());
                if(util::exists(ss->s)){
                    grp_mut_t* pgst = new grp_mut_t();
                    pgst->sample = smp;
                    pgst->amplicon = amp;
                    auto iter = amp2grp.find(amp);
                    if(iter == amp2grp.end()){
                        fprintf(stderr, "Warning: group of %s was not provided, default will be used\n", amp.c_str());
                        pgst->group = "default";
                    }else{
                        pgst->group = iter->second;
                    }
                    util::LineReader* lr = new util::LineReader(ss->s);
                    if(!lr->getline(tmpnl)){
                        fprintf(stderr, "Error parsing position line of %s\n", ss->s);
                    }else{
                        util::split(tmpnl, pgst->pos, "\t");
                    }
                    if(!lr->getline(tmpnl)){
                        fprintf(stderr, "Error parsing frequency line of %s\n", ss->s);
                    }else{
                        util::split(tmpnl, vstr, "\t");
                        for(auto& pval: vstr){
                            pgst->freq.push_back(atof(pval.c_str()));
                        }
                    }
                    if(pgst->pos.empty() || pgst->freq.empty() || (pgst->pos.size() != pgst->freq.size())){
                        fprintf(stderr, "Error parsing file of %s\n", ss->s);
                        delete pgst;
                    }else{
                        afreq[i].push_back(pgst);
                        if(pgst->pos.size() > maxp4focfreq) maxp4focfreq = pgst->pos.size();
                    }
                }
            }
        }
    }
    for(auto& freql: afreq){
        std::sort(freql.begin(), freql.end(), grp_mut_sorter());
    }
    // output
    ss->l = ts->l = 0;
    // html begin
    ksprintf(ss, "<!DOCTYPE html>\n");
    ksprintf(ss, "<html>\n");
    writeHTMLTableHeader(ss, "mutfreq_foc_summary");
    // body begin
    ksprintf(ss, "<body>\n");
    // table begin
    ksprintf(ss, "<table id='mutfreq_foc_summary_table'>\n");
    ksprintf(ss, "<thead>\n");
    ksprintf(ss, "<tr>\n");
    ksprintf(ss, "<th>Sample</th>\n"); ksprintf(ts, "Sample\t");
    ksprintf(ss, "<th>Amplicon</th>\n"); ksprintf(ts, "Amplicon\t");
    ksprintf(ss, "<th>Group</th>\n"); ksprintf(ts, "Group\t");
    ksprintf(ss, "<th>Type</th>\n"); ksprintf(ts, "Type\t");
    for(size_t i = 0; i < maxp4focfreq; ++i){
        ksprintf(ss, "<th>Val%ld</th>\n", i+1);
        ksprintf(ts, "Val%ld\t", i+1);
    }
    ts->s[ts->l-1] = '\n';
    ksprintf(ss, "</tr>\n");
    ksprintf(ss, "</thead>\n");
    ksprintf(ss, "<tbody>\n");
    grp_mut_t* ptsnv = new grp_mut_t();
    for(size_t i = 0; i < afreq.size(); ++i){
        for(size_t j = 0; j < afreq[i].size(); ++j){
            if(j == 0){
                afreq[i][j]->rec2th(ss, nfreq[i], maxp4focfreq);
                afreq[i][j]->rec2td(ss, nfreq[i], maxp4focfreq);
                afreq[i][j]->rec2tsvh(ts, nfreq[i], maxp4focfreq);
                afreq[i][j]->rec2tsvd(ts, nfreq[i], maxp4focfreq);
                ptsnv->sample = "avg_of_" + afreq[i][j]->group;
                ptsnv->amplicon = "avg_of_" + afreq[i][j]->group;
                ptsnv->group = afreq[i][j]->group;
                ptsnv->freq = afreq[i][j]->freq;
                ptsnv->nsmp = 1;
            }
            if(j > 0){
                if(afreq[i][j]->group == afreq[i][j-1]->group){
                    afreq[i][j]->rec2td(ss, nfreq[i], maxp4focfreq);
                    afreq[i][j]->rec2tsvd(ts, nfreq[i], maxp4focfreq);
                    for(size_t k = 0; k < afreq[i][j]->freq.size(); ++k){
                        ptsnv->freq[k] += afreq[i][j]->freq[k];
                    }
                    ++ptsnv->nsmp;
                }else{
                    for(auto& e: ptsnv->freq) e /= ptsnv->nsmp;
                    ptsnv->rec2td(ss, nfreq[i], maxp4focfreq);
                    ptsnv->rec2tsvd(ts, nfreq[i], maxp4focfreq);
                    ptsnv->nsmp = 0;
                    afreq[i][j]->rec2th(ss, nfreq[i], maxp4focfreq);
                    afreq[i][j]->rec2td(ss, nfreq[i], maxp4focfreq);
                    afreq[i][j]->rec2tsvh(ts, nfreq[i], maxp4focfreq);
                    afreq[i][j]->rec2tsvd(ts, nfreq[i], maxp4focfreq);
                    ptsnv->sample = "avg_of_" + afreq[i][j]->group;
                    ptsnv->amplicon = "avg_of_" + afreq[i][j]->group;
                    ptsnv->group = afreq[i][j]->group;
                    ptsnv->freq = afreq[i][j]->freq;
                    ptsnv->nsmp = 1;
                }
            }
        }
        if(ptsnv->nsmp){
            ptsnv->rec2td(ss, nfreq[i], maxp4focfreq);
            ptsnv->rec2tsvd(ts, nfreq[i], maxp4focfreq);
        }
    }
    delete ptsnv; ptsnv = NULL;
    ksprintf(ss, "</tbody>\n");
    // table end
    ksprintf(ss, "</table>\n");
    hmo.printExportButtons(ss, "mutfreq_foc_summary_table", "mutfreq.focus.summary");
    // body end
    ksprintf(ss, "</body>\n");
    // html end
    ksprintf(ss, "</html>\n");
    ps->l = 0;
    ksprintf(ps, "%s/%s/%s", outdir, rptd.snvdir.c_str(), rptd.mutfreq.c_str());
    util::writestr(ss->s, ss->l, ps->s);
    ps->l = 0;
    ksprintf(ps, "%s/%s/%s", outdir, rptd.snvdir.c_str(), rptd.tsv_mutfreq.c_str());
    util::writestr(ts->s, ts->l, ps->s);
}

void rpt_opt_t::gen_biasf(){
    // parsing mutfreq
    std::vector<std::string> amps;
    ss->l = 0;
    std::string tmpnl;
    std::vector<std::string> vstr, tvstr;
    for(auto& smp: samples){
        ss->l = 0;
        ksprintf(ss, "%s/%s/hapcnt", indir, smp.c_str());
        amps.clear();
        util::listDir(ss->s, amps);
        for(auto& amp: amps){
            for(size_t i = 0; i < nbias.size(); ++i){
                ss->l = 0;
                ksprintf(ss, "%s/%s/hapcnt/%s/%s", indir, smp.c_str(), amp.c_str(), outbias[i].c_str());
                if(util::exists(ss->s)){
                    grp_mut_t* pgst = new grp_mut_t();
                    pgst->sample = smp;
                    pgst->amplicon = amp;
                    auto iter = amp2grp.find(amp);
                    if(iter == amp2grp.end()){
                        fprintf(stderr, "Warning: group of %s was not provided, default will be used\n", amp.c_str());
                        pgst->group = "default";
                    }else{
                        pgst->group = iter->second;
                    }
                    util::LineReader* lr = new util::LineReader(ss->s);
                    if(!lr->getline(tmpnl)){
                        fprintf(stderr, "Error parsing position line of %s\n", ss->s);
                    }else{
                        util::split(tmpnl, tvstr, "\t");
                        pgst->pos.clear();
                        for(size_t bb = 1; bb < tvstr.size(); ++bb) pgst->pos.push_back(tvstr[bb]);
                    }
                    pgst->freq.clear();
                    while(lr->getline(tmpnl)){
                        util::split(tmpnl, tvstr, "\t");
                        for(size_t bb = 1; bb < tvstr.size(); ++bb) pgst->freq.push_back(atof(tvstr[bb].c_str()));
                        pgst->row.push_back(tvstr[0]);
                    }
                    if(pgst->pos.empty() || pgst->freq.empty() || pgst->row.empty() || (pgst->pos.size() * pgst->row.size() != pgst->freq.size())){
                        fprintf(stderr, "Error parsing file of %s\n", ss->s);
                        delete pgst;
                    }else{
                        abias[i].push_back(pgst);
                        if(pgst->pos.size() > maxp4bias) maxp4bias = pgst->pos.size();
                    }
                }
            }
        }
    }
    for(auto& bsl: abias){
        std::sort(bsl.begin(), bsl.end(), grp_mut_sorter());
    }
    // output
    ss->l = ts->l = 0;
    // html begin
    ksprintf(ss, "<!DOCTYPE html>\n");
    ksprintf(ss, "<html>\n");
    writeHTMLTableHeader(ss, "mut_bias_summary");
    // body begin
    ksprintf(ss, "<body>\n");
    // table begin
    ksprintf(ss, "<table id='mut_bias_summary_table'>\n");
    ksprintf(ss, "<thead>\n");
    ksprintf(ss, "<tr>\n");
    ksprintf(ss, "<th>Sample</th>\n"); ksprintf(ts, "Sample\t");
    ksprintf(ss, "<th>Amplicon</th>\n"); ksprintf(ts, "Amplicon\t");
    ksprintf(ss, "<th>Group</th>\n"); ksprintf(ts, "Group\t");
    ksprintf(ss, "<th>Type</th>\n"); ksprintf(ts, "Type\t");
    ksprintf(ss, "<th>Ref\\Alt</th>\n"); ksprintf(ts, "Ref\\Alt\t");
    for(size_t i = 0; i < maxp4bias; ++i){
        ksprintf(ss, "<th>Val%ld</th>\n", i+1);
        ksprintf(ts, "Val%ld\t", i+1);
    }
    ts->s[ts->l-1] = '\n';
    ksprintf(ss, "</tr>\n");
    ksprintf(ss, "</thead>\n");
    ksprintf(ss, "<tbody>\n");
    grp_mut_t* ptsnv = new grp_mut_t();
    for(size_t i = 0; i < abias.size(); ++i){
        for(size_t j = 0; j < abias[i].size(); ++j){
            if(j == 0){
                abias[i][j]->recm2th(ss, nbias[i], "Ref2Alt", maxp4bias);
                abias[i][j]->recm2td(ss, nbias[i], maxp4bias);
                abias[i][j]->recm2tsvh(ts, nbias[i], "Ref2Alt", maxp4bias);
                abias[i][j]->recm2tsvd(ts, nbias[i], maxp4bias);
                ptsnv->sample = "avg_of_" + abias[i][j]->group;
                ptsnv->amplicon = "avg_of_" + abias[i][j]->group;
                ptsnv->group = abias[i][j]->group;
                ptsnv->freq = abias[i][j]->freq;
                ptsnv->row = abias[i][j]->row;
                ptsnv->nsmp = 1;
            }
            if(j > 0){
                if(abias[i][j]->group == abias[i][j-1]->group){
                    abias[i][j]->recm2td(ss, nbias[i], maxp4bias);
                    abias[i][j]->recm2tsvd(ts, nbias[i], maxp4bias);
                    for(size_t k = 0; k < abias[i][j]->freq.size(); ++k){
                        ptsnv->freq[k] += abias[i][j]->freq[k];
                    }
                    ++ptsnv->nsmp;
                }else{
                    for(auto& e: ptsnv->freq) e /= ptsnv->nsmp;
                    ptsnv->recm2td(ss, nbias[i], maxp4bias);
                    ptsnv->recm2tsvd(ts, nbias[i], maxp4bias);
                    ptsnv->nsmp = 0;
                    abias[i][j]->recm2th(ss, nbias[i], "Ref2Alt", maxp4bias);
                    abias[i][j]->recm2td(ss, nbias[i], maxp4bias);
                    abias[i][j]->recm2tsvh(ts, nbias[i], "Ref2Alt", maxp4bias);
                    abias[i][j]->recm2tsvd(ts, nbias[i], maxp4bias);
                    ptsnv->sample = "avg_of_" + abias[i][j]->group;
                    ptsnv->amplicon = "avg_of_" + abias[i][j]->group;
                    ptsnv->group = abias[i][j]->group;
                    ptsnv->freq = abias[i][j]->freq;
                    ptsnv->row =abias[i][j]->row;
                    ptsnv->nsmp = 1;
                }
            }
        }
        if(ptsnv->nsmp){
            ptsnv->recm2td(ss, nbias[i], maxp4bias);
            ptsnv->recm2tsvd(ts, nbias[i], maxp4bias);
        }
    }
    delete ptsnv; ptsnv = NULL;
    ksprintf(ss, "</tbody>\n");
    // table end
    ksprintf(ss, "</table>\n");
    hmo.printExportButtons(ss, "mut_bias_summary_table", "mutbias.summary");
    // body end
    ksprintf(ss, "</body>\n");
    // html end
    ksprintf(ss, "</html>\n");
    ps->l = 0;
    ksprintf(ps, "%s/%s/%s", outdir, rptd.snvdir.c_str(), rptd.sbiasfile.c_str());
    util::writestr(ss->s, ss->l, ps->s);
    ps->l = 0;
    ksprintf(ps, "%s/%s/%s", outdir, rptd.snvdir.c_str(), rptd.tsv_sbiasfile.c_str());
    util::writestr(ts->s, ts->l, ps->s);
}

void rpt_opt_t::gen_fmndist(){
    // parsing mutcnt 8 reads
    std::vector<std::string> amps;
    ss->l = 0;
    std::string tmpnl;
    std::vector<std::string> vstr;
    for(auto& smp: samples){
        ss->l = 0;
        ksprintf(ss, "%s/%s/hapcnt", indir, smp.c_str());
        amps.clear();
        util::listDir(ss->s, amps);
        for(auto& amp: amps){
            for(size_t i = 0; i < nrcnt.size(); ++i){
                ss->l = 0;
                ksprintf(ss, "%s/%s/hapcnt/%s/%s", indir, smp.c_str(), amp.c_str(), outrcnt[i].c_str());
                if(util::exists(ss->s)){
                    grp_mut_t* pgst = new grp_mut_t();
                    pgst->sample = smp;
                    pgst->amplicon = amp;
                    auto iter = amp2grp.find(amp);
                    if(iter == amp2grp.end()){
                        fprintf(stderr, "Warning: group of %s was not provided, default will be used\n", amp.c_str());
                        pgst->group = "default";
                    }else{
                        pgst->group = iter->second;
                    }
                    util::LineReader* lr = new util::LineReader(ss->s);
                    if(!lr->getline(tmpnl)){
                        fprintf(stderr, "Error parsing position line of %s\n", ss->s);
                    }else{
                        util::split(tmpnl, pgst->pos, "\t");
                    }
                    if(!lr->getline(tmpnl)){
                        fprintf(stderr, "Error parsing dist count line of %s\n", ss->s);
                    }else{
                        util::split(tmpnl, vstr, "\t");
                        for(auto& pval: vstr){
                            pgst->freq.push_back(atof(pval.c_str()));
                        }
                    }
                    if(pgst->pos.empty() || pgst->freq.empty() || (pgst->pos.size() != pgst->freq.size())){
                        fprintf(stderr, "Error parsing file of %s\n", ss->s);
                        delete pgst;
                    }else{
                        // resize to maxmdistl
                        if(maxmdistl < pgst->pos.size()){
                            for(size_t i = maxmdistl+1; i < pgst->freq.size(); ++i){
                                pgst->freq[maxmdistl] += pgst->freq[i];
                            }
                            pgst->freq.resize(maxmdistl+1);
                            pgst->pos.resize(maxmdistl+1);
                        }
                        arnct[i].push_back(pgst);
                        if(pgst->pos.size() > maxp4focrcnt) maxp4focrcnt = pgst->pos.size();
                    }
                }
            }
        }
    }
    for(auto& ncntl: arnct){
        std::sort(ncntl.begin(), ncntl.end(), grp_mut_sorter());
    }
    // output
    ss->l = ts->l = 0;
    // html begin
    ksprintf(ss, "<!DOCTYPE html>\n");
    ksprintf(ss, "<html>\n");
    writeHTMLTableHeader(ss, "mutdist_foc_summary");
    // body begin
    ksprintf(ss, "<body>\n");
    // table begin
    ksprintf(ss, "<table id='mutdist_foc_summary_table'>\n");
    ksprintf(ss, "<thead>\n");
    ksprintf(ss, "<tr>\n");
    ksprintf(ss, "<th>Sample</th>\n"); ksprintf(ts, "Sample\t");
    ksprintf(ss, "<th>Amplicon</th>\n"); ksprintf(ts, "Amplicon\t");
    ksprintf(ss, "<th>Group</th>\n"); ksprintf(ts, "Group\t");
    ksprintf(ss, "<th>Type</th>\n"); ksprintf(ts, "Type\t");
    for(size_t i = 0; i < maxp4focrcnt; ++i){
        ksprintf(ss, "<th>Val%ld</th>\n", i+1);
        ksprintf(ts, "Val%ld\t", i+1);
    }
    ts->s[ts->l-1] = '\n';
    ksprintf(ss, "</tr>\n");
    ksprintf(ss, "</thead>\n");
    ksprintf(ss, "<tbody>\n");
    grp_mut_t* ptsnv = new grp_mut_t();
    for(size_t i = 0; i < arnct.size(); ++i){
        for(size_t j = 0; j < arnct[i].size(); ++j){
            if(j == 0){
                arnct[i][j]->rec2th(ss, nrcnt[i], maxp4focrcnt);
                arnct[i][j]->rec2td(ss, nrcnt[i], maxp4focrcnt);
                arnct[i][j]->rec2tsvh(ts, nrcnt[i], maxp4focrcnt);
                arnct[i][j]->rec2tsvd(ts, nrcnt[i], maxp4focrcnt);
                ptsnv->sample = "avg_of_" + arnct[i][j]->group;
                ptsnv->amplicon = "avg_of_" + arnct[i][j]->group;
                ptsnv->group = arnct[i][j]->group;
                ptsnv->freq = arnct[i][j]->freq;
                ptsnv->nsmp = 1;
            }
            if(j > 0){
                if(arnct[i][j]->group == arnct[i][j-1]->group){
                    arnct[i][j]->rec2td(ss, nrcnt[i], maxp4focrcnt);
                    arnct[i][j]->rec2tsvd(ts, nrcnt[i], maxp4focrcnt);
                    for(size_t k = 0; k < arnct[i][j]->freq.size(); ++k){
                        ptsnv->freq[k] += arnct[i][j]->freq[k];
                    }
                    ++ptsnv->nsmp;
                }else{
                    for(auto& e: ptsnv->freq) e /= ptsnv->nsmp;
                    ptsnv->rec2td(ss, nrcnt[i], maxp4focrcnt);
                    ptsnv->rec2tsvd(ts, nrcnt[i], maxp4focrcnt);
                    ptsnv->nsmp = 0;
                    arnct[i][j]->rec2th(ss, nrcnt[i], maxp4focrcnt);
                    arnct[i][j]->rec2td(ss, nrcnt[i], maxp4focrcnt);
                    arnct[i][j]->rec2tsvh(ts, nrcnt[i], maxp4focrcnt);
                    arnct[i][j]->rec2tsvd(ts, nrcnt[i], maxp4focrcnt);
                    ptsnv->sample = "avg_of_" + arnct[i][j]->group;
                    ptsnv->amplicon = "avg_of_" + arnct[i][j]->group;
                    ptsnv->group = arnct[i][j]->group;
                    ptsnv->freq = arnct[i][j]->freq;
                    ptsnv->nsmp = 1;
                }
            }
        }
        if(ptsnv->nsmp){
            ptsnv->rec2td(ss, nrcnt[i], maxp4focrcnt);
            ptsnv->rec2tsvd(ts, nrcnt[i], maxp4focrcnt);
        }
    }
    delete ptsnv; ptsnv = NULL;
    ksprintf(ss, "</tbody>\n");
    // table end
    ksprintf(ss, "</table>\n");
    hmo.printExportButtons(ss, "mutdist_foc_summary_table", "mutdist.focus.summary");
    // body end
    ksprintf(ss, "</body>\n");
    // html end
    ksprintf(ss, "</html>\n");
    ps->l = 0;
    ksprintf(ps, "%s/%s/%s", outdir, rptd.snvdir.c_str(), rptd.mutn8rdist.c_str());
    util::writestr(ss->s, ss->l, ps->s);
    ps->l = 0;
    ksprintf(ps, "%s/%s/%s", outdir, rptd.snvdir.c_str(), rptd.tsv_mutn8rdist.c_str());
    util::writestr(ts->s, ts->l, ps->s);
}

void rpt_opt_t::gen_sumqcf(){
    std::vector<std::string> fsFile;
    for(auto& s: samples){
        ss->l = 0;
        if(diffamps) ksprintf(ss, "%s/%s/%s", indir, s.c_str(), GETOOLS_SPLITCMBTSV);
        else ksprintf(ss, "%s/%s/%s", indir, s.c_str(), GETOOLS_SPLITTSV);
        fsFile.push_back(std::string(ss->s));
    }
    ss->l = ts->l = 0;
    // html begin
    ksprintf(ss, "<!DOCTYPE html>\n");
    ksprintf(ss, "<html>\n");
    writeHTMLTableHeader(ss, "qc_summary");
    // body begin
    ksprintf(ss, "<body>\n");
    // table begin
    ksprintf(ss, "<table id='summary_table'>\n");
    ksprintf(ss, "<thead>\n");
    ksprintf(ss, "<tr>\n");
    ksprintf(ss, "<th>Samples</th>\n"); ksprintf(ts, "Samples\t");
    std::ifstream fr(fsFile[0]);
    std::string tmpline;
    std::getline(fr, tmpline);
    std::vector<std::string> vstr;
    util::split(tmpline, vstr, "\t");
    for(auto& e: vstr) ksprintf(ts, "%s\t", e.c_str());
    ts->s[ts->l-1] = '\n';
    tmpline = util::replace(tmpline, "<", "&lt;");
    tmpline = util::replace(tmpline, ">", "&gt;");
    util::split(tmpline, vstr, "\t");
    for(auto& e: vstr) ksprintf(ss, "<th>%s</th>\n", e.c_str());
    ksprintf(ss, "</tr>\n");
    ksprintf(ss, "</thead>\n");
    ksprintf(ss, "<tbody>\n");
    fr.close();
    for(size_t i = 0; i < fsFile.size(); ++i){
        fr.open(fsFile[i].c_str());
        std::getline(fr, tmpline);
        std::getline(fr, tmpline);
        ksprintf(ss, "<tr>\n");
        ksprintf(ss, "<td>%s</td>\n", samples[i].c_str()); ksprintf(ts, "%s\t", samples[i].c_str());
        util::split(tmpline, vstr, "\t");
        for(auto& e: vstr){
            ksprintf(ss, "<td>%s</td>\n", e.c_str());
            ksprintf(ts, "%s\t", e.c_str());
        }
        ts->s[ts->l-1] = '\n';
        fr.close();
    }
    ksprintf(ss, "</tbody>\n");
    // table end
    ksprintf(ss, "</table>\n");
    hmo.printExportButtons(ss, "summary_table", "qc.summary");
    // body end
    ksprintf(ss, "</body>\n");
    // html end
    ksprintf(ss, "</html>\n");
    ps->l = 0;
    ksprintf(ps, "%s/%s/%s", outdir, rptd.qcdir.c_str(), rptd.qcfile.c_str());
    util::writestr(ss->s, ss->l, ps->s);
    ps->l = 0;
    ksprintf(ps, "%s/%s/%s", outdir, rptd.qcdir.c_str(), rptd.tsv_qcfile.c_str());
    util::writestr(ts->s, ts->l, ps->s);
}

void rpt_opt_t::gen_sumedf(){
    std::vector<std::string> fsFile;
    for(auto& s: samples){
        ss->l = 0;
        ksprintf(ss, "%s/%s/%s", indir, s.c_str(), GETOOLS_SUBLIBTSV);
        fsFile.push_back(std::string(ss->s));
    }
    ss->l = ts->l = 0;
    // html begin
    ksprintf(ss, "<!DOCTYPE html>\n");
    ksprintf(ss, "<html>\n");
    writeHTMLTableHeader(ss, "ed_summary");
    // body begin
    ksprintf(ss, "<body>\n");
    // table begin
    ksprintf(ss, "<table id='summary_table'>\n");
    ksprintf(ss, "<thead>\n");
    ksprintf(ss, "<tr>\n");
    ksprintf(ss, "<th>Samples</th>\n"); ksprintf(ts, "Samples\t");
    std::ifstream fr(fsFile[0]);
    std::string tmpline;
    std::getline(fr, tmpline);
    std::vector<std::string> vstr;
    util::split(tmpline, vstr, "\t");
    for(auto& e: vstr) ksprintf(ts, "%s\t", e.c_str());
    ts->s[ts->l-1] = '\n';
    tmpline = util::replace(tmpline, "<", "&lt;");
    tmpline = util::replace(tmpline, ">", "&gt;");
    util::split(tmpline, vstr, "\t");
    for(auto& e: vstr) ksprintf(ss, "<th>%s</th>\n", e.c_str());
    ksprintf(ss, "</tr>\n");
    ksprintf(ss, "</thead>\n");
    ksprintf(ss, "<tbody>\n");
    fr.close();
    for(size_t i = 0; i < fsFile.size(); ++i){
        fr.open(fsFile[i].c_str());
        std::getline(fr, tmpline);
        while(std::getline(fr, tmpline)){
            ksprintf(ss, "<tr>\n");
            ksprintf(ss, "<td>%s</td>\n", samples[i].c_str()); ksprintf(ts, "%s\t", samples[i].c_str());
            util::split(tmpline, vstr, "\t");
            for(auto& e: vstr){
                ksprintf(ss, "<td>%s</td>\n", e.c_str());
                ksprintf(ts, "%s\t", e.c_str());
            }
            ts->s[ts->l-1] = '\n';
        }
        fr.close();
    }
    ksprintf(ss, "</tbody>\n");
    // table end
    ksprintf(ss, "</table>\n");
    hmo.printExportButtons(ss, "summary_table", "ed.summary");
    // body end
    ksprintf(ss, "</body>\n");
    // html end
    ksprintf(ss, "</html>\n");
    ps->l = 0;
    ksprintf(ps, "%s/%s/%s", outdir, rptd.qcdir.c_str(), rptd.edfile.c_str());
    util::writestr(ss->s, ss->l, ps->s);
    ps->l = 0;
    ksprintf(ps, "%s/%s/%s", outdir, rptd.qcdir.c_str(), rptd.tsv_edfile.c_str());
    util::writestr(ts->s, ts->l, ps->s);
}

void rpt_opt_t::gen_qcf(){
    for(auto& s : samples){
        ss->l = 0;
        ps->l = 0;
        ksprintf(ps, "%s/%s/%s", indir, s.c_str(), GETOOLS_SUBLIBTSV);
        tsv2HTMLTable(ps->s, ss);
        ps->l = 0;
        ksprintf(ps, "%s/%s/%s.html", outdir, rptd.sumdir.c_str(), s.c_str());
        util::writestr(ss->s, ss->l, ps->s);
    }
}

void rpt_opt_t::cp_rpf(){
    if(gen4web) return;
    ss->l = 0;
    for(auto& s: samples){
        ss->l = 0;
        ksprintf(ss, "cp -rf %s/%s/html %s/%s/%s/;", indir, s.c_str(), outdir, rptd.rptdir.c_str(), s.c_str());
        ksprintf(ss, "cp -f %s/%s/report.html %s/%s/%s/index.html;", indir, s.c_str(), outdir, rptd.rptdir.c_str(), s.c_str());
        ts->l = 0;
        ksprintf(ts, "%s/%s/%s", indir, s.c_str(), SCANA_AMP_ALLELE_IN_SC_TSV);
        if(util::exists(ts->s)){
            ksprintf(ss, "cp -f %s %s/%s/%s/%s;", ts->s, outdir, rptd.rptdir.c_str(), s.c_str(), SCANA_AMP_ALLELE_IN_SC_TSV);
        }
        ts->l = 0;
        ksprintf(ts, "%s/%s/%s", indir, s.c_str(), SCANA_AMP_CO_EDIT_TSV);
        if(util::exists(ts->s)){
            ksprintf(ss, "cp -f %s %s/%s/%s/%s;", ts->s, outdir, rptd.rptdir.c_str(), s.c_str(), SCANA_AMP_CO_EDIT_TSV);
        }
        ts->l = 0;
        ksprintf(ts, indir, s.c_str(), SCANA_CELL_FIND_JSON);
        if(util::exists(ts->s)){
            ksprintf(ss, "cp -f %s %s/%s/%s/%s;", ts->s, outdir, rptd.rptdir.c_str(), s.c_str(), SCANA_CELL_FIND_JSON);
        }
        system(ss->s);
    }
}

void rpt_opt_t::reportFiles(){
    if(consamp) gen_consf();
    if(haprpt){
        gen_hap();
        gen_snv();
        gen_aac();
        if(extsnv.size()) gen_fsnv();
        gen_fmutf();
        gen_fmndist();
        gen_biasf();
    }
    gen_sumqcf();
    gen_sumedf();
    gen_qcf();
    cp_rpf();
}

void rpt_opt_t::report(){
    init();
    reportHTML();
    reportFiles();
}

void rpt_opt_t::writeHTMLTableHeader(kstring_t* s, std::string title){
    ksprintf(s, "<head>\n");
    ksprintf(ks, "<meta http-equiv=\"Content-Type\" content=\"text/html; charset=UTF-8\">\n");
    ksprintf(s, "<title>%s at %s</title>\n", title.c_str(), util::currentTime().c_str());
    hmo.printJsAndCSSRes(s);
    ksprintf(s, "<style>\n");
    ksprintf(s, "table {\n");
    ksprintf(s, "font-family: arial, sans-serif;\n");
    ksprintf(s, "border-collapse: collapse;\n");
    ksprintf(s, "width: 100%%;\n");
    ksprintf(s, "}\n");
    ksprintf(s, "td, th {\n");
    ksprintf(s, "border: 1px solid #dddddd;\n");
    ksprintf(s, "text-align: left;\n");
    ksprintf(s, "padding: 8px;\n");
    ksprintf(s, "}\n");
    ksprintf(s, "tr:nth-child(even) {\n");
    ksprintf(s, "background-color: #dddddd;\n");
    ksprintf(s, "}\n");
    ksprintf(s, "</style>\n");
    ksprintf(s, "</head>\n");
}

void rpt_opt_t::tsv2HTMLTable(char* f, kstring_t* s){
    ksprintf(s, "<!DOCTYPE html>\n");
    ksprintf(s, "<html>\n");
    writeHTMLTableHeader(s, "sub_lib_info");
    ksprintf(s, "<body>\n");
    ksprintf(ss, "<table id='summary_table'>\n");
    ksprintf(ss, "<thead>\n");
    ksprintf(s, "<tr>\n");
    util::LineReader lr(f);
    std::string tmpline;
    std::vector<std::string> vstr;
    lr.getline(tmpline);
    tmpline = util::replace(tmpline, "<", "&lt;");
    tmpline = util::replace(tmpline, ">", "&gt;");
    util::split(tmpline, vstr, "\t");
    for(auto& e: vstr) ksprintf(s, "<th>%s</th>\n", e.c_str());
    ksprintf(s, "</tr>\n");
    ksprintf(ss, "</thead>\n");
    ksprintf(ss, "<tbody>\n");
    while(lr.getline(tmpline)){
        ksprintf(s, "<tr>\n");
        util::split(tmpline, vstr, "\t");
        for(auto& e: vstr) ksprintf(s, "<th>%s</th>\n", e.c_str());
        ksprintf(s, "</tr>\n");
    }
    ksprintf(ss, "</tbody>\n");
    ksprintf(s, "</table>\n");
    hmo.printExportButtons(ss, "summary_table", "sub.lib.info");
    ksprintf(s, "</body>\n");
    ksprintf(s, "</html>\n");
}

void rpt_usage(rpt_opt_t* opt, char* arg0){
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: %s %s [options]\n", PACKAGE_NAME, arg0);
    fprintf(stderr, "\nInput/output options:\n");
    fprintf(stderr, "         -i STR  getools analysis samples result directory\n");
    fprintf(stderr, "         -o STR  output directory of report files\n");
    fprintf(stderr, "         -s FILE sample list to report(one name per row)\n");
    fprintf(stderr, "         -l STR  sample names to report(seperated by comma)\n");
    fprintf(stderr, "         -t STR  report title\n");
    fprintf(stderr, "         -j STR  json configure file of js/css\n");
    fprintf(stderr, "         -r INT  row number of table\n");
    fprintf(stderr, "         -d      samples are analysised by different config files if set\n");
    fprintf(stderr, "\nEdit consistence/pattern options:\n");
    fprintf(stderr, "         -a STR  amplicon to do consistence analysis\n");
    fprintf(stderr, "         -n INT  number of top shared variants to compare [%d]\n", opt->na4cons);
    fprintf(stderr, "         -w INT  html topN shared variants graph width [%d]\n", opt->topnW);
    fprintf(stderr, "         -z INT  dimension reduction method(%d:%s,%d:%s,%d:%s) [%d]\n",
                                        DR_METHOD_NONE, DR_METHOD2STR_ARR[DR_METHOD_NONE],
                                        DR_METHOD_PCA, DR_METHOD2STR_ARR[DR_METHOD_PCA],
                                        DR_METHOD_UMAP, DR_METHOD2STR_ARR[DR_METHOD_UMAP],
                                        opt->drmethod);
    fprintf(stderr, "         -u      output all unique variants instead of topN if set\n");
    fprintf(stderr, "         -x      show sample names in dimension reduction graph if set\n");
    fprintf(stderr, "         -y      calculate genotype of each amplicon if set\n");
    fprintf(stderr, "\nHaplotype SNV analysis options:\n");
    fprintf(stderr, "         -g STR  amplicon group info list\n");
    fprintf(stderr, "         -e STR  extra snv to focus on, format(REF:ALT,REF:ALT2,...)\n");
    fprintf(stderr, "         -m INT  max mut count in one read to enrich into dist [%ld]\n", opt->maxmdistl);
    fprintf(stderr, "\n");
}

int rpt_main(int argc, char** argv){
    rpt_opt_t opt;
    if(argc == 1){
        rpt_usage(&opt, argv[0]);
        return 0;
    }
    int c = -1;
    while((c = getopt(argc, argv, "i:o:s:l:t:j:r:a:n:w:z:g:e:m:duxyh")) >= 0){
        switch(c){
            case 'i': opt.indir = optarg; break;
            case 'o': opt.outdir = optarg; break;
            case 's': opt.splist = optarg; break;
            case 'l': opt.spstr = optarg; break;
            case 't': opt.title = optarg; break;
            case 'j': opt.jscfg = optarg; break;
            case 'r': opt.nrow = atoi(optarg); break;
            case 'a': opt.consamp = optarg; break;
            case 'n': opt.na4cons = atoi(optarg); break;
            case 'w': opt.topnW = atoi(optarg); break;
            case 'z': opt.drmethod = atoi(optarg); break;
            case 'g': opt.ampgrpf = optarg; break;
            case 'e': opt.extsnv = optarg; break;
            case 'm': opt.maxmdistl = atoi(optarg); break;
            case 'd': opt.diffamps = true; break;
            case 'u': opt.uniqa = true; break;
            case 'x': opt.showt = true; break;
            case 'y': opt.calg = true; break;
            case 'h': rpt_usage(&opt, argv[0]); return 0; break;
            default: break;
        }
    }
    if(opt.valid()){
        opt.report();
    }
    return 0;
}
