#include "pscmp.h"


pscmp_opt_t::pscmp_opt_t(){
}

pscmp_opt_t::~pscmp_opt_t(){
    if(ps){
        if(ps->s) free(ps->s);
        free(ps); ps = NULL;
    }
    if(ts){
        if(ts->s) free(ts->s);
        free(ts); ts = NULL;
    }
    if(vnuc){
        for(int i = 0; i < plen; ++i) free(vnuc[i]);
        free(vnuc); vnuc = NULL;
    }
}


bool pscmp_opt_t::valid(){
    if(!util::exists(indir)){
        fprintf(stderr, "input directory does not exist\n");
        return false;
    }
    if(splist.empty()){
        fprintf(stderr, "sample list must be provided\n");
        return false;
    }
    if(!util::exists(splist)){
        fprintf(stderr, "sample list does not exist\n");
        return false;
    }
    return true;
}

void pscmp_opt_t::init(){
    if(!util::exists(outdir)) util::makedir(outdir);
    outhtml = util::joinpath(outdir, outhtml);
    outpsqc = util::joinpath(outdir, n_outpsqc);
    outpscnt = util::joinpath(outdir, n_outpscnt);
    outlogps = util::joinpath(outdir, n_outlogps);
    util::makeListFromFileByLine(splist, samples);
    ps = (kstring_t*)calloc(1, sizeof(kstring_t));
    ts = (kstring_t*)calloc(1, sizeof(kstring_t));
    hmo.init(hrjsn);
}

void pscmp_opt_t::parse_qc(){
    std::string line;
    std::vector<std::string> vstr;
    std::set<std::string> actgs;
    contig_pscnt_qc_t* pqc = NULL;
    for(auto& s: samples){
        ps->l = 0;
        ksprintf(ps, "%s/%s/%s", indir.c_str(), s.c_str(), PSCNT_RDC_TSV);
        if(!(util::exists(ps->s))){
            fprintf(stderr, "QC file not found:%s\n", ps->s);
            exit(EXIT_FAILURE);
        }
        util::LineReader lr(ps->s);
        lr.getline(line);
        while(lr.getline(line)){
            util::split(line, vstr, "\t");
            if(actgs.find(vstr[0]) != actgs.end()){
                fprintf(stderr, "No duplicated contig names[%s] allowed in all samples\n", vstr[0].c_str());
                exit(EXIT_FAILURE);
            }else{
                actgs.insert(vstr[0]);
            }
            contigs.push_back(vstr[0]);
            auto iter = map4ctgqc.find(vstr[0]);
            if(iter == map4ctgqc.end()){
                pqc = new contig_pscnt_qc_t();
                map4ctgqc[vstr[0]] = pqc;
            }else{
                pqc = iter->second;
            }
            pqc->total_reads = atol(vstr[1].c_str());
            pqc->with_anno_pat = atol(vstr[2].c_str());
            pqc->with_other_pat = atol(vstr[3].c_str());
            pqc->no_pat_found = atol(vstr[4].c_str());
        }
    }
    if(cnorm.size()){
        auto iter = map4ctgqc.find(cnorm);
        nttr = iter->second->with_anno_pat + iter->second->with_other_pat;
    }
}

void pscmp_opt_t::parse_pscnt(){
    std::string line, lctg;
    pscnt_info_t* psc = NULL;
    contig_pscnt_t* pct = NULL;
    contig_pscnt_t* pnct = NULL;
    std::vector<std::string> vstr;
    for(auto& s: samples){
        ps->l = 0;
        ksprintf(ps, "%s/%s/%s", indir.c_str(), s.c_str(), PSCNT_PSC_TSV);
        if(!(util::exists(ps->s))){
            fprintf(stderr, "pscnt file not found:%s\n", ps->s);
            exit(EXIT_FAILURE);
        }
        util::LineReader lr(ps->s);
        lr.getline(line);
        while(lr.getline(line)){
            util::split(line, vstr, "\t");
            if(lctg != vstr[0]){
                pct = new contig_pscnt_t();
                if(vstr[0] == cnorm){
                    pct->isnorm = true;
                    pnct = pct;
                }
                map4pscnt[vstr[0]] = pct;
                lctg = vstr[0];
            }
            pct->ttpr += atol(vstr[2].c_str());
            auto iter = pct->pscnt.find(vstr[1]);
            if(iter == pct->pscnt.end()){
                psc = new pscnt_info_t();
                pct->pscnt[vstr[1]] = psc;
            }else{
                psc = iter->second;
            }
            if(psc->anno.empty()) psc->anno = vstr[4];
            psc->count[0] = atol(vstr[2].c_str());
            psc->freq[0] = atof(vstr[3].c_str());
        }
    }
    // update freq use only pattern seq count as sample space
    for(auto& p1: map4pscnt){
        for(auto& p2: p1.second->pscnt){
            p2.second->freq[0] = (double)p2.second->count[0]/(double)p1.second->ttpr;
        }
    }
    // add control freq back if there is any
    if(pnct){
        for(auto& p1: map4pscnt){
            for(auto& p2: p1.second->pscnt){
                auto iter = pnct->pscnt.find(p2.first);
                if(iter != pnct->pscnt.end()){
                    p2.second->count[1] = iter->second->count[0];
                    p2.second->freq[1] = iter->second->freq[0];
                }
            }
        }
    }
    pcpnorm = pnct;
}

void pscmp_opt_t::parse_freq(){
    std::string line, lctg;
    std::vector<std::string> vstr;
    contig_nuc_freq_t* pnuc = NULL;
    int bidx = 0;
    for(auto& s: samples){
        ps->l = 0;
        ksprintf(ps, "%s/%s/%s", indir.c_str(), s.c_str(), PSCNT_NXC_TSV);
        if(!(util::exists(ps->s))){
            fprintf(stderr, "psnxc file not found:%s\n", ps->s);
            exit(EXIT_FAILURE);
        }
        util::LineReader lr(ps->s);
        lr.getline(line);
        while(lr.getline(line)){
            util::split(line, vstr, "\t");
            plen = vstr.size()-2;
            if(lctg != vstr[0]){
                pnuc = new contig_nuc_freq_t();
                pnuc->cfreq = (double**)malloc(plen * sizeof(double*));
                pnuc->len = plen;
                for(int i = 0; i < plen; ++i) pnuc->cfreq[i] = (double*)calloc(5, sizeof(double));
                map4psnuc[vstr[0]] = pnuc;
                lctg = vstr[0];
            }
            bidx = nuc_to_3bit[(int)vstr[1][0]];
            for(int i = 0; i < plen; ++i) pnuc->cfreq[i][bidx] = atof(vstr[2+i].c_str());
        }
    }
    if(cnorm.size()){
        auto iter = map4psnuc.find(cnorm);
        nfreq = iter->second->cfreq;
    }
}


void pscmp_opt_t::parse_in(){
    parse_qc();
    parse_freq();
    parse_pscnt();
}

void pscmp_opt_t::genrpt(){
    kstring_t* s = (kstring_t*)calloc(1, sizeof(kstring_t));
    html2head(s);
    html2body(s);
    html2foot(s);
    FILE* fp = fopen(outhtml.c_str(), "w");
    fwrite(s->s, sizeof(char), s->l, fp);
    fclose(fp);
}

void pscmp_opt_t::html2head(kstring_t* s){
    hmo.html2head(s, "pscmp");
}

void pscmp_opt_t::html2body(kstring_t* s){
    allqc2html(s);
    allps2html(s);
    allnuc2html(s);
}

void pscmp_opt_t::html2foot(kstring_t* s){
    hmo.html2endbody(s);
    hmo.printJssave2svg(s);
    hmo.html2footonly(s);
}

void pscmp_opt_t::allqc2html(kstring_t* s){
    // sec beg
    ksprintf(s, "<div class='section_div'>\n");
    ksprintf(s, "<div class='section_title' onclick=showOrHide('psqc')><a name='summary'>QC Summary</a></div>\n");
    ksprintf(s, "<div id='psqc'>\n");
    ksprintf(s, "<table class='summary_table' id='qc_summary_table'>\n");
    ts->l = 0;
    contig_pscnt_qc_t::html2thead(s);
    contig_pscnt_qc_t::tsv2head(ts);
    ksprintf(s, "<tdoby>");
    for(auto& pqc: map4ctgqc){
        pqc.second->html2tbody(s, pqc.first);
        pqc.second->tsv2body(ts, pqc.first);
    }
    ksprintf(s, "</tdoby>");
    ksprintf(s, "</table>");
    ksprintf(s, "<div>right click to download <a href=\"./%s\">tsv</a></div>", n_outpsqc.c_str());
    // sec end
    ksprintf(s, "</div>\n");
    ksprintf(s, "</div>\n");
    // write to tsv
    FILE* fp = fopen(outpsqc.c_str(), "w");
    fwrite(ts->s, sizeof(char), ts->l, fp);
    fclose(fp);
}

void pscmp_opt_t::allps2html(kstring_t* s){
    // sec beg
    ksprintf(s, "<div class='section_div'>\n");
    ksprintf(s, "<div class='section_title' onclick=showOrHide('pscnt')><a name='summary'>Pattern Seq Count</a></div>\n");
    ksprintf(s, "<div id='pscnt'>\n");
    ksprintf(s, "<div>right click to download <a href=\"./%s\">tsv</a></div>", n_outpscnt.c_str());
    ts->l = 0;
    pscnt_info_t::tsv2head(ts);
    for(auto& psc: map4pscnt){
        for(auto& pxx: psc.second->pscnt){
            pxx.second->tsv2body(ts, psc.first, pxx.first);
        }
    }
    // sec end
    ksprintf(s, "</div>\n");
    ksprintf(s, "</div>\n");
    // write to tsv
    FILE* fp = fopen(outpscnt.c_str(), "w");
    fwrite(ts->s, sizeof(char), ts->l, fp);
    fclose(fp);
}

void pscmp_opt_t::allnuc2html(kstring_t* s){
    // sec beg
    ksprintf(s, "<div class='section_div'>\n");
    ksprintf(s, "<div class='section_title' onclick=showOrHide('psnuc')><a name='summary'>Pattern Base Dist</a></div>\n");
    ksprintf(s, "<div id='psnuc'>\n");
    if(cnorm.size()) onenuc2html(s, cnorm);
    // sub sec of each contig beg
    for(auto& ctg: contigs){
        if(ctg != cnorm){
            onenuc2html(s, ctg);
        }
    }
    // sub sec of each contig end
    // sec end
    ksprintf(s, "</div>\n");
    ksprintf(s, "</div>\n");
}

void pscmp_opt_t::onenuc2html(kstring_t* s, const std::string& n){
    std::string sid = n + "_psnuc";
    ksprintf(s, "<div class='section_title_lite' onclick=showOrHide('%s')>Nucleotide Frequency and Possible Logo of %s</div>\n", sid.c_str(), n.c_str());
    ksprintf(s, "<div id='%s'>\n", sid.c_str());
    auto iter = map4psnuc.find(n);
    contig_nuc_freq_t* pnuc = iter->second;
    // dist
    nuc2freq(s, n, pnuc);
    // logo
    if(nfreq) nuc2logo4exist2(s, n);
    if(nfreq) nuc2logo4exist1(s, n);
    ksprintf(s, "</div>\n");
}

void pscmp_opt_t::nuc2freq(kstring_t* s, const std::string& n, contig_nuc_freq_t* f){
    const char* colors[5] = {"rgba(128,128,0,1.0)", "rgba(0,255,0,1.0)",
                             "rgba(0,0,255,1.0)", "rgba(128,0,128,1.0)", "rgba(20,20,20,1.0)"}; // ACGTN
    std::string subsect = "Nucleotide frequency along pattern of " + n;
    std::string divName = util::replace(subsect, " ", "_"); 
    std::string title = "'nucleotide frequency along pattern of " + n + "'";
    ksprintf(s, "<div class='section_title_offset_left'><a title='click to hide/show' onclick=showOrHide('%s')>%s</a></div>\n",
                divName.c_str(), subsect.c_str());
    ksprintf(s, "<div id='%s'>\n", divName.c_str());
    ksprintf(s, "<div class='figure' id='plot_%s' style=\"margin-left:20px;\"></div>\n", divName.c_str());
    ksprintf(s, "</div>\n");
    ksprintf(s, "\n<script type=\"text/javascript\">\n");
    // case
    for(int i = 0; i <= 4; ++i){
        ksprintf(s, "var case_nuc_%c = {\n", bit3_to_nuc[i]);
        ksprintf(s, "  x: [");
        for(int j = 0; j < plen; ++j) ksprintf(s, "%d,", j+1);
        ksprintf(s, "],\n");
        ksprintf(s, "  y: [");
        for(int j = 0; j < plen; ++j) ksprintf(s, "%lf,", f->cfreq[j][i]);
        ksprintf(s, "],\n");
        ksprintf(s, "  name: '%c',\n", bit3_to_nuc[i]);
        ksprintf(s, "  type: 'bar',\n");
        ksprintf(s, "  domain: {row: 0,},\n");
        ksprintf(s, "  marker: {color: '%s'},\n", colors[i]);
        ksprintf(s, "  xaxis: 'x1',\n");
        ksprintf(s, "  yaxis: 'y1',\n");
        ksprintf(s, "};\n");
    }
    // fold change
    if(nfreq){
        for(int i = 0; i <= 4; ++i){
            ksprintf(s, "var nuc_fold_%c = {\n", bit3_to_nuc[i]);
            ksprintf(s, "  x: [");
            for(int j = 0; j < plen; ++j) ksprintf(s, "%d,", j+1);
            ksprintf(s, "],\n");
            ksprintf(s, "  y: [");
            for(int j = 0; j < plen; ++j){
                if(f->cfreq[j][i] > 0 && nfreq[j][i] > 0){
                    ksprintf(s, "%lf,", log2(f->cfreq[j][i]/nfreq[j][i]));
                }else{
                    ksprintf(s, ".0,");
                }
            }
            ksprintf(s, "],\n");
            ksprintf(s, "  name: '%c',\n", bit3_to_nuc[i]);
            ksprintf(s, "  mode: 'markers',\n");
            ksprintf(s, "  type: 'scatter',\n");
            ksprintf(s, "  domain: {row: 1,},\n");
            ksprintf(s, "  marker: {color: '%s'},\n", colors[i]);
            ksprintf(s, "  xaxis: 'x2',\n");
            ksprintf(s, "  yaxis: 'y2',\n");
            ksprintf(s, "  showlegend: false,\n");
            ksprintf(s, "};\n");
        }
    }
    ksprintf(s, "var layout={\n");
    ksprintf(s, "  grid: {rows: %d, columns: 1, pattern: 'independent', roworder: 'top to bottom'},\n", nfreq ? 2 : 1);
    ksprintf(s, "  height: 600,\n");
    ksprintf(s, "  barmode: 'stack',\n");
    ksprintf(s, "  bargap: 0.01,\n");
    ksprintf(s, "  xaxis:{\n");
    ksprintf(s, "    title: 'positions',\n");
    ksprintf(s, "    dtick: 1,\n");
    ksprintf(s, "  },\n");
    ksprintf(s, "  yaxis:{\n");
    ksprintf(s, "    title: 'frequency',\n");
    ksprintf(s, "    range: [0, 1],\n");
    ksprintf(s, "  },\n");
    ksprintf(s, "  xaxis2:{\n");
    ksprintf(s, "    title: 'positions',\n");
    ksprintf(s, "    dtick: 1,\n");
    ksprintf(s, "  },\n");
    ksprintf(s, "  yaxis2:{\n");
    ksprintf(s, "    title: 'log2(case/norm)',\n");
    ksprintf(s, "    range: [%lf, %lf],\n", fcymax, fcymin);
    ksprintf(s, "  },\n");
    ksprintf(s, "};\n");

    ksprintf(s, "var data = [");
    for(int i = 0; i <= 4; ++i){
        ksprintf(s, "case_nuc_%c,", bit3_to_nuc[i]);
    }
    if(nfreq){
        for(int i = 0; i <= 4; ++i){
            ksprintf(s, "nuc_fold_%c,", bit3_to_nuc[i]);
        }
    }
    ksprintf(s, "];\n");

    ksprintf(s, "var config = {\n");
    ksprintf(s, "  toImageButtonOptions: {\n");
    ksprintf(s, "    format: 'svg',\n");
    ksprintf(s, "     filename: 'nuc_dist_of_%s',\n", n.c_str());
    ksprintf(s, "     height: %d,\n", 600);
    ksprintf(s, "     width: %d,\n", 600);
    ksprintf(s, "     scale: 1,\n");
    ksprintf(s, "  }\n");
    ksprintf(s, "};\n");
    
    ksprintf(s, "Plotly.newPlot('plot_%s', data, layout, config);\n", divName.c_str());
    ksprintf(s, "</script>\n");
}

void pscmp_opt_t::nuc2logo4exist2(kstring_t* s, const std::string& n){
    // get seqs
    auto iter = map4pscnt.find(n);
    if(!vnuc){
        vnuc = (int64_t**)calloc(plen, sizeof(int64_t*));
        for(int i = 0; i < plen; ++i) vnuc[i] = (int64_t*)calloc(4, sizeof(int64_t));
    }else{
        for(int i = 0; i < plen; ++i) memset(vnuc[i], 0, 4*sizeof(int64_t));
    }
    double fcase = .0, fnorm = .0, rate = .0, flog2 = .0;
    int64_t ttvalr = 0;
    for(auto& ps: iter->second->pscnt){
        if(ps.second->count[0] > 0 && ps.second->count[1] > 0){
            fcase = (double)ps.second->count[0]/(double)iter->second->ttpr;
            fnorm = (double)ps.second->count[1]/(double)nttr;
            rate = fcase/fnorm;
            flog2 = log2(rate);
            if(flog2 < maxfcl){
                for(size_t j = 0; j < ps.first.size(); ++j){
                    vnuc[j][nuc_to_3bit[(int)ps.first[j]]] += ps.second->count[0];
                    ++ttvalr;
                }
            }
        }
    }
    // plot
    std::string subsect = "Nucleotide logo along pattern of " + n + " with " + std::to_string(ttvalr) + " seqs exist in case and normal";
    std::string divName = util::replace(subsect, " ", "_");
    std::string title = "'nucleotide logo along pattern of" + n + "'";
    ksprintf(s, "<div class='section_title_offset_left'><a title='click to hide/show' onclick=showOrHide('%s')>%s</a></div>\n",
                divName.c_str(), subsect.c_str());
    ksprintf(s, "<div id='%s'>\n", divName.c_str());
    ksprintf(s, "<div class='figure' id='plot_%s' style=\"width:300px;height:300px;margin-left:20px;\"></div>\n", divName.c_str());
    ksprintf(s, "<div><button style=\"margin-left:20px;\" class='btn btn-primary' type='submit' onclick=\"save2svg('plot_%s', '%s', 300, 300)\">Click To Save SVG</button></div>", divName.c_str(), divName.c_str());
    ksprintf(s, "</div>\n");
    ksprintf(s, "\n<script type=\"text/javascript\">\n");
    ksprintf(s, "const nuc_fold2_change_of_%s = [", n.c_str());
    for(int j = 0; j < plen; ++j){
        ksprintf(s, "[");
        int64_t ttbr = 0;
        for(int i = 0; i < 4; ++i) ttbr += vnuc[j][i];
        for(int i = 0; i < 4; ++i){
            if(ttbr) ksprintf(s, "%lf,", (double)vnuc[j][i]/(double)ttbr);
            else ksprintf(s, ".0,");
        }
        ksprintf(s, "],");
    }
    ksprintf(s, "];\n");
    ksprintf(s, "logojs.embedDNALogo(document.getElementById(\"plot_%s\"), { ppm: nuc_fold2_change_of_%s });", divName.c_str(), n.c_str());
    ksprintf(s, "</script>\n");
}

void pscmp_opt_t::nuc2logo4exist1(kstring_t* s, const std::string& n){
    // get seqs
    auto iter = map4pscnt.find(n);
    if(!vnuc){
        vnuc = (int64_t**)calloc(plen, sizeof(int64_t*));
        for(int i = 0; i < plen; ++i) vnuc[i] = (int64_t*)calloc(4, sizeof(int64_t));
    }else{
        for(int i = 0; i < plen; ++i) memset(vnuc[i], 0, 4*sizeof(int64_t));
    }
    int64_t ttvalr = 0;
    for(auto& ps: pcpnorm->pscnt){
        if(iter->second->pscnt.find(ps.first) == iter->second->pscnt.end()){
            for(size_t j = 0; j < ps.first.size(); ++j){
                vnuc[j][nuc_to_3bit[(int)ps.first[j]]] += ps.second->count[0];
            }
            ++ttvalr;
        }
    }
    // plot
    std::string subsect = "Nucleotide logo along pattern of " + n + " with " + std::to_string(ttvalr) + " seqs exist only in normal";
    std::string divName = util::replace(subsect, " ", "_");
    std::string title = "'nucleotide logo along pattern of " + n + "'";
    ksprintf(s, "<div class='section_title_offset_left'><a title='click to hide/show' onclick=showOrHide('%s')>%s</a></div>\n",
                divName.c_str(), subsect.c_str());
    ksprintf(s, "<div id='%s'>\n", divName.c_str());
    ksprintf(s, "<div class='figure' id='plot_%s' style=\"width:300px;height:300px;margin-left:20px;\"></div>\n", divName.c_str());
    ksprintf(s, "<div><button style=\"margin-left:20px;\" class='btn btn-primary' type='submit' onclick=\"save2svg('plot_%s',, '%s', 300, 300)\">Click To Save SVG</button></div>", divName.c_str(), divName.c_str());
    ksprintf(s, "</div>\n");
    ksprintf(s, "\n<script type=\"text/javascript\">\n");
    ksprintf(s, "const nuc_fold2_change_of1_%s = [", n.c_str());
    for(int j = 0; j < plen; ++j){
        ksprintf(s, "[");
        int64_t ttbr = 0;
        for(int i = 0; i < 4; ++i) ttbr += vnuc[j][i];
        for(int i = 0; i < 4; ++i){
            if(ttbr) ksprintf(s, "%lf,", (double)vnuc[j][i]/(double)ttbr);
            else ksprintf(s, ".0,");
        }
        ksprintf(s, "],");
    }
    ksprintf(s, "];\n");
    ksprintf(s, "logojs.embedDNALogo(document.getElementById(\"plot_%s\"), { ppm: nuc_fold2_change_of1_%s });", divName.c_str(), n.c_str());
    ksprintf(s, "</script>\n");
}
void pscmp_usage(pscmp_opt_t* opt, char* arg0){
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: %s %s [options]\n\n", PACKAGE_NAME, arg0);
    fprintf(stderr, "Options: -i FILE   input directory of pscnt analysis\n");
    fprintf(stderr, "         -s FILE   sample list file\n");
    fprintf(stderr, "         -f FLOAT  max log2 fold change for logo plot [%f]\n", opt->maxfcl);
    fprintf(stderr, "         -m FLOAT  min fold change limit to plot [%f]\n", opt->fcymin);
    fprintf(stderr, "         -M FLOAT  max fold change limit to plot [%f]\n", opt->fcymax);
    fprintf(stderr, "         -n STR    normal or control sample\n");
    fprintf(stderr, "         -o FILE   output directory [%s]\n", opt->outdir.c_str());
    fprintf(stderr, "         -9 FILE   js/css cdn configure json file\n");
    fprintf(stderr, "\n");
}

int pscmp_main(int argc, char** argv){
    pscmp_opt_t opt;
    if(argc == 1){
        pscmp_usage(&opt, argv[0]);
        return 0;
    }
    int c = 0;
    while((c = getopt(argc, argv, "i:s:f:m:M:n:o:9:h")) >= 0){
        switch(c){
            case 'i': opt.indir = optarg; break;
            case 's': opt.splist = optarg; break;
            case 'f': opt.maxfcl = atof(optarg); break;
            case 'm': opt.fcymin = atof(optarg); break;
            case 'M': opt.fcymax = atof(optarg); break;
            case 'n': opt.cnorm = optarg; break;
            case 'o': opt.outdir = optarg; break;
            case '9': opt.hrjsn = optarg; break;
            case 'h': pscmp_usage(&opt, argv[0]); return 0; break;
            default: break;
        }
    }
    if(opt.valid()){
        opt.init();
        opt.parse_in();
        opt.genrpt();
        return 0;
    }else{
        return 1;
    }
}
