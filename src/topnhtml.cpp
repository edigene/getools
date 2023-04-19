#include "htmlutil.h"
#include "gedetect.h"
#include "geplot.h"

void gep_opt_t::topn2html(kstring_t* s){
    std::string subsect = "Top" + std::to_string(topn) + " allele sequences";
    std::string divName = util::replace(subsect, " ", "_");
    std::string title = "'top " + std::to_string(topn) + " allele sequences'";
    ksprintf(s, "<div class='subsection_title'><a title='click to hide/show' onclick=showOrHide('%s')>%s</a></div>\n",
                divName.c_str(), subsect.c_str());
    ksprintf(s, "<div id='%s'>\n", divName.c_str());
    ksprintf(s, "<div class='sub_section_tips'>Alignment information of each sequence will be shown on mouse over.<br>tn: top number<br>fs: frameshift<br>sc: sequence count<br>af: frequencey<br>cs: alignment cigar<br>aa: accumulated frequence from top1 sequence to this one<br>cn: number of singele cells with this sequence</div>\n");
    ksprintf(s, "<style type=\"text/css\">\n");
    ksprintf(s, ".topnfig {width:1000px;height:%ldpx;}\n", 200 + hvars.size()*15);
    ksprintf(s, "</style>\n");
    ksprintf(s, "<div class='topnfig' id='plot_%s'></div>\n", divName.c_str());
    ksprintf(s, "</div>\n");
    ksprintf(s, "\n<script type=\"text/javascript\">\n");

    std::string xtvs = "[", xtts = "[";
    for(int i = 0; i < rlen; ++ i){
        if(abs(i - sgrbeg) % 5 == 0){
            xtvs.append(std::to_string(i) + ",");
            if(i - sgrbeg > 0) xtts.append("'+" + std::to_string(i-sgrbeg) + "',");
            else xtts.append("'" + std::to_string(i-sgrbeg) + "',");
        }
    }
    xtvs.append("]");
    xtts.append("]");

    std::string cscvs = "var colorscaleValue = [\n";
    cscvs.append("[0.00,'#EA4335'], [0.20,'#4285F4'], [0.40,'#FBBC05'], [0.60,'#34A853'], [0.80, '#9933FF'], [1.00,'#FFFFFF'],]\n"); // ACGTN-(I)
    std::string xvars = "[", yvars = "[", svars = "[", ttvars = "[", tvvars = "[", colvars = "[", tmpcol, tmphov, bmks;
    kstring_t* seq = NULL;
    kstring_t* cseq = NULL;
    ttvars.append("'Input Amplicon Sequence',");
    double accu_af = .0;
    for(int j = hvars.size() - 1; j >= 0; --j){
        tvvars.append(std::to_string(j) + ",");
        seq = hvars[hvars.size()-1-j]->seq;
        cseq = hvars[j]->seq;
        tmpcol = "  ["; tmphov = " [";
        bmks.clear();
        for(int i = 0; i < rlen; ++i){
            int k = hvars.size()-1-j;
            if(i == 0 && k >= 1){
                accu_af += hvars[k]->af;
                bmks.append("tn:" + std::to_string(k) + "<br>fs:" + std::to_string(hvars[k]->fscnt) + "<br>");
                bmks.append("sc:" + std::to_string(hvars[k]->cc) + "<br>");
                bmks.append("af:" + std::to_string(hvars[k]->af*100) + "%<br>");
                bmks.append("cs:" + std::string(hvars[k]->cigar->s) + "<br>");
                bmks.append("aa:" + std::to_string(accu_af*100) + "%<br>");
                bmks.append("cn:" + std::to_string(hvars[k]->cn));
                ttvars.append("'" + bmks + "',");
            }
            xvars.append(std::to_string(i) + ",");
            yvars.append(std::to_string(j) + ",");
            svars.append("'" + std::string(1, seq->s[i]) + "',");
            switch(cseq->s[i]){
                case 'A':
                    tmpcol.append(std::to_string(0.00) + ",");
                    break;
                case 'C':
                    tmpcol.append(std::to_string(0.20) + ",");
                    break;
                case 'G':
                    tmpcol.append(std::to_string(0.40) + ",");
                    break;
                case 'T':
                    tmpcol.append(std::to_string(0.60) + ",");
                    break;
                case 'I':
                    tmpcol.append(std::to_string(0.80) + ",");
                    break;
                default:
                    tmpcol.append(std::to_string(1.00) + ",");
                    break;
            }
        }
        tmpcol.append("],\n");
        tmphov.append("],\n");
        colvars.append(tmpcol);
    }
    xvars.append("]"); yvars.append("]"); svars.append("]"); tvvars.append("]"); ttvars.append("]"), colvars.append("]");

    std::string jsnstr;
    jsnstr.append(cscvs);
    jsnstr.append("var ncsa = {\n");
    jsnstr.append("  x: " + xvars + ",\n");
    jsnstr.append("  y: " + yvars + ",\n");
    jsnstr.append("  mode: 'text',\n");
    jsnstr.append("  text: " + svars + ",\n");
    jsnstr.append("  type: 'scattergl',\n");
    jsnstr.append("  hoverinfo: 'none',\n");
    jsnstr.append("  textfont: {\n");
    jsnstr.append("    size: 10,\n");
    jsnstr.append("  },\n");
    jsnstr.append("}\n");

    jsnstr.append("var ncol = {\n");
    jsnstr.append("  z: " + colvars + ",\n");
    jsnstr.append("  type: 'heatmap',\n");
    jsnstr.append("  hoverinfo: 'y',\n");
    jsnstr.append("  colorscale: colorscaleValue,\n");
    jsnstr.append("  showscale: false,\n");
    jsnstr.append("}\n");

    jsnstr.append("var layout = {\n");
    jsnstr.append("  title: " + title + ",\n");
    jsnstr.append("  width: 1000,\n");
    jsnstr.append("  height: " + std::to_string(200+15 * hvars.size()) + ",\n");
    jsnstr.append("  xaxis: {\n");
    jsnstr.append("    range: [" + std::to_string(hrbeg-0.5) + "," + std::to_string(hrend+0.5) + "],\n");
    jsnstr.append("    showline: false,\n");
    jsnstr.append("    zeroline: false,\n");
    jsnstr.append("    showgrid: false,\n");
    jsnstr.append("    tickmode: 'array',\n");
    jsnstr.append("    tickvals: " + xtvs + ",\n");
    jsnstr.append("    ticktext: " + xtts + ",\n");
    jsnstr.append("    title: {\n");
    jsnstr.append("      text: 'relative position to sgRNA beg',\n");
    jsnstr.append("      standoff: 2,\n");
    jsnstr.append("    },\n");
    jsnstr.append("  },\n");
    jsnstr.append("  yaxis: {\n");
    jsnstr.append("    range: [-0.5, " + std::to_string(hvars.size()+2) + "],\n");
    jsnstr.append("    showline: false,\n");
    jsnstr.append("    zeroline: false,\n");
    jsnstr.append("    showgrid: false,\n");
    jsnstr.append("    visible: true,\n");
    jsnstr.append("    ticks: '',\n");
    jsnstr.append("    standoff: 4,\n");
    jsnstr.append("    automargin: true,\n");
    jsnstr.append("    showticklabels: false,\n");
    jsnstr.append("    tickmode: 'array',\n");
    jsnstr.append("    tickvals: " + tvvars + ",\n");
    jsnstr.append("    ticktext: " + ttvars + ",\n");
    jsnstr.append("  },\n");
    jsnstr.append("}\n");

    jsnstr.append("var config = {\n");
    jsnstr.append("  toImageButtonOptions: {\n");
    jsnstr.append("    format: 'svg',\n");
    jsnstr.append("     filename: 'top" + std::to_string(topn) + "_events_of_" + std::string(name) + "',\n");
    jsnstr.append("     height: " + std::to_string(200+hvars.size()*15) + ",\n");
    jsnstr.append("     width: 1000,\n");
    jsnstr.append("     scale: 1,\n");
    jsnstr.append("  }\n");
    jsnstr.append("};\n");


    jsnstr.append("var data = [ncsa, ncol]\n");

    jsnstr.append("Plotly.newPlot('plot_" + divName + "', data, layout, config);\n");
    ksprintf(s, "%s", jsnstr.c_str());
    ksprintf(s, "</script>\n");
}
