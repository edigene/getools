#include "gedetect.h"
#include "htmlutil.h"

void GEDetector::reportHTML(kstring_t* s){
    summary();
    // compute pos relative to sgRNA beg
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
    // summary graph
    {
        std::string subsect = "Edit efficience computation";
        std::string divName = util::replace(subsect, " ", "_");
        std::string title = "'edited rate(" + std::to_string(edieff * 100) + "%), ";
        title.append("mutated rate(" + std::to_string(muteff * 100) + "%)'");
        ksprintf(s, "<div class='subsection_title'><a title='click to hide/show' onclick=showOrHide('%s')>%s</a></div>\n", divName.c_str(), subsect.c_str());
        ksprintf(s, "<div id='%s'>\n", divName.c_str());
        ksprintf(s, "<div class='sub_section_tips'>Percentage of each type of sequences will be shown on mouse over.</div>\n");
        ksprintf(s, "<div class='figure' id='plot_%s'></div>\n", divName.c_str());
        ksprintf(s, "</div>\n");
        ksprintf(s, "\n<script type=\"text/javascript\">\n");

        std::string jsnstr;
        jsnstr.append("var qdpie = {\n");
        jsnstr.append("  values: [");
        jsnstr.append(std::to_string(dropcnt) + ",");
        jsnstr.append(std::to_string(totcnt) + "],\n");
        jsnstr.append("  labels: ['DroppedSeq', 'QualifiedSeq'],\n");
        jsnstr.append("  title: 'All Assigned Sequences',\n");
        jsnstr.append("  type: 'pie',\n");
        jsnstr.append("  textinfo: 'label',\n");
        jsnstr.append("  textposition: 'inside',\n");
        jsnstr.append("  hoverinfo: 'label+value+percent',\n");
        jsnstr.append("  insidetextorientation: 'radial',\n");
        jsnstr.append("  domain: {column: 0},\n");
        jsnstr.append("};\n");

        jsnstr.append("var edpie = {\n");
        jsnstr.append("  values: [");
        jsnstr.append(std::to_string(refcnt) + ",");
        jsnstr.append(std::to_string(edicnt) + ",");
        jsnstr.append(std::to_string(othcnt) + "],\n");
        jsnstr.append("  labels: ['Refseq', 'EditedSeq', 'OtherSeq'],\n");
        jsnstr.append("  title: 'All Qualified Sequences',\n");
        jsnstr.append("  type: 'pie',\n");
        jsnstr.append("  textinfo: 'label',\n");
        jsnstr.append("  textposition: 'inside',\n");
        jsnstr.append("  hoverinfo: 'label+value+percent',\n");
        jsnstr.append("  insidetextorientation: 'radial',\n");
        jsnstr.append("  domain: {column: 1},\n");
        jsnstr.append("};\n");

        jsnstr.append("var droppie = {\n");
        jsnstr.append("  values: [");
        for(int i = 0; i < TOTAL_DROP_REASON; ++i){
            jsnstr.append(std::to_string(droprcnt[i]) + ",");
        }
        jsnstr.append("],\n");
        jsnstr.append("  labels: [");
        for(int i = 0; i < TOTAL_DROP_REASON; ++i){
            jsnstr.append("'" + std::string(drop_reasons_str[i]) + "',");
        }
        jsnstr.append("],\n");
        jsnstr.append("  title: 'All Dropped Sequences',\n");
        jsnstr.append("  type: 'pie',\n");
        jsnstr.append("  textinfo: 'label',\n");
        jsnstr.append("  textposition: 'inside',\n");
        jsnstr.append("  hoverinfo: 'label+value+percent',\n");
        jsnstr.append("  insidetextorientation: 'radial',\n");
        jsnstr.append("  domain: {column: 2},\n");
        jsnstr.append("};\n");


        jsnstr.append("var layout = {\n");
        jsnstr.append("  title:" + title + ",\n");
        jsnstr.append("  grid: {rows: 1, columns: 3},\n");
        jsnstr.append("  height: 600,\n");
        jsnstr.append("  showlegend: false,\n");
        jsnstr.append("};\n");

        jsnstr.append("var data = [qdpie, edpie, droppie];\n");

        jsnstr.append("var config = {\n");
        jsnstr.append("  toImageButtonOptions: {\n");
        jsnstr.append("    format: 'svg',\n");
        jsnstr.append("     filename: '" + divName + "_of_" + std::string(name) + "',\n");
        jsnstr.append("     height: " + std::to_string(opt->hmo.figh) + ",\n");
        jsnstr.append("     width: " + std::to_string(opt->hmo.figw) + ",\n");
        jsnstr.append("     scale: 1,\n");
        jsnstr.append("  }\n");
        jsnstr.append("};\n");

        jsnstr.append("Plotly.newPlot('plot_" + divName + "', data, layout, config);\n");
        ksprintf(s, "%s", jsnstr.c_str());
        ksprintf(s, "</script>\n");
    }
    // cluster graph
    {
        if(kmcjs->s) kputsn(kmcjs->s, kmcjs->l, s);
    }
    // nucleotide distribution bar plot
    {
        std::string colors[5] = {"rgba(128,128,0,1.0)", "rgba(0,255,0,1.0)",
                                 "rgba(0,0,255,1.0)", "rgba(128,0,128,1.0)", "rgba(20,20,20,1.0)"}; // ACGTN
        std::string subsect = "SNV along amplicon reference sequence";
        std::string divName = util::replace(subsect, " ", "_");
        std::string title = "'snv along amplicon reference'";
        ksprintf(s, "<div class='subsection_title'><a title='click to hide/show' onclick=showOrHide('%s')>%s</a></div>\n",
                    divName.c_str(), subsect.c_str());
        ksprintf(s, "<div id='%s'>\n", divName.c_str());
        ksprintf(s, "<div class='sub_section_tips'>Mismatch information and position relative to sgrRNA beg/end will be shown on mouse over, slide left and right to see more.</div>\n");
        ksprintf(s, "<div class='figure' id='plot_%s'></div>\n", divName.c_str());
        ksprintf(s, "</div>\n");
        ksprintf(s, "\n<script type=\"text/javascript\">\n");
        
        std::string jsnstr;
        jsnstr.append("var refbar = {\n");
        jsnstr.append("  x: [");
        for(int i = 0; i < rlen; ++i) jsnstr.append(std::to_string(i) + ",");
        jsnstr.append("],\n");
        jsnstr.append("  y: [");
        for(int i = 0; i < rlen; ++i) jsnstr.append(std::to_string(refcov[i]) + ",");
        jsnstr.append("],\n");
        jsnstr.append("  name: 'Reference',\n");
        jsnstr.append("  text: [");
        std::string postr;
        for(int i = 0; i < rlen; ++i){
            int64_t ttbi = nccnt[0][i]+nccnt[1][i]+nccnt[2][i]+nccnt[3][i]+nccnt[4][i];
            double afb = 0;
            int b = 0;
            std::string bmk, sgr;
            for(b = 0; b < 5; ++b){
                afb = (double)nccnt[b][i]/(double)ttbi * 100;
                bmk.append(std::string(1, bit3_to_nuc[b]));
                bmk.append("(" + std::to_string(nccnt[b][i]) + ",");
                bmk.append(std::to_string(afb) + "%)<br>");
            }
            sgr.append("sgrBeg");
            if(i < sgrbeg) sgr.append(std::to_string(i-sgrbeg));
            else if(i > sgrbeg) sgr.append("+" + std::to_string(i-sgrbeg));
            sgr.append(",sgrEnd");
            if(i < sgrend) sgr.append(std::to_string(i-sgrend));
            else if(i > sgrend) sgr.append("+" + std::to_string(i-sgrend));
            postr.append("'" + sgr + "',");
            bmk.append(sgr);
            jsnstr.append("'" + bmk + "',");
        }
        jsnstr.append("],\n");
        jsnstr.append("  textposition: 'none',\n");
        jsnstr.append("  type: 'bar',\n");
        jsnstr.append("  domain: {row: 0},\n");
        jsnstr.append("  marker: {color: [");
        for(int i = 0; i < rlen; ++i) jsnstr.append("'" + colors[nuc_to_3bit[(int)ref[i]]] + "',");
        jsnstr.append("]},\n");
        jsnstr.append("  yaxis: {\n");
        jsnstr.append("    title: 'Coverage'\n");
        jsnstr.append("  },\n");
        jsnstr.append("};\n");

        jsnstr.append("var misbar = {\n");
        jsnstr.append("  x: [");
        for(int i = 0; i < rlen; ++i) jsnstr.append(std::to_string(i) + ",");
        jsnstr.append("],\n");
        jsnstr.append("  y: [");
        for(int i = 0; i < rlen; ++i) jsnstr.append(std::to_string(altcov[i]) + ",");
        jsnstr.append("],\n");
        jsnstr.append("  name: 'Mismatches',\n");
        jsnstr.append("  type: 'bar',\n");
        jsnstr.append("  domain: {row: 0},\n");
        jsnstr.append("  marker: {color: 'rgb(160, 160, 160)'},\n");
        jsnstr.append("  xaxis: 'x1',\n");
        jsnstr.append("  yaxis: 'y1',\n");
        jsnstr.append("};\n");

        double maxmr = .0;
        jsnstr.append("var mislin = {\n");
        jsnstr.append("  x: [");
        for(int i = 0; i < rlen; ++i) jsnstr.append(std::to_string(i) + ",");
        jsnstr.append("],\n");
        jsnstr.append("  y: [");
        for(int i = 0; i < rlen; ++i){
            jsnstr.append(std::to_string(1-refrate[i]) + ",");
            if(1-refrate[i] > maxmr) maxmr = 1-refrate[i];
        }
        jsnstr.append("],\n");
        jsnstr.append("  name: 'MismatchRate',\n");
        jsnstr.append("  mode: 'lines',\n");
        jsnstr.append("  domain: {row: 1},\n");
        jsnstr.append("  marker: {color: 'rgb(255, 102, 255)'},\n");
        jsnstr.append("  text: [" + postr + "],\n");
        jsnstr.append("  xaxis: 'x2',\n");
        jsnstr.append("  yaxis: 'y2',\n");
        jsnstr.append("};\n");

        maxmr = MIN(maxmr+MIN(0.2*maxmr, 0.1), 1.0);

        for(size_t i = 0; i < donorsnv.size(); ++i){
            if(donorsnv[i] > 0){
                jsnstr.append("var dsim" + std::to_string(i) + " = {\n");
                jsnstr.append("  x: [" + std::to_string(i) + "," + std::to_string(i) + "],\n");
                jsnstr.append("  y: [0," + std::to_string(maxmr) + "],\n");
                jsnstr.append("  name: 'desired:" + std::string(1, ref[i]) + "2" + std::string(1, kswge_int2nt[donorsnv[i]]) + "',\n");
                jsnstr.append("  mode: 'lines',\n");
                jsnstr.append("  line: {\n");
                jsnstr.append("    dash: 'dot',\n");
                jsnstr.append("    color: 'green',\n");
                jsnstr.append("  },\n");
                jsnstr.append("  hoverinfo: 'name',\n");
                jsnstr.append("  xaxis: 'x2',\n");
                jsnstr.append("  yaxis: 'y2',\n");
                jsnstr.append("};\n");
            }
        }
        jsnstr.append("var data = [refbar, misbar, mislin");
        for(size_t i = 0; i < donorsnv.size(); ++i){
            if(donorsnv[i] > 0) jsnstr.append(", dsim" + std::to_string(i));
        }
        jsnstr.append("];\n");

        int prrbeg = sgrbeg - 25;
        int prrend = sgrend + 25;
        if(prrend - prrbeg != 100){
            int extprl = (100 - (prrend - prrbeg)) / 2;
            prrbeg = MAX(0, prrbeg - extprl);
            prrend = MIN(prrend + extprl, rlen-1);
        }
        jsnstr.append("var layout={\n");
        jsnstr.append("  title:" + title + ",\n");
        jsnstr.append("  grid: {rows: 2, columns: 1, pattern: 'independent', roworder: 'top to bottom'},\n");
        jsnstr.append("  height: 600,\n");
        jsnstr.append("  bargap: 0.05,\n");
        jsnstr.append("  barmode: 'stack',\n");
        jsnstr.append("  showlegend: false,\n");
        jsnstr.append("  xaxis: {\n");
        jsnstr.append("    range: [" + std::to_string(prrbeg) + "," + std::to_string(prrend) + "],\n");
        jsnstr.append("    tickfont: {size: 10,},\n");
        jsnstr.append("    tickangle: 0,\n");
        jsnstr.append("    tickmode: 'array',\n");
        jsnstr.append("    tickvals: [");
        for(int i = 0; i < rlen; ++i) jsnstr.append(std::to_string(i) + ",");
        jsnstr.append("],\n");
        jsnstr.append("    ticktext: [");
        for(int i = 0; i < rlen; ++i) jsnstr.append("'" + std::string(1, ref[i]) + "',");
        jsnstr.append("],\n");
        jsnstr.append("    title: {\n");
        jsnstr.append("      text: 'amplican reference sequence',\n");
        jsnstr.append("      standoff: 2,\n");
        jsnstr.append("    },\n");
        jsnstr.append("  },\n");
        jsnstr.append("  yaxis: {\n");
        jsnstr.append("    title: {\n");
        jsnstr.append("      text: 'coverage sequence counts',\n");
        jsnstr.append("      standoff: 2,\n");
        jsnstr.append("    },\n");
        jsnstr.append("  },\n");
        jsnstr.append("  xaxis2: {\n");
        jsnstr.append("    range: [" + std::to_string(prrbeg) + "," + std::to_string(prrend) + "],\n");
        jsnstr.append("    tickmode: 'array',\n");
        jsnstr.append("    tickvals: " + xtvs + ",\n");
        jsnstr.append("    ticktext: " + xtts + ",\n");
        jsnstr.append("    title: {\n");
        jsnstr.append("      text: 'relative position to sgRNA beg',\n");
        jsnstr.append("      standoff: 2,\n");
        jsnstr.append("    },\n");
        jsnstr.append("  },\n");
        jsnstr.append("  yaxis2: {\n");
        jsnstr.append("    title: {\n");
        jsnstr.append("      text: 'mismatch rate',\n");
        jsnstr.append("      standoff: 2,\n");
        jsnstr.append("    },\n");
        jsnstr.append("  },\n");
        jsnstr.append("};\n");

        jsnstr.append("var config = {\n");
        jsnstr.append("  toImageButtonOptions: {\n");
        jsnstr.append("    format: 'svg',\n");
        jsnstr.append("     filename: '" + divName + "_of_" + std::string(name) + "',\n");
        jsnstr.append("     height: 600,\n");
        jsnstr.append("     width: " + std::to_string(opt->hmo.figw) + ",\n");
        jsnstr.append("     scale: 1,\n");
        jsnstr.append("  }\n");
        jsnstr.append("};\n");

        jsnstr.append("Plotly.newPlot('plot_" + divName + "', data, layout, config);\n");
        ksprintf(s, "%s", jsnstr.c_str());
        ksprintf(s, "</script>\n");
    }
    // variant count barplot of seqs
    {
        std::string subsect = "Variant sequence counts";
        std::string divName = util::replace(subsect, " ", "_");
        std::string title = "'variant sequence counts'";
        ksprintf(s, "<div class='subsection_title'><a title='click to hide/show' onclick=showOrHide('%s')>%s</a></div>\n",
                    divName.c_str(), subsect.c_str());
        ksprintf(s, "<div id='%s'>\n", divName.c_str());
        ksprintf(s, "<div class='sub_section_tips'>Sequence counts of each type/length of variants  will be shown on mouse over.</div>\n");
        ksprintf(s, "<div class='figure' id='plot_%s'></div>\n", divName.c_str());
        ksprintf(s, "</div>\n");
        ksprintf(s, "\n<script type=\"text/javascript\">\n");

        std::string jsnstr;
        jsnstr.append("var totbar = {\n");
        jsnstr.append("  x: [");
        for(int i = 0; i <= tot_cnt_dist->vcmax; ++i){
            jsnstr.append(std::to_string(i) + ",");
        }
        jsnstr.append("],\n");
        jsnstr.append("  y: [");
        for(int i = 0; i <= tot_cnt_dist->vcmax; ++i){
            jsnstr.append(std::to_string(tot_cnt_dist->vcdist_array[i]) + ",");
        }
        jsnstr.append("],\n");
        jsnstr.append("  type: 'bar',\n");
        jsnstr.append("  name: 'AllVar ',\n");
        jsnstr.append("};\n");

        jsnstr.append("var insbar = {\n");
        jsnstr.append("  x: [");
        for(int i = 0; i <= tot_cnt_dist->vcmax; ++i){
            jsnstr.append(std::to_string(i) + ",");
        }
        jsnstr.append("],\n");
        jsnstr.append("  y: [");
        for(int i = 0; i <= tot_cnt_dist->vcmax; ++i){
            if(i <= ins_cnt_dist->vcmax) jsnstr.append(std::to_string(ins_cnt_dist->vcdist_array[i]) + ",");
            else jsnstr.append("0,");
        }
        jsnstr.append("],\n");
        jsnstr.append("  type: 'bar',\n");
        jsnstr.append("  name: 'Insertion ',\n");
        jsnstr.append("};\n");

        jsnstr.append("var delbar = {\n");
        jsnstr.append("  x: [");
        for(int i = 0; i <= tot_cnt_dist->vcmax; ++i){
            jsnstr.append(std::to_string(i) + ",");
        }
        jsnstr.append("],\n");
        jsnstr.append("  y: [");
        for(int i = 0; i <= tot_cnt_dist->vcmax; ++i){
            if(i <= del_cnt_dist->vcmax) jsnstr.append(std::to_string(del_cnt_dist->vcdist_array[i]) + ",");
            else jsnstr.append("0,");
        }
        jsnstr.append("],\n");
        jsnstr.append("  type: 'bar',\n");
        jsnstr.append("  name: 'Deletion ',\n");
        jsnstr.append("};\n");

        jsnstr.append("var dinbar = {\n");
        jsnstr.append("  x: [");
        for(int i = 0; i <= tot_cnt_dist->vcmax; ++i){
            jsnstr.append(std::to_string(i) + ",");
        }
        jsnstr.append("],\n");
        jsnstr.append("  y: [");
        for(int i = 0; i <= tot_cnt_dist->vcmax; ++i){
            if(i <= din_cnt_dist->vcmax) jsnstr.append(std::to_string(din_cnt_dist->vcdist_array[i]) + ",");
            else jsnstr.append("0,");
        }
        jsnstr.append("],\n");
        jsnstr.append("  type: 'bar',\n");
        jsnstr.append("  name: 'Delins ',\n");
        jsnstr.append("};\n");

        jsnstr.append("var snvbar = {\n");
        jsnstr.append("  x: [");
        for(int i = 0; i <= tot_cnt_dist->vcmax; ++i){
            jsnstr.append(std::to_string(i) + ",");
        }
        jsnstr.append("],\n");
        jsnstr.append("  y: [");
        for(int i = 0; i <= tot_cnt_dist->vcmax; ++i){
            if(i <= snv_cnt_dist->vcmax) jsnstr.append(std::to_string(snv_cnt_dist->vcdist_array[i]) + ",");
            else jsnstr.append("0,");
        }
        jsnstr.append("],\n");
        jsnstr.append("  type: 'bar',\n");
        jsnstr.append("  name: 'SNV ',\n");
        jsnstr.append("};\n");

        jsnstr.append("var data = [totbar, insbar, delbar, dinbar, snvbar];\n");

        jsnstr.append("var layout = {\n");
        jsnstr.append("  title: 'variant counts of sequences',\n");
        jsnstr.append("  xaxis: {\n");
        jsnstr.append("    title: {\n");
        jsnstr.append("      text: 'variant counts of one sequence',\n");
        jsnstr.append("      standoff: 2,\n");
        jsnstr.append("    },\n");
        jsnstr.append("  },\n");
        jsnstr.append("  yaxis: {\n");
        jsnstr.append("    title: {\n");
        jsnstr.append("      text: 'sequence counts',\n");
        jsnstr.append("      standoff: 2,\n");
        jsnstr.append("    },\n");
        jsnstr.append("    automargin: true,\n");
        jsnstr.append("    tickmode: 'array',\n");
        jsnstr.append("    tickvals: [\n");
        int nstep = 20;
        int npcnt = 5;
        int64_t step = (totcnt + 5)/nstep;
        for(int i = 0; i < nstep; ++i){
            jsnstr.append(std::to_string(step * i) + ",");
        }
        jsnstr.append(std::to_string(totcnt) + ",");
        jsnstr.append("],\n");
        jsnstr.append("   ticktext: [\n");
        for(int i = 0; i <= nstep; ++i){
            jsnstr.append("'" +  std::to_string(step * i) + "(" + std::to_string(i*npcnt)+"%)',");
        }
        jsnstr.append("],\n");
        jsnstr.append(" },\n");
        jsnstr.append("}\n");

        jsnstr.append("var config = {\n");
        jsnstr.append("  toImageButtonOptions: {\n");
        jsnstr.append("    format: 'svg',\n");
        jsnstr.append("     filename: '" + divName + "_of_" + std::string(name) + "',\n");
        jsnstr.append("     height: " + std::to_string(opt->hmo.figh) + ",\n");
        jsnstr.append("     width: " + std::to_string(opt->hmo.figw) + ",\n");
        jsnstr.append("     scale: 1,\n");
        jsnstr.append("  }\n");
        jsnstr.append("};\n");

        jsnstr.append("Plotly.newPlot('plot_" + divName + "', data, layout, config);\n");
        ksprintf(s, "%s", jsnstr.c_str());
        ksprintf(s, "</script>\n");
    }
    // ins, del, delins length count barplot
    {
        std::string subsect = "Variant length counts";
        std::string divName = util::replace(subsect, " ", "_");
        std::string title = "'variant length counts'";
        ksprintf(s, "<div class='subsection_title'><a title='click to hide/show' onclick=showOrHide('%s')>%s</a></div>\n",
                    divName.c_str(), subsect.c_str());
        ksprintf(s, "<div id='%s'>\n", divName.c_str());
        ksprintf(s, "<div class='sub_section_tips'>Variant length sequence counts of each type of variants will be shown on mouse over.</div>\n");
        ksprintf(s, "<div class='figure' id='plot_%s'></div>\n", divName.c_str());
        ksprintf(s, "</div>\n");
        ksprintf(s, "\n<script type=\"text/javascript\">\n");

        int64_t maxlen = ins_len_dist->vcmax;
        if(del_len_dist->vcmax > maxlen) maxlen = del_len_dist->vcmax;
        if(din_len_dist->vcmax > maxlen) maxlen = din_len_dist->vcmax;

        std::string jsnstr;
        jsnstr.append("var insbar = {\n");
        jsnstr.append("  x: [");
        for(int i = 0; i <= maxlen; ++i){
            jsnstr.append(std::to_string(i) + ",");
        }
        jsnstr.append("],\n");
        jsnstr.append("  y: [");
        for(int i = 0; i <= maxlen; ++i){
            if(i <= ins_len_dist->vcmax) jsnstr.append(std::to_string(ins_len_dist->vcdist_array[i]) + ",");
            else jsnstr.append("0,");
        }
        jsnstr.append("],\n");
        jsnstr.append("  type: 'bar',\n");
        jsnstr.append("  name: 'Insertion ',\n");
        jsnstr.append("};\n");

        jsnstr.append("var delbar = {\n");
        jsnstr.append("  x: [");
        for(int i = 0; i <= maxlen; ++i){
            jsnstr.append(std::to_string(i) + ",");
        }
        jsnstr.append("],\n");
        jsnstr.append("  y: [");
        for(int i = 0; i <= maxlen; ++i){
            if(i <= del_len_dist->vcmax) jsnstr.append(std::to_string(del_len_dist->vcdist_array[i]) + ",");
            else jsnstr.append("0,");
        }
        jsnstr.append("],\n");
        jsnstr.append("  type: 'bar',\n");
        jsnstr.append("  name: 'Deletion ',\n");
        jsnstr.append("};\n");

        jsnstr.append("var dinbar = {\n");
        jsnstr.append("  x: [");
        for(int i = 0; i <= maxlen; ++i){
            jsnstr.append(std::to_string(i) + ",");
        }
        jsnstr.append("],\n");
        jsnstr.append("  y: [");
        for(int i = 0; i <= maxlen; ++i){
            if(i <= din_len_dist->vcmax) jsnstr.append(std::to_string(din_len_dist->vcdist_array[i]) + ",");
            else jsnstr.append("0,");
        }
        jsnstr.append("],\n");
        jsnstr.append("  type: 'bar',\n");
        jsnstr.append("  name: 'Delins ',\n");
        jsnstr.append("};\n");

        jsnstr.append("var data = [insbar, delbar, dinbar];\n");

        jsnstr.append("var layout = {\n");
        jsnstr.append("  title: 'variant lengths of sequences',\n");
        jsnstr.append("  xaxis: {\n");
        jsnstr.append("    title: {\n");
        jsnstr.append("      text: 'variant length',\n");
        jsnstr.append("      standoff: 2,\n");
        jsnstr.append("    },\n");
        jsnstr.append("  },\n");
        jsnstr.append("  yaxis: {\n");
        jsnstr.append("    title: {\n");
        jsnstr.append("      text: 'sequence counts',\n");
        jsnstr.append("      standoff: 2,\n");
        jsnstr.append("    },\n");
        jsnstr.append("    automargin: true,\n");
        jsnstr.append("    tickmode: 'array',\n");
        jsnstr.append("    tickvals: [\n");
        int nstep = 20;
        int npcnt = 5;
        int64_t step = (totcnt + 5)/nstep;
        for(int i = 0; i < nstep; ++i){
            jsnstr.append(std::to_string(step * i) + ",");
        }
        jsnstr.append(std::to_string(totcnt) + ",");
        jsnstr.append("],\n");
        jsnstr.append("   ticktext: [\n");
        for(int i = 0; i <= nstep; ++i){
            jsnstr.append("'" +  std::to_string(step * i) + "(" + std::to_string(i*npcnt)+"%)',");
        }
        jsnstr.append("],\n");
        jsnstr.append(" },\n");
        jsnstr.append("}\n");

        jsnstr.append("var config = {\n");
        jsnstr.append("  toImageButtonOptions: {\n");
        jsnstr.append("    format: 'svg',\n");
        jsnstr.append("     filename: '" + divName + "_of_" + std::string(name) + "',\n");
        jsnstr.append("     height: " + std::to_string(opt->hmo.figh) + ",\n");
        jsnstr.append("     width: " + std::to_string(opt->hmo.figw) + ",\n");
        jsnstr.append("     scale: 1,\n");
        jsnstr.append("  }\n");
        jsnstr.append("};\n");

        jsnstr.append("Plotly.newPlot('plot_" + divName + "', data, layout, config);\n");
        ksprintf(s, "%s", jsnstr.c_str());
        ksprintf(s, "</script>\n");
    }
    // ins, del, delins, snv pos-wise lineplot
    {
        std::string subsect = "Variant counts along amplicon reference sequence";
        std::string divName = util::replace(subsect, " ", "_");
        std::string title = "'amplicon position variant counts'";
        ksprintf(s, "<div class='subsection_title'><a title='click to hide/show' onclick=showOrHide('%s')>%s</a></div>\n",
                    divName.c_str(), subsect.c_str());
        ksprintf(s, "<div id='%s'>\n", divName.c_str());
        ksprintf(s, "<div class='sub_section_tips'>Variant counts of each type on each position will be shown on mouse over.</div>\n");
        ksprintf(s, "<div class='figure' id='plot_%s'></div>\n", divName.c_str());
        ksprintf(s, "</div>\n");
        ksprintf(s, "\n<script type=\"text/javascript\">\n");

        std::string jsnstr;
        jsnstr.append("var totline = {\n");
        jsnstr.append("  x: [");
        for(int i = 0; i < rlen; ++i){
            jsnstr.append(std::to_string(i) + ",");
        }
        jsnstr.append("],\n");
        jsnstr.append("  y: [");
        for(int i = 0; i < rlen; ++i){
            int totvar = ins_pos_dist->vcdist_array[i] + del_pos_dist->vcdist_array[i] + 
                         din_pos_dist->vcdist_array[i] + snv_pos_dist->vcdist_array[i];
            jsnstr.append(std::to_string(totvar) + ",");
        }
        jsnstr.append("],\n");
        jsnstr.append("  type: 'line',\n");
        jsnstr.append("  name: 'AllVar ',\n");
        jsnstr.append("};\n");

        jsnstr.append("var insline = {\n");
        jsnstr.append("  x: [");
        for(int i = 0; i < rlen; ++i){
            jsnstr.append(std::to_string(i) + ",");
        }
        jsnstr.append("],\n");
        jsnstr.append("  y: [");
        for(int i = 0; i < rlen; ++i){
            if(i <= ins_pos_dist->vcmax) jsnstr.append(std::to_string(ins_pos_dist->vcdist_array[i]) + ",");
            else jsnstr.append("0,");
        }
        jsnstr.append("],\n");
        jsnstr.append("  type: 'line',\n");
        jsnstr.append("  name: 'Insertion ',\n");
        jsnstr.append("};\n");

        jsnstr.append("var delline = {\n");
        jsnstr.append("  x: [");
        for(int i = 0; i < rlen; ++i){
            jsnstr.append(std::to_string(i) + ",");
        }
        jsnstr.append("],\n");
        jsnstr.append("  y: [");
        for(int i = 0; i < rlen; ++i){
            if(i <= del_pos_dist->vcmax) jsnstr.append(std::to_string(del_pos_dist->vcdist_array[i]) + ",");
            else jsnstr.append("0,");
        }
        jsnstr.append("],\n");
        jsnstr.append("  type: 'line',\n");
        jsnstr.append("  name: 'Deletion ',\n");
        jsnstr.append("};\n");

        jsnstr.append("var dinline = {\n");
        jsnstr.append("  x: [");
        for(int i = 0; i < rlen; ++i){
            jsnstr.append(std::to_string(i) + ",");
        }
        jsnstr.append("],\n");
        jsnstr.append("  y: [");
        for(int i = 0; i < rlen; ++i){
            if(i <= din_pos_dist->vcmax) jsnstr.append(std::to_string(din_pos_dist->vcdist_array[i]) + ",");
            else jsnstr.append("0,");
        }
        jsnstr.append("],\n");
        jsnstr.append("  type: 'line',\n");
        jsnstr.append("  name: 'Delins ',\n");
        jsnstr.append("};\n");

        jsnstr.append("var snvline = {\n");
        jsnstr.append("  x: [");
        for(int i = 0; i < rlen; ++i){
            jsnstr.append(std::to_string(i) + ",");
        }
        jsnstr.append("],\n");
        jsnstr.append("  y: [");
        for(int i = 0; i < rlen; ++i){
            if(i <= snv_pos_dist->vcmax) jsnstr.append(std::to_string(snv_pos_dist->vcdist_array[i]) + ",");
            else jsnstr.append("0,");
        }
        jsnstr.append("],\n");
        jsnstr.append("  type: 'line',\n");
        jsnstr.append("  name: 'SNV ',\n");
        jsnstr.append("};\n");

        jsnstr.append("var data = [totline, insline, delline, dinline, snvline];\n");

        jsnstr.append("var layout = {\n");
        jsnstr.append("  title: 'variant counts along amplicon',\n");
        jsnstr.append("  xaxis: {\n");
        jsnstr.append("    tickmode: 'array',\n");
        jsnstr.append("    tickvals: " + xtvs + ",\n");
        jsnstr.append("    ticktext: " + xtts + ",\n");
        jsnstr.append("    title: {\n");
        jsnstr.append("      text: 'relative position to sgRNA beg',\n");
        jsnstr.append("      standoff: 2,\n");
        jsnstr.append("    },\n");
        jsnstr.append("  },\n");
        jsnstr.append("  yaxis: {\n");
        jsnstr.append("    title: {\n");
        jsnstr.append("      text: 'variant counts',\n");
        jsnstr.append("      standoff: 2,\n");
        jsnstr.append("    },\n");
        jsnstr.append("    automargin: true,\n");
        jsnstr.append("    tickmode: 'array',\n");
        jsnstr.append("    tickvals: [\n");
        int nstep = 20;
        int npcnt = 5;
        int64_t step = (totcnt + 5)/nstep;
        for(int i = 0; i < nstep; ++i){
            jsnstr.append(std::to_string(step * i) + ",");
        }
        jsnstr.append(std::to_string(totcnt) + ",");
        jsnstr.append("],\n");
        jsnstr.append("   ticktext: [\n");
        for(int i = 0; i <= nstep; ++i){
            jsnstr.append("'" +  std::to_string(step * i) + "(" + std::to_string(i*npcnt)+"%)',");
        }
        jsnstr.append("],\n");
        jsnstr.append(" },\n");
        jsnstr.append("}\n");

        jsnstr.append("var config = {\n");
        jsnstr.append("  toImageButtonOptions: {\n");
        jsnstr.append("    format: 'svg',\n");
        jsnstr.append("     filename: '" + divName + "_of_" + std::string(name) + "',\n");
        jsnstr.append("     height: " + std::to_string(opt->hmo.figh) + ",\n");
        jsnstr.append("     width: " + std::to_string(opt->hmo.figw) + ",\n");
        jsnstr.append("     scale: 1,\n");
        jsnstr.append("  }\n");
        jsnstr.append("};\n");

        jsnstr.append("Plotly.newPlot('plot_" + divName + "', data, layout, config);\n");
        ksprintf(s, "%s", jsnstr.c_str());
        ksprintf(s, "</script>\n");
    }
    // ins, del, delins mean length barplot
    {
        std::string subsect = "Mean variant length along amplicon reference sequence";
        std::string divName = util::replace(subsect, " ", "_");
        std::string title = "'amplicon position mean variant length counts'";
        ksprintf(s, "<div class='subsection_title'><a title='click to hide/show' onclick=showOrHide('%s')>%s</a></div>\n",
                    divName.c_str(), subsect.c_str());
        ksprintf(s, "<div id='%s'>\n", divName.c_str());
        ksprintf(s, "<div class='sub_section_tips'>Mean variant length of each type of variants on amplicon reference positions will be shown on mouse over.</div>\n");
        ksprintf(s, "<div class='figure' id='plot_%s'></div>\n", divName.c_str());
        ksprintf(s, "</div>\n");
        ksprintf(s, "\n<script type=\"text/javascript\">\n");

        double mvlen = 0.0;
        std::string jsnstr;
        
        jsnstr.append("var insline = {\n");
        jsnstr.append("  x: [");
        for(int i = 0; i < rlen; ++i){
            jsnstr.append(std::to_string(i) + ",");
        }
        jsnstr.append("],\n");
        jsnstr.append("  y: [");
        for(int i = 0; i < rlen; ++i){
            mvlen =  0;
            if(ins_pos_dist->vcdist_array[i]){
                mvlen = (double)ins_lps_dist->vcdist_array[i]/(double)ins_pos_dist->vcdist_array[i];
            }
            jsnstr.append(std::to_string(mvlen) + ",");
        }
        jsnstr.append("],\n");
        jsnstr.append("  type: 'line',\n");
        jsnstr.append("  name: 'Insertion ',\n");
        jsnstr.append("};\n");

        jsnstr.append("var delline = {\n");
        jsnstr.append("  x: [");
        for(int i = 0; i < rlen; ++i){
            jsnstr.append(std::to_string(i) + ",");
        }
        jsnstr.append("],\n");
        jsnstr.append("  y: [");
        for(int i = 0; i < rlen; ++i){
            mvlen =  0;
            if(del_pos_dist->vcdist_array[i]){
                mvlen = (double)del_lps_dist->vcdist_array[i]/(double)del_pos_dist->vcdist_array[i];
            }
            jsnstr.append(std::to_string(mvlen) + ",");
        }
        jsnstr.append("],\n");
        jsnstr.append("  type: 'line',\n");
        jsnstr.append("  name: 'Deletion ',\n");
        jsnstr.append("};\n");

        jsnstr.append("var dinline = {\n");
        jsnstr.append("  x: [");
        for(int i = 0; i < rlen; ++i){
            jsnstr.append(std::to_string(i) + ",");
        }
        jsnstr.append("],\n");
        jsnstr.append("  y: [");
        for(int i = 0; i < rlen; ++i){
            mvlen =  0;
            if(din_pos_dist->vcdist_array[i]){
                mvlen = (double)din_lps_dist->vcdist_array[i]/(double)din_pos_dist->vcdist_array[i];
            }
            jsnstr.append(std::to_string(mvlen) + ",");
        }
        jsnstr.append("],\n");
        jsnstr.append("  type: 'line',\n");
        jsnstr.append("  name: 'Delins ',\n");
        jsnstr.append("};\n");

        jsnstr.append("var data = [insline, delline, dinline];\n");

        jsnstr.append("var layout = {\n");
        jsnstr.append("  title: 'variant mean length along amplicon',\n");
        jsnstr.append("  xaxis: {\n");
        jsnstr.append("    tickmode: 'array',\n");
        jsnstr.append("    tickvals: " + xtvs + ",\n");
        jsnstr.append("    ticktext: " + xtts + ",\n");
        jsnstr.append("    title: {\n");
        jsnstr.append("      text: 'relative position to sgRNA beg',\n");
        jsnstr.append("      standoff: 2,\n");
        jsnstr.append("    },\n");
        jsnstr.append("  },\n");
        jsnstr.append("  yaxis: {\n");
        jsnstr.append("    title: {\n");
        jsnstr.append("      text: 'mean variant length',\n");
        jsnstr.append("      standoff: 2,\n");
        jsnstr.append("    },\n");
        jsnstr.append("    automargin: true,\n");
        jsnstr.append("  },\n");
        jsnstr.append("};\n");

        jsnstr.append("var config = {\n");
        jsnstr.append("  toImageButtonOptions: {\n");
        jsnstr.append("    format: 'svg',\n");
        jsnstr.append("     filename: '" + divName + "_of_" + std::string(name) + "',\n");
        jsnstr.append("     height: " + std::to_string(opt->hmo.figh) + ",\n");
        jsnstr.append("     width: " + std::to_string(opt->hmo.figw) + ",\n");
        jsnstr.append("     scale: 1,\n");
        jsnstr.append("  }\n");
        jsnstr.append("};\n");


        jsnstr.append("Plotly.newPlot('plot_" + divName + "', data, layout, config);\n");
        ksprintf(s, "%s", jsnstr.c_str());
        ksprintf(s, "</script>\n");
    }
}
