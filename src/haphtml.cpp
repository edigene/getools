#include "hapcnt.h"

void hap_opt_t::printHeader(kstring_t* s){
    hmo.init(jscdn);
    // beg <html>
    ksprintf(s, "<html>\n");
    // beg <head>
    ksprintf(s, "<head>\n");
    // beg <meta>
    ksprintf(s, "<meta http-equiv=\"content-type\" content=\"text/html;charset=utf-8\" />\n");
    // end <meta>
    // beg <title>
    ksprintf(s, "<title>haplotype snv report at %s </title>\n", util::currentTime().c_str());
    // end <title>
    hmo.printJsAndCSSRes(s);
    // beg <script>
    ksprintf(s, "<script type=\"text/javascript\">\n");
    ksprintf(s, "    function showOrHide(divname) {\n");
    ksprintf(s, "        div = document.getElementById(divname);\n");
    ksprintf(s, "        if(div.style.display == 'none')\n");
    ksprintf(s, "            div.style.display = 'block';\n");
    ksprintf(s, "        else\n");
    ksprintf(s, "            div.style.display = 'none';\n");
    ksprintf(s, "    }\n");
    ksprintf(s, "    function toggle(targetid){\n");
    ksprintf(s, "       if (document.getElementById){\n");
    ksprintf(s, "           target=document.getElementById(targetid)\n");
    ksprintf(s, "               if (target.style.display=='table-row'){\n");
    ksprintf(s, "                   target.style.display='none';\n");
    ksprintf(s, "               } else {\n");
    ksprintf(s, "                   target.style.display='table-row';\n");
    ksprintf(s, "               }\n");
    ksprintf(s, "       }\n");
    ksprintf(s, "   }\n");
    ksprintf(s, "   function toggle_target_list(targetid){\n");
    ksprintf(s, "       if (document.getElementById){\n");
    ksprintf(s, "           target=document.getElementById(targetid);\n");
    ksprintf(s, "               if (target.style.display=='block'){\n");
    ksprintf(s, "                   target.style.display='none';\n");
    ksprintf(s, "                   document.getElementById('target_view_btn').value='view';\n");
    ksprintf(s, "               } else {\n");
    ksprintf(s, "                   document.getElementById('target_view_btn').value='hide';\n");
    ksprintf(s, "                   target.style.display='block';\n");
    ksprintf(s, "               }\n");
    ksprintf(s, "       }\n");
    ksprintf(s, "   }\n");
    ksprintf(s, "</script>\n");
    // end <script>
    // beg <style>
    ksprintf(s, "<style type=\"text/css\">\n");
    ksprintf(s, "td {border:1px solid #dddddd;padding:5px;font-size:12px;}\n");
    ksprintf(s, "table {border:1px solid #999999;padding:2x;border-collapse:collapse; width:1000px}\n");
    ksprintf(s, ".col1 {width:300px; font-weight:bold;}\n");
    ksprintf(s, "img {padding:30px;}\n");
    ksprintf(s, "#menu {font-family:Consolas, 'Liberation Mono', Menlo, Courier, monospace;}\n");
    ksprintf(s, "#menu a {color:#0366d6; font-size:18px;font-weight:600;line-height:28px;text-decoration:none;");
    ksprintf(s, "font-family:-apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol'}\n");
    ksprintf(s, "a:visited {color: #999999}\n");
    ksprintf(s, ".alignleft {text-align:left;}\n");
    ksprintf(s, ".alignright {text-align:right;}\n");
    ksprintf(s, ".figure {width:%dpx;height:%dpx;}\n", hmo.figw, hmo.figh);
    ksprintf(s, ".header {color:#ffffff;padding:1px;height:20px;background:#000000;}\n");
    ksprintf(s, ".section_title {color:#ffffff;font-size:20px;padding:5px;text-align:left;background:#000099; margin-top:10px;}\n");
    ksprintf(s, ".subsection_title {font-size:16px;padding:5px;margin-top:10px;text-align:left;color:#000099}\n");
    ksprintf(s, "#container {text-align:center;padding:3px 3px 3px 10px;font-family:Arail,'Liberation Mono', Menlo, Courier, monospace;}\n");
    ksprintf(s, ".menu_item {text-align:left;padding-top:5px;font-size:18px;}\n");
    ksprintf(s, ".highlight {text-align:left;padding-top:30px;padding-bottom:30px;font-size:20px;line-height:35px;}\n");
    ksprintf(s, "#helper {text-align:left;border:1px dotted #fafafa;color:#777777;font-size:12px;}\n");
    ksprintf(s, "#footer {text-align:left;padding:15px;color:#ffffff;font-size:10px;background:#000099;");
    ksprintf(s, "font-family:Arail,'Liberation Mono', Menlo, Courier, monospace;}\n");
    ksprintf(s, ".sub_section_tips {color:#999999;font-size:10px;padding-left:5px;padding-bottom:3px;}\n");
    ksprintf(s, "</style>\n");
    // end <style>
    ksprintf(s, "</head>\n");
    // end <head>
    // beg <body>
    ksprintf(s, "<body><div id='container'>\n");
}


void hap_opt_t::printFooter(kstring_t* s){
    // end <body><container...
    ksprintf(s, "</div>\n</body>\n");
    // begin footer
    ksprintf(s, "<div id='footer'> ");
    ksprintf(s, "<p>command: getools hapcnt </p>");
    ksprintf(s, "</div>\n");
    // end footer
    // end html
    ksprintf(s, "</html>\n");
}

void hap_opt_t::printBody(kstring_t* s){
    ksprintf(s, "<div class='section_div'>\n");
    ksprintf(s, "<div class='subsection_title'><a title='click to hide/show' onclick=showOrHide('hapgtn')>Top%d HapSNV</a></div>\n", gtn);
    ksprintf(s, "<div id='hapgtn'>\n");
    ksprintf(s, "<div class='figure' id='hapgtn'></div>\n");
    ksprintf(s, "</div>\n");
    ksprintf(s, "\n<script type=\"text/javascript\">\n");
    
    ksprintf(s, "data = [\n");
    hap_rec_t* hr = NULL;
    float maxf = .0, ttf = .0;
    for(int i = 0; i < gtn && i < int(hapl.size()); ++i){
        hr = hapl[i];
        ksprintf(s, "{\n");
        ksprintf(s, "  type: 'scatter',\n");
        ksprintf(s, "  mode: 'lines+markers',\n");
        ksprintf(s, "  name: '%s->%s',\n", hr->refs.c_str(), hr->muts.c_str());
        ksprintf(s, "  text: [");
        for(size_t xp = 0; xp < hr->pos.size(); ++xp) ksprintf(s, "'%c->%c',", hr->refs[xp], hr->muts[xp]);
        ksprintf(s, "],\n");
        ksprintf(s, "  x: [");
        for(size_t xp = 0; xp < hr->pos.size(); ++xp) ksprintf(s, "%d,", hr->pos[xp]+1);
        ksprintf(s, "],\n");
        ksprintf(s, "  y: [");
        for(size_t xp = 0; xp < hr->pos.size(); ++xp){
            ttf = 100*(double)hr->cnt/(double)ttr;
            if(ttf > maxf) maxf = ttf;
            ksprintf(s, "%f,", ttf);
        }
        ksprintf(s, "],\n");
        ksprintf(s, "},\n");
    }
    ksprintf(s, "]\n");
    
    ksprintf(s, "layout = {\n");
    ksprintf(s, "  xaxis: {");
    ksprintf(s, "    title: {text: 'Reference Positions(1based)',},\n");
    ksprintf(s, "    range: [%d, %d],\n", beg, end);
    ksprintf(s, "  },\n");
    ksprintf(s, "  yaxis: {");
    ksprintf(s, "    title: {text: 'Haplotype Frequency(%%)',},\n");
    ksprintf(s, "    range: [%f, %f],\n", .0, maxf + 5);
    ksprintf(s, "  },\n");
    ksprintf(s, "  title: 'top%d hap snv',\n", gtn);
    ksprintf(s, "  hovermode:'closest',");
    ksprintf(s, "}\n");

    ksprintf(s, "var config = {\n");
    ksprintf(s, "  toImageButtonOptions: {\n");
    ksprintf(s, "    format: 'svg',\n");
    ksprintf(s, "     filename: 'topn.hap.snv',\n");
    ksprintf(s, "     height: 600,\n");
    ksprintf(s, "     scale: 1,\n");
    ksprintf(s, "  }\n");
    ksprintf(s, "};\n");

    ksprintf(s, "Plotly.newPlot('hapgtn', data, layout, config)\n");
    ksprintf(s, "</script>\n");
    ksprintf(s, "</div>\n");
}


void hap_opt_t::hap2html(){
    kstring_t ks = {0, 0, 0};
    printHeader(&ks);
    printBody(&ks);
    printFooter(&ks);
    FILE* fp = fopen(outhtml.c_str(), "w");
    fwrite(ks.s, ks.l, sizeof(char), fp);
    fclose(fp);
}
