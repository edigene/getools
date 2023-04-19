#include "splitresult.h"

SplitResult::SplitResult(Options *opt){
    mOptions = opt;
    mDropReads = 0;
    mSplitFailedReads = 0;
    mQCFailedReads = 0;
    mTotalReads = 0;
    mDropCount8Reasons.resize(MATCH_FAIL_RCNT, 0);
    mBaseCount.resize(mOptions->samples.size(), 0);
    mSplitCount.resize(mOptions->samples.size(), 0);
    mSplitRead1.resize(mOptions->samples.size(), NULL);
    mSplitRead2.resize(mOptions->samples.size(), NULL);
    mQCRead1.resize(mOptions->samples.size(), NULL);
    mQCRead2.resize(mOptions->samples.size(), NULL);
    for(size_t i = 0; i < mOptions->samples.size(); ++i){
        mSplitRead1[i] = (kstring_t*)calloc(1, sizeof(kstring_t));
        mSplitRead2[i] = (kstring_t*)calloc(1, sizeof(kstring_t));
        mQCRead1[i] = new QcStat(mOptions);
        mQCRead2[i] = new QcStat(mOptions);
    }
    mFRCount.resize(mOptions->samples.size(), 0);
    mRFCount.resize(mOptions->samples.size(), 0);
    mMMCount.resize(mOptions->maxmm * 2 + 1);
    mTotalFR = 0;
    mTotalRF = 0;
    mSummarized = false;
    if(opt->in2.size()) mIsPE = true;
}

SplitResult::~SplitResult(){
    for(size_t i = 0; i < mOptions->samples.size(); ++i){
        if(mSplitRead1[i]->m) free(mSplitRead1[i]->s);
        if(mSplitRead2[i]->m) free(mSplitRead2[i]->s);
        if(mQCRead1[i]) delete mQCRead1[i];
        if(mQCRead2[i]) delete mQCRead2[i];
    }
}

void SplitResult::init(){
    for(size_t i = 0; i < mOptions->samples.size(); ++i){
        mSplitRead1[i]->l = 0;
        mSplitRead2[i]->l = 0;
    }
}

void SplitResult::addDropRead(krec1_t* r1, krec1_t* r2, bool isQCFail){
    if(!mOptions->dosccal){
        r1->comment.l = 0;
        if(r2) r2->comment.l = 0;
    }
    ++mDropReads;
    ++mTotalReads;
    if(r1->dwhy != MATCH_FAIL_UNKNOWN) ++mDropCount8Reasons[r1->dwhy];
    else if(r2->dwhy != MATCH_FAIL_UNKNOWN) ++mDropCount8Reasons[r2->dwhy];
    else ++mDropCount8Reasons[r1->dwhy];
    if(isQCFail) ++mQCFailedReads;
    else ++mSplitFailedReads;
    r1->sample = mOptions->dropidx;
    mQCRead1[mOptions->dropidx]->statRead(r1);
    if(r2){
        mQCRead2[mOptions->dropidx]->statRead(r2);
        r2->sample = mOptions->dropidx;
    }
    if(mOptions->outspl){
        kputc('@', mSplitRead1[mOptions->dropidx]);
        kputs(MATCH_FAIL_REASONS[r1->dwhy], mSplitRead1[mOptions->dropidx]);
        kputc('_', mSplitRead1[mOptions->dropidx]);
        kputsn(r1->name.s, r1->name.l, mSplitRead1[mOptions->dropidx]);
        if(mOptions->dosccal && r1->comment.l){
            kputc(' ', mSplitRead1[mOptions->dropidx]);
            kputsn(r1->comment.s, r1->comment.l, mSplitRead1[mOptions->dropidx]);
        }
        kputc('\n', mSplitRead1[mOptions->dropidx]);
        kputsn(r1->seq.s, r1->seq.l, mSplitRead1[mOptions->dropidx]);
        kputsn("\n+\n", 3, mSplitRead1[mOptions->dropidx]);
        kputsn(r1->qual.s, r1->qual.l, mSplitRead1[mOptions->dropidx]);
        kputc('\n', mSplitRead1[mOptions->dropidx]);
        if(r2){
            kputc('@', mSplitRead2[mOptions->dropidx]);
            kputs(MATCH_FAIL_REASONS[r2->dwhy], mSplitRead2[mOptions->dropidx]);
            kputc('_', mSplitRead2[mOptions->dropidx]);
            kputsn(r2->name.s, r2->name.l, mSplitRead2[mOptions->dropidx]);
            if(mOptions->dosccal && r2->comment.l){
                kputc(' ', mSplitRead2[mOptions->dropidx]);
                kputsn(r2->comment.s, r2->comment.l, mSplitRead2[mOptions->dropidx]);
            }
            kputc('\n', mSplitRead2[mOptions->dropidx]);
            kputsn(r2->seq.s, r2->seq.l, mSplitRead2[mOptions->dropidx]);
            kputsn("\n+\n", 3, mSplitRead2[mOptions->dropidx]);
            kputsn(r2->qual.s, r2->qual.l, mSplitRead2[mOptions->dropidx]);
            kputc('\n', mSplitRead2[mOptions->dropidx]);
        }
    }
}

void SplitResult::addSplitRead(krec1_t* r1, krec1_t* r2, int m){
    if(!mOptions->dosccal){
        r1->comment.l = 0;
        if(r2) r2->comment.l = 0;
    }
    mMMCount[m] += 1;
    mSplitCount[r1->sample] += 1;
    ++mTotalReads;
    mBaseCount[r1->sample] += r1->seq.l-r1->off;
    mQCRead1[r1->sample]->statRead(r1, r1->off);
    if(r2){
        mBaseCount[r2->sample] += r2->seq.l-r2->off;
        mQCRead2[r2->sample]->statRead(r2, r2->off);
    }
    if(r2){
        if(r1->strand)  mFRCount[r1->sample] += 1;
        else{
            mRFCount[r1->sample] += 1;
            if(mOptions->adjfr) std::swap(r1, r2);
        }
    }
    if(mOptions->outspl){
        if(mOptions->usesem){// only store the one matched in r1
            krec1_t* kr = r1;
            if(r2->match) kr = r2;
            kputc('@', mSplitRead1[kr->sample]);
            kputsn(kr->name.s, kr->name.l, mSplitRead1[kr->sample]);
            if(mOptions->dosccal && r1->comment.l){// keep all comments of read1/2
                kputc(' ', mSplitRead1[kr->sample]);
                kputsn(r1->comment.s, r1->comment.l, mSplitRead1[kr->sample]);
            }
            kputc('\n', mSplitRead1[kr->sample]);
            kputsn(kr->seq.s+kr->off, kr->seq.l-kr->off, mSplitRead1[kr->sample]);
            ksprintf(mSplitRead1[kr->sample], "\n%c\n", KSTRANDS[kr->strand]);
            kputsn(kr->qual.s+kr->off, kr->qual.l-kr->off, mSplitRead1[kr->sample]);
            kputc('\n', mSplitRead1[kr->sample]);
        }else{
            kputc('@', mSplitRead1[r1->sample]);
            kputsn(r1->name.s, r1->name.l, mSplitRead1[r1->sample]);
            if(mOptions->dosccal && r1->comment.l){
                kputc(' ', mSplitRead1[r1->sample]);
                kputsn(r1->comment.s, r1->comment.l, mSplitRead1[r1->sample]);
            }
            kputc('\n', mSplitRead1[r1->sample]);
            kputsn(r1->seq.s+r1->off, r1->seq.l-r1->off, mSplitRead1[r1->sample]);
            if(mOptions->ostrand){
                ksprintf(mSplitRead1[r1->sample], "\n%c\n", KSTRANDS[r1->strand]);
            }else{
                kputsn("\n+\n", 3, mSplitRead1[r1->sample]);
            }
            kputsn(r1->qual.s+r1->off, r1->qual.l-r1->off, mSplitRead1[r1->sample]);
            kputc('\n', mSplitRead1[r1->sample]);
            if(r2){
                kputc('@', mSplitRead2[r2->sample]);
                kputsn(r2->name.s, r2->name.l, mSplitRead2[r2->sample]);
                if(mOptions->dosccal && r2->comment.l){
                    kputc(' ', mSplitRead2[r2->sample]);
                    kputsn(r2->comment.s, r2->comment.l, mSplitRead2[r2->sample]);
                }
                kputc('\n', mSplitRead2[r2->sample]);
                kputsn(r2->seq.s+r2->off, r2->seq.l-r2->off, mSplitRead2[r2->sample]);
                if(mOptions->ostrand){
                    ksprintf(mSplitRead2[r2->sample], "\n%c\n", KSTRANDS[r2->strand]);
                }else{
                    kputsn("\n+\n", 3, mSplitRead2[r2->sample]);
                }
                kputsn(r2->qual.s+r2->off, r2->qual.l-r2->off, mSplitRead2[r2->sample]);
                kputc('\n', mSplitRead2[r2->sample]);
            }
        }
    }
}

SplitResult* SplitResult::mergeResult(const std::vector<SplitResult*>& vresults){
    if(vresults.empty()) return NULL;
    SplitResult* mr = new SplitResult(vresults[0]->mOptions);
    // merge general
    for(uint32_t i =  0; i < vresults.size(); ++i){
        mr->mDropReads += vresults[i]->mDropReads;
        mr->mQCFailedReads += vresults[i]->mQCFailedReads;
        mr->mSplitFailedReads += vresults[i]->mSplitFailedReads;
        mr->mTotalReads += vresults[i]->mTotalReads;
        for(uint32_t s = 0; s < mr->mOptions->samples.size(); ++s){
            mr->mBaseCount[s] += vresults[i]->mBaseCount[s];
            mr->mSplitCount[s] += vresults[i]->mSplitCount[s];
            mr->mFRCount[s] += vresults[i]->mFRCount[s];
            mr->mRFCount[s] += vresults[i]->mRFCount[s];
        }
        for(uint32_t m = 0; m < mr->mMMCount.size(); ++m){
            mr->mMMCount[m] += vresults[i]->mMMCount[m];
        }
        for(int j = 0; j < MATCH_FAIL_RCNT; ++j){
            mr->mDropCount8Reasons[j] += vresults[i]->mDropCount8Reasons[j];
        }
    }
    // merge qc
    std::vector<QcStat*> qs1, qs2;
    for(uint32_t s = 0; s < mr->mOptions->samples.size(); ++s){
        qs1.clear(); qs2.clear();
        for(uint32_t i =  0; i < vresults.size(); ++i){
            qs1.push_back(vresults[i]->mQCRead1[s]);
            qs2.push_back(vresults[i]->mQCRead2[s]);
        }
        delete mr->mQCRead1[s]; 
        delete mr->mQCRead2[s];
        mr->mQCRead1[s] = QcStat::merge(qs1);
        if(mr->mOptions->in2.size()) mr->mQCRead2[s] = QcStat::merge(qs2);
        else mr->mQCRead2[s] = NULL;
    }
    return mr;
}

void SplitResult::summary(){
    if(mSummarized) return;
    for(uint32_t s = 0; s < mOptions->samples.size(); ++s){
        mTotalFR += mFRCount[s];
        mTotalRF += mRFCount[s];
    }
    mGotRate = 1 - (double)mDropReads/(double)mTotalReads;
    mSummarized = true;
}

void SplitResult::reportHTML(kstring_t* s){
    if(!mSummarized) summary();
    // mm and fr/rf and drop reasons
    {
        std::string subsect = "Split summary";
        std::string divName = util::replace(subsect, " ", "_");
        std::string title = "'split got rate (" + std::to_string(mGotRate*100) + "%)'";
        ksprintf(s, "<div class='subsection_title'><a title='click to hide/show' onclick=showOrHide('%s')>%s</a></div>\n", divName.c_str(), subsect.c_str());
        ksprintf(s, "<div id='%s'>\n", divName.c_str());
        ksprintf(s, "<div class='sub_section_tips'>Percentage of each type will be shown on mouse over.</div>\n");
        ksprintf(s, "<div class='figure' id='plot_%s'></div>\n", divName.c_str());
        ksprintf(s, "</div>\n");
        ksprintf(s, "\n<script type=\"text/javascript\">\n");

        std::string jsnstr;
        jsnstr.append("var frpie = {\n");
        jsnstr.append("  values: [");
        jsnstr.append(std::to_string(mTotalFR) + ",");
        jsnstr.append(std::to_string(mTotalRF) + "],\n");
        jsnstr.append("  labels: ['FR', 'RF'],\n");
        jsnstr.append("  type: 'pie',\n");
        jsnstr.append("  title: 'Split Strand Count',\n");
        jsnstr.append("  textinfo: 'label',\n");
        jsnstr.append("  textposition: 'inside',\n");
        jsnstr.append("  hoverinfo: 'label+value+percent',\n");
        jsnstr.append("  insidetextorientation: 'radial',\n");
        jsnstr.append("  domain: {column: 0, row: 0},\n");
        jsnstr.append("};\n");

        jsnstr.append("var mmpie = {\n");
        jsnstr.append("  values: [");
        for(size_t i = 0; i < mMMCount.size(); ++i){
            jsnstr.append(std::to_string(mMMCount[i]) + ",");
        }
        jsnstr.append("],\n");
        jsnstr.append("  labels: [");
        for(size_t i = 0; i < mMMCount.size(); ++i){
            jsnstr.append("'Mismatch" + std::to_string(i) + "',");
        }
        jsnstr.append("],\n");
        jsnstr.append("  type: 'pie',\n");
        jsnstr.append("  title: 'Split Mismatch Count',\n");
        jsnstr.append("  textinfo: 'label',\n");
        jsnstr.append("  textposition: 'inside',\n");
        jsnstr.append("  hoverinfo: 'label+value+percent',\n");
        jsnstr.append("  insidetextorientation: 'radial',\n");
        jsnstr.append("  domain: {column: 1, row: 0},\n");
        jsnstr.append("};\n");

        jsnstr.append("var droppie = {\n");
        jsnstr.append("  values: [");
        for(int i = 0; i < MATCH_FAIL_RCNT; ++i){
            jsnstr.append(std::to_string(mDropCount8Reasons[i]) + ",");
        }
        jsnstr.append("],\n");
        jsnstr.append("  labels: [");
        for(size_t i = 0; i < MATCH_FAIL_RCNT; ++i){
            jsnstr.append("'" + std::string(MATCH_FAIL_REASONS[i]) + "',");
        }
        jsnstr.append("],\n");
        jsnstr.append("  type: 'pie',\n");
        jsnstr.append("  title: 'Split Failure Count',\n");
        jsnstr.append("  textinfo: 'label',\n");
        jsnstr.append("  textposition: 'inside',\n");
        jsnstr.append("  hoverinfo: 'label+value+percent',\n");
        jsnstr.append("  insidetextorientation: 'radial',\n");
        jsnstr.append("  domain: {column: 2, row: 0},\n");
        jsnstr.append("};\n");

        jsnstr.append("var layout = {\n");
        jsnstr.append("  title:" + title + ",\n");
        jsnstr.append("  grid: {rows: 1, columns: 3},\n");
        jsnstr.append("  showlegend: false,\n");
        jsnstr.append("};\n");

        jsnstr.append("var data = [frpie, mmpie, droppie];\n");

        jsnstr.append("var config = {\n");
        jsnstr.append("  toImageButtonOptions: {\n");
        jsnstr.append("    format: 'svg',\n");
        jsnstr.append("     filename: '" + divName + "',\n");
        jsnstr.append("     height: " + std::to_string(mOptions->hmo.figh) + ",\n");
        jsnstr.append("     width: " + std::to_string(mOptions->hmo.figw) + ",\n");
        jsnstr.append("     scale: 1,\n");
        jsnstr.append("  }\n");
        jsnstr.append("};\n");

        jsnstr.append("Plotly.newPlot('plot_" + divName + "', data, layout, config);\n");
        ksprintf(s, "%s", jsnstr.c_str());
        ksprintf(s, "</script>\n");
    }
    // amplican-wise split stat
    {
        std::string subsect = "Read assigned to each amplicons";
        std::string divName = util::replace(subsect, " ", "_");
        std::string title = "'reads assigned to each amplicons'";
        if(mIsPE) title = "'read pairs assigned to each amplicons'";
        ksprintf(s, "<div class='subsection_title'><a title='click to hide/show' onclick=showOrHide('%s')>%s</a></div>\n", divName.c_str(), subsect.c_str());
        ksprintf(s, "<div id='%s'>\n", divName.c_str());
        ksprintf(s, "<div class='sub_section_tips'>Value of each sample will be shown on mouse over.</div>\n");
        ksprintf(s, "<div class='figure' id='plot_%s'></div>\n", divName.c_str());
        ksprintf(s, "</div>\n");
        ksprintf(s, "\n<script type=\"text/javascript\">\n");

        std::string jsnstr;

        jsnstr.append("var splbar = {\n");
        jsnstr.append("  x: [");
        for(size_t i = 0; i < mSplitCount.size()-1; ++i){
            jsnstr.append("'<a href=\"" + mOptions->hmo.subdir + "/" + mOptions->samples[i] + ".html\">" + mOptions->samples[i] + "</a>',");
        }
        jsnstr.append("],\n");
        jsnstr.append("  y: [");
        for(size_t i = 0; i < mSplitCount.size()-1; ++i){
            jsnstr.append(std::to_string(mSplitCount[i]) + ",");
        }
        jsnstr.append("],\n");
        jsnstr.append("  type: 'bar',\n");
        jsnstr.append("  name: 'Split Count',\n");
        jsnstr.append("  hovertext: [");
        for(size_t i = 0; i < mSplitCount.size()-1; ++i){
            jsnstr.append("'" + std::to_string((double)mSplitCount[i]/(double)mTotalReads*100) + "%',");
        }
        jsnstr.append("],\n");
        jsnstr.append("}\n");

        jsnstr.append("var layout = {\n");
        jsnstr.append("  title:" + title + ",\n");
        jsnstr.append("  showlegend: false,\n");
        jsnstr.append("  xaxis:{\n");
        jsnstr.append("    automargin: true,\n");
        jsnstr.append("    dtick: 1,\n");
        jsnstr.append("    title: {\n");
        jsnstr.append("      text: 'amplicon',\n");
        jsnstr.append("      standoff: 20,\n");
        jsnstr.append("    },\n");
        jsnstr.append("},\n");
        jsnstr.append("  yaxis:{\n");
        jsnstr.append("    automargin: true,\n");
        jsnstr.append("    title: {\n");
        if(mIsPE) jsnstr.append("      text: 'read pairs',\n");
        else jsnstr.append("      text: 'reads',\n");
        jsnstr.append("      standoff: 20,\n");
        jsnstr.append("    },\n");
        jsnstr.append("  },\n");
        jsnstr.append("};\n");

        jsnstr.append("var data = [splbar];\n");

        jsnstr.append("var config = {\n");
        jsnstr.append("  toImageButtonOptions: {\n");
        jsnstr.append("    format: 'svg',\n");
        jsnstr.append("     filename: '" + divName + "',\n");
        jsnstr.append("     height: " + std::to_string(mOptions->hmo.figh) + ",\n");
        jsnstr.append("     width: " + std::to_string(mOptions->hmo.figw) + ",\n");
        jsnstr.append("     scale: 1,\n");
        jsnstr.append("  }\n");
        jsnstr.append("};\n");

        jsnstr.append("Plotly.newPlot('plot_" + divName + "', data, layout, config);\n");
        ksprintf(s, "%s", jsnstr.c_str());
        ksprintf(s, "</script>\n");
    }
}

void SplitResult::reportJSON(kstring_t* s, const char* dh, const char* dm){
    if(!mSummarized) summary();
    // Summary
    ksprintf(s, "%s\"Summary\": {\n", dh);
    ksprintf(s, "%s%s\"ReadsInput\": %lld,\n", dh, dm, mTotalReads);
    ksprintf(s, "%s%s\"ReadsGot\": %lld,\n", dh, dm, mTotalReads - mDropReads);
    ksprintf(s, "%s%s\"ReadsDrop\": %lld,\n", dh, dm, mDropReads);
    ksprintf(s, "%s%s\"ReadsFailQC\": %lld,\n", dh, dm, mQCFailedReads);
    ksprintf(s, "%s%s\"ReadsFailSplit\": %lld,\n", dh, dm, mSplitFailedReads);
    ksprintf(s, "%s%s\"DropRate\": %.2lf,\n", dh, dm, 1 - mGotRate);
    ksprintf(s, "%s%s\"GotRate\": %.2lf,\n", dh, dm, mGotRate);
    ksprintf(s, "%s%s\"FRGot\": %lld,\n", dh, dm, mTotalFR);
    ksprintf(s, "%s%s\"RFGot\": %lld,\n", dh, dm, mTotalRF);
    for(size_t i = 0; i < mMMCount.size(); ++i){
        ksprintf(s, "%s%s\"Mismatch%lu\":%lld,\n", dh, dm, i, mMMCount[i]);
    }
    s->l -= 2;
    ksprintf(s, "\n%s},\n", dh);
    // Fail Reasons
    ksprintf(s, "%s\"DropReasons\": {\n", dh);
    for(int i = 0 ; i < MATCH_FAIL_RCNT; ++i){
        ksprintf(s, "%s%s\"%s\":%lld,\n", dh, dm, MATCH_FAIL_REASONS[i],  mDropCount8Reasons[i]);
    }
    s->l -= 2;
    ksprintf(s, "\n%s},\n", dh);
    // LibraryReadsGot
    ksprintf(s, "%s\"LibraryReadsGot\": {\n", dh);
    int msn = mOptions->samples.size() - 1;
    uint64_t ttbase = 0;
    for(int i = 0; i < msn; ++i){
        ksprintf(s, "%s%s\"%s\": %lld,\n", dh, dm, mOptions->samples[i].c_str(), mSplitCount[i]);
        ttbase += mBaseCount[i];
    }
    s->l -= 2;
    ksprintf(s, "\n%s},\n", dh);
    // LibraryBaseGot
    ksprintf(s, "%s\"LibraryBaseGot\": {\n", dh);
    for(int i = 0; i < msn; ++i){
        ksprintf(s, "%s%s\"%s\": \"%lld(%.2lf%%)\",\n", 
                    dh, dm,
                    mOptions->samples[i].c_str(), 
                    mBaseCount[i], (double)mBaseCount[i]/ttbase * 100);
    }
    s->l -= 2;
    ksprintf(s, "\n%s},\n", dh);
    // LibraryStrandStat
    ksprintf(s, "%s\"LibraryStrand\": {\n", dh);
    for(int i = 0; i < msn; ++i){
        ksprintf(s, "%s%s\"%s\": \"FR:%lld,RF:%lld\",\n",
                    dh, dm,
                    mOptions->samples[i].c_str(), 
                    mFRCount[i], mRFCount[i]);
    }
    s->l -= 2;
    ksprintf(s, "\n%s}", dh);
}

void SplitResult::tsvHead(kstring_t* s, bool m){
    // summary
    ksprintf(s, "%s\t%s\t%s\t%s\t%s\t", "ReadsGot", "ReadsDrop", "GotRate", "FRRate", "RFRate");
    // mismatch
    for(size_t i = 0; i < mMMCount.size(); ++i){
        ksprintf(s, "%s%lu\t", "Mismatch", i);
    }
    // split details
    if(m) ksprintf(s, "%s", "SplitCount");
    else{
        int msn = mOptions->samples.size() - 1;
        for(int i = 0; i < msn; ++i){
            ksprintf(s, "%s\t", mOptions->samples[i].c_str());
        }
        s->l -= 1;// drop tailing \t
    }
}

void SplitResult::tsvBody(kstring_t* s, bool m){
    if(!mSummarized) summary();
    ksprintf(s, "%lld\t%lld\t%.2lf\t%.2lf\t%.2lf\t", mTotalReads-mDropReads, mDropReads, double(mTotalReads - mDropReads)/mTotalReads, double(mTotalFR)/mTotalReads, double(mTotalRF)/mTotalReads);
    for(size_t i = 0; i < mMMCount.size(); ++i){
        ksprintf(s, "%f\t", (double)mMMCount[i]/(double)mTotalReads);
    }
    int msn = mOptions->samples.size() - 1;
    if(m){
        for(int i = 0; i < msn; ++i){
            ksprintf(s, "%s(%lld);", mOptions->samples[i].c_str(), mSplitCount[i]);
        }
    }else{
        for(int i = 0; i < msn; ++i){
            ksprintf(s, "%lld\t", mSplitCount[i]);
        }
    }
    s->l -= 1;
}

