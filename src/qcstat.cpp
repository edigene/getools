#include "qcstat.h"

QcStat::QcStat(Options* opt, size_t bufLen){
    mOpt = opt;
    mBuflen = bufLen;
    mReads = mBases = 0;
    mMinReadLen = 0xffff;
    mMaxReadLen = 0;
    mMeanReadLen = 0;
    mGCPercent = .0;
    mMinQual = 127;
    mMaxQual = 33;
    mQ20Chr = '5';
    mQ30Chr = '?';
    mQualBase = 33;
    mLowQChr = mOpt->lowq.minBaseQual + '!';
    if(mOpt->isq64){
        mQ20Chr = 'T';
        mQ30Chr = '^';
        mQualBase = 64;
        mLowQChr = mOpt->lowq.minBaseQual + '@';
    }
    mMeanQual = 0;
    mQ20Total = mQ30Total = 0;
    mLengthSum = 0;
    mLowQual = 0;
    mQualSum = 0;
    mNReads = 0;
    mSummarized = false;
    allocateRes();
}

QcStat::~QcStat(){
    for(int i = 0; i < 5; ++i){
        free(mCycleQ20Bases[i]); mCycleQ20Bases[i] = NULL;
        free(mCycleQ30Bases[i]); mCycleQ30Bases[i] = NULL;
        free(mCycleBaseContents[i]); mCycleBaseContents[i] = NULL;
        free(mCycleBaseQuality[i]); mCycleBaseQuality[i] = NULL;
    }
    free(mCycleTotalBase); mCycleTotalBase = NULL;
    free(mCycleTotalQuality); mCycleTotalQuality = NULL;
    for(auto iter = mContentCurves.begin(); iter != mContentCurves.end(); ++iter){
        free(iter->second); iter->second = NULL;
    }
    for(auto iter = mQualityCurves.begin(); iter != mQualityCurves.end(); ++iter){
        free(iter->second); iter->second = NULL;
    }
}

void QcStat::allocateRes(){
    for(int i = 0; i < 5; ++i){
        mQ20Bases[i] = 0;
        mQ30Bases[i] = 0;
        mBaseContents[i] = 0;
        mCycleQ20Bases[i] = (size_t*)calloc(mBuflen, sizeof(size_t));
        mCycleQ30Bases[i] = (size_t*)calloc(mBuflen, sizeof(size_t));
        mCycleBaseContents[i] = (size_t*)calloc(mBuflen, sizeof(size_t));
        mCycleBaseQuality[i] = (size_t*)calloc(mBuflen, sizeof(size_t));
    }
    mCycleTotalBase = (size_t*)calloc(mBuflen, sizeof(size_t));
    mCycleTotalQuality = (size_t*)calloc(mBuflen, sizeof(size_t));
}

void QcStat::reallocRes(size_t bufLen){
    if(bufLen <= mBuflen) return;
    size_t obfl = mBuflen;
    mBuflen = bufLen;
    for(int i = 0; i < 5; ++i){
        mCycleQ20Bases[i] = (size_t*)realloc(mCycleQ20Bases[i], mBuflen * sizeof(size_t));
        memset(mCycleQ20Bases[i]+obfl, 0, (mBuflen-obfl)*sizeof(size_t));
        mCycleQ30Bases[i] = (size_t*)realloc(mCycleQ30Bases[i], mBuflen * sizeof(size_t));
        memset(mCycleQ30Bases[i]+obfl, 0, (mBuflen-obfl)*sizeof(size_t));
        mCycleBaseContents[i] = (size_t*)realloc(mCycleBaseContents[i], mBuflen * sizeof(size_t));
        memset(mCycleBaseContents[i]+obfl, 0, (mBuflen-obfl)*sizeof(size_t));
        mCycleBaseQuality[i] = (size_t*)realloc(mCycleBaseQuality[i], mBuflen * sizeof(size_t));
        memset(mCycleBaseQuality[i]+obfl, 0, (mBuflen-obfl)*sizeof(size_t));
    }
    mCycleTotalBase = (size_t*)realloc(mCycleTotalBase, mBuflen * sizeof(size_t));
    memset(mCycleTotalBase+obfl, 0, (mBuflen-obfl)*sizeof(size_t));
    mCycleTotalQuality = (size_t*)realloc(mCycleTotalQuality, mBuflen * sizeof(size_t));
    memset(mCycleTotalQuality+obfl, 0, (mBuflen-obfl)*sizeof(size_t));
}

bool QcStat::statRead(const krec1_t* r, int off){
    mGoodRead = true;
    if(r->seq.l > mBuflen) reallocRes(r->seq.l);
    ++mReads;
    mBases += r->seq.l - off;
    mLengthSum += r->seq.l - off;
    if(r->seq.l - off < mMinReadLen) mMinReadLen = r->seq.l;
    if(r->seq.l - off > mMaxReadLen) mMaxReadLen = r->seq.l;
    int ncnt = 0, lowqcnt = 0, bidx = 0;
    for(size_t i = off; i < r->seq.l; ++i){
        bidx = nuc_to_3bit[(int)r->seq.s[i]];
        mMinQual = std::min(r->qual.s[i] - mQualBase, mMinQual);
        mMaxQual = std::max(r->qual.s[i] - mQualBase, mMaxQual);
        if(r->qual.s[i] > mQ30Chr){
            ++mCycleQ30Bases[bidx][i];
            ++mCycleQ20Bases[bidx][i];
            ++mQ20Total;
            ++mQ30Total;
        }else if(r->qual.s[i] > mQ20Chr){
            ++mCycleQ20Bases[bidx][i];
            ++mQ20Total;
        }
        if(r->qual.s[i] < mLowQChr) ++lowqcnt;
        ncnt += (bidx == 4);
        ++mCycleBaseContents[bidx][i];
        mCycleBaseQuality[bidx][i] += r->qual.s[i] - mQualBase;
        ++mCycleTotalBase[i];
        mCycleTotalQuality[i] += r->qual.s[i] - mQualBase;
        ++mBaseContents[bidx];
    }
    // quality and N check
    if(ncnt > mOpt->lowq.maxN || lowqcnt > mOpt->lowq.maxLowQualBase || lowqcnt >  mOpt->lowq.maxLowQualFrac * r->seq.l){
        ++mLowQual;
        mGoodRead = false;
    }
    if(ncnt > mOpt->lowq.maxN) ++mNReads;
    // length check
    if(r->seq.l < mOpt->minlen) mGoodRead = false;
    // leading N check
    // for(size_t i = 0; i < mOpt->minlen; ++i){
    //     if(nuc_to_3bit[(int)r->seq.s[i]] == 4){
    //         mGoodRead = false;
    //         break;
    //     }
    // }
    return mGoodRead;
}

QcStat* QcStat::merge(const std::vector<QcStat*>& qcs){
    if(qcs.size() == 0) return NULL;
    size_t maxbuff = 0;
    for(size_t i = 0; i < qcs.size(); ++i){
        if(maxbuff < qcs[i]->mBuflen) maxbuff = qcs[i]->mBuflen;
    }
    QcStat* q = new QcStat(qcs[0]->mOpt, maxbuff);
    for(size_t i = 0; i < qcs.size(); ++i){
        q->mReads += qcs[i]->mReads;
        q->mBases += qcs[i]->mBases;
        q->mNReads += qcs[i]->mNReads;
        if(q->mMinReadLen > qcs[i]->mMinReadLen) q->mMinReadLen = qcs[i]->mMinReadLen;
        if(q->mMaxReadLen < qcs[i]->mMaxReadLen) q->mMaxReadLen = qcs[i]->mMaxReadLen;
        if(q->mMinQual > qcs[i]->mMinQual) q->mMinQual = qcs[i]->mMinQual;
        if(q->mMaxQual < qcs[i]->mMaxQual) q->mMaxQual = qcs[i]->mMaxQual;
        q->mLowQual += qcs[i]->mLowQual;
        q->mQ20Total += qcs[i]->mQ20Total;
        q->mQ30Total += qcs[i]->mQ30Total;
        q->mLengthSum += qcs[i]->mLengthSum;
        for(size_t j = 0; j < 5; ++j){
            q->mQ20Bases[j] += qcs[i]->mQ20Bases[j];
            q->mQ30Bases[j] += qcs[i]->mQ30Bases[j];
            q->mBaseContents[j] += qcs[i]->mBaseContents[j];
            for(size_t c = 0; c < std::min(q->mBuflen, qcs[i]->mBuflen); ++c){
                q->mCycleQ20Bases[j][c] += qcs[i]->mCycleQ20Bases[j][c];
                q->mCycleQ30Bases[j][c] += qcs[i]->mCycleQ30Bases[j][c];
                q->mCycleBaseContents[j][c] += qcs[i]->mCycleBaseContents[j][c];
                q->mCycleBaseQuality[j][c] += qcs[i]->mCycleBaseQuality[j][c];
            }
        }
        for(size_t c = 0; c < std::min(q->mBuflen, qcs[i]->mBuflen); ++c){
            q->mCycleTotalBase[c] += qcs[i]->mCycleTotalBase[c];
            q->mCycleTotalQuality[c] += qcs[i]->mCycleTotalQuality[c];
        }
    }
    return q;
}

void QcStat::summary(){
    if(mSummarized) return;
    // get cycle and total bases
    mCycle = mBuflen;
    for(mCycle = mBuflen; mCycle >= 1; --mCycle){
        if(mCycleTotalBase[mCycle-1]) break;
    }
    for(size_t i = 0; i < mCycle; ++i){
        mQualSum += mCycleTotalQuality[i];
    }
    if(mBases > 0){
        mMeanQual = (double)mQualSum/(double)mBases;
        mGCPercent = (double)(mBaseContents[nuc_to_3bit[(int)'G']] + mBaseContents[nuc_to_3bit[(int)'C']])/(double)mBases;
    }
    cycleCurve();
    if(mReads > 0) mMeanReadLen = (double)mLengthSum/(double)mReads;

    mSummarized = true;
}
void QcStat::cycleCurve(){
    // mean qual curve
    double *meanQualCurve = (double*)calloc(mCycle, sizeof(double));
    for(size_t i = 0; i < mCycle; ++i){
        if(mCycleTotalBase[i]) meanQualCurve[i] = (double)mCycleTotalQuality[i]/(double)mCycleTotalBase[i];
        else meanQualCurve[i] = 0;
    }
    mQualityCurves["Mean"] = meanQualCurve;
    // qualty curves and base contents curves
    for(size_t b = 0; b < 5; ++b){
        double* qualCurve = (double*)calloc(mCycle, sizeof(double));
        double* contentCurve = (double*)calloc(mCycle, sizeof(double));
        for(size_t c = 0; c < mCycle; ++c){
            if(mCycleBaseContents[b][c] == 0){
                qualCurve[c] = 0;
                contentCurve[c] = 0;
            }else{
                qualCurve[c] = (double)mCycleBaseQuality[b][c]/(double)mCycleBaseContents[b][c];
            }
            if(mCycleTotalBase[c] == 0){
                contentCurve[c] = 0;
            }else{
                contentCurve[c] = (double)mCycleBaseContents[b][c]/(double)mCycleTotalBase[c];
            }
        }
        mQualityCurves[std::string(1, bit3_to_nuc[b])] = qualCurve;
        mContentCurves[std::string(1, bit3_to_nuc[b])] = contentCurve;
    }
    // GC content
    double* gcContentCurve = (double*)calloc(mCycle, sizeof(double));
    int gidx = nuc_to_3bit[(int)'G'];
    int cidx = nuc_to_3bit[(int)'C'];
    for(size_t c = 0; c < mCycle; ++c){
        if(mCycleTotalBase[c] == 0){
            gcContentCurve[c] = .0;
        }else{
            gcContentCurve[c] = (double)(mCycleBaseContents[gidx][c] + mCycleBaseContents[cidx][c])/mCycleTotalBase[c];
        }
    }
    mContentCurves["GC"] = gcContentCurve;
}

void QcStat::reportJSON(kstring_t* s, const char* dh, const char* dm, int r){
    if(!mSummarized) summary();
    ksprintf(s, "%s\"Read%d\": {\n", dh, r);
    ksprintf(s, "%s%s\"TotalReads\": %lu,\n", dh, dm, mReads);
    ksprintf(s, "%s%s\"TotalBases\": %lu,\n", dh, dm, mBases);
    ksprintf(s, "%s%s\"MinReadLen\": %lu,\n", dh, dm, mMinReadLen);
    ksprintf(s, "%s%s\"MaxReadLen\": %lu,\n", dh, dm, mMaxReadLen);
    ksprintf(s, "%s%s\"MeanReadLen\": %lu,\n", dh, dm, mMeanReadLen);
    ksprintf(s, "%s%s\"MinBaseQual\": %d,\n", dh, dm, mMinQual);
    ksprintf(s, "%s%s\"MaxBaseQual\": %d,\n", dh, dm, mMaxQual);
    ksprintf(s, "%s%s\"MeanBaseQual\": %f,\n", dh, dm, mMeanQual);
    ksprintf(s, "%s%s\"Q20Bases\": %lu,\n", dh, dm, mQ20Total);
    ksprintf(s, "%s%s\"Q30Bases\": %lu,\n", dh, dm, mQ30Total);
    ksprintf(s, "%s%s\"Q20Rate\": %f,\n", dh, dm, mBases == 0 ? 0 : mQ20Total/(double)mBases);
    ksprintf(s, "%s%s\"Q30Rate\": %f,\n", dh, dm, mBases == 0 ? 0 : mQ30Total/(double)mBases);
    ksprintf(s, "%s%s\"LowQualReads\": %lu,\n", dh, dm, mLowQual);
    ksprintf(s, "%s%s\"NRichReads\": %lu,\n", dh, dm, mNReads);
    ksprintf(s, "%s%s\"GCContent\": %f,\n", dh, dm, mGCPercent);
    // QualityCurve
    if(mBases > 0 && mCycle > 0){
        ksprintf(s, "%s%s\"QualityCurve\": {", dh, dm);
        for(auto iter = mQualityCurves.begin(); iter != mQualityCurves.end(); ++iter){
            ksprintf(s, "\n%s%s%s\"%s\": [", dh, dm, dm, iter->first.c_str());
            for(size_t c = 0; c < mCycle-1; ++c) ksprintf(s, "%f,", iter->second[c]);
            ksprintf(s, "%f],", iter->second[mCycle-1]);
        }
        s->l -= 1; // trim invalid ','
        ksprintf(s, "\n%s%s},\n", dh, dm);
        // ContentCurve
        ksprintf(s, "%s%s\"ContentCurve\": {", dh, dm);
        for(auto iter = mContentCurves.begin(); iter != mContentCurves.end(); ++iter){
            ksprintf(s, "\n%s%s%s\"%s\": [", dh, dm, dm, iter->first.c_str());
            for(size_t c = 0; c < mCycle-1; ++c) ksprintf(s, "%f,", iter->second[c]);
            ksprintf(s, "%f],", iter->second[mCycle-1]);
        }
        s->l -= 1; // trim invalid ','
        ksprintf(s, "\n%s%s}\n", dh, dm);
    }else{
        s->l -= 2;
        ksprintf(s, "\n");
    }
    ksprintf(s, "%s}", dh);
}

void QcStat::reportHTML(kstring_t* s, int r){
    summary();
    // quality
    {
        std::string subsect = "Read" + std::to_string(r) + " base quality";
        std::string divName = util::replace(subsect, " ", "_");
        std::string title = "";
        ksprintf(s, "<div class='subsection_title'><a title='click to hide/show' onclick=showOrHide('%s')>%s</a></div>\n", divName.c_str(), subsect.c_str());
        ksprintf(s, "<div id='%s'>\n", divName.c_str());
        ksprintf(s, "<div class='sub_section_tips'>Value of each position will be shown on mouse over.</div>\n");
        ksprintf(s, "<div class='figure' id='plot_%s'></div>\n", divName.c_str());
        ksprintf(s, "</div>\n");

        std::string alphabets[5] = {"A", "T", "C", "G", "Mean"};
        std::string colors[5] = {"rgba(128,128,0,1.0)", "rgba(128,0,128,1.0)", "rgba(0,255,0,1.0)", "rgba(0,0,255,1.0)", "rgba(20,20,20,1.0)"};
        ksprintf(s, "\n<script type=\"text/javascript\">\n");
        std::string jsnstr = "var data = [\n";
        int64_t* x = (int64_t*)calloc(mCycle, sizeof(int64_t));
        int64_t total = 0;
        for(size_t i = 0; i < mCycle; ++i){
            x[total] = i + 1;
            ++total;
        }
        std::string xcrds = util::join(x, total, ",");
        free(x); x = NULL;
        for(int b = 0; b < 5; ++b){
            std::string base = std::string(1, bit3_to_nuc[b]);
            jsnstr.append("  {\n");
            jsnstr.append("    x: [" + xcrds + "],\n");
            jsnstr.append("    y: [" + util::join(mQualityCurves[base], total, ",") + "],\n");
            jsnstr.append("    name: '" + base + "',\n");
            jsnstr.append("    mode: 'markers',\n");
            jsnstr.append("    type: 'scatter',\n");
            jsnstr.append("    marker: {color:'" + colors[b] + "', size:6},\n");
            jsnstr.append("  },\n");
        }
        jsnstr.append("];\n");
        jsnstr.append("var layout = {\n");
        jsnstr.append("  title:'" + title + "',\n");
        jsnstr.append("  xaxis: {\n");
        jsnstr.append("    title: 'position',\n");
        jsnstr.append("  },\n");
        jsnstr.append("  yaxis: {\n");
        jsnstr.append("    title: 'quality',\n");
        jsnstr.append("  },\n");
        jsnstr.append("};\n");

        jsnstr.append("var config = {\n");
        jsnstr.append("  toImageButtonOptions: {\n");
        jsnstr.append("    format: 'svg',\n");
        jsnstr.append("     filename: '" + divName + "',\n");
        jsnstr.append("     height: " + std::to_string(mOpt->hmo.figh) + ",\n");
        jsnstr.append("     width: " + std::to_string(mOpt->hmo.figw) + ",\n");
        jsnstr.append("     scale: 1,\n");
        jsnstr.append("  }\n");
        jsnstr.append("};\n");

        jsnstr.append("Plotly.newPlot('plot_" + divName + "', data, layout, config);\n");
        ksprintf(s, "%s", jsnstr.c_str());
        ksprintf(s, "</script>\n");
    }
    // content
    {
        std::string subsect = "Read" + std::to_string(r) + " base contents";
        std::string divName = util::replace(subsect, " ", "_");
        std::string title = "";
        ksprintf(s, "<div class='subsection_title'><a title='click to hide/show' onclick=showOrHide('%s')>%s</a></div>\n", divName.c_str(), subsect.c_str());
        ksprintf(s, "<div id='%s'>\n", divName.c_str());
        ksprintf(s, "<div class='sub_section_tips'>Value of each position will be shown on mouse over.</div>\n");
        ksprintf(s, "<div class='figure' id='plot_%s'></div>\n", divName.c_str());
        ksprintf(s, "</div>\n");

        std::string alphabets[6] = {"A", "T", "C", "G", "N", "GC"};
        std::string colors[6] = {"rgba(128,128,0,1.0)", "rgba(128,0,128,1.0)", "rgba(0,255,0,1.0)", 
                                 "rgba(0,0,255,1.0)", "rgba(20,20,20,1.0)", "rgba(20,20,20,1.0)"};
        ksprintf(s, "\n<script type=\"text/javascript\">\n");
        std::string jsnstr = "var data = [\n";
        int64_t* x = (int64_t*)calloc(mCycle, sizeof(int64_t));
        int64_t total = 0;
        for(size_t i = 0; i < mCycle; ++i){
            x[total] = i + 1;
            ++total;
        }
        std::string xcrds = util::join(x, total, ",");
        free(x); x = NULL;
        for(int b = 0; b < 6; ++b){
            std::string base = alphabets[b];
            jsnstr.append("  {\n");
            jsnstr.append("    x: [" + xcrds + "],\n");
            jsnstr.append("    y: [" + util::join(mContentCurves[base], total, ",") + "],\n");
            jsnstr.append("    name: '" + base + "',\n");
            jsnstr.append("    mode: 'markers',\n");
            jsnstr.append("    type: 'scatter',\n");
            jsnstr.append("    marker: {color:'" + colors[b] + "', size:6}\n");
            jsnstr.append("  },\n");
        }
        jsnstr.append("];\n");

        jsnstr.append("var layout = {\n");
        jsnstr.append("  title:'" + title + "',\n");
        jsnstr.append("  xaxis: {\n");
        jsnstr.append("    title: 'position',\n");
        jsnstr.append("  },\n");
        jsnstr.append("  yaxis: {\n");
        jsnstr.append("    title: 'base content ratios',\n");
        jsnstr.append("  },\n");
        jsnstr.append("};\n");

        jsnstr.append("var config = {\n");
        jsnstr.append("  toImageButtonOptions: {\n");
        jsnstr.append("    format: 'svg',\n");
        jsnstr.append("     filename: '" + divName + "',\n");
        jsnstr.append("     height: " + std::to_string(mOpt->hmo.figh) + ",\n");
        jsnstr.append("     width: " + std::to_string(mOpt->hmo.figw) + ",\n");
        jsnstr.append("     scale: 1,\n");
        jsnstr.append("  }\n");
        jsnstr.append("};\n");

        jsnstr.append("Plotly.newPlot('plot_" + divName + "', data, layout, config);\n");
        ksprintf(s, "%s", jsnstr.c_str());
        ksprintf(s, "</script>\n");
    }
}

void QcStat::tsvHead(kstring_t* s){
    // bases summary
    ksprintf(s, "%s\t%s\t%s\t", "TotalBases", "Q20Bases", "Q30Bases");
    // reads summary
    ksprintf(s, "%s\t%s\t%s\t%s\t", "TotalReads", "HighQualReads", "LowQualReads", "NRichReads");
    // quality
    ksprintf(s, "%s\t%s\t%s\t%s\t", "MeanBaseQual", "<Q20Rate", "Q20Rate", "Q30Rate");
    // base freq
    ksprintf(s, "%s\t%s\t%s\t%s\t%s", "A", "C", "G", "T", "N");
}

void QcStat::tsvBody(kstring_t* s, QcStat* r2){
    if(!mSummarized) summary();
    if(r2 && (!r2->mSummarized)) r2->summary();
    if(r2){
        ksprintf(s, "%lu\t%lu\t%lu\t", mBases + r2->mBases, mQ20Total + r2->mQ20Total, mQ30Total + r2->mQ30Total);
        ksprintf(s, "%lu\t%lu\t%lu\t%lu\t", mReads + r2->mReads, mReads + r2->mReads - mLowQual - r2->mLowQual, mLowQual + r2->mLowQual, mNReads + r2->mNReads);
        if(mBases + r2->mBases == 0){
            ksprintf(s, "0\t0\t0\t0\t");
            for(int i = 0; i < 5; ++i) ksprintf(s, "0\t");
        }else{
            ksprintf(s, "%lf\t%lf\t%lf\t%lf\t", (mMeanQual * mBases + r2->mMeanQual * r2->mBases)/(mBases+r2->mBases),
                                                1 - (double)(mQ20Total+r2->mQ20Total)/(double)(mBases+r2->mBases),
                                                (double)(mQ20Total+r2->mQ20Total)/(double)(mBases+r2->mBases),
                                                (double)(mQ30Total+r2->mQ30Total)/(double)(mBases+r2->mBases));
            for(int i = 0; i < 5; ++i){
                ksprintf(s, "%lf\t", (double)(mBaseContents[i] + r2->mBaseContents[i])/(double)(mBases+r2->mBases));
            }
        }
    }else{
        ksprintf(s, "%lu\t%lu\t%lu\t", mBases, mQ20Total, mQ30Total);
        ksprintf(s, "%lu\t%lu\t%lu\t%lu\t", mReads, mReads - mLowQual, mLowQual, mNReads);
        if(mBases == 0){
            ksprintf(s, "0\t0\t0\t0\t");
            for(int i = 0; i < 5; ++i) ksprintf(s, "0\t");
        }else{
            ksprintf(s, "%lf\t%lf\t%lf\t%lf\t", mMeanQual,
                                                1 - (double)mQ20Total/(double)mBases,
                                                (double)mQ20Total/(double)mBases,
                                                (double)mQ30Total/(double)mBases);
            for(int i = 0; i < 5; ++i){
                ksprintf(s, "%lf\t", (double)(mBaseContents[i])/(double)(mBases));
            }
        }
    }
    s->l -= 1;
}
