#include "htmlutil.h"
#include "processor.h"

void Processor::printHeader(kstring_t* s){
    mOptions->hmo.html2head(s, "getools");
}

void Processor::reportOneHTML(int i){
    kstring_t* s = (kstring_t*)calloc(1, sizeof(kstring_t));
    printHeader(s);
    QcStat* qs1 = mSplitResults->mQCRead1[i];
    QcStat* qs2 = mSplitResults->mQCRead2[i];
    PairMerger* merge = mPeMergers[i];
    DeReper* derep = mDeRepers[i];
    GEDetector* gedect = mGEDetectors[i];
    kstring_t* tops = mTopnHtml[i];
    // report title
     ksprintf(s, "<h1 style='text-align:left;'><a href='https://github.com/wulj2/getools.git' target='_blank' style='color:#000099;text-decoration:none;'>Gene edit report of %s</a></h1>\n", mOptions->samples[i].c_str());
    // Summary beg
    ksprintf(s, "<div class='section_div'>\n");
    ksprintf(s, "<div class='section_title' onclick=showOrHide('summary')><a name='summary'>Summary</a></div>\n");
    ksprintf(s, "<div id='summary'>\n");
    // Summary -> General beg
    ksprintf(s, "<div class='subsection_title' onclick=showOrHide('general')>General</div>\n");
    ksprintf(s, "<div id='general'>\n");
    ksprintf(s, "<table class='summary_table'>\n");
    {
        // getools version
        std::string verinfo = mOptions->soft->ver + " (<a href='https://github.com/wulj2/getools.git'>https://github.com/wulj2/getools.git)</a>"; 
        htmlutil::outTableRow(s, "getools version:",  verinfo.c_str());
        // sequencing
        std::string cycinfo = "";
        if(qs2){
            cycinfo = "paired end";
            cycinfo.append("(" + std::to_string(qs1->mCycle) + " cycles + " + std::to_string(qs2->mCycle) + " cycles)");
        }else{
            cycinfo = "single end(" + std::to_string(qs1->mCycle) + " cycles)";
        }
        htmlutil::outTableRow(s, "sequencing:", cycinfo.c_str());
        // mean length
        std::string meanls = std::to_string(qs1->mMeanReadLen) + "bp";
        if(qs2) meanls.append(", " + std::to_string(qs2->mMeanReadLen) + "bp");
        htmlutil::outTableRow(s, "mean read length:", meanls.c_str());
    }
    ksprintf(s, "</table>\n");
    ksprintf(s, "</div>\n");
    // Summary -> General end
    // Summary -> Fastq QC beg
    ksprintf(s, "<div class='subsection_title' onclick=showOrHide('qc')>Fastq QC</div>\n");
    ksprintf(s, "<div id='qc'>\n");
    ksprintf(s, "<table class='summary_table'>\n");
    {
        // total reads
        std::string totrs = htmlutil::formatNumber(qs1->mReads);
        if(qs2) totrs.append(", " + htmlutil::formatNumber(qs2->mReads));
        htmlutil::outTableRow(s, "total reads:", totrs.c_str());
        // total bases
        std::string totbs = htmlutil::formatNumber(qs1->mBases);
        if(qs2) totbs.append(", " + htmlutil::formatNumber(qs2->mBases));
        htmlutil::outTableRow(s, "total bases:",  totbs.c_str());
        // Q20 bases
        std::string tq20 = htmlutil::formatNumber(qs1->mQ20Total);
        tq20.append("(" + std::to_string((double)qs1->mQ20Total/(double)qs1->mBases*100) + "%)");
        if(qs2){
            tq20.append(", " + htmlutil::formatNumber(qs2->mQ20Total));
            tq20.append("(" + std::to_string((double)qs2->mQ20Total/(double)qs2->mBases*100) + "%)");
        }
        htmlutil::outTableRow(s, "Q20 bases:",  tq20.c_str());
        // Q30 bases
        std::string tq30 = htmlutil::formatNumber(qs1->mQ30Total);
        tq30.append("(" + std::to_string((double)qs1->mQ30Total/(double)qs1->mBases*100) + "%)");
        if(qs2){
            tq30.append(", " + htmlutil::formatNumber(qs2->mQ30Total));
            tq30.append("(" + std::to_string((double)qs2->mQ30Total/(double)qs2->mBases*100) + "%)");
        }
        htmlutil::outTableRow(s, "Q30 bases:",  tq30.c_str());
        // GC contents
        std::string gcs = std::to_string(qs1->mGCPercent*100) + "%";
        if(qs2) gcs.append(", " + std::to_string(qs2->mGCPercent*100) + "%");
        htmlutil::outTableRow(s, "GC Content:", gcs.c_str());
        // Reads with Adapter
        if(qs2 && merge && merge->total){
            std::string rwas = std::to_string(merge->peovhang);
            rwas.append("(" + std::to_string((double)merge->peovhang/(double)(merge->total)*100) + "%)");
            htmlutil::outTableRow(s, "pairs with adapters:", rwas.c_str());
        }
    }
    ksprintf(s, "</table>\n");
    ksprintf(s, "</div>\n");
    // Summary -> Fastq QC end
    // Summary -> Computation beg
    ksprintf(s, "<div class='subsection_title' onclick=showOrHide('cal')>Computation</div>\n");
    ksprintf(s, "<div id='cal'>\n");
    ksprintf(s, "<table class='summary_table'>\n");
    {
        // total seqs
        std::string tts = htmlutil::formatNumber(qs1->mReads);
        htmlutil::outTableRow(s, "total sequences:", tts.c_str());
        // merged seqs
        if(qs2){
            std::string mgs = htmlutil::formatNumber(merge->merged);
            mgs.append(" (" + std::to_string((double)merge->merged/(double)merge->total*100) + "%)");
            htmlutil::outTableRow(s, "merged sequences:", mgs.c_str());
        }
        // uniq seqs
        std::string uqs = htmlutil::formatNumber(derep->clusters);
        uqs.append(" (" + std::to_string(100*(1-derep->duprate)) + "%)");
        htmlutil::outTableRow(s, "unique sequences:", uqs.c_str());
        // drop count
        std::string dps = htmlutil::formatNumber(gedect->dropcnt);
        dps.append(" (" + std::to_string((double)gedect->dropcnt/(double)derep->sequencecount*100) + "%)");
        htmlutil::outTableRow(s, "dropped sequences:", dps.c_str());
        // ref count
        std::string rfs = htmlutil::formatNumber(gedect->refcnt);
        rfs.append(" (" + std::to_string((double)gedect->refcnt/(double)gedect->totcnt*100) + "%)");
        htmlutil::outTableRow(s, "reference sequences:", rfs.c_str());
        // other count
        std::string ots = htmlutil::formatNumber(gedect->othcnt);
        ots.append(" (" + std::to_string((double)gedect->othcnt/(double)gedect->totcnt*100) + "%)");
        htmlutil::outTableRow(s, "other sequences:", ots.c_str());
        // edit count
        std::string eds = htmlutil::formatNumber(gedect->edicnt);
        eds.append(" (" + std::to_string(gedect->edieff*100) + "%)");
        htmlutil::outTableRow(s, "allele sequences:", eds.c_str());
        // donor count
        if(mOptions->recstat){
            std::string rcs = htmlutil::formatNumber(gedect->reccnt);
            rcs.append(" (" + std::to_string(gedect->receff*100) + "%)");
            htmlutil::outTableRow(s, "total HDR sequences:", rcs.c_str());
            rcs = htmlutil::formatNumber(gedect->recent);
            rcs.append(" (" + std::to_string(gedect->reeeff*100) + "%)");
            htmlutil::outTableRow(s, "exact matched HDR sequences:", rcs.c_str());
            rcs = htmlutil::formatNumber(gedect->recpct);
            rcs.append(" (" + std::to_string(gedect->recpef*100) + "%)");
            htmlutil::outTableRow(s, "any HDR variant hit sequences:", rcs.c_str());
        }
    }
    ksprintf(s, "</table>\n");
    ksprintf(s, "</div>\n");
    // Summary -> Computation end
    ksprintf(s, "</div>\n");
    ksprintf(s, "</div>\n");
    // Summary end
    // Fastq QC beg
    ksprintf(s, "<div class='section_div'>\n");
    ksprintf(s, "<div class='section_title' onclick=showOrHide('qcplot')><a name='fqcsummary'>Fastq QC</a></div>\n");
    ksprintf(s, "<div id='qcplot'>\n");
    if(qs1) qs1->reportHTML(s, 1);
    if(qs2) qs2->reportHTML(s, 2);
    ksprintf(s, "</div>\n");
    ksprintf(s, "</div>\n");
    // Fastq QC end
    // Merge Stat beg
    if(qs2){
        ksprintf(s, "<div class='section_div'>\n");
        ksprintf(s, "<div class='section_title' onclick=showOrHide('merge')><a name='rmsummary'>Reads merge statistics</a></div>\n");
        ksprintf(s, "<div id='merge'>\n");
        merge->reportHTML(s);
        ksprintf(s, "</div>\n");
        ksprintf(s, "</div>\n");
    }
    // Merge Stat end
    // Derep Stat beg
    ksprintf(s, "<div class='section_div'>\n");
    ksprintf(s, "<div class='section_title' onclick=showOrHide('drplt')><a name='sdsummary'>Sequence deduplication statistics</a></div>\n");
    ksprintf(s, "<div id='drplt'>\n");
    derep->reportHTML(s);
    ksprintf(s, "</div>\n");
    ksprintf(s, "</div>\n");
    // Derep Stat end
    // Edit Stat beg
    ksprintf(s, "<div class='section_div'>\n");
    ksprintf(s, "<div class='section_title' onclick=showOrHide('edplt')><a name='gecsummary'>Gene edit computation</a></div>\n");
    ksprintf(s, "<div id='edplt'>\n");
    gedect->reportHTML(s);
    ksprintf(s, "%s", tops->s);
    ksprintf(s, "</div>\n");
    ksprintf(s, "</div>\n");
    // Edit Stat end
    printFooter(s);
    // write to file
    FILE* fp = fopen(mOutHTMLs[i].c_str(), "w");
    fwrite(s->s, sizeof(char), s->l, fp);
    fclose(fp);
    // release res
    free(s->s); free(s);
}

void Processor::reportQCOneHTML(int i){
    kstring_t* s = (kstring_t*)calloc(1, sizeof(kstring_t));
    printHeader(s);
    QcStat* qs1 = mSplitResults->mQCRead1[i];
    QcStat* qs2 = mSplitResults->mQCRead2[i];
    // report title
     ksprintf(s, "<h1 style='text-align:left;'><a href='https://github.com/wulj2/getools.git' target='_blank' style='color:#000099;text-decoration:none;'>Library split report of %s</a></h1>\n", mOptions->samples[i].c_str());
    // Summary beg
    ksprintf(s, "<div class='section_div'>\n");
    ksprintf(s, "<div class='section_title' onclick=showOrHide('summary')><a name='summary'>Summary</a></div>\n");
    ksprintf(s, "<div id='summary'>\n");
    // Summary -> General beg
    ksprintf(s, "<div class='subsection_title' onclick=showOrHide('general')>General</div>\n");
    ksprintf(s, "<div id='general'>\n");
    ksprintf(s, "<table class='summary_table'>\n");
    {
        // getools version
        std::string verinfo = mOptions->soft->ver + " (<a href='https://github.com/wulj2/getools.git'>https://github.com/wulj2/getools.git)</a>"; 
        htmlutil::outTableRow(s, "getools version:",  verinfo.c_str());
        // sequencing
        std::string cycinfo = "";
        if(qs2){
            cycinfo = "paired end";
            cycinfo.append("(" + std::to_string(qs1->mCycle) + " cycles + " + std::to_string(qs2->mCycle) + " cycles)");
        }else{
            cycinfo = "single end(" + std::to_string(qs1->mCycle) + " cycles)";
        }
        htmlutil::outTableRow(s, "sequencing:", cycinfo.c_str());
        // mean length
        std::string meanls = std::to_string(qs1->mMeanReadLen) + "bp";
        if(qs2) meanls.append(", " + std::to_string(qs2->mMeanReadLen) + "bp");
        htmlutil::outTableRow(s, "mean read length:", meanls.c_str());
    }
    ksprintf(s, "</table>\n");
    ksprintf(s, "</div>\n");
    // Summary -> General end
    // Summary -> Fastq QC beg
    ksprintf(s, "<div class='subsection_title' onclick=showOrHide('qc')>Fastq QC</div>\n");
    ksprintf(s, "<div id='qc'>\n");
    ksprintf(s, "<table class='summary_table'>\n");
    {
        // total reads
        std::string totrs = htmlutil::formatNumber(qs1->mReads);
        if(qs2) totrs.append(", " + htmlutil::formatNumber(qs2->mReads));
        htmlutil::outTableRow(s, "total reads:", totrs.c_str());
        // total bases
        std::string totbs = htmlutil::formatNumber(qs1->mBases);
        if(qs2) totbs.append(", " + htmlutil::formatNumber(qs2->mBases));
        htmlutil::outTableRow(s, "total bases:",  totbs.c_str());
        // Q20 bases
        std::string tq20 = htmlutil::formatNumber(qs1->mQ20Total);
        tq20.append("(" + std::to_string((double)qs1->mQ20Total/(double)qs1->mBases*100) + "%)");
        if(qs2){
            tq20.append(", " + htmlutil::formatNumber(qs2->mQ20Total));
            tq20.append("(" + std::to_string((double)qs2->mQ20Total/(double)qs2->mBases*100) + "%)");
        }
        htmlutil::outTableRow(s, "Q20 bases:",  tq20.c_str());
        // Q30 bases
        std::string tq30 = htmlutil::formatNumber(qs1->mQ30Total);
        tq30.append("(" + std::to_string((double)qs1->mQ30Total/(double)qs1->mBases*100) + "%)");
        if(qs2){
            tq30.append(", " + htmlutil::formatNumber(qs2->mQ30Total));
            tq30.append("(" + std::to_string((double)qs2->mQ30Total/(double)qs2->mBases*100) + "%)");
        }
        htmlutil::outTableRow(s, "Q30 bases:",  tq30.c_str());
        // GC contents
        std::string gcs = std::to_string(qs1->mGCPercent*100) + "%";
        if(qs2) gcs.append(", " + std::to_string(qs2->mGCPercent*100) + "%");
        htmlutil::outTableRow(s, "GC Content:", gcs.c_str());
    }
    ksprintf(s, "</table>\n");
    ksprintf(s, "</div>\n");
    // Summary -> Fastq QC end
    ksprintf(s, "</div>\n");
    ksprintf(s, "</div>\n");
    // Summary end
    // Fastq QC beg
    ksprintf(s, "<div class='section_div'>\n");
    ksprintf(s, "<div class='section_title' onclick=showOrHide('qcplot')><a name='fqcsummary'>Fastq QC</a></div>\n");
    ksprintf(s, "<div id='qcplot'>\n");
    if(qs1) qs1->reportHTML(s, 1);
    if(qs2) qs2->reportHTML(s, 2);
    ksprintf(s, "</div>\n");
    ksprintf(s, "</div>\n");
    ksprintf(s, "</div>\n");
    // Fastq QC end
    printFooter(s);
    // write to file
    FILE* fp = fopen(mOutHTMLs[i].c_str(), "w");
    fwrite(s->s, sizeof(char), s->l, fp);
    fclose(fp);
    // release res
    free(s->s); free(s);
}

void Processor::reportAllHTML(){
    kstring_t* s = (kstring_t*)calloc(1, sizeof(kstring_t));
    printHeader(s);
    QcStat* qs1 = mQcStats1;
    QcStat* qs2 = mQcStats2;
    SplitResult* spr = mSplitResults;
    // report title
    ksprintf(s, "<h1 style='text-align:left;'><a href='https://github.com/wulj2/getools.git' target='_blank' style='color:#000099;text-decoration:none;'>Gene edit report by getools</a></h1>\n");
    // Summary beg
    ksprintf(s, "<div class='section_div'>\n");
    ksprintf(s, "<div class='section_title' onclick=showOrHide('summary')><a name='summary'>Summary</a></div>\n");
    ksprintf(s, "<div id='summary'>\n");
    // Summary -> General beg
    ksprintf(s, "<div class='subsection_title' onclick=showOrHide('general')>General</div>\n");
    ksprintf(s, "<div id='general'>\n");
    ksprintf(s, "<table class='summary_table'>\n");
    {
        // getools version
        std::string verinfo = mOptions->soft->ver + " (<a href='https://github.com/wulj2/getools.git'>https://github.com/wulj2/getools.git)</a>"; 
        htmlutil::outTableRow(s, "getools version:",  verinfo.c_str());
        // sequencing
        std::string cycinfo = "";
        if(qs2){
            cycinfo = "paired end";
            cycinfo.append("(" + std::to_string(qs1->mCycle) + " cycles + " + std::to_string(qs2->mCycle) + " cycles)");
        }else{
            cycinfo = "single end(" + std::to_string(qs1->mCycle) + " cycles)";
        }
        htmlutil::outTableRow(s, "sequencing:", cycinfo.c_str());
        // mean length
        std::string meanls = std::to_string(qs1->mMeanReadLen) + "bp";
        if(qs2) meanls.append(", " + std::to_string(qs2->mMeanReadLen) + "bp");
        htmlutil::outTableRow(s, "mean read length:", meanls.c_str());
    }
    ksprintf(s, "</table>\n");
    ksprintf(s, "</div>\n");
    // Summary -> General end
    // Summary -> Fastq QC beg
    ksprintf(s, "<div class='subsection_title' onclick=showOrHide('qc')>Fastq QC</div>\n");
    ksprintf(s, "<div id='qc'>\n");
    ksprintf(s, "<table class='summary_table'>\n");
    {
        // total reads
        std::string totrs = htmlutil::formatNumber(qs1->mReads);
        if(qs2) totrs.append(", " + htmlutil::formatNumber(qs2->mReads));
        htmlutil::outTableRow(s, "total reads:", totrs.c_str());
        // total bases
        std::string totbs = htmlutil::formatNumber(qs1->mBases);
        if(qs2) totbs.append(", " + htmlutil::formatNumber(qs2->mBases));
        htmlutil::outTableRow(s, "total bases:",  totbs.c_str());
        // Q20 bases
        std::string tq20 = htmlutil::formatNumber(qs1->mQ20Total);
        tq20.append("(" + std::to_string((double)qs1->mQ20Total/(double)qs1->mBases*100) + "%)");
        if(qs2){
            tq20.append(", " + htmlutil::formatNumber(qs2->mQ20Total));
            tq20.append("(" + std::to_string((double)qs2->mQ20Total/(double)qs2->mBases*100) + "%)");
        }
        htmlutil::outTableRow(s, "Q20 bases:",  tq20.c_str());
        // Q30 bases
        std::string tq30 = htmlutil::formatNumber(qs1->mQ30Total);
        tq30.append("(" + std::to_string((double)qs1->mQ30Total/(double)qs1->mBases*100) + "%)");
        if(qs2){
            tq30.append(", " + htmlutil::formatNumber(qs2->mQ30Total));
            tq30.append("(" + std::to_string((double)qs2->mQ30Total/(double)qs2->mBases*100) + "%)");
        }
        htmlutil::outTableRow(s, "Q30 bases:",  tq30.c_str());
        // GC contents
        std::string gcs = std::to_string(qs1->mGCPercent*100) + "%";
        if(qs2) gcs.append(", " + std::to_string(qs2->mGCPercent*100) + "%");
        htmlutil::outTableRow(s, "GC Content:", gcs.c_str());
    }
    ksprintf(s, "</table>\n");
    ksprintf(s, "</div>\n");
    // Summary -> Fastq QC end
    // Split Stat beg
    ksprintf(s, "<div class='subsection_title' onclick=showOrHide('gsplrpt')>Pooled library split</div>\n");
    ksprintf(s, "<div id='gsplrpt'>\n");
    ksprintf(s, "<table class='summary_table'>\n");
    {
        std::string prefix = "Reads ";
        if(qs2) prefix = "Read pairs ";
        // Reads Input
        std::string inrs = htmlutil::formatNumber(spr->mTotalReads);
        std::string key = prefix + "to split:";
        htmlutil::outTableRow(s, key, inrs);
        // Reads Drop
        std::string rdrop = htmlutil::formatNumber(spr->mDropReads);
        rdrop.append("(QCFail:" + std::to_string((double)(spr->mQCFailedReads/(double)spr->mDropReads*100)) + "%,");
        rdrop.append(" SplitFail:" + std::to_string((double)(spr->mSplitFailedReads/(double)spr->mDropReads*100)) + "%)");
        key = prefix + "drop:";
        htmlutil::outTableRow(s, key, rdrop);
        // Reads Got
        std::string rgot = htmlutil::formatNumber(spr->mTotalReads-spr->mDropReads);
        rgot.append("(" + std::to_string((double)(spr->mTotalReads - spr->mDropReads)/(double)spr->mTotalReads*100) + "%)");
        key = prefix + "got:";
        htmlutil::outTableRow(s, key, rgot);
        // FR stat
        if(qs2){
            std::string frgv = "FR: " + htmlutil::formatNumber(spr->mTotalFR);
            frgv.append("(" + std::to_string((double)spr->mTotalFR/(double)spr->mTotalReads*100) + "%), ");
            frgv.append("RF: " + htmlutil::formatNumber(spr->mTotalRF));
            frgv.append("(" + std::to_string((double)spr->mTotalRF/(double)spr->mTotalReads*100) + "%)");
            key = prefix + " strands:";
            htmlutil::outTableRow(s, key, frgv);
        }
    }
    ksprintf(s, "</table>\n");
    ksprintf(s, "</div>\n");
    // Single cell beg
    if(mOptions->dosccal){
        ksprintf(s, "<div class='subsection_title' onclick=showOrHide('scrpt')>Single cell analysis</div>\n");
        ksprintf(s, "<div id='scrpt'>\n");
        ksprintf(s, "<table class='summary_table'>\n");
        {
            // raw cells
            htmlutil::outTableRow(s, "Raw cells", mScAna->rawc);
            // valid cells
            htmlutil::outTableRow(s, "Valid cells", mScAna->caffrac);
            // raw amplicons
            htmlutil::outTableRow(s, "Amplicons in panel", mScAna->rawa);
            // amplicons passed depth
            htmlutil::outTableRow(s, "QC passed amplicons", mScAna->aafth4mon);
            // allele across cell
            ksprintf(s, "<tr><td class='col1'><a href='./%s'>Allele table</a></td>", SCANA_AMP_ALLELE_IN_SC_TSV);
            ksprintf(s, "<td class='col2'>Top%d alleles of amplicons in valid cells</td></tr>", mOptions->edo.topn);
            // co-edit
            ksprintf(s, "<tr><td class='col1'><a href='./%s'>Co-edit table</a></td>", SCANA_AMP_CO_EDIT_TSV);
            ksprintf(s, "<td class='col2'>Amplicons overall edit status across valid cells</td></tr>");
            // ana-details
            ksprintf(s, "<tr><td class='col1'><a href='./%s'>Cell finder process</a></td>", SCANA_CELL_FIND_JSON);
            ksprintf(s, "<td class='col2'>Valid cell computation procedure</td></tr>");
        }
        ksprintf(s, "</table>\n");
        ksprintf(s, "</div>\n");
    }
    // Single cell end
    // Summary -> Split end
    ksprintf(s, "</div>\n");
    ksprintf(s, "</div>\n");
    // Summary end
    // Fastq QC beg
    ksprintf(s, "<div class='section_div'>\n");
    ksprintf(s, "<div class='section_title' onclick=showOrHide('qcplot')><a name='pfqsummary'>Pooled Fastq QC</a></div>\n");
    ksprintf(s, "<div id='qcplot'>\n");
    if(qs1) qs1->reportHTML(s, 1);
    if(qs2) qs2->reportHTML(s, 2);
    ksprintf(s, "</div>\n");
    ksprintf(s, "</div>\n");
    // Split QC beg
    ksprintf(s, "<div class='section_div'>\n");
    ksprintf(s, "<div class='section_title' onclick=showOrHide('splrpt')><a name='lsrsummary'>Library split result statistics</a></div>\n");
    ksprintf(s, "<div id='splrpt'>\n");
    mSplitResults->reportHTML(s);
    ksprintf(s, "</div>\n");
    ksprintf(s, "</div>\n");
    // Split QC end
    // Reads Ratio stat beg
    {
        ksprintf(s, "<div class='section_div'>\n");
        ksprintf(s, "<div class='section_title' onclick=showOrHide('rddist')><a name='rdssummary'>Reads distribution of each amplicon</a></div>\n");
        ksprintf(s, "<div id='rddist'>\n");
        std::string subsect = "Reads distribution of each amplicon";
        std::string divName = util::replace(subsect, " ", "_");
        std::string title = "'reads distribution of each amplicons'";
        ksprintf(s, "<div class='subsection_title'><a title='click to hide/show' onclick=showOrHide('%s')>%s</a></div>\n", divName.c_str(), subsect.c_str());
        ksprintf(s, "<div id='%s'>\n", divName.c_str());
        ksprintf(s, "<div class='sub_section_tips'>Value of each sample will be shown on mouse over,<br>Click each amplicon name to detailed sequnce distribution</div>\n");
        ksprintf(s, "<div class='figure' id='plot_%s'></div>\n", divName.c_str());
        ksprintf(s, "</div></div></div>\n");
        ksprintf(s, "\n<script type=\"text/javascript\">\n");

        std::string sse = "Edit efficience computation";
        std::string ssd = "plot_" + util::replace(sse, " ", "_");

        std::string jsnstr;
        jsnstr.append("var dropbar = {\n");
        jsnstr.append("  x: [");
        for(size_t i = 0; i < mOptions->samples.size()-1; ++i){
            jsnstr.append("'<a href=\"html/" + mOptions->samples[i] + ".html#" + ssd + "\">" + mOptions->samples[i] + "</a>',");
        }
        jsnstr.append("],\n");
        jsnstr.append("  y: [");
        for(size_t i = 0; i < mOptions->samples.size()-1; ++i){
            if(mGEDetectors[i]->allcnt){
                jsnstr.append(std::to_string((double)mGEDetectors[i]->dropcnt/(double)mGEDetectors[i]->allcnt*100) + ",");
            }else{
                jsnstr.append("0,");
            }
        }
        jsnstr.append("],\n");
        jsnstr.append("  type: 'bar',\n");
        jsnstr.append("  name: 'Dropped',\n");
        jsnstr.append("};\n");

        jsnstr.append("var refbar = {\n");
        jsnstr.append("  x: [");
        for(size_t i = 0; i < mOptions->samples.size()-1; ++i){
            jsnstr.append("'<a href=\"html/" + mOptions->samples[i] + ".html#" + ssd + "\">" + mOptions->samples[i] + "</a>',");
        }
        jsnstr.append("],\n");
        jsnstr.append("  y: [");
        for(size_t i = 0; i < mOptions->samples.size()-1; ++i){
            if(mGEDetectors[i]->allcnt){
                jsnstr.append(std::to_string((double)mGEDetectors[i]->refcnt/(double)mGEDetectors[i]->allcnt*100) + ",");
            }else{
                jsnstr.append("0,");
            }
        }
        jsnstr.append("],\n");
        jsnstr.append("  type: 'bar',\n");
        jsnstr.append("  name: 'RefType',\n");
        jsnstr.append("};\n");

        jsnstr.append("var otherbar = {\n");
        jsnstr.append("  x: [");
        for(size_t i = 0; i < mOptions->samples.size()-1; ++i){
            jsnstr.append("'<a href=\"html/" + mOptions->samples[i] + ".html#" + ssd + "\">" + mOptions->samples[i] + "</a>',");
        }
        jsnstr.append("],\n");
        jsnstr.append("  y: [");
        for(size_t i = 0; i < mOptions->samples.size()-1; ++i){
            if(mGEDetectors[i]->allcnt){
                jsnstr.append(std::to_string((double)mGEDetectors[i]->othcnt/(double)mGEDetectors[i]->allcnt*100) + ",");
            }else{
                jsnstr.append("0,");
            }
        }
        jsnstr.append("],\n");
        jsnstr.append("  type: 'bar',\n");
        jsnstr.append("  name: 'Others',\n");
        jsnstr.append("};\n");

        jsnstr.append("var edcbar = {\n");
        jsnstr.append("  x: [");
        for(size_t i = 0; i < mOptions->samples.size()-1; ++i){
            jsnstr.append("'<a href=\"html/" + mOptions->samples[i] + ".html#" + ssd + "\">" + mOptions->samples[i] + "</a>',");
        }
        jsnstr.append("],\n");
        jsnstr.append("  y: [");
        for(size_t i = 0; i < mOptions->samples.size()-1; ++i){
            if(mGEDetectors[i]->allcnt){
                jsnstr.append(std::to_string((double)mGEDetectors[i]->edicnt/(double)mGEDetectors[i]->allcnt*100) + ",");
            }else{
                jsnstr.append("0,");
            }
        }
        jsnstr.append("],\n");
        jsnstr.append("  type: 'bar',\n");
        jsnstr.append("  name: 'Edited',\n");
        jsnstr.append("};\n");


        jsnstr.append("var layout = {\n");
        jsnstr.append("  title:" + title + ",\n");
        jsnstr.append("  showlegend: true,\n");
        jsnstr.append("  barmode: 'stack',\n");
        jsnstr.append("  xaxis:{\n");
        jsnstr.append("    automargin: true,\n");
        jsnstr.append("    dtick: 1,\n");
        jsnstr.append("    title: {\n");
        jsnstr.append("      text: 'amplicons',\n");
        jsnstr.append("      standoff: 20,\n");
        jsnstr.append("    },\n");
        jsnstr.append("  },\n");
        jsnstr.append("  yaxis:{\n");
        jsnstr.append("    automargin: true,\n");
        jsnstr.append("    title: {\n");
        jsnstr.append("      text: 'percents(%)',\n");
        jsnstr.append("      standoff: 20,\n");
        jsnstr.append("    },\n");
        jsnstr.append("  },\n");
        jsnstr.append("};\n");

        jsnstr.append("var data = [dropbar, refbar, otherbar, edcbar];\n");

        jsnstr.append("var config = {\n");
        jsnstr.append("  toImageButtonOptions: {\n");
        jsnstr.append("    format: 'svg',\n"); jsnstr.append("     filename: '" + divName + "',\n");
        jsnstr.append("     height: " + std::to_string(mOptions->hmo.figh) + ",\n");
        jsnstr.append("     width: " + std::to_string(mOptions->hmo.figw) + ",\n");
        jsnstr.append("     scale: 1,\n");
        jsnstr.append("  }\n");
        jsnstr.append("};\n");

        jsnstr.append("Plotly.newPlot('plot_" + divName + "', data, layout, config);\n");
        ksprintf(s, "%s", jsnstr.c_str());
        ksprintf(s, "</script>\n");
    }
    // Reads Ratio stat end
    // Edit efficience beg
    {
        ksprintf(s, "<div class='section_div'>\n");
        ksprintf(s, "<div class='section_title' onclick=showOrHide('edrpt')><a name='effsummary'>Edit efficience of each amplicon</a></div>\n");
        ksprintf(s, "<div id='edrpt'>\n");
        std::string subsect = "Edit efficience of each amplicons";
        std::string divName = util::replace(subsect, " ", "_");
        std::string title = "'edit efficience of each amplicons'";
        ksprintf(s, "<div class='subsection_title'><a title='click to hide/show' onclick=showOrHide('%s')>%s</a></div>\n", divName.c_str(), subsect.c_str());
        ksprintf(s, "<div id='%s'>\n", divName.c_str());
        ksprintf(s, "<div class='sub_section_tips'>Value of each sample will be shown on mouse over.</div>\n");
        ksprintf(s, "<div class='figure' id='plot_%s'></div>\n", divName.c_str());
        ksprintf(s, "</div>\n");
        ksprintf(s, "\n<script type=\"text/javascript\">\n");
        std::string jsnstr;

        jsnstr.append("var edibar = {\n");
        jsnstr.append("  x: [");
        for(size_t i = 0; i < mOptions->samples.size()-1; ++i){
            jsnstr.append("'<a href=\"html/" + mOptions->samples[i] + ".html\">" + mOptions->samples[i] + "</a>',");
        }
        jsnstr.append("],\n");
        jsnstr.append("  y: [");
        for(size_t i = 0; i < mOptions->samples.size()-1; ++i){
            jsnstr.append(std::to_string(mGEDetectors[i]->edieff) + ",");
        }
        jsnstr.append("],\n");
        jsnstr.append("  type: 'bar',\n");
        jsnstr.append("  name: 'Edit Efficience',\n");
        jsnstr.append("  hovertext: [");
        for(size_t i = 0; i < mOptions->samples.size()-1; ++i){
            jsnstr.append("'total: " + std::to_string(mGEDetectors[i]->totcnt) + "<br>");
            jsnstr.append("edited: " + std::to_string(mGEDetectors[i]->edicnt));
            jsnstr.append("(" + std::to_string((double)(mGEDetectors[i]->edicnt)/(double)mGEDetectors[i]->totcnt*100) + "%)");
            jsnstr.append("<br>");
            jsnstr.append("reftype: " + std::to_string(mGEDetectors[i]->refcnt));
            jsnstr.append("(" + std::to_string((double)(mGEDetectors[i]->refcnt)/(double)mGEDetectors[i]->totcnt*100) + "%)");
            jsnstr.append("<br>");
            jsnstr.append("others: " + std::to_string(mGEDetectors[i]->othcnt));
            jsnstr.append("(" + std::to_string((double)(mGEDetectors[i]->othcnt)/(double)mGEDetectors[i]->totcnt*100) + "%)");
            jsnstr.append("',");
        }
        jsnstr.append("],\n");
        jsnstr.append("  hoverinfo: 'text',\n");
        jsnstr.append("}\n");

        jsnstr.append("var layout = {\n");
        jsnstr.append("  title:" + title + ",\n");
        jsnstr.append("  showlegend: false,\n");
        jsnstr.append("  xaxis:{\n");
        jsnstr.append("    automargin: true,\n");
        jsnstr.append("    dtick: 1,\n");
        jsnstr.append("    title: {\n");
        jsnstr.append("      text: 'amplicons',\n");
        jsnstr.append("      standoff: 20,\n");
        jsnstr.append("    },\n");
        jsnstr.append("  },\n");
        jsnstr.append("  yaxis:{\n");
        jsnstr.append("    automargin: true,\n");
        jsnstr.append("    title: {\n");
        jsnstr.append("      text: 'edit rate',\n");
        jsnstr.append("      standoff: 20,\n");
        jsnstr.append("    },\n");
        jsnstr.append("  },\n");
        jsnstr.append("};\n");

        jsnstr.append("var data = [edibar];\n");

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
    // Edit efficience ends
    // Detailed report begs
    {
        ksprintf(s, "<div class='subsection_title' onclick=showOrHide('edit_eff_table')>Links to detailed report of each amplicons</div>\n");
        ksprintf(s, "<div id='edit_eff_table'>\n");
        ksprintf(s, "<div class='sub_section_tips'>Click amplicon name to go to its detailed report,<br>Click edited rate to topN allele sequences,<br>Click mutated rate to variants counts along amplicon,<br>Amplicon with qualified sequences less than %d will be masked in red color<br>Amplicon with qualified sequences less than %d will be masked in yellow color</div>\n", mOptions->edo.minqs1k, mOptions->edo.minqs10k);
        ksprintf(s, "<table class='summary_table' id='eff_table_xx'>\n");
        ksprintf(s, "<thead>\n");
        ksprintf(s, "<tr><th class='col1'>Amplicon</th>");
        ksprintf(s, "<th class='col2'>Edited Rate</th>");
        if(mOptions->recstat){
            ksprintf(s, "<th class='col2'>Total Recomb Rate</th>");
            ksprintf(s, "<th class='col2'>Exact Recomb Rate</th>");
            ksprintf(s, "<th class='col2'>Recomb Hit Rate</th>");
        }
        ksprintf(s, "<th class='col2'>Overall Mutated Rate</th>");
        ksprintf(s, "<th class='col2'>Qualfied Sequnces</th>");
        ksprintf(s, "<th class='col2'>Unqualified Rate</th></tr>\n");
        ksprintf(s, "</thead>\n");
        ksprintf(s, "<tbody>\n");

        std::string subsect = "Variant counts along amplicon reference sequence";
        std::string vcid = util::replace(subsect, " ", "_");
        subsect = "Top" + std::to_string(mOptions->edo.topn) + " allele sequences";
        std::string tnid = util::replace(subsect, " ", "_");
        std::string sse = "Edit efficience computation";
        std::string ssd = "plot_" + util::replace(sse, " ", "_");
        std::vector<std::string> bgcolor = {"red", "yellow", "white"};
        int bgidx = 2;
        for(size_t i = 0; i < mOptions->samples.size()-1; ++i){
            bgidx = 2;
            if(mGEDetectors[i]->totcnt < mOptions->edo.minqs1k) bgidx = 0;
            else if(mGEDetectors[i]->totcnt < mOptions->edo.minqs10k) bgidx = 1;
            std::string linkk = "<a href=\"" + mOptions->hmo.subdir + "/" + mOptions->samples[i] + ".html\">" + mOptions->samples[i] + "</a>";
            std::string linkv = "<a href=\"" + mOptions->hmo.subdir + "/" + mOptions->samples[i] + ".html#" + tnid + "\">" + std::to_string(mGEDetectors[i]->edieff*100) + "%</a>";
            std::string recs = std::to_string(mGEDetectors[i]->receff*100) + "%";
            std::string rees = std::to_string(mGEDetectors[i]->reeeff*100) + "%";
            std::string reas = std::to_string(mGEDetectors[i]->recpef*100) + "%";
            std::string linkr = "<a href=\"" + mOptions->hmo.subdir + "/" + mOptions->samples[i] + ".html#" + vcid + "\">" + std::to_string(mGEDetectors[i]->muteff*100) + "%</a>";
            std::string linkq = "<a href=\"" + mOptions->hmo.subdir + "/" + mOptions->samples[i] + ".html#" + ssd + "\">" + std::to_string(mGEDetectors[i]->totcnt) + "</a>";
            std::string linkd = "<a href=\"" + mOptions->hmo.subdir + "/" + mOptions->samples[i] + ".html#" + ssd + "\">" + std::to_string((double)mGEDetectors[i]->dropcnt/(double)mGEDetectors[i]->allcnt*100) + "%</a>";
            ksprintf(s, "<tr><td bgcolor='%s' class='col1'>%s</td>", bgcolor[bgidx].c_str(), linkk.c_str());
            ksprintf(s, "<td bgcolor='%s' class='col2'>%s</td>", bgcolor[bgidx].c_str(), linkv.c_str());
            if(mOptions->recstat){
                ksprintf(s, "<td bgcolor='%s' class='col2'>%s</td>", bgcolor[bgidx].c_str(), recs.c_str());
                ksprintf(s, "<td bgcolor='%s' class='col2'>%s</td>", bgcolor[bgidx].c_str(), rees.c_str());
                ksprintf(s, "<td bgcolor='%s' class='col2'>%s</td>", bgcolor[bgidx].c_str(), reas.c_str());
            }
            ksprintf(s, "<td bgcolor='%s' class='col2'>%s</td>", bgcolor[bgidx].c_str(), linkr.c_str());
            ksprintf(s, "<td bgcolor='%s' class='col2'>%s</td>", bgcolor[bgidx].c_str(), linkq.c_str());
            ksprintf(s, "<td bgcolor='%s' class='col2'>%s</td></tr>", bgcolor[bgidx].c_str(), linkd.c_str());
        }
        ksprintf(s, "</tbody>\n");
        ksprintf(s, "</table>\n");
        mOptions->hmo.printExportButtons(s, "eff_table_xx", "gedit.eff.summary");
        ksprintf(s, "</div></div></div>\n");

    }
    // Detailed report ends
    printFooter(s);
    // write to file
    FILE* fp = fopen(mOptions->hmr.c_str(), "w");
    fwrite(s->s, sizeof(char), s->l, fp);
    fclose(fp);
    // release res
    free(s->s); free(s);
}

void Processor::reportQCAllHTML(){
    kstring_t* s = (kstring_t*)calloc(1, sizeof(kstring_t));
    printHeader(s);
    QcStat* qs1 = mQcStats1;
    QcStat* qs2 = mQcStats2;
    SplitResult* spr = mSplitResults;
    // report title
    ksprintf(s, "<h1 style='text-align:left;'><a href='https://github.com/wulj2/getools.git' target='_blank' style='color:#000099;text-decoration:none;'>Library split report by getools</a></h1>\n");
    // Summary beg
    ksprintf(s, "<div class='section_div'>\n");
    ksprintf(s, "<div class='section_title' onclick=showOrHide('summary')><a name='summary'>Summary</a></div>\n");
    ksprintf(s, "<div id='summary'>\n");
    // Summary -> General beg
    ksprintf(s, "<div class='subsection_title' onclick=showOrHide('general')>General</div>\n");
    ksprintf(s, "<div id='general'>\n");
    ksprintf(s, "<table class='summary_table'>\n");
    {
        // getools version
        std::string verinfo = mOptions->soft->ver + " (<a href='https://github.com/wulj2/getools.git'>https://github.com/wulj2/getools.git)</a>"; 
        htmlutil::outTableRow(s, "getools version:",  verinfo.c_str());
        // sequencing
        std::string cycinfo = "";
        if(qs2){
            cycinfo = "paired end";
            cycinfo.append("(" + std::to_string(qs1->mCycle) + " cycles + " + std::to_string(qs2->mCycle) + " cycles)");
        }else{
            cycinfo = "single end(" + std::to_string(qs1->mCycle) + " cycles)";
        }
        htmlutil::outTableRow(s, "sequencing:", cycinfo.c_str());
        // mean length
        std::string meanls = std::to_string(qs1->mMeanReadLen) + "bp";
        if(qs2) meanls.append(", " + std::to_string(qs2->mMeanReadLen) + "bp");
        htmlutil::outTableRow(s, "mean read length:", meanls.c_str());
    }
    ksprintf(s, "</table>\n");
    ksprintf(s, "</div>\n");
    // Summary -> General end
    // Summary -> Fastq QC beg
    ksprintf(s, "<div class='subsection_title' onclick=showOrHide('qc')>Fastq QC</div>\n");
    ksprintf(s, "<div id='qc'>\n");
    ksprintf(s, "<table class='summary_table'>\n");
    {
        // total reads
        std::string totrs = htmlutil::formatNumber(qs1->mReads);
        if(qs2) totrs.append(", " + htmlutil::formatNumber(qs2->mReads));
        htmlutil::outTableRow(s, "total reads:", totrs.c_str());
        // total bases
        std::string totbs = htmlutil::formatNumber(qs1->mBases);
        if(qs2) totbs.append(", " + htmlutil::formatNumber(qs2->mBases));
        htmlutil::outTableRow(s, "total bases:",  totbs.c_str());
        // Q20 bases
        std::string tq20 = htmlutil::formatNumber(qs1->mQ20Total);
        tq20.append("(" + std::to_string((double)qs1->mQ20Total/(double)qs1->mBases*100) + "%)");
        if(qs2){
            tq20.append(", " + htmlutil::formatNumber(qs2->mQ20Total));
            tq20.append("(" + std::to_string((double)qs2->mQ20Total/(double)qs2->mBases*100) + "%)");
        }
        htmlutil::outTableRow(s, "Q20 bases:",  tq20.c_str());
        // Q30 bases
        std::string tq30 = htmlutil::formatNumber(qs1->mQ30Total);
        tq30.append("(" + std::to_string((double)qs1->mQ30Total/(double)qs1->mBases*100) + "%)");
        if(qs2){
            tq30.append(", " + htmlutil::formatNumber(qs2->mQ30Total));
            tq30.append("(" + std::to_string((double)qs2->mQ30Total/(double)qs2->mBases*100) + "%)");
        }
        htmlutil::outTableRow(s, "Q30 bases:",  tq30.c_str());
        // GC contents
        std::string gcs = std::to_string(qs1->mGCPercent*100) + "%";
        if(qs2) gcs.append(", " + std::to_string(qs2->mGCPercent*100) + "%");
        htmlutil::outTableRow(s, "GC Content:", gcs.c_str());
    }
    ksprintf(s, "</table>\n");
    ksprintf(s, "</div>\n");
    // Summary -> Fastq QC end
    // Split Stat beg
    ksprintf(s, "<div class='subsection_title' onclick=showOrHide('gsplrpt')>Pooled library split</div>\n");
    ksprintf(s, "<div id='gsplrpt'>\n");
    ksprintf(s, "<table class='summary_table'>\n");
    {
        std::string prefix = "Reads ";
        if(qs2) prefix = "Read pairs ";
        // Reads Input
        std::string inrs = htmlutil::formatNumber(spr->mTotalReads);
        std::string key = prefix + "to split:";
        htmlutil::outTableRow(s, key, inrs);
        // Reads Got
        std::string rgot = htmlutil::formatNumber(spr->mTotalReads-spr->mDropReads);
        rgot.append("(" + std::to_string((double)(spr->mTotalReads - spr->mDropReads)/(double)spr->mTotalReads*100) + "%)");
        key = prefix + "got:";
        htmlutil::outTableRow(s, key, rgot);
        // FR stat
        if(qs2){
            std::string frgv = "FR: " + htmlutil::formatNumber(spr->mTotalFR);
            frgv.append("(" + std::to_string((double)spr->mTotalFR/(double)spr->mTotalReads*100) + "%),");
            frgv.append("RF: " + htmlutil::formatNumber(spr->mTotalRF));
            frgv.append("(" + std::to_string((double)spr->mTotalRF/(double)spr->mTotalReads*100) + "%)");
            key = prefix + " strands:";
            htmlutil::outTableRow(s, key, frgv);
        }
    }
    ksprintf(s, "</table>\n");
    ksprintf(s, "</div>\n");
    // Summary -> Split end
    ksprintf(s, "</div>\n");
    ksprintf(s, "</div>\n");
    // Summary end
    // Fastq QC beg
    ksprintf(s, "<div class='section_div'>\n");
    ksprintf(s, "<div class='section_title' onclick=showOrHide('qcplot')><a name='pqcsummary'>Pooled Fastq QC</a></div>\n");
    ksprintf(s, "<div id='qcplot'>\n");
    if(qs1) qs1->reportHTML(s, 1);
    if(qs2) qs2->reportHTML(s, 2);
    ksprintf(s, "</div>\n");
    ksprintf(s, "</div>\n");
    // Split QC beg
    ksprintf(s, "<div class='section_div'>\n");
    ksprintf(s, "<div class='section_title' onclick=showOrHide('splrpt')><a name='lsrsummary'>Library split result statistics</a></div>\n");
    ksprintf(s, "<div id='splrpt'>\n");
    mSplitResults->reportHTML(s);
    // Split QC end
    // Detailed report begs
    ksprintf(s, "<div class='subsection_title' onclick=showOrHide('assign_rate_table')>Links to detailed report of each amplicons</div>\n");
    ksprintf(s, "<div id='assign_rate_table'>\n");
    ksprintf(s, "<div class='sub_section_tips'>Click to go to its detailed report</div>\n");
    ksprintf(s, "<table class='summary_table' id='split_table_xx'>\n");
    ksprintf(s, "<thead>\n");
    ksprintf(s, "<tr><th class='col1'>Amplicon</th>");
    ksprintf(s, "<th class='col2'>Number</th>");
    ksprintf(s, "<th class='col2'>Percents</th>\n</tr>\n</thead>");
    ksprintf(s, "<tdoby>\n");
    for(size_t i = 0; i < mOptions->samples.size()-1; ++i){
        std::string linkv = "<a href=\"" + mOptions->hmo.subdir + "/" + mOptions->samples[i] + ".html\">" + std::to_string(mSplitResults->mSplitCount[i]/(double)mSplitResults->mTotalReads * 100) + "%</a>";
        ksprintf(s, "<tr><td class='col2'>%s</td>", mOptions->samples[i].c_str());
        ksprintf(s, "<td class='col2'>%lld</td>", mSplitResults->mSplitCount[i]);
        ksprintf(s, "<td class='col2'>%s</td>", linkv.c_str());
        ksprintf(s, "</tr>\n");
    }
    ksprintf(s, "</tbody></table>\n");
    mOptions->hmo.printExportButtons(s, "split_table_xx", "gesplit.summary");
    ksprintf(s, "</div>\n");
    ksprintf(s, "</div>\n");
    ksprintf(s, "</div>\n");
    // Detailed report ends
    printFooter(s);
    // write to file
    FILE* fp = fopen(mOptions->hmr.c_str(), "w");
    fwrite(s->s, sizeof(char), s->l, fp);
    fclose(fp);
    // release res
    free(s->s); free(s);
}

void Processor::printFooter(kstring_t* s){
    // end <body><container...
    ksprintf(s, "</div>\n</body>\n");
    // begin footer
    ksprintf(s, "<div id='footer'> ");
    ksprintf(s, "<p>command: %s %s</p>", mOptions->soft->nam.c_str(), mOptions->soft->cmd.c_str());
    ksprintf(s, "<p>version: %s %s</p>", mOptions->soft->nam.c_str(), mOptions->soft->ver.c_str()); 
    ksprintf(s, "<p>workdir: %s</p>", mOptions->soft->cwd.c_str());
    ksprintf(s, "<p>runtime: %s</p>", mOptions->soft->getExecutionTime().c_str());
    ksprintf(s, "<p>endtime: %s</p>", util::currentTime().c_str());
    ksprintf(s, "</div>\n");
    // end footer
    // end html
    ksprintf(s, "</html>\n");
}
