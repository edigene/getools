#include "derep.h"

int DeReper::seqcmp(char* a, char* b, int n){
    if(n <= 0) return 0;
    char* p = a;
    char* q = b;
    while((n-- > 0) && (nuc_to_3bit[(int)(*p)] == nuc_to_3bit[(int)(*q)])){
        if((n == 0) || (*p == 0) || (*q == 0)) break;
        p++; q++;
    }
    return nuc_to_3bit[(int)(*p)] - nuc_to_3bit[(int)(*q)];
}

void DeReper::reverse_complement(char* rc, char* seq, int64_t len){
    for(int64_t i = 0; i < len; ++i){
        rc[i] = nuc_to_cmp[(int)seq[len-1-i]];
    }
    rc[len] = 0;
}

void DeReper::string_normalize(char* normalized, char* s, unsigned int len){
    char *p = s;
    char *q = normalized;
    for(unsigned int i = 0; i < len; ++i){
        *q++ = nuc_to_norm[(int)(*p++)];
    }
    *q = 0;
}

void DeReper::rehash(){
    // double size
    uint64_t new_hashtablesize = 2 * hashtablesize;
    uint64_t new_hash_mask = new_hashtablesize - 1;
    bucket_t* new_hashtable = (bucket_t*)malloc(sizeof(bucket_t)*new_hashtablesize);
    memset(new_hashtable, 0, sizeof(bucket_t)*new_hashtablesize);
    // rehash
    for(uint64_t i = 0; i < hashtablesize; ++i){
        bucket_t* old_bp = hashtable + i;
        if(old_bp->size){
            uint64_t k = old_bp->hash & new_hash_mask;
            while(new_hashtable[k].size) k = (k+1) & new_hash_mask;
            bucket_t* new_bp = new_hashtable + k;
            *new_bp = *old_bp;
        }
    }
    // update
    free(hashtable);
    hashtable = new_hashtable;
    alloc_clusters *= 2;
    hashtablesize = new_hashtablesize;
    hash_mask = new_hash_mask;
}

bool DeReper::derep_rec(const krec1_t* r){
    return derep_rec(r->seq.s, r->seq.l, r->qual.s, false, NULL, 0);
}

bool DeReper::derep_rec(char* seq, int seqlen, char* qual, bool qrev,  char* extseq, int extlen, bool appext){
    // filter
    if(seqlen < opt_minseqlength){
        ++discarded_short;
        return true;
    }
    if(seqlen > opt_maxseqlength){
        ++discarded_long;
        return true;
    }
    // stat
    nucleotidecount += seqlen;
    if(seqlen > longest) longest = seqlen;
    if(seqlen < shortest) shortest = seqlen;
    // realloc if necessary
    if(seqlen > alloc_seqlen){
        alloc_seqlen = seqlen;
        seq_up = (char*)realloc(seq_up, (alloc_seqlen+1)*sizeof(char));
        rc_seq_up = (char*)realloc(rc_seq_up, (alloc_seqlen+1)*sizeof(char));
    }
    // rehash if necessary
    if(clusters + 1 > alloc_clusters) rehash();
    // do dedup
    // normalize sequence: uppercase and replace U by T
    string_normalize(seq_up, seq, seqlen);
    // reverse complement if necessary
    if(opt_strand > 1) reverse_complement(rc_seq_up, seq_up, seqlen);
    // find free bucket or bucket for identical sequence
    uint64_t hash = CityHash64(seq_up, seqlen);
    uint64_t hash_ext = 0;
    if((!appext) && extseq && extlen){
        hash_ext = CityHash64(extseq, extlen);
        hash ^= hash_ext;
    }
    uint64_t j = hash & hash_mask;
    bucket_t *bp = hashtable + j;
    while((bp->size) && 
          ((hash != bp->hash) || 
           (seqcmp(seq_up, bp->seq, seqlen)) || 
           ((!appext) && extseq && extlen && strcmp(extseq, bp->ext)))){
        j = (j+1) & hash_mask;
        bp = hashtable + j;
    }
    if((opt_strand > 1) && !bp->size){ // no match on plus strand, check minus strand if needed
        uint64_t rc_hash = CityHash64(rc_seq_up, seqlen);
        if((!appext) && extseq && extlen) rc_hash ^= hash_ext;
        uint64_t k = rc_hash & hash_mask;
        bucket_t* rc_bp = hashtable + k;
        while((rc_bp->size) &&
              ((rc_hash != rc_bp->hash) ||
               (seqcmp(rc_seq_up, rc_bp->seq, seqlen)) || 
               ((!appext) && extseq && extlen && strcmp(extseq, bp->ext)))){
            k = (k+1) & hash_mask;
            rc_bp = hashtable + k;
        }
        if(rc_bp->size){
            bp = rc_bp;
            j = k;
            qrev = (!qrev);
        }
    }
    if(bp->size){ // dup
        bp->seqno_last = sequencecount;
        if(qual && bp->size < MAX_QUAL_DEP){
            if(qrev){
                for(int i = 0; i < seqlen; ++i) bp->qual[seqlen-1-i] += qual[i];
            }else{
                for(int i = 0; i < seqlen; ++i) bp->qual[i] += qual[i];
            }
        }
        bp->size += 1;
        ++sequencecount;
        if(appext && extseq && extlen){
            if(!bp->emap) bp->emap = kh_init(es);
            khiter_t k = kh_get(es, bp->emap, extseq);
            if(k == kh_end(bp->emap)){
                int r;
                k = kh_put(es, bp->emap, strndup(extseq, extlen), &r);
                kh_val(bp->emap, k) = 1;
            }else{
                kh_val(bp->emap, k) += 1;
            }
        }
        return true;
    }else{ // newly
        bp->size = 1;
        bp->hash = hash;
        bp->seqno_first = sequencecount;
        bp->seqno_last = sequencecount;
        bp->seq = strndup(seq, seqlen);
        if(qual){
            bp->qual = (uint32_t*)malloc(seqlen*sizeof(uint32_t));
            if(qrev){
                for(int i = 0; i < seqlen; ++i) bp->qual[seqlen-1-i] = qual[i];
            }else{
                for(int i = 0; i < seqlen; ++i) bp->qual[i] = qual[i];
            }
        }
        if(extseq && extlen){
            if(appext){
                if(!bp->emap) bp->emap = kh_init(es);
                int r;
                khiter_t k = kh_get(es, bp->emap, extseq);
                k = kh_put(es, bp->emap, strndup(extseq, extlen), &r);
                kh_val(bp->emap, k) = 1;
            }else{
                bp->ext = strndup(extseq, extlen);
            }
        }

        bp->len = seqlen;
        bp->id = clusters;
        ++clusters;
        ++sequencecount;
        return false;
    }
}

void DeReper::summary(){
    if(summarized) return;
    if(sequencecount > 0) duprate = 1-(double)clusters/(double)sequencecount;
    maxidx = 0;
    for(uint64_t i = 0; i < hashtablesize; ++i){
        if(hashtable[i].size > maxidx) maxidx = hashtable[i].size;
    }
    dupdist = (int64_t*)calloc(maxidx+1, sizeof(int64_t));
    for(uint64_t i = 0; i < hashtablesize; ++i){
        if(hashtable[i].size) dupdist[hashtable[i].size] += hashtable[i].size;
    }
    minidx = 0;
    modidx = 0;
    int64_t old_sum = 0;
    int64_t acc_sum = old_sum = 0;
    int64_t quat_arr[4] = {0};
    for(int i = 1; i <= 4; ++i) quat_arr[i-1] = sequencecount * (0.25 * i);
    for(int64_t i = 0; i <= maxidx; ++i){
        if(dupdist[i]){
            if(minidx == 0) minidx = i;
            if(dupdist[i] > dupdist[modidx]) modidx = i;
            old_sum = acc_sum;
            acc_sum += dupdist[i];
            for(int q = 0; q < 4; ++q){
                if(quat_arr[q] > old_sum && quat_arr[q] <= acc_sum){
                    quat_dup[q] = i;
                }
            }
        }
    }
    summarized = true;
}

void DeReper::reportJSON(kstring_t* s, const char* dh, const char* dm){
    summary();
    ksprintf(s, "%s\"DerepSummary\": {\n", dh);
    ksprintf(s, "%s%s\"TotalSeqs\": %lld,\n", dh, dm, sequencecount);
    ksprintf(s, "%s%s\"UniqSeqs\": %lld,\n", dh, dm, clusters);
    ksprintf(s, "%s%s\"DupRate\": %lf,\n", dh, dm, duprate);
    // dup dist
    {
        ksprintf(s, "%s%s\"DupCountMax\": %lld,\n", dh, dm, maxidx);
        ksprintf(s, "%s%s\"DupCountMin\": %lld,\n", dh, dm, minidx);
        ksprintf(s, "%s%s\"DupCountMod\": %lld,\n", dh, dm, modidx);
        ksprintf(s, "%s%s\"DupCountQuartile\": [", dh, dm);
        for(int i = 0; i < 4; ++i){
            ksprintf(s, "%lld,", quat_dup[i]);
        }
        s->s[s->l-1] = ']';
        ksprintf(s, ",\n");
        ksprintf(s, "%s%s\"DupCountDist\": [", dh, dm);
        for(int i = 0; i <= maxidx; ++i){
            ksprintf(s, "%lld,", dupdist[i]);
        }
        if(s->s[s->l-1] == ',') s->s[s->l-1] = ']';
        else kputc(']', s);
    }
    ksprintf(s, "\n%s}", dh);
}

void DeReper::reportHTML(kstring_t* s){
    summary();
    std::string subsect = "Deduplication of sequences";
    std::string divName = util::replace(subsect, " ", "_");
    std::string title = "'duplication rate (" + std::to_string(duprate*100) + "%)'";
    ksprintf(s, "<div class='subsection_title'><a title='click to hide/show' onclick=showOrHide('%s')>%s</a></div>\n", divName.c_str(), subsect.c_str());
    ksprintf(s, "<div id='%s'>\n", divName.c_str());
    ksprintf(s, "<div class='sub_section_tips'>Percentage of each repeat level will be shown on mouse over.</div>\n");
    ksprintf(s, "<div class='figure' id='plot_%s'></div>\n", divName.c_str());
    ksprintf(s, "</div>\n");
    ksprintf(s, "\n<script type=\"text/javascript\">\n");

    std::string jsnstr;
    jsnstr.append("var repie = {\n");
    jsnstr.append("  values: [");
    for(int i = 0; i <= maxidx; ++i){
        if(dupdist[i]) jsnstr.append(std::to_string(dupdist[i]) + ",");
    }
    jsnstr.append("],\n");
    jsnstr.append("  labels: [");
    for(int i = 0; i <= maxidx; ++i){
        if(dupdist[i]) jsnstr.append("'rep#" + std::to_string(i) + "',");
    }
    jsnstr.append("],\n");
    jsnstr.append("  type: 'pie',\n");
    jsnstr.append("  textinfo: 'label',\n");
    jsnstr.append("  textposition: 'inside',\n");
    jsnstr.append("  hoverinfo: 'label+percent',\n");
    jsnstr.append("  insidetextorientation: 'radial',\n");
    jsnstr.append("};\n");

    jsnstr.append("var data = [repie];\n");
    jsnstr.append("var layout = {\n");
    jsnstr.append("  title:" + title + ",\n");
    jsnstr.append("  showlegend: false,\n");
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

void DeReper::tsvHead(kstring_t* s){
    ksprintf(s, "%s\t%s\t%s", "TotalSeqs", "UniqueSeqs", "DupRate");
}

void DeReper::tsvBody(kstring_t* s){
    summary();
    ksprintf(s, "%lld\t%lld\t%lf", sequencecount, clusters, duprate);
}

void derep_usage(const char* arg0){
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage: %s %s [options]\n\n", PACKAGE_NAME, arg0);
    fprintf(stderr, "Options: -i FILE input fastq file\n");
    fprintf(stderr, "         -o FILE output tsv file\n");
    fprintf(stderr, "         -r      reverse complement sequences counted as duplication\n");
    fprintf(stderr, "\n");
}

int derep_main(int argc, char** argv){
    if(argc == 1){
        derep_usage(argv[0]);
        return 0;
    }
    const char* infq = NULL;
    const char* outf = NULL;
    bool rev_derep = false;
    int c = -1;
    while((c = getopt(argc, argv, "i:o:urh")) >= 0){
        switch(c){
            case 'i': infq = strdup(optarg); break;
            case 'o': outf = strdup(optarg); break;
            case 'r': rev_derep = true; break;
            case 'h': derep_usage(argv[0]); return 0; break;
            default: break;
        }
    }
    if(infq == NULL || outf == NULL){
        fprintf(stderr, "input/output file must be provided\n");
        return 1;
    }
    DeReper *dr = new DeReper();
    dr->opt_strand = (rev_derep ? 2 : 1);
    gzFile fp = gzopen(infq, "r");
    kseq1_t* ks = kseq1_init(fp);
    krec1_t* kr = krec1_init();
    while(kseq1_read(ks, kr) >= 0) dr->derep_rec(kr);
    fprintf(stderr, "total seqs: %llu\n", dr->sequencecount);
    fprintf(stderr, "uniq  seqs: %llu\n", dr->clusters);
    FILE* fout = fopen(outf, "w");
    bucket_t::out_bucket_head(fout);
    for(uint64_t i = 0; i < dr->hashtablesize; ++i){
        dr->hashtable[i].out_bucket_body(fout);
    }
    delete dr;
    fclose(fout);
    gzclose(fp);
    kseq1_destroy(ks);
    krec1_destroy(kr);
    return 0;
}
