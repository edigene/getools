#include "mergepairs.h"

int PairMerger::get_qual(char q){
    int qual = q - opt_fastq_ascii;
    if(qual < opt_fastq_qmin || qual > opt_fastq_qmax) return -1;
    else return qual;
}

double PairMerger::q_to_p(int q){
  int x = q - opt_fastq_ascii;
  if(x < 2) return 0.75;
  else return exp10(-x/10.0);
}

void PairMerger::precompute_qual(){
  /* Precompute tables of scores etc */
  for(int x = 33; x <= 126; x++){
      double px = q_to_p(x);
      q2p[x] = px;
      for(int y = 33; y <= 126; y++){
          double py = q_to_p(y);
          double p, q;
          /* Quality score equations from Edgar & Flyvbjerg (2015) */
          /* Match */
          p = px * py / 3.0 / (1.0 - px - py + 4.0 * px * py / 3.0);
          q = round(-10.0 * log10(p));
          q = MIN(q, opt_fastq_qmaxout);
          q = MAX(q, opt_fastq_qminout);
          merge_qual_same[x][y] = opt_fastq_ascii + q;
          /* Mismatch, x is highest quality */
          p = px * (1.0 - py / 3.0) / (px + py - 4.0 * px * py / 3.0);
          q = round(-10.0 * log10(p));
          q = MIN(q, opt_fastq_qmaxout);
          q = MAX(q, opt_fastq_qminout);
          merge_qual_diff[x][y] = opt_fastq_ascii + q;
          // observed match,
          // p = probability that they truly are identical,
          // given error probabilites of px and py, resp.
          // Given two initially identical aligned bases, and
          // the error probabilities px and py,
          // what is the probability of observing a match (or a mismatch)?
          p = 1.0 - px - py + px * py * 4.0 / 3.0;
          match_score[x][y] = log2(p/0.25);
          // Use a minimum mismatch penalty
          mism_score[x][y] = MIN(log2((1.0-p)/0.75), merge_mismatchmax);
        }
    }
}

void PairMerger::merge_sym(char* sym, char* qual, char fwd_sym, char rev_sym, char fwd_qual, char rev_qual){
    if(rev_sym == 'N'){
        *sym = fwd_sym;
        *qual = fwd_qual;
    }else if(fwd_sym == 'N'){
        *sym = rev_sym;
        *qual = rev_qual;
    }else if (fwd_sym == rev_sym){
        /* agreement */
        *sym = fwd_sym;
        *qual = merge_qual_same[(unsigned)fwd_qual][(unsigned)rev_qual];
    }else{
        /* disagreement */
        if(fwd_qual > rev_qual){
            *sym = fwd_sym;
            *qual = merge_qual_diff[(unsigned)fwd_qual][(unsigned)rev_qual];
        }else{
            *sym = rev_sym;
            *qual = merge_qual_diff[(unsigned)rev_qual][(unsigned)fwd_qual];
        }
    }
}

void PairMerger::keep(){
    ++merged;
    sum_fragment_length += ip->merged_length;
    sum_squared_fragment_length += ip->merged_length * ip->merged_length;
    sum_ee_merged += ip->ee_merged;
    sum_ee_fwd += ip->ee_fwd;
    sum_ee_rev += ip->ee_rev;
    sum_errors_fwd += ip->fwd_errors;
    sum_errors_rev += ip->rev_errors;
    ++merge_ops_cnt[ip->reason];
    ++merge_len_cnt[ip->merged_length];
    if(ip->fwd_length > ip->merged_length || ip->rev_length > ip->merged_length) ++peovhang;
}

void PairMerger::discard(){
    ++notmerged;
    ++merge_ops_cnt[ip->reason];
}

void PairMerger::merge(){
    /* length of 5' overhang of the forward sequence not merged with the reverse sequence */
    int64_t fwd_5prime_overhang = ip->fwd_trunc > ip->offset ? ip->fwd_trunc - ip->offset : 0;
    ip->ee_merged = 0.0;
    ip->ee_fwd = 0.0;
    ip->ee_rev = 0.0;
    ip->fwd_errors = 0;
    ip->rev_errors = 0;
    char sym, qual;
    char fwd_sym, fwd_qual, rev_sym, rev_qual;
    int64_t fwd_pos, rev_pos, merged_pos;
    double ee;
    merged_pos = 0;
    // 5' overhang in forward sequence
    fwd_pos = 0;
    while(fwd_pos < fwd_5prime_overhang){
        sym = ip->fwd_sequence[fwd_pos];
        qual = ip->fwd_quality[fwd_pos];
        ip->merged_sequence[merged_pos] = sym;
        ip->merged_quality[merged_pos] = qual;
        ee = q2p[(unsigned)qual];
        ip->ee_merged += ee;
        ip->ee_fwd += ee;
        ++fwd_pos;
        ++merged_pos;
    }

    // Merged region
    int64_t rev_3prime_overhang = ip->offset > ip->fwd_trunc ? ip->offset - ip->fwd_trunc : 0;
    rev_pos = ip->rev_trunc - 1 - rev_3prime_overhang;
    while((fwd_pos < ip->fwd_trunc) && (rev_pos >= 0)){
        fwd_sym = ip->fwd_sequence[fwd_pos];
        rev_sym = nuc_to_cmp[(int)(ip->rev_sequence[rev_pos])];
        fwd_qual = ip->fwd_quality[fwd_pos];
        rev_qual = ip->rev_quality[rev_pos];
        merge_sym(&sym, &qual, fwd_sym, rev_sym, fwd_qual, rev_qual);
        if(sym != fwd_sym) ip->fwd_errors++;
        if(sym != rev_sym) ip->rev_errors++;
        ip->merged_sequence[merged_pos] = sym;
        ip->merged_quality[merged_pos] = qual;
        ip->ee_merged += q2p[(unsigned)qual];
        ip->ee_fwd += q2p[(unsigned)fwd_qual];
        ip->ee_rev += q2p[(unsigned)rev_qual];
        ++fwd_pos;
        --rev_pos;
        ++merged_pos;
    }
    
    // 5' overhang in reverse sequence
    while(rev_pos >= 0){
        sym = nuc_to_cmp[(int)(ip->rev_sequence[rev_pos])];
        qual = ip->rev_quality[rev_pos];
        ip->merged_sequence[merged_pos] = sym;
        ip->merged_quality[merged_pos] = qual;
        merged_pos++;
        ee = q2p[(unsigned)qual];
        ip->ee_merged += ee;
        ip->ee_rev += ee;
        --rev_pos;
    }
    int64_t mergelen = merged_pos;
    ip->merged_length = mergelen;
    ip->merged_sequence[mergelen] = 0;
    ip->merged_quality[mergelen] = 0;
    if(ip->ee_merged <= opt_fastq_maxee){
        ip->reason = okbykmers;
        ip->merged = 1;
    }else{
        ip->reason = maxee;
    }
}

int64_t PairMerger::optimize(){
  /* ungapped alignment in each diagonal */
  int64_t tl = ip->fwd_trunc + ip->rev_trunc;
  int64_t i1 = 1;
  int64_t i2 = tl - 1;
  double best_score = 0.0;
  int64_t best_i = 0;
  int64_t best_diffs = 0;
  int hits = 0;
  int kmers = 0;
  int diags[tl];
  memset(diags, 0, tl * sizeof(int));
  kh_insert_kmers(kmerhash, k, ip->fwd_sequence, ip->fwd_trunc);
  kh_find_diagonals(kmerhash, k, ip->rev_sequence, ip->rev_trunc, diags);
  for(int64_t i = i1; i <= i2; i++){
      int diag = tl - i;
      int diagcount = diags[diag];
      if(diagcount >= merge_mindiagcount){// possible overlap region
          kmers = 1;
          if(diagcount >= merge_minrepeatdiagcount){// possible repeat region
              hits++;
              if(hits > 1) break;// really repeat region found
          }
          /* for each interesting diagonal */
          int64_t fwd_3prime_overhang = i > ip->rev_trunc ? i - ip->rev_trunc : 0; // fwd 3' unmapped out of rev
          int64_t rev_3prime_overhang = i > ip->fwd_trunc ? i - ip->fwd_trunc : 0; // rev 3' unmapped out of fwd
          int64_t overlap = i - fwd_3prime_overhang - rev_3prime_overhang;
          int64_t fwd_pos_start = ip->fwd_trunc - fwd_3prime_overhang - 1;
          int64_t rev_pos_start = ip->rev_trunc - rev_3prime_overhang - overlap;
          int64_t fwd_pos = fwd_pos_start;
          int64_t rev_pos = rev_pos_start;
          double score = 0.0;
          int64_t diffs = 0;
          double score_high = 0.0;
          double dropmax = 0.0;
          for(int64_t j = 0; j < overlap; ++j){
              /* for each pair of bases in the overlap */
              char fwd_sym = ip->fwd_sequence[fwd_pos];
              char rev_sym = nuc_to_cmp[(int)(ip->rev_sequence[rev_pos])];
              unsigned int fwd_qual = ip->fwd_quality[fwd_pos];
              unsigned int rev_qual = ip->rev_quality[rev_pos];
              fwd_pos--;
              rev_pos++;
              if(fwd_sym == rev_sym){
                  score += match_score[fwd_qual][rev_qual];
                  if(score > score_high) score_high = score;
              }else{
                  score += mism_score[fwd_qual][rev_qual];
                  diffs++;
                  if(score < score_high - dropmax) dropmax = score_high - score;
              }
          }
          if(dropmax >= merge_dropmax) score = 0.0;
          if(score > best_score){
              best_score = score;
              best_i = i;
              best_diffs = diffs;
            }
        }
  }
  if(hits > 1){
      ip->reason = repeat;
      return 0;
  }
  if((!opt_fastq_allowmergestagger) && (best_i > ip->fwd_trunc)){
      ip->reason = staggered;
      return 0;
  }
  if(best_diffs > opt_fastq_maxdiffs){
      ip->reason = maxdiffs;
      return 0;
  }
  if((100.0 * best_diffs / best_i) > opt_fastq_maxdiffpct){
      ip->reason = maxdiffpct;
      return 0;
  }
  if(kmers == 0){
      ip->reason = nokmers;
      return 0;
  }
  if(best_score < merge_minscore){
      ip->reason = minscore;
      return 0;
  }
  if(best_i < opt_fastq_minovlen){
      ip->reason = minovlen;
      return 0;
  }
  int mergelen = ip->fwd_trunc + ip->rev_trunc - best_i;
  if(mergelen < opt_fastq_minmergelen){
      ip->reason = minmergelen;
      return 0;
  }
  if(mergelen > opt_fastq_maxmergelen){
      ip->reason = maxmergelen;
      return 0;
  }
  return best_i;
}

void PairMerger::process(){
  ip->merged = 0;
  bool skip = 0;
  /* check length */
  if((ip->fwd_length < opt_fastq_minlen) || (ip->rev_length < opt_fastq_minlen)){
      ip->reason = minlen;
      skip = 1;
  }
  if((ip->fwd_length > opt_fastq_maxlen) || (ip->rev_length > opt_fastq_maxlen)){
      ip->reason = maxlen;
      skip = 1;
  }
  /* truncate sequences by quality */
  int64_t fwd_trunc = ip->fwd_length;
  if(!skip && opt_trunc_by_qual){
      for(int64_t i = 0; i < ip->fwd_length; i++)
          if(get_qual(ip->fwd_quality[i]) <= opt_fastq_truncqual){
              fwd_trunc = i;
              break;
          }
      if(fwd_trunc < opt_fastq_minlen){
          ip->reason = minlen;
          skip = 1;
      }
  }
  ip->fwd_trunc = fwd_trunc;
  int64_t rev_trunc = ip->rev_length;
  if(!skip && opt_trunc_by_qual){
      for(int64_t i = 0; i < ip->rev_length; i++)
          if(get_qual(ip->rev_quality[i]) <= opt_fastq_truncqual){
              rev_trunc = i;
              break;
          }
      if(rev_trunc < opt_fastq_minlen){
          ip->reason = minlen;
          skip = 1;
      }
  }
  ip->rev_trunc = rev_trunc;
  /* count n's */
  /* replace quality of N's by zero */
  if(!skip && opt_filter_nbase){
      int64_t fwd_ncount = 0;
      for(int64_t i = 0; i < fwd_trunc; i++){
          if(ip->fwd_sequence[i] == 'N'){
              ip->fwd_quality[i] = opt_fastq_ascii;
              fwd_ncount++;
          }
      }
      if(fwd_ncount > opt_fastq_maxns){
          ip->reason = maxns;
          skip = 1;
      }
  }
  if(!skip && opt_filter_nbase){
      int64_t rev_ncount = 0;
      for(int64_t i = 0; i < rev_trunc; i++){
          if(ip->rev_sequence[i] == 'N'){
              ip->rev_quality[i] = opt_fastq_ascii;
              rev_ncount++;
          }
      }
      if(rev_ncount > opt_fastq_maxns){
          ip->reason = maxns;
          skip = 1;
      }
  }
  ip->offset = 0;
  if(!skip) ip->offset = optimize();
  if(ip->offset > 0) merge();
  if(!ip->merged) merge_by_olp();
  keep_or_discard();
}

void PairMerger::parse_pair(krec1_t* fwd, krec1_t* rev){
    ip->fwd_header = fwd->name.s;
    ip->rev_header = rev->name.s;
    ip->fwd_length = fwd->seq.l-fwd->off;
    ip->rev_length = rev->seq.l-rev->off;
    sum_read_length += ip->fwd_length + ip->rev_length;
    ip->fwd_sequence = fwd->seq.s+fwd->off;
    ip->rev_sequence = rev->seq.s+rev->off;
    ip->fwd_quality = fwd->qual.s+fwd->off;
    ip->rev_quality = rev->qual.s+rev->off;
    int64_t merged_seq_needed = ip->fwd_length + ip->rev_length + 1;
    if(merged_seq_needed > ip->merged_seq_alloc){
        ip->merged_seq_alloc = merged_seq_needed;
        ip->merged_sequence = (char*)realloc(ip->merged_sequence, merged_seq_needed*sizeof(char));
        ip->merged_quality = (char*)realloc(ip->merged_quality, merged_seq_needed*sizeof(char));
    }
    if(merged_seq_needed > merge_len_alloc){
        int64_t omlcnt = merge_len_alloc;
        merge_len_alloc = merged_seq_needed;
        if(merge_len_cnt){
            merge_len_cnt = (uint64_t*)realloc(merge_len_cnt, merge_len_alloc * sizeof(uint64_t));
            memset(merge_len_cnt + omlcnt, 0, (merge_len_alloc-omlcnt) * sizeof(uint64_t));
        }else merge_len_cnt = (uint64_t*)calloc(merge_len_alloc, sizeof(uint64_t));
    }
    ip->merged_sequence[0] = 0;
    ip->merged_quality[0] = 0;
    ip->merged = 0;
    ++total;
}

void PairMerger::keep_or_discard(){
    if(ip->merged) keep();
    else discard();
}

void PairMerger::merge(krec1_t *fwd, krec1_t *rev){
    parse_pair(fwd, rev);
    process();
}

void PairMerger::merge(const char* sf, const char* sr, const char* qf, const char* qr){
    krec1_t *fwd = krec1_init();
    krec1_t *rev = krec1_init();
    ksprintf(&fwd->name, "fwd");
    ksprintf(&rev->name, "rev");
    ksprintf(&fwd->seq, "%s", sf);
    ksprintf(&rev->seq, "%s", sr);
    if(qf && qr){
        ksprintf(&fwd->qual, "%s", qf);
        ksprintf(&rev->qual, "%s", qr);
    }else{
        ksprintf(&fwd->qual, "%s", fwd->seq.s);
        ksprintf(&rev->qual, "%s", rev->seq.s);
    }
    fwd->strand = 1;
    rev->strand = 1;
    merge(fwd, rev);
    krec1_destroy(fwd);
    krec1_destroy(rev);
}

void PairMerger::merge(const std::string& sf, const std::string& sr){
    const char* fwd = sf.c_str();
    const char* rev = sr.c_str();
    merge(fwd, rev, NULL, NULL);
}

void PairMerger::merge_by_olp(){
    int offset = 0, diff = 0, ovlen = 0, maxovdiff = 0, olptype = 0, bestoff = 0;
    int maxseqerr = opt_olpm_maxseqerr * MAX(ip->fwd_length, ip->rev_length);
    double mindiffp = 1.0, diffp = 1.0;
    int64_t mindiffc = ip->fwd_length + ip->rev_length;
    int i = 0;
    while(offset <= ip->fwd_length - opt_olpm_minolplen){
        ovlen = MIN(ip->fwd_length - offset, ip->rev_length);
        maxovdiff = MIN(ovlen * opt_olpm_maxdiffpct, maxseqerr);
        i = diff = 0;
        for(i = 0; i < ovlen; ++i){
            if(ip->fwd_sequence[offset+i] != nuc_to_cmp[(int)(ip->rev_sequence[ip->rev_length-1-i])]){
                ++diff;
            }
        }
        if(diff <= maxovdiff){
            diffp = (double)diff/(double)ovlen;
            if(diffp < mindiffp && diff < mindiffc){
                mindiffp = diffp;
                mindiffc = diff;
                bestoff = offset;
            }
            olptype = 1;
        }
        ++offset;
    }
    if(!olptype){
        offset = 0; diffp = 0.0; mindiffp = 1.0; mindiffc = ip->rev_length + ip->fwd_length;
        while(offset <= ip->rev_length - opt_olpm_minolplen){
            ovlen = MIN(ip->fwd_length, ip->rev_length - offset);
            maxovdiff = MIN(ovlen * opt_olpm_maxdiffpct, maxseqerr);
            diff = 0;
            int i = 0;
            for(i = 0; i < ovlen; ++i){
                if(nuc_to_cmp[(int)ip->rev_sequence[ip->rev_length-1-offset-i]] != ip->fwd_sequence[i]){
                    ++diff;
                }
            }
            if(diff <= maxovdiff){
                diffp = (double)diff/(double)ovlen;
                if(diffp < mindiffp && diff < mindiffc){
                    mindiffp = diffp;
                    mindiffc = diff;
                    bestoff = offset;
                }
                olptype = 2; // staggered overlap
            }
            ++offset;
        }
    }
    if(olptype == 0) return;
    ip->ee_merged = 0.0;
    ip->ee_fwd = 0.0;
    ip->ee_rev = 0.0;
    ip->fwd_errors = 0;
    ip->rev_errors = 0;
    char sym, qual;
    char fwd_sym, fwd_qual, rev_sym, rev_qual;
    int64_t fwd_pos, rev_pos, merged_pos;
    double ee;
    merged_pos = 0;
    int32_t fwd_5prime_overhang = 0, rev_3prime_overhang = 0;
    if(olptype == 1){
        fwd_5prime_overhang = bestoff;
        rev_3prime_overhang = 0;
    }else if(olptype == 2){
        fwd_5prime_overhang = 0;
        rev_3prime_overhang = bestoff;
    }
    fwd_pos = 0;
    // 5' overhang in forward sequence
    fwd_pos = 0;
    while(fwd_pos < fwd_5prime_overhang){
        sym = ip->fwd_sequence[fwd_pos];
        qual = ip->fwd_quality[fwd_pos];
        ip->merged_sequence[merged_pos] = sym;
        ip->merged_quality[merged_pos] = qual;
        ee = q2p[(unsigned)qual];
        ip->ee_merged += ee;
        ip->ee_fwd += ee;
        fwd_pos++;
        merged_pos++;
    }
    // Merged region
    rev_pos = ip->rev_trunc - 1 - rev_3prime_overhang;
    while((fwd_pos < ip->fwd_trunc) && (rev_pos >= 0)){
        fwd_sym = ip->fwd_sequence[fwd_pos];
        rev_sym = nuc_to_cmp[(int)(ip->rev_sequence[rev_pos])];
        fwd_qual = ip->fwd_quality[fwd_pos];
        rev_qual = ip->rev_quality[rev_pos];
        merge_sym(&sym, &qual, fwd_sym, rev_sym, fwd_qual, rev_qual);
        if(sym != fwd_sym) ip->fwd_errors++;
        if(sym != rev_sym) ip->rev_errors++;
        ip->merged_sequence[merged_pos] = sym;
        ip->merged_quality[merged_pos] = qual;
        ip->ee_merged += q2p[(unsigned)qual];
        ip->ee_fwd += q2p[(unsigned)fwd_qual];
        ip->ee_rev += q2p[(unsigned)rev_qual];
        fwd_pos++;
        rev_pos--;
        merged_pos++;
    }
    // 5' overhang in reverse sequence
    while(rev_pos >= 0){
        sym = nuc_to_cmp[(int)(ip->rev_sequence[rev_pos])];
        qual = ip->rev_quality[rev_pos];
        ip->merged_sequence[merged_pos] = sym;
        ip->merged_quality[merged_pos] = qual;
        merged_pos++;
        ee = q2p[(unsigned)qual];
        ip->ee_merged += ee;
        ip->ee_rev += ee;
        rev_pos--;
    }
    int64_t mergelen = merged_pos;
    ip->merged_length = mergelen;
    ip->merged_sequence[mergelen] = 0;
    ip->merged_quality[mergelen] = 0;
    if(ip->ee_merged <= opt_fastq_maxee){
        ip->reason = okbyoverlap;
        ip->merged = 1;
    }
}

void PairMerger::summary(){
    if(summarized) return;
    errors_rate_fwd = errors_rate_rev = .0;
    if(merged > 0){
        errors_rate_fwd = (double)sum_errors_fwd/(double)merged;
        errors_rate_rev = (double)sum_errors_rev/(double)merged;
    }
    while(merge_len_alloc > 1 && merge_len_cnt[merge_len_alloc-1] == 0) --merge_len_alloc;
    uint64_t acc_sum = merge_len_cnt[0];
    uint64_t old_sum = acc_sum;
    uint64_t quat_arr[4] = {0};
    for(int i = 1; i <= 4; ++i) quat_arr[i-1] = merged * (0.25 * i);
    for(int64_t i = 0; i < merge_len_alloc; ++i){
        if(merge_len_cnt[i]){
            if(merge_len_min == 0) merge_len_min = i;
            if(i > merge_len_max) merge_len_max = i;
            if(merge_len_cnt[i] > merge_len_cnt[merge_len_mod]) merge_len_mod = i;
            old_sum = acc_sum;
            acc_sum += merge_len_cnt[i];
            for(int q = 0; q < 4; ++q){
                if(quat_arr[q] > old_sum && quat_arr[q] <= acc_sum){
                    merge_len_quant[q] = i;
                }
            }
        }
    }
    summarized = true;
}

void PairMerger::reportJSON(kstring_t* s, const char* dh, const char* dm){
    summary();
    ksprintf(s, "%s\"MergeSummary\": {\n", dh);
    ksprintf(s, "%s%s\"TotalPairs\": %lld,\n", dh, dm, total);
    ksprintf(s, "%s%s\"MergedPirs\": %lld,\n", dh, dm, merged);
    ksprintf(s, "%s%s\"PairsWithAdapters\": %lld,\n", dh, dm, peovhang);
    ksprintf(s, "%s%s\"UnMergedPairs\": %lld,\n", dh, dm, notmerged);
    ksprintf(s, "%s%s\"ForwardMismatchCount\": %llu,\n", dh, dm, sum_errors_fwd);
    ksprintf(s, "%s%s\"ReverseMismatchCount\": %llu,\n", dh, dm, sum_errors_rev);
    ksprintf(s, "%s%s\"ForwardMismatchRate\": %lf,\n", dh, dm, errors_rate_fwd);
    ksprintf(s, "%s%s\"ReverseMismatchRate\": %lf,\n", dh, dm, errors_rate_rev);
    ksprintf(s, "%s%s\"MergedLengthMax\": %lld,\n", dh, dm, merge_len_max);
    ksprintf(s, "%s%s\"MergedLengthMin\": %lld,\n", dh, dm, merge_len_min);
    ksprintf(s, "%s%s\"MergedLengthMod\": %lld,\n", dh, dm, merge_len_mod);
    ksprintf(s, "%s%s\"MergedLenthQuartile\": [", dh, dm);
    for(int i = 0; i < 4; ++i){
        ksprintf(s, "%lld,", merge_len_quant[i]);
    }
    s->s[s->l-1] = ']';
    ksprintf(s, ",\n");
    ksprintf(s, "%s%s\"MergedLenthDist\": [", dh, dm);
    for(int i = 0; i < merge_len_alloc; ++i){
        ksprintf(s, "%lld,", merge_len_cnt[i]);
    }
    if(s->s[s->l-1] == ',') s->s[s->l-1] = ']';
    else kputc(']', s);
    ksprintf(s, "\n%s},\n", dh);

    ksprintf(s, "%s\"MergeStatus\": {", dh);
    for(int i = 1; i < MERGE_REASON_CNT; ++i){
        ksprintf(s, "\n%s%s\"%s\": %llu,", dh, dm, reason_str[i], merge_ops_cnt[i]);
    }
    s->s[s->l-1] = '\n';
    ksprintf(s, "%s}", dh);
}

void PairMerger::reportHTML(kstring_t* s){
    summary();
    std::string subsect = "Reads merge status counts";
    std::string divName = util::replace(subsect, " ", "_");
    std::string title = "'merge rate (" + std::to_string((double)merged/(double)total*100) + "%)'";
    ksprintf(s, "<div class='subsection_title'><a title='click to hide/show' onclick=showOrHide('%s')>%s</a></div>\n", divName.c_str(), subsect.c_str());
    ksprintf(s, "<div id='%s'>\n", divName.c_str());
    ksprintf(s, "<div class='sub_section_tips'>Merge status counts and percentage of each type of reasons will be shown on mouse over.</div>\n");
    ksprintf(s, "<div class='figure' id='plot_%s'></div>\n", divName.c_str());
    ksprintf(s, "</div>\n");
    ksprintf(s, "\n<script type=\"text/javascript\">\n");
    
    std::string jsnstr;
    jsnstr.append("var mbar = {\n");
    jsnstr.append("  x: [");
    for(int i = 1; i < MERGE_REASON_CNT; ++i) jsnstr.append("'" + std::string(reason_str[i]) + "',");
    jsnstr.append("],\n");
    jsnstr.append("  y: [");
    for(int i = 1; i < MERGE_REASON_CNT; ++i) jsnstr.append(std::to_string(merge_ops_cnt[i]) + ",");
    jsnstr.append("],\n");
    jsnstr.append("  text: [");
    for(int i = 1; i < MERGE_REASON_CNT; ++i){
        jsnstr.append("'" + std::to_string(merge_ops_cnt[i]) + "(");
        jsnstr.append(std::to_string((double)merge_ops_cnt[i]/(double)total*100) + "%)',");
    }
    jsnstr.append("],\n");
    jsnstr.append("  type: 'bar',\n");
    jsnstr.append("  marker: {\n");
    jsnstr.append("    color: ['rgb(153,255,153)', 'rgb(153,255,153)'");
    for(int i = 3; i < MERGE_REASON_CNT; ++i) jsnstr.append(", 'rgb(255,51,51)'");
    jsnstr.append("],\n");
    jsnstr.append("  },\n");
    jsnstr.append("};\n");

    jsnstr.append("var layout = {\n");
    jsnstr.append("  title:" + title + ",\n");
    jsnstr.append("  xaxis: {\n");
    jsnstr.append("    tickangle: -45, tickfont: {size: 8}\n");
    jsnstr.append("  },\n");
    jsnstr.append("  yaxis: {\n");
    jsnstr.append("    title: 'sequence count',\n");
    jsnstr.append("    zroline: false,\n");
    jsnstr.append("    gridwidth: 2,\n");
    jsnstr.append("  },\n");
    jsnstr.append("  showlegend: false,\n");
    jsnstr.append("  bargap: 0.05,\n");
    jsnstr.append("};\n");

    jsnstr.append("var data = [mbar];\n");

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

void PairMerger::tsvHead(kstring_t* s){
    ksprintf(s, "%s\t%s\t%s\t%s\t%s", "MergedPairs", "UnmergedPairs", "PairsWithAapters", "ForwardMismatchRate", "ReverseMismatchRate");
}

void PairMerger::tsvBody(kstring_t* s){
    summary();
    ksprintf(s, "%lld\t%lld\t%lld\t%lf\t%lf", merged, notmerged, peovhang, errors_rate_fwd, errors_rate_rev);
}

void mergefq_usage(PairMerger* pm, char* arg0){
    fprintf(stderr, "\nUsage: %s %s [options]\n\n", PACKAGE_NAME, arg0);
    fprintf(stderr, "Options: -f read1 fastq file to overlap\n");
    fprintf(stderr, "         -r read2 fastq file to overlap\n");
    fprintf(stderr, "         -p merged fastq output\n");
    fprintf(stderr, "         -o unmerged read1 output\n");
    fprintf(stderr, "         -O unmerged read2 output\n");
    fprintf(stderr, "         -j output json stat file\n");
    fprintf(stderr, "         -l overlap min length of overlap merge [%lld]\n", pm->opt_olpm_minolplen);
    fprintf(stderr, "         -m overlap max mismatch rate of overlap merge [%f]\n", pm->opt_olpm_maxdiffpct);
    fprintf(stderr, "         -e overall max mismatch rate of overlap merge [%f]\n", pm->opt_olpm_maxseqerr);
    fprintf(stderr, "         -L overlap min length of kmer merge [%lld]\n", pm->opt_fastq_minovlen);
    fprintf(stderr, "         -M overlap max mismatch rate of kmer merge [%f]\n", pm->opt_fastq_maxdiffpct);
    fprintf(stderr, "\n");

}

int mergefq_main(int argc, char** argv){
    PairMerger* pm = new PairMerger();
    int c = 0;
    const char* f1 = NULL;
    const char* f2 = NULL;
    const char* om = NULL;
    const char* o1 = NULL;
    const char* o2 = NULL;
    const char* oj = NULL;
    while((c = getopt(argc, argv, "f:r:p:o:O:j:l:m:e:L:M")) >= 0){
        switch(c){
            case 'f': f1 = strdup(optarg); break;
            case 'r': f2 = strdup(optarg); break;
            case 'p': om = strdup(optarg); break;
            case 'o': o1 = strdup(optarg); break;
            case 'O': o2 = strdup(optarg); break;
            case 'j': oj = strdup(optarg); break;
            case 'l': pm->opt_olpm_minolplen = atoi(optarg); break;
            case 'L': pm->opt_fastq_minovlen = atoi(optarg); break;
            case 'm': pm->opt_olpm_maxdiffpct = atof(optarg); break;
            case 'M': pm->opt_fastq_maxdiffpct = atof(optarg); break;
            case 'e': pm->opt_olpm_maxseqerr = atof(optarg); break;
            default: mergefq_usage(pm, argv[0]); delete pm; return 0;
        }
    }
    if(f1 == NULL || f2 == NULL){
        mergefq_usage(pm, argv[0]); return 0;
    }
    if(om == NULL) om = "merged.fq.gz";
    if(o1 == NULL) o1 = "unmerged.r1.fq.gz";
    if(o2 == NULL) o2 = "unmerged.r2.fq.gz";
    if(oj == NULL) oj = "merge.stat.json";
    krec1_t *fwd = krec1_init();
    krec1_t *rev = krec1_init();
    gzFile fp1 = gzopen(f1, "r");
    gzFile fp2 = gzopen(f2, "r");
    gzFile fpm = gzopen(om, "wb");
    gzFile uf1 = gzopen(o1, "wb");
    gzFile uf2 = gzopen(o2, "wb");
    kseq1_t* ks1 = kseq1_init(fp1);
    kseq1_t* ks2 = kseq1_init(fp2);
    while(kseq1_read(ks1, fwd) >= 0 && kseq1_read(ks2, rev) >= 0){
        fwd->strand = 1;
        rev->strand = 0;
        pm->merge(fwd, rev);
        if(pm->ip->merged){
            gzputc(fpm, '@');
            gzputs(fpm, fwd->name.s);
            gzputc(fpm, '\n');
            gzputs(fpm, pm->ip->merged_sequence);
            gzputs(fpm, "\n+\n");
            gzputs(fpm, pm->ip->merged_quality);
            gzputc(fpm, '\n');
        }else{
            fwd->strand = rev->strand = 1;
            krec1_outgz(fwd, uf1);
            krec1_outgz(rev, uf2);
        }
    }
    gzclose(fp1); gzclose(fp2); gzclose(uf1); gzclose(uf2); gzclose(fpm);
    krec1_destroy(fwd); krec1_destroy(rev);
    kseq1_destroy(ks1); kseq1_destroy(ks2);
    kstring_t ks = {0, 0, 0};
    ksprintf(&ks, "{\n");
    pm->reportJSON(&ks, " ", " ");
    ksprintf(&ks, "}\n");
    FILE* ft = fopen(oj, "w");
    fprintf(ft, "%s", ks.s);
    fclose(ft);
    free(ks.s);
    delete pm;
    return 0;
}

void mergepe_usage(PairMerger* pm, char* arg0){
    fprintf(stderr, "\nUsage: %s %s [options]\n\n", PACKAGE_NAME, arg0);
    fprintf(stderr, "Options: -f forward nucleotide sequence to overlap\n");
    fprintf(stderr, "         -r reverse nucleotide sequence to overlap\n");
    fprintf(stderr, "         -q forward nucleotide sequence quality\n");
    fprintf(stderr, "         -Q reverse nucleotide sequence quality\n");
    fprintf(stderr, "         -l overlap min length of overlap merge [%lld]\n", pm->opt_olpm_minolplen);
    fprintf(stderr, "         -m overlap max mismatch rate of overlap merge [%f]\n", pm->opt_olpm_maxdiffpct);
    fprintf(stderr, "         -e overall max mismatch rate of overlap merge [%f]\n", pm->opt_olpm_maxseqerr);
    fprintf(stderr, "         -L overlap min length of kmer merge [%lld]\n", pm->opt_fastq_minovlen);
    fprintf(stderr, "         -M overlap max mismatch rate of kmer merge [%f]\n", pm->opt_fastq_maxdiffpct);
    fprintf(stderr, "\n");

}

int mergepe_main(int argc, char** argv){
    PairMerger* pm = new PairMerger();
    int c = 0;
    const char* s1 = NULL;
    const char* s2 = NULL;
    const char* q1 = NULL;
    const char* q2 = NULL;
    while((c = getopt(argc, argv, "f:r:l:L:m:M:e:q:Q:")) >= 0){
        switch(c){
            case 'f': s1 = strdup(optarg); break;
            case 'r': s2 = strdup(optarg); break;
            case 'q': q1 = strdup(optarg); break;
            case 'Q': q2 = strdup(optarg); break;
            case 'l': pm->opt_olpm_minolplen = atoi(optarg); break;
            case 'L': pm->opt_fastq_minovlen = atoi(optarg); break;
            case 'm': pm->opt_olpm_maxdiffpct = atof(optarg); break;
            case 'M': pm->opt_fastq_maxdiffpct = atof(optarg); break;
            case 'e': pm->opt_olpm_maxseqerr = atof(optarg); break;
            default: mergepe_usage(pm, argv[0]); delete pm; return 0;
        }
    }
    if(s1 == NULL || s2 == NULL){
        mergepe_usage(pm, argv[0]); return 0;
    }
    fprintf(stderr, "fwd(%2zu): %s\nrev(%2zu): %s\n", strlen(s1), s1, strlen(s2), s2);
    krec1_t *fwd = krec1_init();
    krec1_t *rev = krec1_init();
    ksprintf(&fwd->name, "fwd");
    ksprintf(&rev->name, "rev");
    ksprintf(&fwd->seq, "%s", s1);
    ksprintf(&rev->seq, "%s", s2);
    if(q1 && q2){
        ksprintf(&fwd->qual, "%s", q1);
        ksprintf(&rev->qual, "%s", q2);
    }else{
        ksprintf(&fwd->qual, "%s", fwd->seq.s);
        ksprintf(&rev->qual, "%s", rev->seq.s);
    }
    fwd->strand = 1;
    rev->strand = 1;
    pm->merge(fwd, rev);
    if(pm->ip->merged){
        fprintf(stderr, "merge method: %s\n", reason_str[pm->ip->reason]);
        fprintf(stderr, "merge seqlen: %2lld\n", pm->ip->merged_length);
        fprintf(stderr, "merge seqstr: %s\n",  pm->ip->merged_sequence);
        fprintf(stderr, "merge qualstr: %s\n",  pm->ip->merged_quality);
    }else{
        fprintf(stderr, "fail merge by kmer: %s\n", reason_str[pm->ip->reason]);
    }
    pm->ip->merged = false;
    pm->merge_by_olp();
    if(pm->ip->merged){
        fprintf(stderr, "merge method: %s\n", reason_str[pm->ip->reason]);
        fprintf(stderr, "merge seqlen: %2lld\n", pm->ip->merged_length);
        fprintf(stderr, "merge seqstr: %s\n",  pm->ip->merged_sequence);
        fprintf(stderr, "merge qualstr: %s\n",  pm->ip->merged_quality);
    }else{
        fprintf(stderr, "fail merge by overlap\n");
    }
    krec1_destroy(fwd);
    krec1_destroy(rev);
    delete pm;
    return 0;
}
