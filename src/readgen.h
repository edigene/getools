#ifndef READ_GEN_H
#define READ_GEN_H

#include "readpack.h"
#include <string>
#include <thread>
#include <atomic>
#include <mutex>

// reads generator
struct ReadGenerator{
    PackRepo** repos = NULL; // reads pack repos
    const char* reads[2] = {NULL, NULL}; // reads input file
    gzFile fp[2] = {NULL, NULL}; // read file handler
    kseq1_t *ks[2] = {NULL, NULL}; // read file parser
    uint64_t rc[2] = {0, 0}; // read counts
    int incnt = 0; // input file count
    std::atomic_bool finish[2]; // input finish file count
    std::atomic_int_fast64_t pcnt; // packs passed to other processor
    std::thread** pt = NULL; // threads for read
#ifdef DEBUG
    std::mutex logmtx; // log mtx
#endif

    // constructor
    ReadGenerator(const char* r1, const char* r2 = NULL){
        reads[0] = r1;
        reads[1] = r2;
        incnt = 1;
        if(r2) incnt = 2;
        finish[0] = false;
        finish[1] = true;
        if(incnt == 2) finish[1] = false;
        pcnt = 0;
        rc[0] = rc[1] = 0;
    }

    // destructor
    ~ReadGenerator(){
        for(int i = 0; i < incnt; ++i){
            pt[i]->join();
            delete pt[i];
            delete repos[i];
            repos[i] = NULL;
            gzclose(fp[i]);
            kseq1_destroy(ks[i]);
        }
        if(repos) free(repos);
    }

    // initializer
    void init(int32_t p, int32_t r){
        repos = (PackRepo**)malloc(incnt * sizeof(PackRepo*));
        for(int i = 0; i < incnt; ++i){
            repos[i] = new PackRepo(p);
            repos[i]->init(r);
            fp[i] = gzopen(reads[i], "r");
            ks[i] = kseq1_init(fp[i]);
        }
    }
   
    // fill repo from input file
    void fill(int r){
        bool unfinish = true;
        int64_t pidx = 0;
        while(unfinish){
            for(int i = 0; i < repos[r]->m; ++i){
                if(repos[r]->packs[i]->w){
                    repos[r]->packs[i]->n = 0;
                    int32_t c = 0;
                    while(c < repos[r]->packs[i]->m && kseq1_read(ks[r], repos[r]->packs[i]->reads[c]) >= 0){
                        ++rc[r];
                        ++c;
                    }
                    if(c > 0){
                        repos[r]->packs[i]->n = c;
                        repos[r]->packs[i]->k = pidx++;
                        repos[r]->packs[i]->w = false;
                        repos[r]->packs[i]->r = true;
                        repos[r]->n += 1;
#ifdef DEBUG
                        logmtx.lock();
                        fprintf(stderr, "thread %d generate %d reads in pack %lld\n", r, c, pidx);
                        logmtx.unlock();
#endif
                    }else{
                        repos[r]->packs[i]->w = false;
                        repos[r]->packs[i]->r = false;
                        unfinish = false;
                        break;
                    }
                }
            }
        }
        finish[r] = true;
#ifdef DEBUG
        logmtx.lock();
        fprintf(stderr, "thread %d finish reading all reads: %lld\n", r, rc[r]);
        logmtx.unlock();
#endif
    }
     
    // start generator
    void start(){
        pt = new std::thread*[incnt];
        for(int i = 0; i < incnt; ++i){
            pt[i] = new std::thread(&ReadGenerator::fill, this, i);
        }
    }

    // get read packs, must set rps[i]->w to true and rps[i]->r to false after process packs got for reuse
    bool getPack(ReadPack** rps){
        for(int i = 0; i < incnt; ++i) rps[i] = NULL;
        int found = 0;
        int rfinish = 0;
        while(true){
            for(int r = 0; r < incnt; ++r){
                for(int p = 0; p < repos[r]->m; ++p){
                    if(repos[r]->packs[p]->r && repos[r]->packs[p]->k == pcnt){
                        rps[r] = repos[r]->packs[p];
                    }
                }
            }
            found = 0;
            for(int i = 0; i < incnt; ++i) if(rps[i]) ++found;
            if(found == incnt) break; // got packs
            if(finish[0] && finish[1]){// all finish read
                rfinish = 0;
                for(int r = 0; r < incnt; ++r){
                    if(pcnt >= repos[r]->n) ++rfinish;
                }
                if(rfinish == incnt) break;
                { // buggy result of all finish, but one pack is empty, other is filled
                    for(int r = 0; r < incnt; ++r){
                        for(int p = 0; p < repos[r]->m; ++p){
                            if(repos[r]->packs[p]->r && repos[r]->packs[p]->k == pcnt){
                                rps[r] = repos[r]->packs[p];
                            }
                        }
                    }
                    found = 0;
                    for(int i = 0; i < incnt; ++i) if(rps[i]) ++found;
                    if(found == incnt) break; // got packs
                    rfinish = 0;
                    for(int r = 0; r < incnt; ++r){
                        if(pcnt >= repos[r]->n) ++rfinish;
                    }
                    if(rfinish == incnt) break;
                    else{
                        time_t tt = time(NULL);
                        tm* t = std::localtime(&tt);
                        fprintf(stderr, "[%d-%02d-%02d %02d:%02d:%02d] read1 and read2 have unequal reads, program aborted!!!\n",  t->tm_year + 1900, t->tm_mon + 1, t->tm_mday, t->tm_hour, t->tm_min, t->tm_sec);
                        exit(EXIT_FAILURE);
                    }
                }
            }else if(finish[0] != finish[1]){
                if(finish[0] && reads[1] && rc[0] < rc[1]){
                    time_t tt = time(NULL);
                    tm* t = std::localtime(&tt);
                    fprintf(stderr, "[%d-%02d-%02d %02d:%02d:%02d] read2 has more reads than read1, program aborted!!!\n",  t->tm_year + 1900, t->tm_mon + 1, t->tm_mday, t->tm_hour, t->tm_min, t->tm_sec);
                    exit(EXIT_FAILURE);
                }else if(finish[1] && reads[1] && rc[1] < rc[0]){
                    time_t tt = time(NULL);
                    tm* t = std::localtime(&tt);
                    fprintf(stderr, "[%d-%02d-%02d %02d:%02d:%02d] read1 has more reads than read2, program aborted!!!\n",  t->tm_year + 1900, t->tm_mon + 1, t->tm_mday, t->tm_hour, t->tm_min, t->tm_sec);
                    exit(EXIT_FAILURE);
                }
            }else{// some file not finish read, wait
                usleep(1);
            }
        }
        if(found == incnt){
            ++pcnt;
            return true;
        }else{
            return false;
        }
    }
};

#endif // READ_GEN_H
