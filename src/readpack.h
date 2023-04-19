#ifndef READ_PACH_H
#define READ_PACH_H

#include "krec.h"
#include <atomic>

// structure to store a buncn of reads
struct ReadPack{
    krec1_t** reads = NULL; // reads pointers
    int64_t m = 0; // capacity
    int64_t n; // total reads stored
    int64_t k; // pack index;
    std::atomic_bool w; // writable if true
    std::atomic_bool r; // readable if true

    // constructor
    ReadPack(){
        n = k = 0;
        markw();
    };
    
    // destructor
    ~ReadPack(){
        if(reads){
            for(int i = 0; i < m; ++i){
                if(reads[i]){ krec1_destroy(reads[i]); reads[i] = NULL; }
            }
            free(reads);
        }
    };

    // initializer
    void init(int64_t c){
        m = c;
        n = k = 0;
        reads = (krec1_t**)malloc(m * sizeof(krec1_t*));
        for(int64_t i = 0; i < m; ++i){
            reads[i] = krec1_init();
        }
    }

    // w-only marker
    void markw(){
        r = false;
        w = true;

    }

    // r-only marker
    void markr(){
        w = false;
        r = true;
    }
};

// structure to store bunch of packs
struct PackRepo{
    ReadPack** packs = NULL; // packs
    int32_t m = 0; // capacity
    int32_t n = 0; // total packs processed

    // constructor
    PackRepo(int32_t c){
        m = c;
        packs = (ReadPack**)malloc(m * sizeof(ReadPack*));
        for(int i = 0; i < m; ++i){
            packs[i] = new ReadPack();
        }
    }

    // destructor
    ~PackRepo(){
        if(packs){
            for(int i = 0; i < m; ++i){
                if(packs[i]){
                    delete packs[i];
                    packs[i] = NULL;
                }
            }
            free(packs);
        }
    }
    
    // initializer
    void init(int32_t r){
        for(int i = 0; i < m; ++i){
            packs[i]->init(r);
        }
    }
};

#endif
