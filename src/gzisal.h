#ifndef GZF_ISAL_H
#define GZF_ISAL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <sys/stat.h>
#include <igzip_lib.h>

#ifndef IGZIP_MIN_COM_LEVEL
#define IGZIP_MIN_COM_LEVEL 0
#endif

#ifndef IGZIP_MAX_COM_LEVEL
#define IGZIP_MAX_COM_LEVEL 3
#endif

#ifndef UNIX
#define UNIX 3
#endif

#ifndef IGZIP_INNER_BUF_SIZE
#define IGZIP_INNER_BUF_SIZE (1<<22)
#endif

#ifndef IGZIP_HEADER_SIZE
#define IGZIP_HEADER_SIZE (1<<16)
#endif

#ifndef IGZIP_COMP_LEVEL_DEFAULT
#define IGZIP_COMP_LEVEL_DEFAULT 2
#endif

const int IGZIP_COMP_LEVEL_SIZE_BUF_LEN[4] = {
    ISAL_DEF_LVL0_DEFAULT, ISAL_DEF_LVL1_DEFAULT, 
    ISAL_DEF_LVL2_DEFAULT, ISAL_DEF_LVL3_DEFAULT,
};

typedef struct{
    FILE* _file;
    char* _mode;
    int _is_plain;
    struct isal_gzip_header* _gzip_header;
    struct inflate_state* _gzip_state;
    struct isal_zstream* _gzip_stream;
    uint8_t* _gzip_input_buffer;
    size_t _gzip_input_buffer_size;
    uint8_t* _gzip_output_buffer;
    size_t _gzip_output_buffer_size;
}__igzFile;

typedef __igzFile* igzFile;

extern inline uint32_t get_poxit_filetime(FILE* fp);
extern inline int is_gz(FILE* fp);
extern inline igzFile igzopen(const char* file, const char* mode);
extern inline void igzclose(igzFile igf);
extern inline int igzread(igzFile file, void* buf, unsigned int len);
extern inline int igzwrite(igzFile file, void* buf, unsigned int len);
extern inline int igzputsn(igzFile file, const char* buf, unsigned int len);
extern inline int isal_inflate_set_compress_level(igzFile file, int level);
extern inline int igzeof(igzFile file);

uint32_t get_poxit_filetime(FILE* fp){
    struct stat file_stats;
    fstat(fileno(fp), &file_stats);
    return file_stats.st_mtime;
}

int is_gz(FILE* fp){
    if(!fp) return 0;
    char buf[2];
    int gzip = 0;
    if(fread(buf, 1, 2, fp) == 2){
        if(((int)buf[0] == 0x1f) && ((int)(buf[1]&0xFF) == 0x8b)) gzip = 1;
    }
    fseek(fp, 12, SEEK_SET);
    if(fread(buf, 1, 2, fp) == 2){
        if((int)buf[0] == 0x42 && (int)(buf[1]&0xFF) == 0x43) gzip = 2;
    }
    fseek(fp, 0, SEEK_SET);
    return gzip;
}

igzFile igzopen(const char* file, const char* mode){
    igzFile igf = (igzFile)calloc(1, sizeof(__igzFile));
    igf->_file = fopen(file, mode);
    if(!igf->_file){
        igzclose(igf);
        return NULL;
    }
    igf->_mode = strdup(mode);
    // plain file
    if(mode[0] == 'r'){
        igf->_is_plain = !is_gz(igf->_file);
        if(igf->_is_plain) return igf;
    }
    // gz file
    igf->_gzip_header = (struct isal_gzip_header*)calloc(1, sizeof(struct isal_gzip_header));
    isal_gzip_header_init(igf->_gzip_header);
    // read
    if(mode[0] == 'r'){
        igf->_gzip_state = (struct inflate_state*)calloc(1, sizeof(struct inflate_state));
        igf->_gzip_input_buffer_size = IGZIP_INNER_BUF_SIZE;
        igf->_gzip_input_buffer = (uint8_t*)malloc(igf->_gzip_input_buffer_size * sizeof(uint8_t));
        isal_inflate_init(igf->_gzip_state);
        igf->_gzip_state->crc_flag = ISAL_GZIP_NO_HDR_VER;
        igf->_gzip_state->next_in = igf->_gzip_input_buffer;
        igf->_gzip_state->avail_in = fread(igf->_gzip_state->next_in, 1, igf->_gzip_input_buffer_size, igf->_file);
        int ret = isal_read_gzip_header(igf->_gzip_state, igf->_gzip_header);
        if(ret != ISAL_DECOMP_OK){
            igzclose(igf);
            return NULL;
        }
    }else if(mode[0] == 'w'){
        // write
        igf->_gzip_header->os = UNIX;
        igf->_gzip_header->time = get_poxit_filetime(igf->_file);
        igf->_gzip_header->name = strdup(file) ; 
        igf->_gzip_header->name_buf_len = strlen(igf->_gzip_header->name)+1;
        igf->_gzip_output_buffer_size = IGZIP_INNER_BUF_SIZE;
        igf->_gzip_output_buffer = (uint8_t*)calloc(igf->_gzip_output_buffer_size, sizeof(uint8_t));
        igf->_gzip_stream =(struct isal_zstream*)calloc(1, sizeof(struct isal_zstream));
        isal_deflate_init(igf->_gzip_stream);
        igf->_gzip_stream->avail_in = 0;
        igf->_gzip_stream->flush = NO_FLUSH;
        igf->_gzip_stream->level = IGZIP_COMP_LEVEL_DEFAULT;
        igf->_gzip_stream->level_buf_size = IGZIP_COMP_LEVEL_SIZE_BUF_LEN[igf->_gzip_stream->level];
        igf->_gzip_stream->level_buf = (uint8_t*)calloc(igf->_gzip_stream->level_buf_size, sizeof(uint8_t));
        igf->_gzip_stream->gzip_flag = IGZIP_GZIP_NO_HDR;
        igf->_gzip_stream->avail_out = igf->_gzip_output_buffer_size;
        igf->_gzip_stream->next_out = igf->_gzip_output_buffer;
        int ret = isal_write_gzip_header(igf->_gzip_stream, igf->_gzip_header);
        if(ret != ISAL_DECOMP_OK){
            igzclose(igf);
            return NULL;
        }
    }
    return igf;
}

int isal_inflate_set_compress_level(igzFile file, int level){
    if(!file) return -1;
    if(!file->_mode) return -1;
    if(file->_mode[0] != 'w') return -1;
    if(level < IGZIP_MIN_COM_LEVEL || level > IGZIP_MAX_COM_LEVEL) return -1;
    if(file->_gzip_stream->level != (uint32_t)level){
        file->_gzip_stream->level = level;
        file->_gzip_stream->level_buf_size = IGZIP_COMP_LEVEL_SIZE_BUF_LEN[file->_gzip_stream->level];
        file->_gzip_stream->level_buf = (uint8_t*)realloc(file->_gzip_stream->level_buf, file->_gzip_stream->level_buf_size*sizeof(uint8_t));
    }
    return 0;
}

void igzclose(igzFile igf){
    if(!igf) return;
    if(igf->_mode) free(igf->_mode);
    if(igf->_gzip_stream && igf->_file) igzwrite(igf, NULL, 0);
    if(igf->_gzip_header){
        if(igf->_gzip_header->name) free(igf->_gzip_header->name);
        free(igf->_gzip_header);
    }
    if(igf->_gzip_state) free(igf->_gzip_state);
    if(igf->_gzip_input_buffer) free(igf->_gzip_input_buffer);
    if(igf->_gzip_output_buffer) free(igf->_gzip_output_buffer);
    if(igf->_gzip_stream){
        if(igf->_gzip_stream->level_buf) free(igf->_gzip_stream->level_buf);
        free(igf->_gzip_stream);
    }
    if(igf->_file) fclose(igf->_file);
    free(igf);
}

int igzread(igzFile file, void* buf, unsigned int len){
    int declen = 0;
    if(file->_is_plain){
        if(!feof(file->_file)) declen = fread((uint8_t*)buf, 1, len, file->_file);
        return declen;
    }
    // Start reading in compressed data and decompress
    do{
        if(!file->_gzip_state->avail_in){
            file->_gzip_state->next_in = file->_gzip_input_buffer;
            file->_gzip_state->avail_in = fread(file->_gzip_state->next_in, 1, file->_gzip_input_buffer_size, file->_file);
        }
        file->_gzip_state->next_out = (uint8_t *)buf;
        file->_gzip_state->avail_out = len;
        if(isal_inflate(file->_gzip_state) != ISAL_DECOMP_OK) return -3;
        declen = file->_gzip_state->next_out - (uint8_t*)buf;
        if(declen) return declen;
    }while(file->_gzip_state->block_state != ISAL_BLOCK_FINISH && (!feof(file->_file) || !file->_gzip_state->avail_out));
    // Add the following to look for and decode additional concatenated files
    if(!feof(file->_file) && !file->_gzip_state->avail_in){
        file->_gzip_state->next_in = file->_gzip_input_buffer;
        file->_gzip_state->avail_in = fread(file->_gzip_state->next_in, 1, file->_gzip_input_buffer_size, file->_file);
    }
    while(file->_gzip_state->avail_in > 0 && file->_gzip_state->next_in[0] == 31){
        // Look for magic numbers for gzip header. Follows the gzread() decision whether to treat as trailing junk
        if(file->_gzip_state->avail_in > 1 && file->_gzip_state->next_in[1] != 139) break;
        isal_inflate_reset(file->_gzip_state);
        file->_gzip_state->crc_flag = ISAL_GZIP;
        do{
            if(!file->_gzip_state->avail_in && !feof(file->_file)){
                file->_gzip_state->next_in = file->_gzip_input_buffer;
                file->_gzip_state->avail_in = fread(file->_gzip_state->next_in, 1, file->_gzip_input_buffer_size, file->_file);
            }
            file->_gzip_state->next_out = (uint8_t *)buf;
            file->_gzip_state->avail_out = len;
            if(isal_inflate(file->_gzip_state) != ISAL_DECOMP_OK) return -3;
            declen = file->_gzip_state->next_out - (uint8_t*)buf;
            if(declen) return declen;
        }while(file->_gzip_state->block_state != ISAL_BLOCK_FINISH && (!feof(file->_file) || !file->_gzip_state->avail_out));
        if(!feof(file->_file) && !file->_gzip_state->avail_in){
            file->_gzip_state->next_in = file->_gzip_input_buffer;
            file->_gzip_state->avail_in = fread(file->_gzip_state->next_in, 1, file->_gzip_input_buffer_size, file->_file);
        }
    }
    return declen;
}

int igzwrite(igzFile file, void* buf, unsigned int len){
    file->_gzip_stream->next_in = (uint8_t*)buf;
    file->_gzip_stream->avail_in = len;
    file->_gzip_stream->end_of_stream = !buf;
    int cclen = 0;
    do{
        if(!file->_gzip_stream->next_out){
            file->_gzip_stream->next_out = file->_gzip_output_buffer;
            file->_gzip_stream->avail_out = file->_gzip_output_buffer_size;
        }
        int ret = isal_deflate(file->_gzip_stream);
        if(ret != ISAL_DECOMP_OK) return -3;
        cclen += fwrite(file->_gzip_output_buffer, 1, file->_gzip_stream->next_out-file->_gzip_output_buffer, file->_file);
        file->_gzip_stream->next_out = NULL;
    }while(!file->_gzip_stream->avail_out);
    return cclen;
}


int igzputsn(igzFile file, const char* buf, unsigned int len){
    file->_gzip_stream->next_in = (uint8_t*)buf;
    file->_gzip_stream->avail_in = len;
    file->_gzip_stream->end_of_stream = !buf;
    int cclen = 0;
    do{
        if(!file->_gzip_stream->next_out){
            file->_gzip_stream->next_out = file->_gzip_output_buffer;
            file->_gzip_stream->avail_out = file->_gzip_output_buffer_size;
        }
        int ret = isal_deflate(file->_gzip_stream);
        if(ret != ISAL_DECOMP_OK) return -3;
        cclen += fwrite(file->_gzip_output_buffer, 1, file->_gzip_stream->next_out-file->_gzip_output_buffer, file->_file);
        file->_gzip_stream->next_out = NULL;
    }while(!file->_gzip_stream->avail_out);
    return cclen;
}

int igzeof(igzFile file){
    if(!file) return 0;
    if(file->_mode[0] != 'w' || file->_mode[0] != 'r') return 0;
    return file->_mode[0] == 'r' ? feof(file->_file) : 0;
}

#endif
