#ifndef UTIL_H
#define UTIL_H

#ifdef WINDOWS
#include <direct.h>
#else
#include <unistd.h>
#endif

#include <set>
#include <string>
#include <cerrno>
#include <cstdio>
#include <cctype>
#include <vector>
#include <mutex>
#include <numeric>
#include <fstream>
#include <sstream>
#include <utility>
#include <cstring>
#include <iostream>
#include <algorithm>
#include <functional>
#include <ftw.h>
#include <zlib.h>
#include <math.h>
#include <stdarg.h>
#include <dirent.h>
#include <sys/stat.h>

// useful macros beg
#ifndef exp10
#define exp10(x) (pow(10.0,(x)))
#endif

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif

#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif
// useful macors end

/** +->- and -->+ */
const unsigned char rev_strand[256] = {
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','-','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',

    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','+',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',

    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',

    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N'
};

/** nucleotide to uppercase */
const unsigned char nuc_to_upper[256] = {
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',

    'N','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O',
    'P','Q','R','S','T','U','V','W','X','Y','Z','N','N','N','N','N',
    'N','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O',
    'P','Q','R','S','T','U','V','W','X','Y','Z','N','N','N','N','N',

    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',

    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N'
};

/** nucleotide to lowercase */
const unsigned char nuc_to_lower[256] = {
    'n','n','n','n','n','n','n','n','n','n','n','n','n','n','n','n',
    'n','n','n','n','n','n','n','n','n','n','n','n','n','n','n','n',
    'n','n','n','n','n','n','n','n','n','n','n','n','n','n','n','n',
    'n','n','n','n','n','n','n','n','n','n','n','n','n','n','n','n',

    'n','a','b','c','d','e','f','g','h','i','j','k','l','m','n','o',
    'p','q','r','s','t','u','v','w','x','y','z','n','n','n','n','n',
    'n','a','b','c','d','e','f','g','h','i','j','k','l','m','n','o',
    'p','q','r','s','t','u','v','w','x','y','z','n','n','n','n','n',

    'n','n','n','n','n','n','n','n','n','n','n','n','n','n','n','n',
    'n','n','n','n','n','n','n','n','n','n','n','n','n','n','n','n',
    'n','n','n','n','n','n','n','n','n','n','n','n','n','n','n','n',
    'n','n','n','n','n','n','n','n','n','n','n','n','n','n','n','n',

    'n','n','n','n','n','n','n','n','n','n','n','n','n','n','n','n',
    'n','n','n','n','n','n','n','n','n','n','n','n','n','n','n','n',
    'n','n','n','n','n','n','n','n','n','n','n','n','n','n','n','n',
    'n','n','n','n','n','n','n','n','n','n','n','n','n','n','n','n'
};
/** nucleotide to uppercase and U to T */
const unsigned char nuc_to_norm[256] = {
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',

    'N','A','B','C','D','N','N','G','H','N','N','K','N','M','N','N',
    'N','N','R','S','T','T','V','W','N','Y','N','N','N','N','N','N',
    'N','A','B','C','D','N','N','G','H','N','N','K','N','M','N','N',
    'N','N','R','S','T','T','V','W','N','Y','N','N','N','N','N','N',

    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',

    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N',
    'N','N','N','N','N','N','N','N','N','N','N','N','N','N','N','N'
};

/** nucleotide complement table */
const unsigned char nuc_to_cmp[256] = {
       0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15,
      16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  31,
      32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,
      48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,
      64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
     'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',  91,  92,  93,  94,  95,
      64, 'T', 'v', 'G', 'h', 'e', 'f', 'C', 'd', 'i', 'j', 'm', 'l', 'k', 'N', 'o',
     'p', 'q', 'y', 's', 'A', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127
 };

/** nucleotide to 2bit table, Aa:0, Cc:1, Gg:2, TtUu:3, All others: 0  */
const uint8_t nuc_to_2bit[256] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

/** nucleotide to 3bit table, Aa:0, Cc:1, Gg:2, Tt:3, Nn:4, All others: 4 */
static const uint8_t nuc_to_3bit[256] = {
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
};

/** 3bit int to ACGTN */
const char bit3_to_nuc[5] = {'A', 'C', 'G', 'T', 'N'};

/** ambiguous nucleotide mask, mask everything but AaCcGgTtUu to 1*/
const unsigned int ambig_nuc_mask[256] = {
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
};

/** utility to operate on strings and directories */
namespace util{
    // beg LineReader
    /** Class to hold a file line reader in the format .gz or plain text */
    struct LineReader{
        std::string mFileName; // name of file
        gzFile mGzipFile; // gzFile to store opened mZipped file handler
        FILE* mFile; // FILE pointer to store opened plain file handler
        bool mZipped; // the file is gzipped if true
        char* mBuf; // mBuffer to store a chunk of characters read from mGzipFile or mFile
        int mBufDataLen; // the length of characters read into the mBuffer after the last read
        int mBufUsedLen; // the length of characters already consumed in the mBuffer
        bool mStdinMode; // read from stdin if true
        bool mNoLineBreakAtEnd; // the file has no '\n' as a line break at the last line if true
        int32_t mReadBufSize; // the mBuffer size used to read

        /** Construct a file line reader with filename
         * @param filename Name of the file
         */
        LineReader(const std::string& filename){
            mFileName = filename;
            mGzipFile = NULL;
            mFile = NULL;
            mStdinMode = false;
            mReadBufSize = (1 << 20);
            mBuf = new char[mReadBufSize];
            mBufDataLen = 0;
            mBufUsedLen = 0;
            mNoLineBreakAtEnd = false;
            init();
        }

        /** LineReader Destructor */
        ~LineReader(){
            close();
            delete mBuf;
            mBuf = nullptr;
        }

        /** Tell whether the file is zipped or not
         * @return true if the fasq file is zipped
         */
        inline bool isZipped(){
            return mZipped;
        }

        /** Get the number of Bytes read from(write to) the file in the meantime, store it in bytesRead
         *  Get the total number of Bytes in the file, store it in bytesTotal
         */
        inline void getBytes(size_t& bytesRead, size_t& bytesTotal){
            if(mZipped){
                bytesRead = gzoffset(mGzipFile);
            }else{
                bytesRead = std::ftell(mFile);
            }
            // use another ifstream without affecting the current reader
            std::ifstream is(mFileName);
            is.seekg(0, is.end);
            bytesTotal = is.tellg();
        }

        /** tell whether the file is a zipped file
         * @param filename file name
         * @return true if the file is a zipped file(with suffix ".gz")
         */
        inline static bool isZippedFile(const std::string& filename){
            if(filename.length() < 4){
                return false;
            }
            if(filename.substr(filename.length() - 3) == ".gz"){
                return true;
            }
            return false;
        }

        /** read one line into line from buffer
         * @param line strin to store result
         * @return true if read successful
         */
        inline bool getline(std::string& line){
            if(mBufUsedLen >= mBufDataLen && eof()){
                return false;
            }
            if(mZipped && mGzipFile == NULL){
                return false;
            }
            line = getlineFromBuffer();
            return true;
        }

        private:
        /** get just one line from the mBuf, update the mBuf if needed
         */
        inline std::string getlineFromBuffer(){
            int start = mBufUsedLen;
            int end = start;
            // look for '\r' or '\n' until the end of mBuf
            while(end < mBufDataLen){
                if(mBuf[end] != '\r' && mBuf[end] != '\n'){
                    ++end;
                }else{
                    break;
                }
            }
            // if '\r' or '\n' found(this line well contained in this mBuf)
            // or this is the last mBuf of file
            if(end < mBufDataLen || mBufDataLen < mReadBufSize){
                int len = end - start;
                std::string line(mBuf + start, len);
                ++end;
                if(end < mBufDataLen - 1 && mBuf[end - 1] == '\r' && mBuf[end] == '\n'){
                    ++end;
                }
                mBufUsedLen = end;
                return line;
            }
            // if '\r' or '\n' not found && this is not the last mBuf of file
            // then this line is not contained in this mBuf, we should read new mBuf
            std::string str(mBuf + start, mBufDataLen - start);
            while(true){
                readToBuf();
                start = 0;
                end = 0;
                // look for '\r' or '\n' until the end of mBuf
                while(end < mBufDataLen){
                    if(mBuf[end] != '\r' && mBuf[end] != '\n'){
                        ++end;
                    }else{
                        break;
                    }
                }
                // if '\r' or '\n' found(this line well contained in this mBuf)
                // or this is the last mBuf of file
                if(end < mBufDataLen || mBufDataLen < mReadBufSize){
                    int len = end - start;
                    str.append(mBuf + start, len);
                    ++end;
                    if(end < mBufDataLen - 1 && mBuf[end - 1] == '\r' && mBuf[end] == '\n'){
                        ++end;
                    }
                    mBufUsedLen = end;
                    return str;
                }
                // if '\r' or '\n' not found && this is not the last mBuf of file
                // then this line is not contained in this mBuf, we should read new mBuf
                str.append(mBuf + start, mBufDataLen);
            }
            return std::string();
        }

        /** Tell whether the file has no '\n' as a line break at the last line
         * @return true if the file has no '\n' as a line break at the last line
         */
        inline bool hasNoLineBreakAtEnd(){
            return mNoLineBreakAtEnd;
        }

        /** Tell whether the LineReader has reach the endof file
         * @return true if eof reached
         */
        inline bool eof(){
            if(mZipped){
                return gzeof(mGzipFile);
            }else{
                return std::feof(mFile);
            }
        }

        /** initialize the LineReader:
         * 1, open file and store file handler into mGzipFile or mFile or read from stdin
         * 2, set the starting position for the next read on compressed file stream file to the beginning of file
         * 3, update the file format mZipped
         * 4, call readToBuf() to try to fill the mBuf from first reading
         */
        inline void init(){
            if(isZippedFile(mFileName)){
                mGzipFile = gzopen(mFileName.c_str(), "r");
                mZipped = true;
                gzrewind(mGzipFile);
            }else{
                if(mFileName == "/dev/stdin"){
                    mFile = stdin;
                }else{
                    mFile = std::fopen(mFileName.c_str(), "rb");
                }
                if(mFile == NULL){
                    std::cerr << "Failed to open file: " <<  mFileName << std::endl;
                    std::exit(1);
                }
                mZipped = false;
            }
            readToBuf();
        }

        /** close the LineReader:
         * 1, close file handler
         * 2, set file handler to NULL
         */
        inline void close(){
            if(mZipped && mGzipFile){
                gzclose(mGzipFile);
                mGzipFile = NULL;
            }else if(mFile){
                std::fclose(mFile);
                mFile = NULL;
            }else{
                return;
            }
        }

        /** trim \n, \r or \r\n in the tail of the line */
        inline void clearLineBreaks(char* line){
            int len = strlen(line);
            if(len >= 2){
                if(line[len - 1] == '\n' || line[len - 1] == '\r'){
                    line[len - 1] = '\0';
                    if(line[len - 2] == '\r'){
                        line[len - 2] = '\0';
                    }
                }
            }
        }

        /** try a reading action to fill the mBuf,
         * update mBufDataLen to Bytes read during this reading action
         * reset mBufUsedLen to zero
         * if read the last line(mBuf is not filled), update mNoLineBreakAtEnd
         */
        inline void readToBuf(){
             if(isZipped()){
                 mBufDataLen = gzread(mGzipFile, mBuf, mReadBufSize);
                 if(mBufDataLen == -1){
                     std::cerr << "Error to read gzip file" << std::endl;
                     std::exit(1);
                 }
             }else{
                 mBufDataLen = std::fread(mBuf, 1, mReadBufSize, mFile);
             }
             mBufUsedLen = 0;
             if(mBufDataLen < mReadBufSize){
                 if(mBuf[mBufDataLen - 1] != '\n'){
                     mNoLineBreakAtEnd = true;
                 }
             }
         }
    };
    // end LineReader
    
    /** whether a string starts with some substring
     * @param str whole string
     * @param pre substring
     * @return true if str starts with pre
     */
    inline bool startsWith(const std::string& str, const std::string& pre){
        if(str.length() < pre.length()){
            return false;
        }else{
            return std::equal(pre.begin(), pre.end(), str.begin());
        }
    }

    /** whether a string ends with some substring
     * @param str whole string
     * @param suf substring
     * @return true if str ends with suf
     */
    inline bool endsWith(const std::string& str, const std::string& suf){
        if(str.length() < suf.length()){
            return false;
        }else{
            return std::equal(suf.rbegin(), suf.rend(), str.rbegin());
        }
    }

    /** get rid of the leading and ending white space characters of a string
     * @param str string to be stripped in both ends
     * @param pat string to be stripped from back
     * @return a string with white spaces stripped in both ends
     */
    inline std::string strip(const std::string& str, const std::string& pat = " \t\n\v\f\r"){
        std::string::size_type ipos = str.find_first_not_of(pat);
        if(ipos == std::string::npos){
            return "";
        }
        std::string::size_type epos = str.find_last_not_of(pat);
        if(epos == ipos){
            return str.substr(ipos);
        }else{
            return str.substr(ipos, epos - ipos + 1);
        }
    }

    /** get rid of the left leading white space characters of a string
     * @param str string to be stripped from front
     * @param pat string to be stripped from back
     * @return a string with left leading white spaces stripped
     */
    inline std::string lstrip(const std::string& str, const std::string& pat = " \t\n\v\f\r"){
        std::string::size_type pos = str.find_first_not_of(pat);
        if(pos == std::string::npos){
            return "";
        }
        return str.substr(pos);
    }

    /** get rid of the trailling white space characters of a string
     * @param str string to be stripped from back
     * @param pat string to be stripped from back
     * @return a string with right ending white spaces stripped
     */
    inline std::string rstrip(const std::string& str, const std::string& pat = " \t\n\v\f\r"){
        std::string::size_type pos = str.find_last_not_of(pat);
        if(pos == std::string::npos){
           return "";
        }
        return str.substr(0, pos + 1);
    } 

    /** remove prefix from a string
     * @param str string to be remove prefix
     * @param rms prefix to be removed
     */
    inline std::string lrmstr(const std::string& str, const std::string& rms){
        if(startsWith(str, rms)){
            return str.substr(rms.size());
        }
        return str;
    }

    /** remove suffix from a string
     * @param str string to be remove suffix
     * @param rms suffix to be removed
     */
    inline std::string rrmstr(const std::string& str, const std::string& rms){
        if(endsWith(str, rms)){
            return str.substr(0, str.size() - rms.size());
        }
        return str;
    }

    /** split a string by predefined seperator into a vector
     * @param str string
     * @param vec vector to store the split results
     * @param sep seperators, can contain a series of seperators
     */
    inline void split(const std::string& str, std::vector<std::string>& vec, std::string sep, bool append=false){
        if(!append) vec.clear(); // clear output vector
        if(str.empty()) return; // empty string, return directly
        std::string::size_type las, cur;
        las = cur = str.find_first_of(sep);
        if(las == std::string::npos){ // only one field without any sep
            vec.push_back(str);
            return;
        }
        if(las == 0) vec.push_back(""); // begin with sep, so there is a blank str before
        if(las != 0) vec.push_back(str.substr(0, las)); // begin with a non-sep field
        while(las != std::string::npos){
            if(las + 1 == std::string::npos){
                vec.push_back(""); // ends with sep, so there is a blank str after
                return;
            }
            cur = str.find_first_of(sep, las + 1);
            if(cur == std::string::npos){
                vec.push_back(str.substr(las + 1)); // last field
            }else{
                if(cur - las == 1){
                    vec.push_back(""); // consecutive sep
                }else{
                    vec.push_back(str.substr(las + 1, cur - las - 1));
                }
            }
            las = cur;
        }
    }

    /** join a list of strings by an seperator
     * @param vec vector to store the split results
     * @param ret joined string
     * @param sep seperator used to join strings
     */
    inline void join(const std::vector<std::string>& vec, std::string& ret, const std::string& sep = " "){
        if(vec.empty()) return;
        if(vec.size() == 1){
            ret = vec[0];
            return;
        }
        for(uint32_t i = 0; i < vec.size() - 1; ++i){
            ret.append(vec[i] + sep);
        }
        ret.append(vec[vec.size() - 1]);
    }
    
    /** join a list of strings by an seperator
     * @param vec vector to store the split results
     * @param sep seperator used to join strings
     * @return joined string
     */
    template <typename T>
    inline std::string join(const std::vector<T>& vec, const std::string& sep = " "){
        if(vec.empty()) return "";
        std::stringstream ss;
        ss << vec[0];
        if(vec.size() == 1){
            return ss.str();
        }
        for(uint32_t i = 1; i < vec.size(); ++i){
            ss << sep << vec[i];
        }
        return ss.str();
    }

    /** join a list of strings by an seperator
     * @param vec array to store the split results
     * @param len length of array
     * @param sep seperator used to join strings
     * @return joined string
     */
    template <typename T>
    inline std::string join(const T* vec, int len,  const std::string& sep = " "){
        if(len == 0) return "";
        std::stringstream ss;
        ss << vec[0];
        if(len == 1){
            return ss.str();
        }
        for(int i = 1; i < len; ++i){
            ss << sep << vec[i];
        }
        return ss.str();
    }

    /** replace a substr apearing in a string with another string
     * @param str string 
     * @param pat substr of string to be replaced
     * @param des string to be used to replaced with pat
     * @return a string with each pat replaced by des
     */
    inline std::string replace(const std::string& str, const std::string& pat, const std::string& des){
        std::string ret;
        std::string::size_type las = 0, cur = 0;
        while((cur = str.find(pat, cur)) != std::string::npos){
            ret.append(str.substr(las, cur - las));
            ret.append(des);
            cur += pat.length();
            las = cur;
        }
        if(las != std::string::npos){
            ret.append(str.substr(las));
        }
        return ret;
    }

    /** replace a char apearing in a string with another char
     * @param str string 
     * @param pat char of string to be replace
     * @param des char to be used to replaced with pat
     */
    inline void replace(std::string& str, const char& pat, const char& des){
        for(auto& e: str) if(e == pat) e = des;
    }

    /** get a reverse sequence of str
     * @param str input string sequence
     * @return reverse sequence of str
     */
    inline std::string reverse(const std::string& str){
        std::string ret = str;
        std::reverse(ret.begin(), ret.end());
        return ret;
    }

    /** get the basename of a path string
     * @param path name
     * @param ext to be removed
     * @return basename of path
     */
    inline std::string basename(const std::string& path, const std::string ext = ""){
        std::string fpath = util::replace(path, "~", std::string(std::getenv("HOME")) + "/");
        if(fpath.find_first_of("\t\n\v\f\r") != std::string::npos){
            return "";
        }
        std::string::size_type pos1 = fpath.find_last_of("/\\");
        if(pos1 == std::string::npos){
            return rrmstr(fpath, ext);
        }
        std::string::size_type pos2 = fpath.find_last_not_of("/\\");
        if(pos2 == fpath.size() - 1){
            return rrmstr(fpath.substr(pos1 + 1, pos2 - pos1), ext);
        }
        std::string::size_type pos3 = fpath.find_last_of("/\\", pos2);
        if(pos3 == std::string::npos){
            return rrmstr(fpath.substr(0, pos2 + 1), ext);
        }
        return rrmstr(fpath.substr(pos3 + 1, pos2 - pos3), ext);
    }

    /** get the dirname of a path string
     * @param path name
     * @return dirname of path
     */
    inline std::string dirname(const std::string& path){
        std::string fpath = util::replace(path, "~", std::string(std::getenv("HOME")) + "/");
        std::string::size_type pos = fpath.find_last_of("/\\");
        if(pos == std::string::npos){
#ifdef _WIN32
            return ".\\";
#else
            return "./";
#endif
        }
        if(pos == fpath.size() - 1){
            std::string::size_type pos1 = fpath.find_last_not_of("/\\");
            if(pos1 == std::string::npos){
#ifdef _WIN32
                return "\\";
#else
                return "/";
#endif
            }
            std::string::size_type pos2 = fpath.find_last_of("/\\", pos1);
            if(pos2 == std::string::npos){
#ifdef _WIN32
                return ".\\";
#else
                return "./";
#endif
            }else{
                return fpath.substr(0, pos2 + 1);
            }
        }else{
            return fpath.substr(0, pos + 1);
        }
    }

    /** get absolute path of a path
     * @param path string of path
     * @return absolute path string
     */
    inline std::string abspath(const std::string& path){
        std::string newPath = util::strip(path);
        if(newPath.length() == 0){
            return "";
        }
        char cpath[FILENAME_MAX];
#ifdef _WIN32
        _realpath(newPath.c_str(), cpath);
#else
        realpath(newPath.c_str(), cpath);
#endif
        return cpath;
    }

    /** get current working directory
     * @return current working directory
     */
    inline std::string cwd(void){
        char cpath[FILENAME_MAX];
#ifdef _WIN32
        _getcwd(cpath, FILENAME_MAX);
#else
        getcwd(cpath, FILENAME_MAX);
#endif
        return cpath;
    }

    /** join dirname and basename into a path
     * @param dirname dirname of the path
     * @param basename basename of the path
     * @return full path
     */
    inline std::string joinpath(const std::string& dirname, const std::string& basename){
#ifdef _WIN32
        return dirname + "\\" + basename;
#else
        return dirname + "/" + basename;
#endif
    }

    /** check a string is a regular file or not
     * @param path string of a file/directory
     * @return true if path is an existing regular file
     */
    inline bool isfile(const std::string& path){
#ifdef _WIN32
        struct _stat info;
        if(_stat(path.c_str(), &info) != 0){
            return false;
        }
        return (info.st_mode & _S_IFREG);
#else
        struct stat info;
        if(stat(path.c_str(), &info) !=0){
            return false;
        }
        return (info.st_mode & S_IFREG);
#endif
    }

    /** check a string is a directory or not
     * @param path string of a file/directory
     * @return true if path is an existing path
     */
    inline bool isdir(const std::string& path){
#ifdef _WIN32
        struct _stat info;
        if(_stat(path.c_str(), &info) != 0){
            return false;
        }
        return (info.st_mode & _S_IFDIR);
#else
        struct stat info;
        if(stat(path.c_str(), &info) != 0){
            return false;
        }
        return (info.st_mode & S_IFDIR);
#endif
    }
    
    /** check a file/directory exists or not
     * @param path string of a file/directory
     * @return true if path exists
     */
    inline bool exists(const std::string& path){
#ifdef _WIN32
        struct _stat s;
        return (_stat(path.c_str(), &s) == 0);
#else
        struct stat s;
        return (stat(path.c_str(), &s) == 0);
#endif
    }

    /** exit and print string to std::cerr
     * @param msg string to print to std::cerr
     * @param f log file fp
     */
    inline void errorExit(const std::string& msg, FILE* f = stderr){
        time_t tt = time(NULL);
        tm* t = std::localtime(&tt);
        char date[60] = {0};
        std::snprintf(date, 60, "[%d-%02d-%02d %02d:%02d:%02d] ",
                t->tm_year + 1900, t->tm_mon + 1, t->tm_mday,
                t->tm_hour, t->tm_min, t->tm_sec);
        fprintf(f, "%s%s\n", date, msg.c_str());
        fflush(f);
        exit(-1);
    }
    
    /** check a file existence status, if not exists, exit 
     * @param path string of file/directory
     */
    inline void validFile(const std::string& path){
        if(util::isdir(path)){
            util::errorExit("this is not a file path!");
        }
        if(!util::isfile(path)){
            util::errorExit("file does not exist");
        }
    }

    /** make directories recursively
     * @param path path of directory to be created
     * @return true if make directories successfully
     */
    inline bool makedir(const std::string& path){
        std::string fpath = util::replace(path, "~", std::string(std::getenv("HOME")) + "/");
#ifdef _WIN32
        int ret = _mkdir(fpath.c_str());
#else
        mode_t mode = 0755;
        int ret = mkdir(fpath.c_str(), mode);
#endif
        if(ret == 0){
            return true;
        }
        switch(errno){
            case ENOENT:
                {
                    if(!util::makedir(util::dirname(fpath))){
                        return false;
                    }
                }
#ifdef _WIN32
                return 0 == _mkdir(fpath.c_str());
#else
                return 0 == mkdir(fpath.c_str(), mode);
#endif
            case EEXIST:
                return util::isdir(fpath);
            default:
                return false;
        }
    }

    /** touch a file
     * @param path file path to touch
     */
    inline void touch(const std::string& path){
        if(util::exists(path)) return;
        std::string dirpath = util::dirname(path);
        if(!util::exists(dirpath)) util::makedir(dirpath);
        std::ofstream fw(path);
        fw.close();
    }

    /** remove non-alpha characters from a string
     * @param str string to be filtered
     * @return a string without non-alpha characters
     */
    inline std::string getAlpha(const std::string& str){
        std::string ret;
        std::copy_if(str.cbegin(), str.cend(), std::back_inserter(ret), (int(*)(int))std::isalpha);
        return ret;
    }

    /** remove invalid sequence characters from a string
     * @param str string to be filtered
     * @param upper convert result to uppercase if true
     */
    inline void getValid(std::string& str, bool upper = false){
        uint16_t total = 0;
        for(uint16_t i = 0; i < str.size(); ++i){
            if(std::isalpha(str[i]) || str[i] == '-' || str[i] == '*'){
                str[total++] = (upper ? std::toupper(str[i]) : str[i]);
            }
        }
        str.resize(total);
    }

    /** make a string each character uppercased
     * @param str string to be uppercased
     */
    inline void str2upper(std::string& str){
        std::transform(str.begin(), str.end(), str.begin(), (int (*)(int))std::toupper);
    }

    inline void nuc2upper(std::string& str){
        for(size_t i = 0; i < str.size(); ++i) str[i] = nuc_to_upper[(int)str[i]];
    }

    /** make a string each character lowercased
     * @param str string to be lowercased
     */
    inline void str2lower(std::string& str){
        std::transform(str.begin(), str.end(), str.begin(), (int (*)(int))std::tolower);
    }
    inline void nuc2lower(std::string& str){
        for(size_t i = 0; i < str.size(); ++i) str[i] = nuc_to_lower[(int)str[i]];
    }

    /** get hamming distance of two strings
     * @param str1 string 1
     * @param str2 string 2
     * @return hamming distance of string 1 and string 2
     */
    inline int hamming(const std::string& str1, const std::string& str2){
        int diff = std::abs((int)(str1.size() - str2.size()));
        for(uint16_t i = 0; i < std::min(str1.size(), str2.size()); ++i){
            diff += (str1[i] == str2[i] ? 0 : 1);
        }
        return diff;
    }

    /** convert number to 33 based score character
     * @param num number of score
     * @return 33 based score character
     */
    inline char num2qual(int num){
        if(num > 127 - 33){
            num = 127 - 33;
        }
        if(num < 0){
            num = 0;
        }
        char c = num + 33;
        return c;
    }

    /** convert a quality uint8_t to phred33 score string
     * @param a pointer to uint8_t likely array
     * @param l length of array
     * @return string representation of phred33 based score
     */
    template<typename T>
    inline std::string qual2str(T *a, uint16_t l){
        std::string qstr(l, '\0');
        for(uint16_t i = 0; i < l; ++i){
            qstr[i] = (char)(33 + a[i]);
        }
        return qstr;
    }

    /** get complement base of a nucleotide base
     * @param base nucleotide base character
     * @return the uppercased complementary nucleotide base
     */
    inline char complement(const char& base){
        return nuc_to_cmp[(int)base];
    }

    /** get reverse completement sequence of a nucleotide sequence
     * @param seq a nucleotide sequence
     * @return the reverse completement sequence of seq
     */
    inline std::string revComp2NewSeq(const std::string& seq){
        std::string retSeq(seq.length(), '\0');
        for(int32_t i = retSeq.length() - 1; i >= 0; --i){
            retSeq[i] = nuc_to_cmp[(int)seq[seq.length() - 1 - i]];
        }
        return retSeq;
    }
    
    /** reverse completement an sequence of a nucleotide sequence
     * @param seq a nucleotide sequence to be reverse complemented
     */
    inline void revComp2OriSeq(std::string& seq){
        int32_t i = seq.size() - 1;
        int32_t j = 0;
        while(j < i){
            seq[i] = nuc_to_cmp[(int)seq[i]];
            seq[j] = nuc_to_cmp[(int)seq[j]];
            seq[i] = seq[i] ^ seq[j];
            seq[j] = seq[i] ^ seq[j];
            seq[i] = seq[i] ^ seq[j];
            ++j;
            --i;
        }
        if(i == j) seq[i] = nuc_to_cmp[(int)seq[i]];
    }

    /** reverse complement an sequence of a nucleotide sequence
     * @param seq a nucleotide sequence to be reverse complemented
     * @param len seq length
     */
    inline void revComp2OriSeq(char* seq, int len){
        int i = len - 1;
        int j = 0;
        while(j < i){
              seq[i] = nuc_to_cmp[(int)seq[i]];
              seq[j] = nuc_to_cmp[(int)seq[j]];
              seq[i] = seq[i] ^ seq[j];
              seq[j] = seq[i] ^ seq[j];
              seq[i] = seq[i] ^ seq[j];
              ++j;
              --i;
        }
        if(i == j) seq[i] = nuc_to_cmp[(int)seq[i]];
    }

    /** reverse complement an sequence of a nucleotide sequence
     * @param seq a nucleotide sequence to be reverse complemented
     * @param len seq length
     * @return new seq rev comp
     */
    inline char* revComp2NewSeq(char* seq, int len){
        char* ret = (char*)malloc(len*sizeof(char));
        for(int j = len-1; j >= 0; --j){
            ret[len-1-j] = nuc_to_cmp[(int)seq[j]];
        }
        return ret;
    }

    template<typename T>
    inline void reverse(T* seq, int len){
        int bpos = 0, epos = len - 1;
        while(bpos < epos){
            seq[bpos] = seq[bpos] ^ seq[epos];
            seq[epos] = seq[bpos] ^ seq[epos];
            seq[bpos] = seq[bpos] ^ seq[epos];
            ++bpos;
            --epos;
        }
    }
    
    /** get forward completment sequene of a nucleotide sequence
     * @param seq a nucleotide sequence
     * @return the forward completement sequence of seq
     */
    inline std::string forwardComplement(const std::string& seq){
        std::string retSeq(seq.length(), '\0');
        for(uint32_t i = 0; i < retSeq.length(); ++i){
            retSeq[i] = complement(seq[i]);
        }
        return retSeq;
    }

    /** write a log message in a thread-safe way
     * @param s log message 
     * @param logmtx reference to a std::mutex object
     * @param f log file
     */
    inline void loginfo(const char* s, std::mutex& logmtx, FILE* f = stderr){
        std::lock_guard<std::mutex> l(logmtx);
        time_t tt = time(NULL);
        tm* t = std::localtime(&tt);
        fprintf(f, "[%d-%02d-%02d %02d:%02d:%02d] %s\n", t->tm_year + 1900, t->tm_mon + 1, t->tm_mday, t->tm_hour, t->tm_min, t->tm_sec, s);
        fflush(f);
    }

    /** write a log message 
     * @param s log message 
     */
    inline void loginfo(const char* s, FILE* f = stderr){
        time_t tt = time(NULL);
        tm* t = std::localtime(&tt);
        fprintf(f, "[%d-%02d-%02d %02d:%02d:%02d] %s\n", t->tm_year + 1900, t->tm_mon + 1, t->tm_mday, t->tm_hour, t->tm_min, t->tm_sec, s);
        fflush(f);
    }

    /** write a log message
     * @param f output destination
     * @param fmt format
     */
    inline void loginfo(FILE* f, const char* fmt, ...){
        va_list ap;
        va_start(ap, fmt);
        time_t tt = time(NULL);
        tm* t = localtime(&tt);
        fprintf(f, "[%d-%02d-%02d %02d:%02d:%02d] ", t->tm_year + 1900, t->tm_mon + 1, t->tm_mday, t->tm_hour, t->tm_min, t->tm_sec);
        vfprintf(f, fmt, ap);
        fprintf(f, "\n");
        va_end(ap);
        fflush(f);
    }

    /** write a log message in thread-safe way
     * @param f output destination
     * @param fmt format
     */
    inline void loginfo(FILE* f, std::mutex& logmtx, const char* fmt, ...){
        std::lock_guard<std::mutex> l(logmtx);
        va_list ap;
        va_start(ap, fmt);
        time_t tt = time(NULL);
        tm* t = localtime(&tt);
        fprintf(f, "[%d-%02d-%02d %02d:%02d:%02d] ", t->tm_year + 1900, t->tm_mon + 1, t->tm_mday, t->tm_hour, t->tm_min, t->tm_sec);
        vfprintf(f, fmt, ap);
        fprintf(f, "\n");
        va_end(ap);
        fflush(f);
    }

    /** write a log message to std::cerr in a thread-safe way
     * @param s log message 
     * @param logmtx reference to a std::mutex object
     * @param f log file
     */
    inline void loginfo(const std::string& s, std::mutex& logmtx, FILE* f = stderr){
        loginfo(s.c_str(), logmtx, f);
    }

    /** write a log message to std::cerr directly
     * @param s log message 
     */
    inline void loginfo(const std::string& s, FILE* f = stderr){
        loginfo(s.c_str(), f);
    }

    /** get current time
     * @return year-mm-dd hh-mm-ss of current time
     */
    inline std::string currentTime(void){
        time_t tt = time(NULL);
        tm* t = std::localtime(&tt);
        char date[60] = {0};
        std::snprintf(date, 60, "%d-%02d-%02d %02d:%02d:%02d",
                t->tm_year + 1900, t->tm_mon + 1, t->tm_mday,
                t->tm_hour, t->tm_min, t->tm_sec);
        return date;
    }

    /** make a list from file by line 
     * @param filename input file
     * @param ret vector to store stripped line
     */ 
    inline void makeListFromFileByLine(const std::string& filename, std::vector<std::string>& ret){
        LineReader fr(filename);
        std::string line;
        while(fr.getline(line)){
            util::strip(line);
            ret.push_back(line);
        }
    }

    /** mak a set fro file by line
     * @param filename input file
     * @param ret set to store stripped line
     */
    template <typename T>
    inline void makeSetFromFileByLine(const std::string& filename, T& ret){
        LineReader fr(filename);
        std::string line;
        while(fr.getline(line)){
            util::strip(line);
            ret.insert(line);
        }
    }

    /** make a list from file by fields
     * @param filename input file
     * @param ret vector to store fields
     * @param field field number, 0 based
     * @param sep seperator
     */
    inline void makeListFromFileByField(const std::string& filename, std::vector<std::string>& ret, int field, std::string sep = "\t"){
        LineReader fr(filename);
        std::string line;
        std::vector<std::string> vstr;
        while(fr.getline(line)){
            split(line, vstr, sep);
            ret.push_back(vstr[field]);
        }
    }
    
    /** make a set from file by fields
     * @param filename input file
     * @param ret set to store fields
     * @param field field number, 0 based
     * @param sep seperator
     */
    inline void makeSetFromFileByField(const std::string& filename, std::set<std::string>& ret, int field, std::string sep = "\t"){
        LineReader fr(filename);
        std::string line;
        std::vector<std::string> vstr;
        while(fr.getline(line)){
            split(line, vstr, sep);
            ret.insert(vstr[field]);
        }
    }

    /** make map from file by line, row1 is key, row2 is value, row1 shold be unique
     * @param filename input file
     * @param ret map to store key value pairs
     * @param sep seperator of fields in line
     */
    template<typename M>
    inline void makeMapPairFromFileByLine(const std::string& filename, M& ret, const std::string& sep = "\t"){
        LineReader fr(filename);
        std::string line;
        std::vector<std::string> vstr;
        while(fr.getline(line)){
            util::split(line, vstr, sep);
            ret[vstr[0]] = vstr[1];
        }
    }

    /** make map from file by line, row1 is key, row2 is value, row1 might have dup
     * @param filename input file
     * @param ret map to store key value pairs
     * @param sep seperator of fields in line
     */
    template<typename M>
    inline void makeMapSetFromFileByLine(const std::string& filename, M& ret, const std::string& sep = "\t"){
        LineReader fr(filename);
        std::string line;
        std::vector<std::string> vstr;
        while(fr.getline(line)){
            util::split(line, vstr, sep);
            ret[vstr[0]].insert(vstr[1]);
        }
    }

    /** test whether an element is in a vector
     * @param v vector
     * @param e element
     * @return true if e in v
     */
    template<typename T>
    inline bool inVector(const std::vector<T>& v, const T& e){
        return std::find(v.cbegin(), v.cend(), e) != v.cend();
    }

    /** convert a vector of strings to integers
     * @param vs vector of strings
     * @param vi vector of integers
     */
    template<typename T>
    inline void strvec2intvec(const std::vector<std::string>& vs, std::vector<T>& vi){
        vi.clear();
        if(vs.empty()) return;
        vi.reserve(vs.size());
        for(auto& s: vs) vi.push_back(atoi(s.c_str()));
    }

    /** convert two vector to a vector of pairs
     * @param v1 vector 1
     * @param v2 vector 2
     * @param vp pair vector
     */
    template<typename T>
    inline void vec2pairvec(const std::vector<T>& v1, const std::vector<T>& v2, std::vector<std::pair<T, T>>& vp){
        vp.clear();
        int maxlen = std::min(v1.size(), v2.size());
        for(int i = 0; i < maxlen; ++i){
            vp.push_back(std::make_pair(v1[i], v2[i]));
        }
    }

    /** list contents in an directory
     * @param path string of path
     * @param vname vector to store contents under path
     * @param nohide do not list hiden files
     */
    inline void listDir(const std::string& path, std::vector<std::string>& vname, bool nohide = true){
        std::string dirpath = util::replace(path, "~", std::string(std::getenv("HOME")) + "/");
        DIR *dir;
        struct dirent *ent;
        if((dir = opendir(dirpath.c_str())) != NULL){
            while((ent = readdir(dir)) != NULL){
                if(nohide && ent->d_name[0] == '.') continue;
                vname.push_back(std::string(ent->d_name));
            }
            closedir(dir);
        }
    }

    /** print an array to std::cout 
     * @param a pointer to an array
     * @param l length of the array
     * @param n array name to show
     */
    template<typename T>
    inline void showArray(const T* a, uint16_t l, const std::string& n){
        std::cout << n << ":";
        for(uint16_t i = 0; i < l; ++i){
            std::cout << a[i] << " ";
        }
        std::cout << std::endl;
    }

    /** count different neighbor pairs in a string
     * @param str a string
     * @return different neighbor pairs in this string
     */
    inline int neighborDiffCount(const std::string& str){
        int diff = 0;
        for(uint32_t i = 0; i < str.length() - 1; ++i){
            if(str[i] != str[i+1]){
                ++diff;
            }
        }
        return diff;
    }
    
    /** get mismatch ratio of two uncleotide sequence
     * @param s1 nucleotide sequence
     * @param s2 nucleotide sequence
     * @return mismatch ratio to the shorter string
     */
    inline float mismatchRatio(const std::string& s1, const std::string& s2){
        uint32_t minLen = std::min(s1.length(), s2.length());
        if(minLen == 0 || s1.find_first_of("ATCG") == std::string::npos || s2.find_first_of("ATCG") == std::string::npos){
            return 1.0;
        }
        if(minLen == 1){
            return s1[0] == s2[0];
        }
        size_t beg = 0;
        for(beg = 0; beg < minLen; ++beg){
            if(s1[beg] == 'N' || s2[beg] == 'N'){
                ++beg;
            }
        }
        size_t end = minLen - 1;
        for(end = minLen - 1; end > beg; --end){
            if(s1[end] == 'N' || s2[end] == 'N'){
                --end;
            }
        }
        int32_t mismatch = 0;
        for(uint32_t i = beg; i <= end; ++i){
            if(s1[i] == 'N' || s2[i] == 'N'){
                continue;
            }
            if(s1[i] != s2[i]){
                ++mismatch;
            }
        }
        return float(mismatch)/(end - beg + 1);
    }

    /** get median of a list of values
     * @param v a vector of values
     * @return median of values in v
     */
    template<typename T>
    inline double median(std::vector<T>& v){
        std::nth_element(v.begin(), v.begin() + v.size()/2, v.end()); 
        return v[v.size()/2];
    }

    /** test whether a string contains only AaCcGgTt (ie. is DNA)
     * @param seq string of characters
     * @return true if seq contains only AaCcGgTt
     */
    inline bool isDNA(const std::string& seq){
        return seq.find_first_not_of("AaTtCcGg") == std::string::npos;
    }

    /** count number of characters in a string
     * @param str full string
     * @param chrs characters to be found in str
     */
    inline int32_t countChrs(const std::string& str, const std::string& chrs){
        int32_t count = 0;
        std::string::size_type pos = 0;
        while((pos = str.find_first_of(chrs, pos)) != std::string::npos){
            ++pos;
            ++count;
        }
        return count;
    }

    /** format an integer to include commas
     * @param n data number to format
     * @return string with formatted number containing commas
     */
    template<typename T>
    inline std::string formatIntWithCommas(const T& n){
        std::string s = std::to_string(n);
        if(s.length() > 3){
            for(int i = s.length() - 3; i > 0; i -=3){
                s.insert(i, ",");
            }
        }
        return s;
    }

    /** divide vector index into even split
     * @param vsize original vector size to split
     * @param nsplit number of splits expected
     * @param vpidx split result vector
     * @param atmostn at most nsplit if true
     * @return actual split got
     */
    inline int32_t divideVecIdx(int32_t vsize, int32_t nsplit, std::vector<std::pair<int32_t, int32_t>>& vpidx, bool atmostn = true){
        vpidx.clear();
        if(!vsize) return 0;
        if(!nsplit){
            vpidx.push_back({0, vsize});
            return 1;
        }
        nsplit = std::min(vsize, nsplit);
        int32_t perSize = ceil((double)vsize / (double)nsplit);
        for(int32_t iidx = 0; iidx < nsplit; ++iidx){
            std::pair<int32_t, int32_t> p;
            p.first = iidx * perSize;
            p.second = (iidx + 1) * perSize;
            if(p.second <= vsize) vpidx.push_back(p);
            else break;
        }
        if(vpidx.size()){
            if(vpidx[vpidx.size()-1].second < vsize){
                if(atmostn){
                    vpidx[vpidx.size()-1].second = vsize;
                }else{
                    std::pair<int32_t, int32_t> p;
                    p.first = vpidx[vpidx.size()-1].second;
                    p.second = vsize;
                    vpidx.push_back(p);
                }
            }
        }else vpidx.push_back({0, vsize});
        return vpidx.size();
    }

    /** compute combination number by Pascal's triangle(https://en.wikipedia.org/wiki/Pascal%27s_triangle)
     * @param n total number of items to select from
     * @param k number of items to get
     * @return number of ways to get k items from total of n
     */
    template<typename T>
    inline T combNum(T n, T k){
        if(n < k) return 0;
        std::vector<T> row(std::max(n + 1, k + 1), 0);
        T i, j;
        row[0] = 1;
        for(i = 1; i <= n; ++i){
            for(j = i; j > 0; --j){
                row[j] += row[j - 1];
            }
        }
        return row[k];
    }

    /** generate combination set
     * @param n total number of items to select from
     * @param k number of items to get
     * @param ret store bit masks of all ways to select k items from total n
     */
    inline void combSet(int n, int k, std::vector<std::vector<int>>& ret){
        std::vector<bool> bitmask(n, 0);
        ret.reserve(combNum(n, k));
        for(int i = 0; i < k; ++i) bitmask[i] = 1;
        do{
            std::vector<int> sel; sel.reserve(k);
            for(int i = 0; i < n; ++i){
                if(bitmask[i]) sel.push_back(i);
            }
            ret.push_back(sel);
        }while(std::prev_permutation(bitmask.begin(), bitmask.end()));
    }

    /** remove en elements in an vector
     * @param vec a vector storing elements of type T
     * @param e element to be removed
     * @param k index of last element removed
     * @param return true if e is removed from vec
     */
    template<typename T>
    inline bool rmElemFromVec(std::vector<T>& vec, const T& e, int32_t& k){
        size_t i = 0, j = 0;
        for(i = 0; i < vec.size(); ++i){
            if(vec[i] != e){
                vec[j++] = vec[i];
            }else{
                k = i;
            }
        }
        vec.resize(j);
        return j < i;
    }

    /** remove ith element of an vector
     * @param vec a vector storing elements of type T
     * @param idx index at which element will be removed
     */
    template<typename T>
    inline void rmIdxFromVec(std::vector<T>& vec, const size_t& idx){
        size_t i = idx, j = idx;
        for(i = idx + 1; i < vec.size(); ++i){
           vec[j++] = vec[i];
        } 
        vec.resize(j);
    }

    /** remove a directory recursivelly, https://stackoverflow.com/a/5467788/2297529 */
    inline int unlink_cb(const char *fpath, const struct stat *, int, struct FTW *){
        int rv = remove(fpath);
        if (rv) perror(fpath);
        return rv;
    }
    inline int rmrf(const char *path){
        return nftw(path, unlink_cb, 64, FTW_DEPTH | FTW_PHYS);
    }

    /** test whether a command is in PATH env
     * @param bin program to test
     * @param env environment to find bin
     * @return path/to/bin
     */
    inline char* which(const char* name, const char* env = NULL){
#ifdef _WIN32
#define WHICH_DELIMITER   ";"
#else
#define WHICH_DELIMITER   ":"
#endif
        if(!env) env = "PATH";
        char* path = strdup(getenv(env));
        char *tok = strtok(path, WHICH_DELIMITER);
        while(tok){
            // path
            int len = strlen(tok) + 2 + strlen(name);
            char *file = (char*)malloc(len*sizeof(char));
            if(!file){ free(path); return NULL;}
            snprintf(file, len, "%s/%s", tok, name);
            // executable
            if(access(file, X_OK) == 0){
                free(path);
                return file;
            }
            // next token
            tok = strtok(NULL, WHICH_DELIMITER);
            free(file);
        }
        free(path);
        return NULL;
    }

    /** get longest single character repeat sequence length in a string
     * @param seq string to check
     * @param len length of seq
     * @return longest single char repeat length in str
     */
    inline int longestRepeatLen(const char* seq, int len){
        if(len == 0) return 0;
        int maxrpl = 0;
        int tmprpl = 0;
        for(int i = 0; i < len; ++i){
            tmprpl = 1;
            while(i+1 < len && seq[i+1] == seq[i]){
                ++i;
                ++tmprpl;
            }
            if(tmprpl > maxrpl) maxrpl = tmprpl;
        }
        return maxrpl;
    }
    
    /** Find the first occurrence of find in s, where the search is limited to the first slen characters of s
     * @param s full str to search into
     * @param find patt to find
     * @param len len of find
     * @param slen search length
     */
    inline char* strnstr( const char *s, const char *find, int len,  int slen){
	char c, sc;
	if((c = *find++) != '\0'){
            do{
                do{
                    if(slen-- < 1 || (sc = *s++) == '\0') return NULL;
                }while (sc != c);
                if(len > slen+1) return NULL;
            }while(strncmp(s, find, len-1) != 0);
	    s--;
	}
	return (char*)s;
    }

    /** write a string of characters to file
     * @param seq seuqnce
     * @param len seq length
     * @para file file name
     */
    inline void writestr(char* seq, int len, char* file){
        FILE* fp = fopen(file, "w");
        fwrite(seq, len, sizeof(char), fp);
        fclose(fp);
    }
}
#endif
