#ifndef _WRITER_H
#define _WRITER_H

#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <zlib.h>
#include "util.h"

// Class to write to gz file or ofstream
class Writer{
        std::string mFilename;  // output filename
        gzFile mGzFile;         // gzFile file handler
        std::ofstream* mStream; // pointer to ofstream
        bool mZipped;           // output file is mZipped or not
        int mCompressLevel;     // compression level for gz file

    public:
        /** Writer constructor
         * @param filename output filename
         * @param compression compression level for gzFile
         */
        Writer(const std::string& filename, const int& compression = 3);
        
        /** Writer constructor
         * @param mStream pointer to ofstream
         */
        Writer(std::ofstream* mStream);
        
        /** Writer constructor
         * @param gzfile gzFile handler
         */
        Writer(gzFile gzfile);
        
        /** Writer destructor
         */
        ~Writer();

        /** Whether the output file is mZipped or not
         * @return true if output filename ends with .gz
         */
        bool isZipped();
        
        /** write a string to file without append \\n
         * @param str string to be written 
         * @return true if successfully written
         */
        bool writeString(const std::string& str);
        
        /** write a string to file with additional \\n appended 
         * @param linestr string to be written
         * @return true if successfully written
         */
        bool writeLine(const std::string& linestr);
        
        /** write a C string to file
         * @param cstr pointer to char(C string)
         * @param size length of the C string to be written
         * @return true if successfully written
         */
        bool write(char* cstr, size_t size);
        
        /** get filename of output file
         * @return output filename
         */
        std::string getFilename();

        /** initialize Writer, detect file format and open file handler
         */
        void init();
        
        /** flush buffer and close file handler
         */
        void close();
};

#endif
