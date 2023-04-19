#ifndef SOFTWARE_H
#define SOFTWARE_H

#include <string>
#include <chrono>
#include <unistd.h>
#include <libgen.h>
#include "timer.h"
#include "common.h"
#include "htslib/kstring.h"

// class to store software environment
struct Software{
    std::string ver; // software version
    std::string cmd; // software execution command
    std::string cwd; // software execution directory
    std::string dur; // software execution time
    std::string cmp; // software update time
    std::string nam; // software name
    std::string sub; // sub module name
    std::string ctt; // contact
    Timer* timer;    // timer to recording execution time

    // constructor
    Software(){
        ver = PACKAGE_VERSION;
        nam = PACKAGE_NAME;
        ctt = PACKAGE_MAINTER_EMAIL;
        cmp = cmp = std::string(__TIME__) + " " + std::string(__DATE__);
        char cpath[FILENAME_MAX];
        getcwd(cpath, FILENAME_MAX);
        cwd = cpath;
        timer = new Timer();
    }

    // destructor
    ~Software(){
        if(timer){
            delete timer;
        }
    }

    // update
    inline void update(int argc, char** argv){
        for(int i = 0; i < argc; ++i){
            cmd.append(argv[i]);
            if(i < argc - 1) cmd.append(" ");
        }
        sub = argv[0];
    }

    // get exe time
    inline void end(){
        dur = timer->toStr();
    }

    /** get string representation of execution time
     * @return string representation of execution time
     */
    inline std::string getExecutionTime(){
        return timer->toStr();
    }

    // report json
    void reportJSON(kstring_t* s, const char* dh, const char* dm){
        ksprintf(s, "%s%s\"SoftwareName\": \"%s\",\n", dh, dm, PACKAGE_NAME);
        ksprintf(s, "%s%s\"SoftwareVersion\": \"%s\",\n", dh, dm,  PACKAGE_VERSION);
        ksprintf(s, "%s%s\"WorkingDirectory\": \"%s\",\n", dh, dm,  cwd.c_str());
        ksprintf(s, "%s%s\"CommandExecuted\": \"%s\",\n", dh, dm, cmd.c_str());
        ksprintf(s, "%s%s\"TimeExhausted\": \"%s\"", dh, dm, getExecutionTime().c_str());
    }
};

#endif
