#ifndef TIMER_H
#define TIMER_H

#include <chrono>
#include <string>
#include <sstream>

// a simple class to represent a timer for time elapse calculation
class Timer{
    public:
        // constructor
        Timer() : beg_(clock_::now()) {}
        
        // reset the timer, set beginning time to current time
        void reset() { beg_ = clock_::now(); }
        
        /** calculate elapsed time from the beginning
         * @return seconds elapsed from the beginning
         */
        double elapsed() const {
            return std::chrono::duration_cast<second_>(clock_::now() - beg_).count(); 
        }

        /** output time elapsed to ostream
         * @param os reference of ostream
         * @param t reference of Timer
         * @return reference of ostream
         */
        friend inline std::ostream& operator<<(std::ostream& os, const Timer& t){
            int64_t seconds = std::chrono::duration_cast<second_>(clock_::now() - t.beg_).count();
            os << seconds / 86400 << "D,";
            os << (seconds % 86400) / 3600 << "H,";
            os << (seconds % 3600) / 60 << "M,";
            os << seconds % 60 << "S";
            return os;
        }
        
        /** get string representation of time elapsed
         * @return human readble time elspse
         */
        std::string toStr() const {
            int64_t seconds = std::chrono::duration_cast<second_>(clock_::now() - beg_).count();
            std::stringstream oss;
            int32_t dd = seconds / 86400;
            int32_t hh = (seconds % 86400) / 3600;
            int32_t mm = (seconds % 3600) / 60;
            int32_t ss = (seconds % 60);
            oss << dd << "D," << hh << "H," << mm << "M," << ss << "S";
            return oss.str();
        }
    
    private:
        // a clock_ type to represent wall clock time
        typedef std::chrono::high_resolution_clock clock_;

        // a second_ type to represent second units
        typedef std::chrono::duration<double, std::ratio<1> > second_;
        
        // beg_ type to represent time point
        std::chrono::time_point<clock_> beg_;
};

#endif
