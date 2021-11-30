#ifndef Date_H_
#define Date_H_

#include <string>
#include <sstream>
#include <iomanip>

class Date {
private:
    int hh_{0};
    int mm_{0};
    int ss_{0};

public:
    Date() {}
    Date(int hh, int mm, int ss) : hh_(hh), mm_(mm), ss_(ss) {}

    ~Date() {}

public:
    void set_hh(int hh) { hh_ = hh; }
    void set_mm(int mm) { mm_ = mm; }
    void set_ss(int ss) { ss_ = ss; }

    int hh() const { return hh_; };
    int mm() const { return mm_; };
    int ss() const { return ss_; };

    int hhmmss() const { return 10000*hh() + 100*mm() + ss(); }

    std::string hhmmss_str() const { 
        std::ostringstream sout;
        sout << std::setfill('0') << std::setw(6) << hhmmss();

        return sout.str();
    }
};


#endif
