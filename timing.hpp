#include <sys/time.h>

class TimeAcc {
  int m_times;
  long long int m_elapsed;
  timeval m_tval;
public:
  TimeAcc(): m_times(), m_elapsed(), m_tval(){}
  void start();
  void stop(bool = true);
  void report(const char*);
  void report2(const char*);
};
