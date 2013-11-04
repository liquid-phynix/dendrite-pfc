#include "timing.hpp"
#include <cstdio>

void TimeAcc::start(){ gettimeofday(&m_tval, NULL); }
void TimeAcc::stop(bool incr){
  timeval end;
  gettimeofday(&end, NULL);
  if(incr) m_times++;
  m_elapsed += (end.tv_sec - m_tval.tv_sec) * 1000 * 1000 + (end.tv_usec - m_tval.tv_usec);
}
void TimeAcc::report(const char* msg){
  int t = m_times == 0 ? 0 : (float) m_elapsed / (float) m_times / 1000.0;
  printf(msg, t);
  m_times = 0;
  m_elapsed = 0;
}
void TimeAcc::report2(const char* msg){
  int t = m_times == 0 ? 0 : (float) m_elapsed / (float) m_times / 1000.0;
  printf(msg, t, m_times);
  m_times = 0;
  m_elapsed = 0;
}
