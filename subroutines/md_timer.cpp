#include <iostream>
#include <fstream>
#include "md_timer.h"

using namespace std;

void initialize_timer (struct MD_Timer * time) {
  time->clock_holder.tv_sec = 0;
  time->clock_holder.tv_usec = 0;
  time->duration.tv_sec = 0;
  time->duration.tv_usec = 0;
}

void start_timer(struct MD_Timer * time) {
  gettimeofday (&time->clock_holder, NULL);
}

void stop_timer(struct MD_Timer * time) {
  struct timeval end_tv;
  gettimeofday (&end_tv, NULL);
  time->duration.tv_sec += (end_tv.tv_sec - time->clock_holder.tv_sec);
  time->duration.tv_usec += (end_tv.tv_usec - time->clock_holder.tv_usec);
}

double timer_duration(const struct MD_Timer time) {
  return time.duration.tv_sec + 1.0e-6 * (double)time.duration.tv_usec;
}

