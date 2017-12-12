#include <time.h>
#include <sys/time.h>

struct MD_Timer {
  struct timeval clock_holder;
  struct timeval duration;
};

void initialize_timer (struct MD_Timer* time);
void start_timer (struct MD_Timer* time);
void stop_timer (struct MD_Timer* time);
double timer_duration(const struct MD_Timer time);
