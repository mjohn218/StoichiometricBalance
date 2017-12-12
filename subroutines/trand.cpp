#include "utility_calls.h"

int trand()
{
  int i=rand();
  while(i==RAND_MAX)
    i=rand();
  return i;
}
