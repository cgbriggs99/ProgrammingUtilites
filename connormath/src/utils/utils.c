 
#include "connormath.h"

int betweenxx(double bound1, double bound2, double value) {
  if(bound1 < bound2) {
    return (value < bound2 && bound1 < value);
  } else {
    return (value > bound2 && bound1 > value);
  }
}

int betweenxi(double bound1, double bound2, double value) {
  if(bound1 < bound2) {
    return (value <= bound2 && bound1 < value);
  } else {
    return (value >= bound2 && bound1 > value);
  }
}

int betweenix(double bound1, double bound2, double value) {
  if(bound1 < bound2) {
    return (value < bound2 && bound1 <= value);
  } else {
    return (value > bound2 && bound1 >= value);
  }
}

int betweenii(double bound1, double bound2, double value) {
  if(bound1 < bound2) {
    return (value <= bound2 && bound1 <= value);
  } else {
    return (value >= bound2 && bound1 >= value);
  }
}

int relnear(double a, double b, double maxdiff) {
  return (2 * fabs(a - b) / (fabs(a) + fabs(b)) <= fabs(maxdiff));
}

int absnear(double a, double b, double maxdiff) {
  return (fabs(a - b) <= fabs(maxdiff));
}
