#include "nlopt-util.h"

// We don't use the time-based stopping cretirion.
double nlopt_seconds() {
    return 0;
}
