#ifndef QML_SRC_CONMIN_H
#define QML_SRC_CONMIN_H

#include "util.h"

F* take_param(K x, I* n, S* err);
F* make_param(K x, F* param, K* r);


// krr() might return NULL or an error object.
// This dummy object represents no error.
extern struct k0 no_error_;
#define no_error (&no_error_)

// Used internally as both a default empty list and as a flag.
extern struct k0 empty_con_;
#define empty_con (&empty_con_)


#endif // QML_SRC_CONMIN_H
