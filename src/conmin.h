#ifndef QML_SRC_CONMIN_H
#define QML_SRC_CONMIN_H

#include "util.h"

F* take_param(K x, I* n, S* err);
F* make_param(K x, F* param, K* r);


struct call_info {
    int arg; // 0 = make_param(start), 1/-1 = base+arg*scalar
    F base;
    K start;
    K error;
} call;

F call_param(struct call_info* info, int sign, K f, F* param);


struct eval_info {
    F first_member; // for CONMAX, makes passing struct as pttbl well-defined
    struct call_info call;
    K fun, con;
    I contyp;
    int con_sign;
};

F eval_param(struct eval_info* info,
             I which, F* param, I n, F* grad, int grad_step, I* contyp);


// krr() might return NULL or an error object.
// This dummy object represents no error.
extern struct k0 no_error_;
#define no_error (&no_error_)

// Used internally as both a default empty list and as a flag.
extern struct k0 empty_con_;
#define empty_con (&empty_con_)


#endif // QML_SRC_CONMIN_H
