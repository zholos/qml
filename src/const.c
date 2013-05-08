#include <float.h>

#include "util.h"

#define QUOTE_(x) #x
#define QUOTE(x) QUOTE_(x)

// Embed information that can be retrieved by what or strings
char ident[] =
    "@(#)qml " QUOTE(QML_VERSION) " " QUOTE(KXVER) "/" QUOTE(KXARCH);

K
qml_const(K x) {
    check_type(qt(x) == -KL,);
    switch (ql(x)) {
    case 0:
        return ks(QUOTE(QML_VERSION));
    case 1:
        return kf(3.1415926535897932384626433832795028842);
    case 2:
        return kf(2.7182818284590452353602874713526624978);
    case 3:
        return kf(DBL_EPSILON);
    }
    return krr("domain");
}
