#include <string.h>

#include "opt.h"

#include "util.h"


static int
get(K x, I i, int flat, const S s, const struct optn* n, union optv* v)
{
    if (!s)
        return 0;
    if (*s == '\0')
        return 1;

    I k;
    for (k = 0; n[k].s; k++)
        if (!strcmp(s, n[k].s))
            goto found;
    return 0;

found:
    if (!n[k].t && flat) {
        v[k].i = 1;
        return 1;
    } else {
        if (!x || !has_n(x) || i >= qn(x))
            return 0;
        if (n[k].t == -KF)
            return 2 * !!item_F(&v[k].f, x, i);
        if (n[k].t == -KI || !n[k].t)
            return 2 * !!item_I(&v[k].i, x, i);
        assert(0);
        return 0;
    }
}


int
take_opt(K x, const struct optn* n, union optv* v) {
    switch (qt(x)) {
    case 0:
        repeat (i, qn(x)) {
            S s = qt(qK(x, i)) == -KS ? qK(x, i)->s : NULL;
            int r = get(x, i+1, 1, s, n, v);
            if (!r)
                return 0;
            i += r-1;
        }
        break;
    case -KS:
        return get(NULL, 0, 1, xs, n, v);
    case KS:
        repeat (i, qn(x))
            if (!get(NULL, 0, 1, xS[i], n, v))
                return 0;
        break;
    case XD:
        if (qt(xx) != KS && qt(xx) != 0)
            return 0;
        repeat (i, qn(xx)) {
            S s = qt(xx) == KS         ? kS(xx)[i]    :
                  qt(qK(xx, i)) == -KS ? qK(xx, i)->s : NULL;
            if (!get(xy, i, 0, s, n, v))
                return 0;
        }
    }
    return 1;
}
