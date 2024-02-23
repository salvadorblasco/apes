"""Generate data for linsearch test."""


def generate_data(beta, stoich, analyticalc, x0, dx):
    def g(lam):
        import numpy as np
        from libeq.fobj import fobj
        from libeq.cpluss import cpluss
        for l in iter(lam):
            x1 = x0+l*dx
            Fast = fobj(cpluss(x1, beta, stoich, full=True), stoich, analyticalc)
            gast = 0.5*np.sum(np.square(Fast), axis=1)
            yield gast

    return g


def generate_data2(beta, stoich, analyticalc, x0, dx):
    def g(lam):
        import numpy as np
        from libeq.fobj import fobj
        from libeq.cpluss import cpluss
        x1 = x0+lam[:, None]*dx
        Fast = fobj(cpluss(x1, beta, stoich, full=True), stoich, analyticalc)
        gast = 0.5*np.sum(np.square(Fast), axis=1)
        return gast

    return g
