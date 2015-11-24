import numpy as np
import sympy
Assignment = sympy.printing.codeprinter.Assignment
import numpy as np
import sympy
Assignment = sympy.printing.codeprinter.Assignment


from pyne import nucname
from pyne import data
from pyne import rxname


t, G = sympy.symbols('t G')
phi = sympy.MatrixSymbol('phi', G, 1)

decay_rxs = ['bminus', 'bplus', 'ec', 'alpha', 'it', 'sf', 'bminus_n']
xs_rxs = ['gamma', 'z_2n', 'z_3n', 'alpha', 'fission', 'proton', 'gamma_1', 'z_2n_1']


def child_decays(nuc):
    childname = nucname.name(nuc)
    expr = 0
    for rx in decay_rxs:
        try:
            parent = rxname.parent(nuc, rx, b'decay')
        except RuntimeError:
            continue
        if data.branch_ratio(parent, nuc) < 1e-16:
            continue
        parname = nucname.name(parent)
        par0, gamma, lambda_par = sympy.symbols('{0}_0 gamma_{0}_{1}_{2} lambda_{0}'.format(parname, childname, rx))
        expr += gamma * sympy.exp(lambda_par * t) * par0
    return expr

def child_xss(nuc):
    expr = 0
    for rx in xs_rxs:
        try:
            parent = rxname.parent(nuc, rx)
        except RuntimeError:
            continue
        parname = nucname.name(parent)
        par0 = sympy.symbols('{0}_0'.format(parname))
        sigma_rx_par = sympy.MatrixSymbol('sigma_{0}_{1}'.format(rx, parname), 1, G)
        expr += sympy.exp((sigma_rx_par*phi)[0] * t) * par0
    return expr

def gennuc(nuc):
    name = nucname.name(nuc)
    nuc0, nuc1 = sympy.symbols('{0}_0 {0}_1'.format(name))
    lambda_nuc = sympy.symbols('lambda_{0}'.format(name))
    sigma_a_nuc = sympy.MatrixSymbol('sigma_a_{0}'.format(name), 1, G)
    rhs = sympy.exp(-((sigma_a_nuc*phi)[0] + lambda_nuc)*t) * nuc0
    rhs += child_decays(nuc)
    rhs += child_xss(nuc)
    eq = Assignment(nuc1, rhs)
    return eq

if __name__ == '__main__':
    data.atomic_mass('U235')
    nucs = set(data.atomic_mass_map.keys())
    for nuc in data.atomic_mass_map:
        nucm = nuc + 1
        if nucname.anum(nuc) == 0 or data.decay_const(nucm) < 1e-16 or data.decay_const(nuc) == data.decay_const(nucm):
            continue
        nucs.add(nucm)
    nucs = [nuc for nuc in nucs if nucname.anum(nuc) > 0 and
                                not np.isnan(data.decay_const(nuc)) and
                                nuc < 200000000]
    nucs.sort()

    system = list(map(gennuc, nucs))

    with open('system.txt', 'w') as f:
        for eq in system:
            f.write(str(eq) + '\n')
