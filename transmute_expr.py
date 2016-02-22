import re
import json
import numpy as np
import sympy
from sympy.codegen.ast import Assignment, CodeBlock

with open('transmute_data.json') as f:
    DATA = json.load(f)


t, G = sympy.symbols('t G')
phi = sympy.MatrixSymbol('phi', G, 1)

decay_rxs = ['bminus', 'bplus', 'ec', 'alpha', 'it', 'sf', 'bminus_n']
xs_rxs = ['gamma', 'z_2n', 'z_3n', 'alpha', 'fission', 'proton', 'gamma_1', 'z_2n_1']

gamma_base = re.compile('gamma_([A-Z][a-z]?\d+)_')

def child_decays(nuc):
    expr = 0
    for rx in decay_rxs:
        r = re.compile(gamma_base + nuc + '_' + rx)
        for key in DATA['symbols']:
            m = r.match(key)
            if m is not None:
                parname = m.group(1)
                gammaname = m.group(0)
                break
        else:
            continue
        gamma = DATA['symbols'][gammaname]
        lambda_par = DATA['symbols']['lambda_' + parname]
        par0 = sympy.symbols('{0}_0'.format(parname))
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

    system = CodeBlock(*map(gennuc, nucs))

    with open('system.txt', 'w') as f:
        for eq in system.args:
            f.write(sympy.srepr(eq) + '\n')

    with open('system-C.txt', 'w') as f:
        f.write(sympy.ccode(system))

    system_cse = system.cse()

    with open('system-cse.txt', 'w') as f:
        for eq in system_cse.args:
            f.write(str(eq) + '\n')

    with open('system-cse-C.txt', 'w') as f:
        f.write(sympy.ccode(system_cse))
