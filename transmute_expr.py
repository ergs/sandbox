import textwrap
import re
import json
import sympy
from sympy.codegen.ast import Assignment, CodeBlock

TEMPLATE = r"""\
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define I (%%I%%)
#define R (9)

double* transmute(double* N0, double t, double phi, double sigma[I][R]);

double* transmute(double* N0, double t, double phi, double sigma[I][R]) {

    double* N1 = malloc(I*sizeof(double));

%%CODE%%

    return(N1);
}

int main() {
    int i;
    double N0[I] = %%N0%%;
    double* N1;
    double sigma[I][R] = %%SIGMA_ARRAY%%;

    N1 = transmute(N0, 10.0, 4e-10, sigma);

    for (i=0; i < I; i++) {
        printf("%d %f\n", i, N1[i]);
    }
    return(0);
}

"""

with open('transmute_data.json') as f:
    DATA = json.load(f)


t= sympy.symbols('t')
# G = 1
# phi = sympy.MatrixSymbol('phi', G, 1)
phi = sympy.symbols('phi')

decay_rxs = ['bminus', 'bplus', 'ec', 'alpha', 'it', 'sf', 'bminus_n']
xs_rxs = ['gamma', 'z_2n', 'z_3n', 'alpha', 'fission', 'proton', 'gamma_1', 'z_2n_1']

gamma_base = '^gamma_([A-Z][a-z]?\d+)_'

def child_decays(nuc):
    symbols = DATA['symbols']
    expr = 0
    for rx in decay_rxs:
        r = re.compile(gamma_base + nuc + '_' + rx + '$')
        for key in symbols:
            m = r.match(key)
            if m is not None:
                parname = m.group(1)
                gammaname = m.group(0)
                break
        else:
            continue
        if parname not in DATA['nucs']:
            continue
        gamma = symbols[gammaname]
        lambda_par = symbols['lambda_' + parname]
        par0 = sympy.symbols('{0}_0'.format(parname))
        if lambda_par >= 0: # Avoid nan
            expr += gamma * sympy.exp(lambda_par * t) * par0
    return expr

def child_xss(nuc):
    rxs = DATA['channels'][nuc]
    terms = []
    for rx in xs_rxs:
        if rx not in rxs:
            continue
        parname = rxs[rx]
        par0 = sympy.symbols('{0}_0'.format(parname))
        # sigma_rx_par = sympy.MatrixSymbol('sigma_{0}_{1}'.format(rx, parname), 1, G)
        sigma_rx_par = sympy.Symbol('sigma_{0}_{1}'.format(rx, parname))
        # expr += sympy.exp((sigma_rx_par*phi)[0] * t) * par0
        terms.append(sympy.exp((sigma_rx_par*phi) * t) * par0)
    return sympy.Add(*terms)

def gennuc(nuc):
    nuc0, nuc1 = sympy.symbols('{0}_0 {0}_1'.format(nuc))
    lambda_nuc = DATA['symbols'].get('lambda_{0}'.format(nuc), sympy.oo)
    # sigma_a_nuc = sympy.MatrixSymbol('sigma_a_{0}'.format(nuc), 1, G)
    sigma_a_nuc = sympy.Symbol('sigma_a_{0}'.format(nuc))
    # rhs = sympy.exp(-((sigma_a_nuc*phi)[0] + lambda_nuc)*t) * nuc0
    if lambda_nuc == sympy.oo:
        rhs = 0
    else:
        rhs = sympy.exp(-((sigma_a_nuc*phi) + lambda_nuc)*t) * nuc0
    rhs += child_decays(nuc)
    rhs += child_xss(nuc)
    eq = Assignment(nuc1, rhs)
    return eq

def sigma_symbol_to_indexed():
    sigma_symbols = [[sympy.Symbol('sigma_{0}_{1}'.format(rx, nuc)) for rx in xs_rxs + ['a']] for nuc in DATA['nucs']]
    mapping = {sigma_symbols[i][j]: sympy.Symbol('sigma[{0}][{1}]'.format(i, j)) for i in range(len(sigma_symbols)) for j in range(9)}

    return mapping

def nuc_symbol_to_indexed():
    nucs = DATA['nucs']
    Symbol = sympy.Symbol
    return {
        **{Symbol('{0}_0'.format(nuc)): Symbol('N0[{0}]'.format(i)) for i,
            nuc in enumerate(nucs)},
        **{Symbol('{0}_1'.format(nuc)): Symbol('N1[{0}]'.format(i)) for i,
            nuc in enumerate(nucs)},
        }

def generate_sigma_array():
    sigma_symbols = [['sigma_{0}_{1}'.format(rx, nuc) for rx in xs_rxs + ['a']] for nuc in DATA['nucs']]

    with open('sigma.json') as f:
        sigma = json.load(f)

    # We don't use all nucs
    used_sigmas = set()
    for i in sigma:
        *_, nuc = i.rpartition('_')
        if nuc in DATA['nucs']:
            used_sigmas.add(i)

    return [[sigma[i] if i in used_sigmas else 0.0 for i in j] for j in sigma_symbols]

if __name__ == '__main__':
    system = CodeBlock(*list(map(gennuc, DATA['nucs'])))

    sigma_symbols = sorted([i.name for i in system.free_symbols if
        i.name.startswith('sigma')])

    sigma_map = sigma_symbol_to_indexed()
    nuc_map = nuc_symbol_to_indexed()
    # nuc_map = {}

    system = system.xreplace({**sigma_map, **nuc_map})

    sigma_array = generate_sigma_array()

    code = sympy.ccode(system)

    generated_code = TEMPLATE
    for val, repl in {
        "I": len(DATA['nucs']),
        "SIGMA_ARRAY": str(sigma_array).replace('[', '{').replace(']', '}'),
        "CODE": textwrap.indent(code, '    '),
        # For testing
        "N0": str([0.]*(len(DATA['nucs']) - 1) + [1.0]).replace('[',
            '{').replace(']', '}'),
    }.items():
        generated_code = generated_code.replace("%%" + val + "%%", str(repl))

    with open("sigma_array.txt", 'w') as f:
        f.write('[' + ',\n'.join(map(str, sigma_array)) + ']\n')

    with open('system.txt', 'w') as f:
        for eq in system.args:
            f.write(str(eq) + '\n')

    with open('system-C.txt', 'w') as f:
        f.write(code)

    with open('transmute.c', 'w') as f:
        f.write(generated_code)

    #system_cse = system.cse()

    #with open('system-cse.txt', 'w') as f:
    #    for eq in system_cse.args:
    #        f.write(str(eq) + '\n')

    #with open('system-cse-C.txt', 'w') as f:
    #    f.write(sympy.ccode(system_cse))
