import textwrap
import json
from functools import lru_cache
from collections import defaultdict

import networkx

import sympy
from sympy.codegen.ast import Assignment, CodeBlock

TEMPLATE = r"""
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define I (%%I%%)
#define R (9)
#define T (%%T%%)
#define PHI (%%PHI%%)

double* transmute(double* N0, double t, double phi, double sigma[I][R]);

double* transmute(double* N0, double t, double phi, double sigma[I][R]) {

    double* N1 = malloc(I*sizeof(double));

%%CODE%%

    return(N1);
}

int main() {
    int i;
    double* N1;
    double N0[I] = %%N0%%;

    double sigma[I][R] = %%SIGMA_ARRAY%%;

    N1 = transmute(N0, T, PHI, sigma);

    for (i=0; i < I; i++) {
        printf("%d %e\n", i, N1[i]);
    }
    return(0);
}

"""

def load_data():
    global DATA
    global SIGMA

    with open('transmute_data.json') as f:
        DATA = json.load(f)

    with open('sigma.json') as f:
        SIGMA = json.load(f)

load_data()

decay_rxs = ['bminus', 'bplus', 'ec', 'alpha', 'it', 'sf', 'bminus_n']
xs_rxs = ['gamma', 'z_2n', 'z_3n', 'alpha', 'fission', 'proton', 'gamma_1', 'z_2n_1']

gamma_base = '^gamma_([A-Z][a-z]?\d+)_'

# Create from -> to nuclide mapping
FROM_TO = defaultdict(lambda: defaultdict(set))
CHAIN_GRAPH = set()

def create_from_to():
    for key in DATA['symbols'].keys():
        if not key.startswith('gamma_'):
            continue
        _, f, t, *_ = key.split('_')
        FROM_TO[f][t].add(key)

    for sig, (val, f, t) in SIGMA.items():
        if t is None or val < 1e-200:
            continue
        FROM_TO[f][t].add(sig)

create_from_to()

def create_chains():
    global CHAINS

    CHAINS = set()
    cutoff_nucs = DATA['nucs']

    for nuc in cutoff_nucs:
        #print(nuc)
        CHAIN_GRAPH.add(('start', nuc))
        CHAIN_GRAPH.add((nuc, 'end'))
        for t in FROM_TO[nuc]:
            CHAIN_GRAPH.add((nuc, t))

    G = networkx.DiGraph(list(CHAIN_GRAPH))
    CHAINS = sorted(networkx.all_simple_paths(G, 'start', 'end', 22), key=lambda c: list(reversed(c)))

    CHAINS = [i[1:][:-1] for i in CHAINS]

create_chains()

print(len(CHAINS))
print(max(CHAINS, key=len))
print(min(CHAINS, key=len))

t = sympy.symbols('t')
# G = 1
# phi = sympy.MatrixSymbol('phi', G, 1)
phi = sympy.symbols('phi')

@lru_cache(1024)
def decay_const(nuc):
    return DATA['symbols']['lambda_' + nuc]

@lru_cache(1024)
def gamma(f, t):
    prefix = 'gamma_{0}_{1}_'.format(f, t)
    possible = FROM_TO[f][t]
    for p in possible:
        if p.startswith(prefix):
            rx = p
            break
    else:
        return 0
    return DATA['symbols'].get(rx, 0)


@lru_cache(1024)
def sigma_rx(f, t):
    possible = FROM_TO[f][t]
    for p in possible:
        if p.startswith('sigma_') and p.endswith(f):
            rx = p
            break
    else:
        return 0
    sigma = SIGMA.get(rx, [0])[0]
    if sigma > 0:
        sigma = sympy.symbols(rx)
    return sigma


@lru_cache(1024)
def sigma_a(nuc):
    sigma_a_name = 'sigma_a_{0}'.format(nuc)
    sig_a = SIGMA.get(sigma_a_name, [0])[0]
    if sig_a > 0:
        sig_a = sympy.symbols(sigma_a_name)
    return sig_a


@lru_cache(1024)
def genexponent(nuc):
    lambda_1 = decay_const(nuc)
    sig_a = sigma_a(nuc)
    try:
        return sympy.exp(-(lambda_1 + sig_a*phi)*t)
    except:
        import pdb; pdb.set_trace()


@lru_cache(1024)
def gentotalbranch(chain):
    terms = []
    for f, t in zip(chain[:-1], chain[1:]):
        lambda_i = decay_const(f)
        gamma_rx = gamma(f, t)
        sig_rx = sigma_rx(f, t)
        term = (gamma_rx * lambda_i) + (sig_rx * phi)
        terms.append(term)
    return sympy.Mul(*terms)


@lru_cache(1024)
def genci(nuc, chain):
    terms = []
    lambda_i = decay_const(nuc)
    sig_a_i = sigma_a(nuc)
    part_i = lambda_i + (sig_a_i * phi)
    for j in chain:
        if j == nuc:
            continue
        lambda_j = decay_const(j)
        sig_a_j = sigma_a(j)
        part_j = lambda_j + sig_a_j * phi
        term = 1 / (part_j - part_i)
        terms.append(term)
    return sympy.Mul(*terms)


@lru_cache(1024)
def genciexp(chain):
    terms = []
    for nuc in chain:
        ci = genci(nuc, chain)
        exp = genexponent(nuc)
        term = ci * exp
        terms.append(term)
    return sympy.Add(*terms)


@lru_cache(1024)
def genchainexpr(chain):
    nuc0 = sympy.symbols('{0}_0'.format(chain[0]))
    if len(chain) == 1:
        return nuc0 * genexponent(chain[0])
    tb = gentotalbranch(chain)
    ce = genciexp(chain)
    return nuc0 * tb * ce


def gennuc(nuc):
    nuc1 = sympy.symbols('{0}_1'.format(nuc))
    terms = []
    NUM = 0
    for chain in CHAINS:
        if chain[-1] != nuc:
            continue
        NUM += 1
        terms.append(genchainexpr(chain))
    #print(NUM, nuc)
    rhs = sympy.Add(*terms)
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

    # We don't use all nucs
    used_sigmas = set()
    for i in SIGMA:
        *_, nuc = i.rpartition('_')
        if nuc in DATA['nucs']:
            used_sigmas.add(i)

    return [[SIGMA[i][0] if i in used_sigmas else 0.0 for i in j] for j in sigma_symbols]

def main():
    NUCS = DATA['nucs']
    nucs = ['K40']
    system = CodeBlock(*list(map(gennuc, nucs)))

    sigma_map = sigma_symbol_to_indexed()
    nuc_map = nuc_symbol_to_indexed()
    # nuc_map = {}

    system = system.xreplace({**sigma_map, **nuc_map})

    sigma_array = generate_sigma_array()

    code = sympy.ccode(system, order='none')

    generated_code = TEMPLATE

    input_data = [0.0]*len(sigma_array)

    input_data[NUCS.index("K39")] = 1.0
    input_time = 81.0

    for val, repl in {
        "I": len(sigma_array),
        "T": input_time,
        "PHI": 4e-10,
        "SIGMA_ARRAY": str(sigma_array).replace('[', '{').replace(']', '}'),
        "CODE": textwrap.indent(code, '    '),
        # For testing
        "N0": str(input_data).replace('[', '{').replace(']', '}'),
    }.items():
        generated_code = generated_code.replace("%%" + val + "%%", str(repl))

    with open("sigma_array.txt", 'w') as f:
        f.write('[' + ',\n'.join(map(str, sigma_array)) + ']\n')

    # with open('system.txt', 'w') as f:
    #     for eq in system.args:
    #         f.write(str(eq) + '\n')

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

if __name__ == '__main__':
    main()
