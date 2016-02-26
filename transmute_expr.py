import textwrap
import re
import json
import itertools
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

with open('sigma.json') as f:
    SIGMA = json.load(f)



decay_rxs = ['bminus', 'bplus', 'ec', 'alpha', 'it', 'sf', 'bminus_n']
xs_rxs = ['gamma', 'z_2n', 'z_3n', 'alpha', 'fission', 'proton', 'gamma_1', 'z_2n_1']

gamma_base = '^gamma_([A-Z][a-z]?\d+)_'

# Create from -> to nuclide mapping
FROM_TO = {}
def add_from_to(f, t, k):
    if f not in FROM_TO:
        FROM_TO[f] = {}
    to_nucs = FROM_TO[f]
    if t not in to_nucs:
        to_nucs[t] = set()
    rxs = to_nucs[t]
    rxs.add(k)

for key in DATA['symbols'].keys():
    if not key.startswith('gamma_'):
        continue
    _, f, t, *_ = key.split('_')
    add_from_to(f, t, key)

for sig, (val, f, t) in SIGMA.items():
    if t is None or val < 1e-200:
        continue
    add_from_to(f, t, sig)


def make_chains(f, curr=()):
    if len(curr) == 0:
        curr = (f,)
    chains = [curr]
    if len(curr) > 21:
        return chains
    for t in FROM_TO.get(f, ()):
        if t in curr:
            continue  # get rid of cycles
        tchain = curr + (t,)
        newchains = make_chains(t, tchain)
        chains.extend(newchains)
    return chains

CHAINS = set()
for nuc in DATA['nucs'][:346]:
    #print(nuc)
    CHAINS.update(make_chains(nuc))
CHAINS = sorted(CHAINS, key=lambda c: c[-1])

t = sympy.symbols('t')
# G = 1
# phi = sympy.MatrixSymbol('phi', G, 1)
phi = sympy.symbols('phi')

def decay_const(nuc):
    return DATA['symbols']['lambda_' + nuc]


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


def sigma_a(nuc):
    sigma_a_name = 'sigma_a_{0}'.format(nuc)
    sig_a = SIGMA.get(sigma_a_name, [0])[0]
    if sig_a > 0:
        sig_a = sympy.symbols(sigma_a_name)
    return sig_a


def genexponent(nuc):
    lambda_1 = decay_const(nuc)
    sig_a = sigma_a(nuc)
    try:
        return sympy.exp(-(lambda_1 + sig_a*phi)*t)
    except:
        import pdb; pdb.set_trace()


def gentotalbranch(chain):
    terms = []
    for f, t in zip(chain[:-1], chain[1:]):
        lambda_i = decay_const(f)
        gamma_rx = gamma(f, t)
        sig_rx = sigma_rx(f, t)
        term = (gamma_rx * lambda_i) + (sig_rx * phi)
        terms.append(term)
    return sympy.Mul(*terms)


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


def genciexp(chain):
    terms = []
    for nuc in chain:
        ci = genci(nuc, chain)
        exp = genexponent(nuc)
        term = ci * exp
        terms.append(term)
    return sympy.Add(*terms)


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
        #terms.append(genchainexpr(chain))
    print(NUM, nuc)
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

if __name__ == '__main__':
    nucs = DATA['nucs'][:346]
    system = CodeBlock(*list(map(gennuc, nucs)))

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
        "I": len(nucs),
        "SIGMA_ARRAY": str(sigma_array).replace('[', '{').replace(']', '}'),
        "CODE": textwrap.indent(code, '    '),
        # For testing
        "N0": str([0.]*(len(nucs) - 1) + [1.0]).replace('[',
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
