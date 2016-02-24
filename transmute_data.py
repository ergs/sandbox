import json
import warnings
import numpy as np
warnings.simplefilter('ignore')

from pyne.utils import toggle_warnings
toggle_warnings()
from pyne import nucname
from pyne import data
from pyne import rxname

DECAY_RXS = ['bminus', 'bplus', 'ec', 'alpha', 'it', 'sf', 'bminus_n']
XS_RXS = ['gamma', 'z_2n', 'z_3n', 'alpha', 'fission', 'proton', 'gamma_1', 
          'z_2n_1']


def add_child_decays(nuc, symbols):
    childname = nucname.name(nuc)
    for rx in DECAY_RXS:
        try:
            parent = rxname.parent(nuc, rx, b'decay')
        except RuntimeError:
            continue
        if data.branch_ratio(parent, nuc) < 1e-16:
            continue
        parname = nucname.name(parent)
        symbols['lambda_' + parname] = data.decay_const(parent)
        gamma = 'gamma_{0}_{1}_{2}'.format(parname, childname, rx)
        symbols[gamma] = data.branch_ratio(parent, nuc)


def add_child_xss(nuc, channels, parents):
    childname = nucname.name(nuc)
    channels[childname] = rxs = {}
    for rx in XS_RXS:
        try:
            parent = rxname.parent(nuc, rx)
        except RuntimeError:
            continue
        if parent not in parents:
            continue
        parname = nucname.name(parent)
        rxs[rx] = parname


def main():
    # get list of nuclides
    data.atomic_mass('U235')
    nucs = set(data.atomic_mass_map.keys())
    for nuc in data.atomic_mass_map:
        nucm = nuc + 1
        if nucname.anum(nuc) == 0 or data.decay_const(nucm) < 1e-16 or \
                            data.decay_const(nuc) == data.decay_const(nucm):
            continue
        nucs.add(nucm)
    nucs = [nuc for nuc in nucs if nucname.anum(nuc) > 0 and
                                not np.isnan(data.decay_const(nuc)) and
                                nuc < 200000000]
    nucs.sort()
    # get symbols
    symbols = {}
    channels = {}
    for nuc in nucs:
        add_child_decays(nuc, symbols)
        add_child_xss(nuc, channels, nucs)
    # print symbols
    d = {'symbols': symbols, 'nucs': list(map(nucname.name, nucs)), 
         'channels': channels}
    s = json.dumps(d, indent=4, sort_keys=True)
    print(s)


if __name__ == '__main__':
    main()
