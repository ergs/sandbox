from itertools import repeat
from sympy import Symbol, symbols, Add, Sum, Mul, simplify, Pow, exp

def variance_prop(expr, consts=()):
    args = expr.args
    if len(args) == 0:
        if isinstance(expr, Symbol) and expr not in consts:
            signame = 'var_' + expr.name
            return symbols(signame)
        else:
            return 0
    var_args = list(map(variance_prop, args, repeat(consts)))
    if isinstance(expr, Add):
        return Add(*var_args)
    elif isinstance(expr, Mul):
        terms = [v/a**2 for a, v in zip(args, var_args)]
        return simplify(expr**2 * Add(*terms))
    elif isinstance(expr, Pow):
        b = args[1]
        v = var_args[0] * (expr * b / args[0])**2
        return simplify(v)
    elif isinstance(expr, exp):
        return simplify(var_args[0] * expr**2)
    else:
        raise RuntimeError("unknown operator")


if __name__ == '__main__':
    x, y, z = symbols('x y z')
    phi, t = consts = symbols('phi t')
    cases = [x + y, x + y + z, 2*x, x*y, 1/x, x/y, exp(x),
             exp(2*x), exp(-x*t)]
    for case in cases:
        print(case, "=>")
        print(variance_prop(case, consts=consts))
        print('~'*10)

