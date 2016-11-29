from sympy import (nsolve, symbols, exp, simplify, chebyshevt_root, Tuple,
    diff, plot, Rational, together, Poly)

from sympy.utilities.decorator import conserve_mpmath_dps

t = symbols('t')

def nsolve_intervals(expr, bounds, division=30, **kwargs):
    """
    Divide bounds into division intervals and nsolve in each one
    """
    roots = []
    L = bounds[1] - bounds[0]
    for i in range(division):
        interval = [bounds[0] + i*L/division, bounds[0] + (i + 1)*L/division]
        try:
            root = nsolve(expr, interval, solver='bisect', **kwargs)
        except ValueError:
            continue
        else:
            roots.append(root)

    return roots

@conserve_mpmath_dps
def CRAM_exp2(loops=2):
    import mpmath
    prec = 128
    mpmath.mp.dps = prec

    epsilon = symbols("epsilon")
    p0, p1, p2, q1, q2 = symbols("p0, p1, p2, q1, q2")
    i = symbols("i")

    r = (p0 + p1*t + p2*t**2)/(1 + q1*t + q2*t**2)
    E = exp(-(-t - 1)/(2*t - 2)) - r
    expr = E + (-1)**i*epsilon
    expr = expr*(1 + q1*t + q2*t**2)
    expr = simplify(expr)

    points = [chebyshevt_root(7, 6 - j) for j in range(1, 7)]
    for iteration in range(loops):
        print("Iteration", iteration)
        system = Tuple(*[expr.subs({i: j, t: points[j]}) for j in range(5)])
        system = system + Tuple(expr.replace(exp, lambda i: 0).subs({i: 5, t: 1}))
        #print(system)
        sol = dict(zip([p0, p1, p2, q1, q2, epsilon], nsolve(system, [p0, p1, p2, q1, q2, epsilon], [1, 1, 1, 1, 1, 0])))
        D = diff(E.subs(sol), t)
        # plot(E.subs(sol), (t, -1, 1))
        # More 9's here means more accuracy, but we can't use 1 because of the singularity
        points = [-1, *nsolve_intervals(D, [-1, 0.99], maxsteps=300), 1]# mpf('0.9999999999999999999999999999999999999999999999999999999999999999999999999999')]
        print(points)
        Evals = [E.evalf(prec, subs={**sol, t: point}) for point in points[:-1]] + [-r.evalf(prec, subs={**sol, t: 1})]
        print(Evals)
        print('max - min', max(Evals) - min(Evals))
        print('epsilon', sol[epsilon])
        assert len(points) == 6

    print(sol)
    sol = {i: Rational(str(sol[i])) for i in sol}
    print(sol)
    n, d = together(r.subs(sol).subs(t, (2*t - 1)/(2*t + 1))).as_numer_denom() # simplify/cancel here will add degree to the numerator and denominator
    rat_func = (Poly(n)/Poly(d).TC())/(Poly(d)/Poly(d).TC())
    return rat_func.evalf(prec)

#D = CRAM_exp2()

rat_func = CRAM_exp2(8)
print(rat_func)
# plot(rat_func - exp(-t), (t, 0, 100))
