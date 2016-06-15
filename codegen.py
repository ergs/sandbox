from sympy import ccode, fcode, dsolve

def codegen(expr, lang, indent='    ', ):
    if lang == 'C':
        code = ccode
    elif lang == 'Fortran':
        code = fcode
    else:
        raise ValueError("Lang must be 'C' or 'Fortran'")

    try:
        sol = dsolve(expr)
    except ValueError:
        # Not an ODE
        return code(expr)

    return ccode(sol.rhs, assign_to=sol.lhs.func.__name__)

if __name__ == '__main__':
    from sympy import symbols, Function
    x = Function('x')
    a, t = symbols('a t')

    ode = x(t).diff(t) - a*x(t)
    print(codegen(ode, 'C'))
