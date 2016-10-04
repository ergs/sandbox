``` python

In [35]: @numba.jit(nopython=True)
def s(d):
    N = 0
    for i in range(1, d):
        N += 1/i**2
    return N
   ....:

In [54]: @numba.jit(nopython=True)
def s2(d):
    N = 0
    for i in range(d-1, 0, -1):
        N += 1/i**2
    return N
   ....:

In [58]: (pi**2/6.).evalf()
Out[58]: 1.64493406684823

In [59]: s2(3000000000)
Out[59]: 1.644934066514893

In [60]: s(3000000000)
Out[60]: 1.644934057834575

```

```

In [1]: var("i")
Out[1]: i

In [2]: summation(1/i**2, (i, 1, n))
Out[2]: harmonic(n, 2)

In [3]: harmonic(3000000000, 2).evalf()
<hangs>
```
