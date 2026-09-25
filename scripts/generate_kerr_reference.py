"""Generate independent reference values for the Kerr metric in Boyer-Lindquist coordinates.

Everything here is derived from the metric components alone, with SymPy: the inverse
metric by matrix inversion, the derivatives by symbolic differentiation, and the
Christoffel symbols from CLAUDE.md §2's definition. Nothing is taken from the TypeScript
implementation, so the fixture is a genuine independent check rather than a snapshot.

Run:  python3 scripts/generate_kerr_reference.py > tests/fixtures/kerr-reference.ts
"""

import sympy as sp

t, r, th, ph, M, a = sp.symbols('t r theta phi M a', real=True)

Sigma = r**2 + a**2 * sp.cos(th)**2
Delta = r**2 - 2*M*r + a**2

g = sp.zeros(4, 4)
g[0, 0] = -(1 - 2*M*r/Sigma)
g[0, 3] = g[3, 0] = -2*M*a*r*sp.sin(th)**2/Sigma
g[1, 1] = Sigma/Delta
g[2, 2] = Sigma
g[3, 3] = (r**2 + a**2 + 2*M*a**2*r*sp.sin(th)**2/Sigma)*sp.sin(th)**2

g_inv = g.inv()
x = [t, r, th, ph]

christoffel = [[[sum(g_inv[mu, s]*(sp.diff(g[s, al], x[be]) + sp.diff(g[s, be], x[al]) - sp.diff(g[al, be], x[s]))
                     for s in range(4))/2
                 for be in range(4)] for al in range(4)] for mu in range(4)]

d_g_inv = [[[sp.diff(g_inv[mu, nu], x[al]) for nu in range(4)] for mu in range(4)] for al in range(4)]

POINTS = [
    (1, sp.Rational(9, 10), sp.Rational(37, 10), sp.Rational(11, 10)),
    (1, sp.Rational(3, 10), 12, sp.Rational(2, 5)),
    (1, sp.Rational(998, 1000), sp.Rational(5, 2), sp.Rational(12, 5)),
    (1, 0, 8, sp.Rational(157, 100)),
    (2, sp.Rational(7, 5), 7, sp.Rational(9, 10)),
    (1, sp.Rational(-3, 5), 6, sp.Rational(3, 2)),
]

DIGITS = 25


def number(expr, sub):
    value = sp.N(expr.subs(sub), DIGITS)
    return sp.printing.str.sstr(value, full_prec=False)


def flat(matrix_like, sub):
    return ', '.join(number(matrix_like[i, j], sub) for i in range(4) for j in range(4))


print('/* GENERATED FILE — do not edit by hand.')
print(' *')
print(' * Reference values for the Kerr metric in Boyer-Lindquist coordinates, produced by')
print(' * scripts/generate_kerr_reference.py with SymPy. The inverse metric is obtained by')
print(' * matrix inversion and the Christoffel symbols from the definition in CLAUDE.md §2,')
print(' * so nothing here shares an algebraic route with the closed forms under test.')
print(' *')
print(f' * Values are exact rationals evaluated to {DIGITS} significant digits.')
print(' */')
print()
print('export interface KerrReferencePoint {')
print('  readonly mass: number;')
print('  readonly spin: number;')
print('  readonly r: number;')
print('  readonly theta: number;')
print('  /** g_mu_nu, index mu * 4 + nu. */')
print('  readonly metric: readonly number[];')
print('  /** g^{mu nu}, index mu * 4 + nu. */')
print('  readonly inverseMetric: readonly number[];')
print('  /** d_alpha g^{mu nu}, index alpha * 16 + mu * 4 + nu. */')
print('  readonly inverseMetricDerivatives: readonly number[];')
print('  /** Gamma^mu_{alpha beta}, index mu * 16 + alpha * 4 + beta. */')
print('  readonly christoffel: readonly number[];')
print('}')
print()
print('export const KERR_REFERENCE_POINTS: readonly KerrReferencePoint[] = [')
for (Mv, av, rv, thv) in POINTS:
    sub = {M: Mv, a: av, r: rv, th: thv}
    print('  {')
    print(f'    mass: {number(sp.Integer(1) * Mv, {})},')
    print(f'    spin: {number(sp.Integer(1) * av, {})},')
    print(f'    r: {number(sp.Integer(1) * rv, {})},')
    print(f'    theta: {number(sp.Integer(1) * thv, {})},')
    print(f'    metric: [{flat(g, sub)}],')
    print(f'    inverseMetric: [{flat(g_inv, sub)}],')
    derivs = ', '.join(number(d_g_inv[al][mu][nu], sub)
                       for al in range(4) for mu in range(4) for nu in range(4))
    print(f'    inverseMetricDerivatives: [{derivs}],')
    syms = ', '.join(number(christoffel[mu][al][be], sub)
                     for mu in range(4) for al in range(4) for be in range(4))
    print(f'    christoffel: [{syms}],')
    print('  },')
print('];')
