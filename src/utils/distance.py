import sympy

def dist_le(lhs, rhs):
    assert (isinstance(lhs, sympy.Expr) and isinstance(rhs, sympy.Expr))
    return sympy.Piecewise( (lhs - rhs,lhs > rhs), (0, True))

def dist_lt(lhs, rhs):
    assert (isinstance(lhs, sympy.Expr) and isinstance(rhs, sympy.Expr))
    return dist_le(lhs + _theta(), rhs)

def dist_ge(lhs, rhs):
    assert (isinstance(lhs, sympy.Expr) and isinstance(rhs, sympy.Expr))
    return dist_le(rhs, lhs)

def dist_gt(lhs, rhs):
    assert (isinstance(lhs, sympy.Expr) and isinstance(rhs, sympy.Expr))
    return dist_lt(rhs, lhs)

def dist_eq(lhs, rhs):
    assert (isinstance(lhs, sympy.Expr) and isinstance(rhs, sympy.Expr))
    return sympy.Abs(lhs - rhs)

def dist_distinct(lhs, rhs):
    assert (isinstance(lhs, sympy.Expr) and isinstance(rhs, sympy.Expr))
    return sympy.Piecewise( (_theta(), sympy.Eq(lhs, rhs)), (0, True))

# self-implemented theta(), minimum subnormal positive number(IEEE754 double format)
def _theta():
    return sympy.Float(2)**(-1074)

def _theta_single():
    return sympy.Float(2)**(-149)