import z3

def is_fpNeg(a):
    return a.decl().kind() == z3.Z3_OP_FPA_NEG

def is_fpDiv(a):
    return a.decl().kind() == z3.Z3_OP_FPA_DIV

def is_fpMul(a):
    return a.decl().kind() == z3.Z3_OP_FPA_MUL

def is_RNE(a):
    return a.decl().kind() == z3.Z3_OP_FPA_RM_NEAREST_TIES_TO_EVEN

def is_fpAdd(a):
    return a.decl().kind() == z3.Z3_OP_FPA_ADD

def is_lt(a):
    return a.decl().kind() == z3.Z3_OP_LT or a.decl().kind() == z3.Z3_OP_FPA_LT

def is_le(a):
    return a.decl().kind() == z3.Z3_OP_LE or a.decl().kind() == z3.Z3_OP_FPA_LE

def is_ge(a):
    return a.decl().kind() == z3.Z3_OP_GE or a.decl().kind() == z3.Z3_OP_FPA_GE

def is_eq(a):
    return a.decl().kind() == z3.Z3_OP_EQ or a.decl().kind() == z3.Z3_OP_FPA_EQ

def is_distinct(a):
    return a.decl().kind()==z3.Z3_OP_DISTINCT # no FPA_DISTINCT

def is_gt(a):
    return a.decl().kind()==z3.Z3_OP_GT or a.decl().kind()==z3.Z3_OP_FPA_GT

def is_true(a):
    return a.decl().kind()==z3.Z3_OP_TRUE

def is_false(a):
    return a.decl().kind()==z3.Z3_OP_FALSE

def is_fpFP(a):
    return a.decl().kind() == z3.Z3_OP_FPA_TO_FP

def is_variable(a):
    return z3.is_const(a) and a.decl().kind() == z3.Z3_OP_UNINTERPRETED

def is_value(a):
    return z3.is_const(a) and a.decl().kind() != z3.Z3_OP_UNINTERPRETED