import collections
import z3
from sort import Sort

def _is_true(a):
    return a.decl().kind()==z3.Z3_OP_TRUE

def replace_all(text, dic):
    for i, j in dic.iteritems():
        text = text.replace(i, j)
    return text

def rename_var(var):
    reps = {':':'_', '@':'_', '|':'_',"#":'_',"!":"_"}
    return replace_all(var,reps)

def _getSort(expr_z3):
    assert isinstance(expr_z3, z3.ExprRef)
    if expr_z3.sort()==z3.Float32():
        return Sort.Float32
    if expr_z3.sort()==z3.Float64():
        return Sort.Float64
    if expr_z3.sort()==z3.RealSort():
        return Sort.Real
    if expr_z3.sort()==z3.IntSort():
        return Sort.Int
    return Sort.UNKNOWN

def var_hash(expr_z3):
    return "_x_"+expr_z3.hash()

def verify_solution(ez, X_star, symbolTable, printModel=False):
    assert isinstance(symbolTable, collections.OrderedDict)
    assert isinstance(ez, z3.ExprRef)
    assert len(symbolTable) == X_star.size
    model = []
    # (sympy.Symbol('x'),
    for (s, val) in zip(symbolTable.items(), X_star):
        var, sort = s[0], s[1]
        var = rename_var(var)
        if sort == Sort.Float32:
            var_z3 = z3.FP(str(var), z3.Float32())
            val_z3 = z3.FPVal(val, z3.Float32())
        elif sort == Sort.Float64:
            var_z3 = z3.FP(str(var), z3.Float64())
            val_z3 = z3.FPVal(val, z3.Float64())
        else:
            raise NotImplementedError("Unexpected type %s" % sort)
        model.append((var_z3, val_z3))
    if printModel:
        print("model: " + str(model))

    return _is_true(z3.simplify(z3.substitute(ez, *model)))