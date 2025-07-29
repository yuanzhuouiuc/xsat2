import collections
import z3
from src.utils.sort import Sort
import src.utils.z3_util as z3_util

def replace_all(text, dic):
    for i, j in dic.items():
        text = text.replace(i, j)
    return text

def rename_var(var):
    reps = {':':'_', '@':'_', '|':'_', "#":'_', "!":"_"}
    return replace_all(var, reps)

def var_name(expr_z3):
    return "_t_" + str(expr_z3.get_id())

def _getSort(expr_z3):
    assert isinstance(expr_z3, z3.ExprRef)
    if expr_z3.sort() == z3.Float32():
        return Sort.Float32
    if expr_z3.sort() == z3.Float64():
        return Sort.Float64
    if expr_z3.sort() == z3.RealSort():
        return Sort.Real
    if expr_z3.sort() == z3.IntSort():
        return Sort.Int
    return Sort.UNKNOWN

def var_hash(expr_z3):
    return "_x_" + expr_z3.hash()

def _collect_vars(all_vars, expr):
    if z3.is_const(expr) and expr.decl().kind() == z3.Z3_OP_UNINTERPRETED:
        all_vars.add(expr)
    for child in expr.children():
        _collect_vars(all_vars, child)

def verify_solution(ez, X_star, symbolTable, printModel = False):
    assert isinstance(symbolTable, collections.OrderedDict)
    assert isinstance(ez, z3.ExprRef)
    assert len(symbolTable) == X_star.size
    model = []
    all_vars = set()
    _collect_vars(all_vars, ez)
    var_dict = {}
    for var in all_vars:
        var_name = rename_var(var.decl().name())
        var_dict[var_name] = var
    for (s, val) in zip(symbolTable.items(), X_star):
        var, sort = s[0], s[1]
        var = rename_var(var)
        if var in var_dict:
            original_var = var_dict[var]
            if sort == Sort.Float32:
                val_z3 = z3.FPVal(val, z3.Float32())
            elif sort == Sort.Float64:
                val_z3 = z3.FPVal(val, z3.Float64())
            else:
                raise NotImplementedError("Unexpected type %s" % sort)
            model.append((original_var, val_z3))
        else:
            print(f"Warning: Variable {var} not found in expression")
    if printModel:
        print("model: " + str(model))

    return z3_util.is_true(z3.simplify(z3.substitute(ez, *model)))