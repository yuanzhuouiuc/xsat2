import collections
import z3
import math
import numpy as np
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
    return "_x_" + str(expr_z3.hash())

def _collect_vars(all_vars, expr):
    if z3.is_const(expr) and expr.decl().kind() == z3.Z3_OP_UNINTERPRETED:
        all_vars.add(expr)
    for child in expr.children():
        _collect_vars(all_vars, child)

def _to_z3_fp(val, sort):
    # Flatten numpy scalars/0-d arrays to Python numbers
    if isinstance(val, np.generic):
        val = val.item()
    elif isinstance(val, np.ndarray):
        val = val.astype(float).item()
    if isinstance(val, float):
        if math.isnan(val):
            return z3.fpNaN(sort)
        if math.isinf(val):
            return z3.fpPlusInfinity(sort) if val > 0 else z3.fpMinusInfinity(sort)
    if sort == z3.Float32():
        return z3.FPVal(np.float32(val).item(), sort)
    elif sort == z3.Float64():
        return z3.FPVal(float(val), sort)
    else:
        raise z3.Z3Exception(f"Unsupported FP sort: {sort}")

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
    num32, num64 = 0, 0
    for (s, val) in zip(symbolTable.items(), X_star):
        var, sort = s[0], s[1]
        var = rename_var(var)
        if var in var_dict:
            original_var = var_dict[var]
            if sort == Sort.Float32:
                # val_z3 = z3.FPVal(val, z3.Float32())
                val_z3 = _to_z3_fp(val, z3.Float32())
                num32 += 1
            elif sort == Sort.Float64:
                # val_z3 = z3.FPVal(val, z3.Float64())
                val_z3 = _to_z3_fp(val, z3.Float64())
                num64 += 1
            else:
                raise NotImplementedError("Unexpected type %s" % sort)
            model.append((original_var, val_z3))
        else:
            print(f"Warning: Variable {var} not found in expression")
    if printModel:
        print("model: " + str(model))
    print("num32: " + str(num32))
    print("num64: " + str(num64))

    return z3_util.is_true(z3.simplify(z3.substitute(ez, *model)))

def z3_verify(ez, X_star, symbolTable, printModel = False):
    # --- 1. Basic assertions and variable mapping (same as verify_solution) ---
    assert isinstance(symbolTable, collections.OrderedDict)
    assert isinstance(ez, z3.ExprRef)
    assert len(symbolTable) == X_star.size
    all_vars = set()
    _collect_vars(all_vars, ez)
    var_dict = {}
    for var in all_vars:
        var_name = rename_var(var.decl().name())
        var_dict[var_name] = var
    # --- 2. Use a Z3 Solver for verification ---
    s = z3.Solver()
    # Add the original problem constraints to the solver
    s.add(ez)
    if printModel:
        print("Model being checked by z3_verify:")
    # Iterate through the proposed solution and add equality constraints
    for (symbol, val) in zip(symbolTable.items(), X_star):
        var_name, sort = symbol[0], symbol[1]
        var_name = rename_var(var_name)
        if var_name in var_dict:
            original_var = var_dict[var_name]
            # Convert the Python/Numpy value to a Z3 value
            if sort == Sort.Float32:
                # val_z3 = z3.FPVal(val, z3.Float32())
                val_z3 = _to_z3_fp(val, z3.Float32())
            elif sort == Sort.Float64:
                # val_z3 = z3.FPVal(val, z3.Float64())
                val_z3 = _to_z3_fp(val, z3.Float64())
            else:
                # This can be extended to handle Sort.Int, Sort.Real, etc.
                raise NotImplementedError(f"Unexpected type {sort}")
            # Add the constraint: variable == proposed_value
            constraint = (original_var == val_z3)
            s.add(constraint)
            if printModel:
                print(f"  {constraint}")
        else:
            print(f"Warning: Variable {var_name} not found in expression")
    res = s.check() == z3.sat
    return s.check() == z3.sat