import sympy
import z3
import numpy as np
import scipy.optimize as op
import argparse
import sys, os
import time
import collections
import subprocess
import multiprocessing as mp
import warnings
import struct
import pickle
import src.utils.z3_util as z3_util
import src.utils.verification as verification
from src.utils.sort import Sort

DEBUG=False

def _print_xsatInfo():
    try:
        logo = open ('logo.txt',"r").read().strip('\n')
        print(logo)
    except:
        pass
    print()
    print("*"*50)
    print("XSat Version 04/04/2016")
    print ("Contributors: Zhoulai Fu and Zhendong Su")
    print("*"*50)

def _get_template():
    """Template that handles both float32 and float64 variables natively"""
    template = """#include <Python.h>
#include "xsat.h"

static PyObject* R(PyObject* self, PyObject *args){
  // Mixed type declarations
  %(var_declarations)s
  if (!PyArg_ParseTuple(args,"%(parse_formats)s", %(var_refs)s))
    return NULL;
  %(x_body)s
  return Py_BuildValue("d",%(x_expr)s);  // Always return double for distance
}

static PyMethodDef methods[] = {
  {"R", R, METH_VARARGS, NULL},
  {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "foo",
  NULL,
  -1,
  methods,
  NULL, NULL, NULL, NULL
};

PyMODINIT_FUNC
PyInit_foo(void)
{
  PyObject* module = PyModule_Create(&moduledef);
  if (module == NULL)
    return NULL;
  PyModule_AddIntConstant(module, "dim", %(x_dim)s);
  return module;
}
"""
    return template

def _get_operand_type(expr, symbolTable, cache):
    """Determine if an operand is float32 or float64"""
    if z3_util.is_variable(expr):
        var_name = verification.rename_var(expr.decl().name())
        if var_name in symbolTable:
            return symbolTable[var_name]
    elif z3_util.is_value(expr):
        if expr.sort() == z3.Float32():
            return Sort.Float32
        elif expr.sort() == z3.Float64():
            return Sort.Float64
    elif expr.get_id() in cache:
        # Check the type of intermediate result
        if expr.sort() == z3.Float32():
            return Sort.Float32
        elif expr.sort() == z3.Float64():
            return Sort.Float64
    return Sort.Float64  # Default to float64

def _get_comparison_function(base_func, lhs_type, rhs_type):
    """Get the appropriate comparison function based on operand types"""
    if lhs_type == Sort.Float32 and rhs_type == Sort.Float32:
        return f"{base_func}_f32"
    elif lhs_type == Sort.Float32 and rhs_type == Sort.Float64:
        return f"{base_func}_mixed_fd"
    elif lhs_type == Sort.Float64 and rhs_type == Sort.Float32:
        return f"{base_func}_mixed_df"
    else:
        return base_func  # Both float64 or default

def _gen(expr_z3, symbolTable, cache, result):
    ###Leaf: var
    if z3_util.is_variable(expr_z3):
        if DEBUG:
            print("-- Branch _is_variable with ", expr_z3)
        symVar = expr_z3.decl().name()
        symVar = verification.rename_var(symVar)
        if z3.is_int(expr_z3):
            symType = Sort.Int
        elif z3.is_fp(expr_z3):
            if expr_z3.sort() == z3.Float64():
                symType = Sort.Float64
            elif expr_z3.sort() == z3.Float32():
                symType = Sort.Float32
            else:
                raise NotImplementedError("Unexpected sort.", expr_z3.sort())
        elif z3.is_real(expr_z3):
            symType = Sort.Float
            warnings.warn("****WARNING****: Real variable '%s' treated as floating point" % symVar)
        else:
            raise NotImplementedError("Unexpected type")
        if (symVar in symbolTable.keys()):
            assert symType == symbolTable[symVar]
        else:
            symbolTable[symVar] = symType
        return symVar
    ###Leaf: val
    if z3_util.is_value(expr_z3):
        if DEBUG:
            print("-- Branch _is_value")
        if z3.is_fp(expr_z3) or z3.is_real(expr_z3):
            if DEBUG:
                print("---- Sub-Branch FP or Real")
            if isinstance(expr_z3, z3.FPNumRef):
                if DEBUG:
                    print("------- Sub-Sub-Branch _is_FPNumRef")
                # Check for special values first
                if expr_z3.isNaN():
                    str_ret = "NAN"
                elif expr_z3.isInf() and expr_z3.decl().kind() == z3.Z3_OP_FPA_PLUS_INF:
                    str_ret = "INFINITY"
                elif expr_z3.isInf() and expr_z3.decl().kind() == z3.Z3_OP_FPA_MINUS_INF:
                    str_ret = "- INFINITY"
                else:
                    # Handle normal values
                    try:
                        str_ret = str(sympy.Float(str(expr_z3), 17))
                    except ValueError:
                        # Handle other edge cases
                        offset = 127 if expr_z3.sort() == z3.Float32() else 1023
                        # Z3 new version needs the offset to be taken into consideration
                        expr_z3_exponent = expr_z3.exponent_as_long() - offset
                        str_ret = str(sympy.Float(
                            (-1) ** float(expr_z3.sign()) * float(str(expr_z3.significand())) * 2 ** (expr_z3_exponent),
                            17))
            else:
                if DEBUG:
                    print("------- Sub-Sub-Branch other than FPNumRef, probably FPRef")
                str_ret = str(sympy.Float(str((expr_z3)), 17))
        elif z3.is_int(expr_z3):
            if DEBUG:
                print("---- Sub-Branch Integer")
            str_ret = str(sympy.Integer(str(expr_z3)))
        elif z3_util.is_true(expr_z3):
            str_ret = "0"
        elif z3_util.is_false(expr_z3):
            str_ret = "1"
        else:
            raise NotImplementedError("[XSat: Coral Benchmarking] type not considered ")
        if expr_z3.sort() == z3.Float32():
            str_ret = str_ret + "f"
        return str_ret

    # if (expr_z3 in cache): return cache[expr_z3]

    # cache will be a set of defined IDs
    # if (verification.var_name(expr_z3) in cache): return cache[expr_z3]

    if (expr_z3.get_id() in cache): return verification.var_name(expr_z3)

    cache.add(expr_z3.get_id())
    # cache[expr_z3]=verification.var_name(expr_z3)

    sort_z3 = expr_z3.decl().kind()

    expr_type = 'double'
    if expr_z3.sort() == z3.FPSort(8, 24): expr_type = 'float'
    ###
    if sort_z3 == z3.Z3_OP_FPA_LE:
        if DEBUG:
            print("-- Branch _is_le")
        lhs = _gen(expr_z3.arg(0), symbolTable, cache, result)
        rhs = _gen(expr_z3.arg(1), symbolTable, cache, result)
        lhs_type = _get_operand_type(expr_z3.arg(0), symbolTable, cache)
        rhs_type = _get_operand_type(expr_z3.arg(1), symbolTable, cache)
        func_name = _get_comparison_function("DLE", lhs_type, rhs_type)
        toAppend = "double %s = %s(%s,%s);" % (
            verification.var_name(expr_z3), func_name, lhs, rhs)
        result.append(toAppend)
        return verification.var_name(expr_z3)

    #########!!!!!!!!!!!! need to do something
    if sort_z3 == z3.Z3_OP_FPA_TO_FP:
        if DEBUG:
            print("-- Branch _is_fpFP")
        assert expr_z3.num_args() == 2
        if not (z3_util.is_RNE(expr_z3.arg(0))):
            warnings.warn("WARNING!!! I expect the first argument of fpFP is RNE, but it is ", expr_z3.arg(0))
        x = _gen(expr_z3.arg(1), symbolTable, cache, result)
        if expr_z3.sort() == z3.FPSort(8, 24):
            toAppend = "float %s = (float)(%s);" % (verification.var_name(expr_z3), x)
        else:
            toAppend = "double %s = (double)(%s);" % (verification.var_name(expr_z3), x)
        result.append(toAppend)
        return verification.var_name(expr_z3)

    if sort_z3 == z3.Z3_OP_FPA_LT:
        if DEBUG:
            print("-- Branch _is_lt")
        lhs = _gen(expr_z3.arg(0), symbolTable, cache, result)
        rhs = _gen(expr_z3.arg(1), symbolTable, cache, result)
        lhs_type = _get_operand_type(expr_z3.arg(0), symbolTable, cache)
        rhs_type = _get_operand_type(expr_z3.arg(1), symbolTable, cache)
        func_name = _get_comparison_function("DLT", lhs_type, rhs_type)
        toAppend = "double %s = %s(%s,%s);" % (
            verification.var_name(expr_z3), func_name, lhs, rhs)
        result.append(toAppend)
        return verification.var_name(expr_z3)

    if z3_util.is_eq(expr_z3):
        if DEBUG:
            print("-- Branch _is_eq")
        lhs = _gen(expr_z3.arg(0), symbolTable, cache, result)
        rhs = _gen(expr_z3.arg(1), symbolTable, cache, result)
        lhs_type = _get_operand_type(expr_z3.arg(0), symbolTable, cache)
        rhs_type = _get_operand_type(expr_z3.arg(1), symbolTable, cache)
        func_name = _get_comparison_function("DEQ", lhs_type, rhs_type)
        toAppend = "double %s = %s(%s,%s);" % (
            verification.var_name(expr_z3), func_name, lhs, rhs)
        result.append(toAppend)
        return verification.var_name(expr_z3)

    if z3_util.is_fpMul(expr_z3):
        if DEBUG:
            print("-- Branch _is_fpMul")
        if not z3_util.is_RNE(expr_z3.arg(0)):
            warnings.warn("WARNING!!! arg(0) is not RNE but is treated as RNE. arg(0) = ", expr_z3.arg(0))
        assert expr_z3.num_args() == 3
        lhs = _gen(expr_z3.arg(1), symbolTable, cache, result)
        rhs = _gen(expr_z3.arg(2), symbolTable, cache, result)
        if expr_type == 'float':
            toAppend = "float %s = (float)(%s) * (float)(%s);" % (verification.var_name(expr_z3), lhs, rhs)
        else:
            toAppend = "double %s = %s * %s;" % (verification.var_name(expr_z3), lhs, rhs)
        result.append(toAppend)
        return verification.var_name(expr_z3)

    if z3_util.is_fpDiv(expr_z3):
        if DEBUG:
            print("-- Branch _is_fpDiv")
        if not z3_util.is_RNE(expr_z3.arg(0)):
            warnings.warn("WARNING!!! arg(0) is not RNE but is treated as RNE. arg(0) = ", expr_z3.arg(0))
        assert expr_z3.num_args() == 3
        lhs = _gen(expr_z3.arg(1), symbolTable, cache, result)
        rhs = _gen(expr_z3.arg(2), symbolTable, cache, result)
        if expr_type == 'float':
            toAppend = "float %s = (float)(%s) / (float)(%s);" % (verification.var_name(expr_z3), lhs, rhs)
        else:
            toAppend = "double %s = %s / %s;" % (verification.var_name(expr_z3), lhs, rhs)
        result.append(toAppend)
        return verification.var_name(expr_z3)

    if z3_util.is_fpAdd(expr_z3):
        if DEBUG:
            print("-- Branch _is_fpAdd")
        if not z3_util.is_RNE(expr_z3.arg(0)):
            warnings.warn("WARNING!!! arg(0) is not RNE but is treated as RNE. arg(0) = ", expr_z3.arg(0))
        assert expr_z3.num_args() == 3
        lhs = _gen(expr_z3.arg(1), symbolTable, cache, result)
        rhs = _gen(expr_z3.arg(2), symbolTable, cache, result)
        if expr_type == 'float':
            toAppend = "float %s = (float)(%s) + (float)(%s);" % (verification.var_name(expr_z3), lhs, rhs)
        else:
            toAppend = "double %s = %s + %s;" % (verification.var_name(expr_z3), lhs, rhs)
        result.append(toAppend)
        return verification.var_name(expr_z3)

    if z3.is_and(expr_z3):
        if DEBUG: print("-- Branch _is_and")
        ##TODO Not sure if symbolTable will be treated in a multi-threaded way
        toAppendExpr = _gen(expr_z3.arg(0), symbolTable, cache, result)
        for i in range(1, expr_z3.num_args()):
            toAppendExpr = 'BAND( %s,%s )' % (toAppendExpr, _gen(expr_z3.arg(i), symbolTable, cache, result))
        toAppend = "double %s = %s; " \
                   % (verification.var_name(expr_z3), \
                      toAppendExpr, \
                      )
        result.append(toAppend)
        return verification.var_name(expr_z3)

    if z3.is_not(expr_z3):
        if DEBUG:
            print("-- Branch _is_not")
        assert expr_z3.num_args() == 1
        if not (expr_z3.arg(0).num_args() == 2):
            warnings.warn("WARNING!!! arg(0) is not RNE but is treated as RNE. arg(0) = ", expr_z3.arg(0))
        op1 = _gen(expr_z3.arg(0).arg(0), symbolTable, cache, result)
        op2 = _gen(expr_z3.arg(0).arg(1), symbolTable, cache, result)
        lhs_type = _get_operand_type(expr_z3.arg(0).arg(0), symbolTable, cache)
        rhs_type = _get_operand_type(expr_z3.arg(0).arg(1), symbolTable, cache)
        if z3_util.is_ge(expr_z3.arg(0)):
            func = _get_comparison_function("DLT", lhs_type, rhs_type)
        elif z3_util.is_gt(expr_z3.arg(0)):
            func = _get_comparison_function("DLE", lhs_type, rhs_type)
        elif z3_util.is_le(expr_z3.arg(0)):
            func = _get_comparison_function("DGT", lhs_type, rhs_type)
        elif z3_util.is_lt(expr_z3.arg(0)):
            func = _get_comparison_function("DGE", lhs_type, rhs_type)
        elif z3_util.is_eq(expr_z3.arg(0)):
            func = _get_comparison_function("DNE", lhs_type, rhs_type)
        elif z3_util.is_distinct(expr_z3.arg(0)):
            func = _get_comparison_function("DEQ", lhs_type, rhs_type)
        else:
            raise NotImplementedError("Not implemented case")
        a = "%s(%s,%s)" % (func, op1, op2)
        toAppend = "double %s = %s;" % (verification.var_name(expr_z3), a)
        result.append(toAppend)
        return verification.var_name(expr_z3)

    if z3_util.is_fpNeg(expr_z3):
        if DEBUG:
            print("-- Branch _is_fpNeg")
        assert expr_z3.num_args() == 1
        op1 = _gen(expr_z3.arg(0), symbolTable, cache, result)
        toAppend = "%s %s =  - %s ;" \
                   % (expr_type, verification.var_name(expr_z3), \
                      op1, \
                      )
        result.append(toAppend)
        return verification.var_name(expr_z3)

    raise NotImplementedError(
        "Not implemented case 002 for expr_z3  =  %s, kind(%s)" % (expr_z3, expr_z3.decl().kind()))

def gen(expr_z3):
    symbolTable = collections.OrderedDict()
    cache = set()
    result = []
    _gen(expr_z3, symbolTable, cache, result)   #########STOHERE
    if len(symbolTable)==0:
        return symbolTable,'int main(){return 0;}'
    # Build variable declarations with correct types
    var_declarations = []
    parse_formats = []
    var_refs = []
    for var_name, var_type in symbolTable.items():
        if var_type == Sort.Float32:
            var_declarations.append(f"float {var_name};")
            parse_formats.append("f")  # 'f' for float in PyArg_ParseTuple
            var_refs.append(f"&{var_name}")
        elif var_type == Sort.Float64:
            var_declarations.append(f"double {var_name};")
            parse_formats.append("d")  # 'd' for double
            var_refs.append(f"&{var_name}")
        else:
            raise NotImplementedError("Unknown types in smt")
    x_expr = verification.var_name(expr_z3)   #the last var
    x_body = '\n  '.join(result)
    x_dim = len(symbolTable)
    return symbolTable, _get_template() % {
        "var_declarations": "\n  ".join(var_declarations),
        "parse_formats": "".join(parse_formats),
        "var_refs": ", ".join(var_refs),
        "x_expr": x_expr,
        "x_dim": x_dim,
        "x_body": x_body
    }

def get_parser():
    parser = argparse.ArgumentParser(prog='XSat')
    parser.add_argument('smt2_file', help='specify the smt2 file to analyze.', type=argparse.FileType('r'))
    parser.add_argument('-v', '--version', action='version', version='%(prog) version 12/18/2015')
    parser.add_argument('--niter', help='niter in basinhopping', action='store', type=int, required=False, default=100)
    parser.add_argument('--nStartOver', help='startOver times', action='store', type=int, required=False, default=2)
    parser.add_argument('--method', help='Local minimization procedure', default='powell',
                        choices=['powell', 'slsqp', 'cg', 'l-bfgs-b', 'cobyla', 'tnc', 'bfgs', 'nelder-mead',
                                 'noop_min']
                        )
    parser.add_argument('--showTime', help='show the time-related info (default: false)', action='store_true',
                        default=False)
    parser.add_argument('--showResult', help='show the basinhopping output (default:false)', action='store_true',
                        default=False)
    parser.add_argument('--stepSize', help='parameter of basinhopping', type=float, default=10.0);
    parser.add_argument('--stepSize_round2', help='parameter of basinhopping', type=float, default=100.0);
    parser.add_argument('--verify', help='verify the model', action='store_true', default=False)
    parser.add_argument('--verify2', help='verify the model (method 2)', action='store_true', default=False)
    parser.add_argument('--showModel', help='show the model as a var->value mapping', action='store_true',
                        default=False)
    parser.add_argument('--showSymbolTable', help='show the symbol table, var->type', action='store_true',
                        default=False)
    parser.add_argument('--showConstraint', help='show the constraint, using the Z3 frontend', action='store_true',
                        default=False)
    parser.add_argument('--showVariableNumber', help='show variable number, using the Z3 frontend', action='store_true',
                        default=False)

    parser.add_argument('--command_compilation', help='the command used to compile the generated foo.c to foo.so',
                        default='clang -O3 -fbracket-depth=2048 -fPIC -I /usr/local/Cellar/python/2.7.9/Frameworks/Python.framework/Versions/2.7/include/python2.7/ %(file)s.c -dynamiclib -o %(file)s.so -L /usr/local/Cellar/python/2.7.9/Frameworks/Python.framework/Versions/Current/lib/ -lpython2.7')
    parser.add_argument('--startPoint', help='start point in a single dimension', action='store', type=float,
                        default=1.0);
    parser.add_argument("--multi", help="multi-processing (default: false)", default=False, action='store_true')
    parser.add_argument("--multiMessage", help="multi-processing message", default=False, action='store_true')
    parser.add_argument("--round2", help="activate round2 when unsat (default: false)", default=False,
                        action='store_true')
    parser.add_argument("--niter_round2", help="niter for round2", action='store', type=int, required=False,
                        default=100)
    parser.add_argument("--suppressWarning", help="Suppress warnings", default=False, action='store_true')
    return parser


if __name__ == "__main__":
    parser = get_parser()
    if len(sys.argv[1:]) == 0:
        _print_xsatInfo()
        parser.print_help()
        parser.exit()
    args = parser.parse_args()
    if args.suppressWarning:
        warnings.filterwarnings("ignore")
    try:
        expr_z3_lis = z3.parse_smt2_string(args.smt2_file.read())
        expr_z3 = z3.And(expr_z3_lis)
        expr_z3 = z3.simplify(expr_z3, arith_lhs=False, hoist_cmul=False)
    except z3.Z3Exception as e:
        print(e)
        sys.stderr.write("[Xsat] The Z3 fornt-end crashes.\n")
    symbolTable, foo_dot_c = gen(expr_z3)
    args.smt2_file.close()

    # dump symbolTable for future verification step (in xsat.py)
    pickle.dump(symbolTable, open("build/foo.symbolTable", "wb"))
    print(foo_dot_c)
