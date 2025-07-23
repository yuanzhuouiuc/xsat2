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
import importlib
import threading

from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent / 'src'))
from src.core.config import SolverConfig
import src.optimization.mcmc as op_mcmc
import src.utils.verification as verification

def str2bool(v: str) -> bool:
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')



def parser():
    parser = argparse.ArgumentParser(prog='Xsat')
    parser.add_argument('-v', '--version', action='version', version='%(prog) version 2.0.0')
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
    parser.add_argument('--stepSize', help='parameter of basinhopping', type=float, default=10.0)
    parser.add_argument('--round2_stepsize', help='parameter of basinhopping', type=float, default=100.0)
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

    # TODO: previously hardcode compile command, need to check it out
    parser.add_argument('--command_compilation', help='the command used to compile the generated foo.c to foo.so',
                        default='gcc -O3 -fbracket-depth=2048 -fPIC -I /usr/local/Cellar/python/2.7.9/Frameworks/Python.framework/Versions/2.7/include/python2.7/ %(file)s.c -dynamiclib -o %(file)s.so -L /usr/local/Cellar/python/2.7.9/Frameworks/Python.framework/Versions/Current/lib/ -lpython2.7')

    parser.add_argument('--startPoint', help='start point in a single dimension', action='store', type=float,
                        default=1.0)
    parser.add_argument('--round2_threshold', help='threshold_low for round2', action='store', type=float, default=1e9)
    parser.add_argument('--round3_threshold', help='threshold  for round3', action='store', type=float, default=1e10)
    parser.add_argument("--multi", help="multi-processing (default: true)", default=True, action='store', type=str2bool)
    # parser.add_argument("--single", help="single processor  (default: true)",default=True,action='store')
    # parser.add_argument("--round2", help="activate round2 when unsat (default: false)",default=False,action='store_true')
    parser.add_argument("--round2_niter", help="niter for round2", action='store', type=int, required=False, default=50)
    parser.add_argument("--round3_niter", help="niter for round3", action='store', type=int, required=False,
                        default=1000)
    parser.add_argument("--round3_stepsize", help="stepsize for round3", action='store', type=float, required=False,
                        default=5.0)
    parser.add_argument("--suppressWarning", help="Suppress warnings", default=False, action='store_true')
    parser.add_argument("--debug", help="debug mode (with verify and showresults, etc.)", default=True,
                        action='store_true')
    parser.add_argument("--printModel", help="print the model", default=False, action='store_true')
    parser.add_argument("--bench", help="benchmarking mode", default=False, action='store_true')
    parser.add_argument("--genOnly", help="generate code only, without deciding satisfiability", default=False,
                        action='store_true')
    return parser

def configure(args):
    config = SolverConfig()
    config.niter = args.niter
    config.method = args.method
    config.step_size = args.stepSize
    config.multi_process = args.multi

    if args.bench:
        args.debug = False
        args.verify = False
        args.verify2 = False
        args.showResult = False
        args.showTime = False
        args.suppressWarning = True
        args.multi = True
    if args.debug:
        args.verify = True
        args.verify2 = True
        args.showResult = True
        args.showTime = True
        args.suppressWarning = False
    return config

def main():
    args = parser().parse_args()

    if args.suppressWarning:
        warnings.filterwarnings("ignore")

    t_start = time.time()
    # use z3 frontend
    with open("XSAT_IN.txt") as f:
        try:
            expr_z3 = z3.simplify(z3.parse_smt2_file(f.read().rstrip()))
        except z3.Z3Exception:
            sys.stderr.write("[Xsat] The Z3 fornt-end fails when verifying the model.\n")
    with open("build/foo.symbolTable", "rb") as f:
        symbolTable = pickle.load(f)
        if len(symbolTable) == 0:
            print("sat")
            sys.exit(0)
    if not args.multi:
        # round1
        t_round1_start = time.time()
        (X_star, R_star) = op_mcmc.mcmc_round1(args)
        t_round1_end = time.time()
        satisfiable_round1 = verification.verify_solution(expr_z3, X_star, symbolTable, printModel=args.printModel)
        if satisfiable_round1 and R_star != 0:
            sys.stderr.write("WARNING!!!!!!!!!!!!!!!! Actually sat.\n")
        elif not satisfiable_round1 and R_star == 0:
            sys.stderr.write("WARNING!!!!!!!!!!!!!!!  Wrong model Maybe unsat!\n")
        else:
            pass

        # round2
        t_round2_start = time.time()
        if not satisfiable_round1:
            X_star, R_star = op_mcmc.mcmc_round2(args, X_star)
        t_round2_end = time.time()

        # round3
        t_round3_start = time.time()
        if R_star > 0 and R_star < args.round3_threshold:
            (X_star, R_star) = op_mcmc.mcmc_round3(args, X_star)
        t_round3_end = time.time()

        # print out results (optional part)
        if args.showResult:
            print("X_star (final)", X_star)
            print("R_star (final)", R_star)
        # print out results (mandatory part)
        if R_star == 0:
            print('sat')
        else:
            print('unsat')
        # verification part is shown at the end of this file.
    else:
        if args.showTime:
            print("[Xsat] ENTERING: main_multi")
        results_pool = []
        result_lock = threading.Lock()
        pool = mp.Pool()
        # execute it quickly, since a lock is set
        def log_result(result):
            X_star, R_star = result
            if args.showTime:
                print(f"[Xsat-multi] ENTERING: {mp.current_process().name} , log_result Minimum =, {R_star}")
            # assert len(results_pool) <= 1
            with result_lock:
                if len(results_pool) == 0:
                    results_pool.append(result)
                else:
                    X_star_pool, R_star_pool = results_pool[0]
                    if R_star < R_star_pool:
                        results_pool[0] = result
                        if R_star_pool == 0:
                            if args.showTime:
                                print("[Xsat-multi] I kill the other process now!!!")
                            pool.terminate()
        for i in range(mp.cpu_count()):
            p = pool.apply_async(op_mcmc.mcmc, args=(args, i,), callback=log_result)
        pool.close()
        pool.join()
        assert len(results_pool) == 1
        (X_star, R_star) = results_pool[0]
        if X_star.ndim == 0: X_star = np.array([X_star[()]])
        if R_star == 0:
            print('sat')
        else:
            print('unsat')
        if args.showResult:
            print(f"X_star (final) {X_star}")
            print(f"R_star (final) {R_star}")
    t_mcmc = time.time()
    if args.verify:
        if args.showTime:
            print("[Xsat] verify X_star with z3 front-end")
        verified = verification.verify_solution(expr_z3, X_star, symbolTable, printModel=args.printModel)
        if verified and R_star != 0:
            sys.stderr.write("WARNING!!!!!!!!!!!!!!!! Actually sat.\n")
        elif not verified and R_star == 0:
            sys.stderr.write("WARNING!!!!!!!!!!!!!!!  Wrong model !\n")
        else:
            pass
    if args.verify2:
        if args.showTime:
            print("[Xsat] verify X_star with build/R_verify")
        sys.path.insert(0, os.path.join(os.getcwd(), "build/R_verify"))
        import foo as foo_verify
        importlib.reload(foo_verify)  # necessary because name 'foo' now still points to foo_square
        verify_res = foo_verify.R(*X_star) if foo_verify.dim == 1 else foo_verify.R(*(X_star))
        if verify_res == 0 and R_star != 0:
            sys.stderr.write("WARNING from verify2 (using include/R_verify/xsat.h) !!!!!!!!!!!!!!!! Actually sat.\n")
        elif verify_res != 0 and R_star == 0:
            sys.stderr.write("WARNING from verify2  (using include/R_verify/xsat.h) !!!!!!!!!!!!!!!  Wrong model ! \n")
        else:
            pass
    t_verify = time.time()
    if args.showSymbolTable:
        print(symbolTable)
    if args.showConstraint:
        print(expr_z3)
    if args.showVariableNumber:
        print("nVar = ", len(symbolTable))
    if args.showTime:
        print("[Xsat] Time elapsed:")
        #        print "  parse_and_compile    : %g seconds " % (t_parse_and_compile-t_start)
        print("  solve (all)  : %g seconds" % (t_mcmc - t_start))
        if not args.multi:
            print("        round1  : %g seconds" % (t_round1_end - t_round1_start))
            print("        round2  : %g seconds" % (t_round2_end - t_round2_start))
            print("        round3  : %g seconds" % (t_round3_end - t_round3_start))
            print("        verification after round1 also takes a little time")
        print("  verify : %g seconds" % (t_verify - t_mcmc))
        print("\n  Total        : %g seconds" % (t_verify - t_start))


if __name__ == "__main__":
    main()
