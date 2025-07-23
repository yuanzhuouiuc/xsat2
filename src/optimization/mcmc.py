import os
import sys
import struct
import importlib
import warnings
import numpy as np
import scipy.optimize as op
import multiprocessing as mp
from typing import Tuple

foundit = mp.Event()
def _callback_global(x, f, accepted):
    #    print("at minimum %.4f  %s" % (f, 'accepted' if (accepted) else 'not accepted'))
    if f == 0 or foundit.is_set():
        foundit.set()
        return True

#to handle an issue due to 'powell': it returns a zero-dimensional array even if the starting point is of one dimension.
def tr_help(X):
        if X.ndim == 0: return np.array([X])
        else: return X

def scales():
    return [(lambda x: x ** 11, lambda x: np.sign(x) * np.abs(x) ** (1.0/11)),
            (lambda x: x ** 17, lambda x: np.sign(x) * np.abs(x) ** (1.0/17)),
            (lambda x: x ** 25, lambda x: np.sign(x) * np.abs(x) ** (1.0/25))]

def noop_min(fun, x0, args, **options):
    return op.OptimizeResult(x=x0, fun=fun(x0), success=True, nfev=1)

def scale(X, i):
    return X ** (2 * i + 1)
    # if i ==0:
    #    return X
    # else:
    #    return X**(4*i+3)

def R_quick(X,i,f):
    return f(* scale(X,i))
#    return foo.R(* (X**(2*i+1)))

def mcmc_bis(i):
    print("*******value of i = ", i)

# little-endian
@np.vectorize
def nth_fp_vectorized(n, x):
    if x < 0: return -nth_fp_vectorized(-n, -x)
    n = int(n)
    s = struct.pack('<d', x)
    i = struct.unpack('<Q', s)[0]
    m = i + n
    # m = n + struct.unpack('!i',struct.pack('!f',x))[0]
    if m < 0:
        sign_bit = 0x8000000000000000
        m = -m
    else:
        sign_bit = 0
    if m >= 0x7ff0000000000000:
        warnings.warn("Value out of range, with n= %g,x=%g,m=%g, process=%g" % (n, x, m, mp.current_process.name))
        m = 0x7ff0000000000000
    #        raise ValueError('out of range')
    bit_pattern = struct.pack('Q', m | sign_bit)
    return struct.unpack('d', bit_pattern)[0]

# def nth_fp_vectorized2(N,X):
#     return np.vectorize(nth_fp)(N,X)

# round1 use R_square/foo.so to quickly converge to a minimum point. This round1 refers to single-processor case.
def mcmc_round1(args):
    sys.path.insert(0, os.path.join(os.getcwd(), "build/R_ulp"))
    # sys.path.insert(0,os.path.join(os.getcwd(),"build/R_square"))
    import foo as foo_square
    importlib.reload(foo_square)
    sp = np.zeros(foo_square.dim) + args.startPoint
    # transfer to obj function
    obj = lambda X: foo_square.R(*X)
    res = op.basinhopping(obj, sp, niter=args.niter, stepsize=args.stepSize, minimizer_kwargs={'method': args.method},
                          callback=_callback_global)
    if args.showResult:
        print("Result round 1 with single processor")
        print(res)
        print()
    return tr_help(res.x), res.fun

# round two uses scales and ulp. The starting point is from round1's result.
# X_star: round1 optimization result
def mcmc_round2(args, X_star: np.ndarray) -> Tuple[np.ndarray, float]:
    sys.path.insert(0, os.path.join(os.getcwd(), "build/R_ulp"))
    import foo as foo_ulp
    importlib.reload(foo_ulp)
    if args.showTime:
        print("[Xsat] round2 with single processor")
    res_round1_in_ulp = foo_ulp.R(*X_star)
    # if round1 is good enough, jump round2
    if res_round1_in_ulp <= args.round2_threshold:
        if args.showResult:
            print("round 2 is dismissed. XSat uses the X_star from round1")
            print(f" ==> ulp distance: {res_round1_in_ulp}\n")
        return (X_star, res_round1_in_ulp)
    # scale
    R_star = res_round1_in_ulp
    for (scale, scale_inv) in scales():
        res = op.basinhopping(lambda x: foo_ulp.R(*scale(x)), scale_inv(X_star), niter=args.round2_niter,
                              minimizer_kwargs={'method': args.method}, callback=_callback_global,
                              stepsize=args.round2_stepsize)
        if res.fun < R_star:
            X_star = scale(tr_help(res.x))
            R_star = res.fun
        if args.showResult:
            print("result (round 2): where scale(0.1) =", scale(0.1))
            print(res)
            print()
        if R_star < args.round2_threshold: break
    return X_star, R_star

def mcmc_round3(args, X_star: np.ndarray) -> Tuple[np.ndarray, float]:
    sys.path.insert(0,os.path.join(os.getcwd(), "build/R_ulp"))
    import foo as foo_ulp
    importlib.reload(foo_ulp)
    if args.showTime:
        print("[Xsat] round3 with single processor")
    obj_near = lambda N: foo_ulp.R(* nth_fp_vectorized(N, X_star))
    res = op.basinhopping(obj_near, np.zeros(foo_ulp.dim), niter=args.round3_niter, minimizer_kwargs={'method':args.method},
                          callback=_callback_global, stepsize=args.round3_stepsize)
    if args.showResult:
        print("result (round 3):")
        print(res)
        print()
    R_star = res.fun
    X_star = tr_help(nth_fp_vectorized(res.x, X_star))
    return X_star, R_star

def mcmc(args, i):
    # sys.path.insert(0,os.path.join(os.getcwd(),"build/R_square"))
    sys.path.insert(0, os.path.join(os.getcwd(), "build/R_ulp"))
    import foo
    importlib.reload(foo)  # necessary because name 'foo' now still points to foo_square
    np.random.seed()
    if args.method == 'noop_min':
        _minimizer_kwargs = dict(method=noop_min)
    else:
        _minimizer_kwargs = dict(method=args.method)
    sp = np.zeros(foo.dim) + args.startPoint + i
    res = op.basinhopping(lambda X: R_quick(X, i, foo.R), sp, niter=args.niter, stepsize=args.stepSize,
                          minimizer_kwargs=_minimizer_kwargs, callback=_callback_global)
    if args.showResult:
        print("result (round 1) with i = ", i, ":")
        print(res)
        print()
    # do some change here. If the first round gives a good/bad enough result, no need for the second.

    X_star = scale(res.x, i)
    R_star = res.fun
    if res.fun != 0 and res.fun < args.round2_threshold:
        # if args.round2:
        if args.showTime:
            print("[Xsat] round2_move")
        sys.path.insert(0, os.path.join(os.getcwd(), "build/R_ulp"))
        import foo
        importlib.reload(foo)
        sp = np.array([res.x + 0]) if res.x.ndim == 0 else res.x
        obj_near = lambda N: foo.R(*nth_fp_vectorized(N, scale(sp, i)))
        # print op.fmin_powell(obj_near,np.zeros(foo.dim))
        res_round2 = op.basinhopping(obj_near, np.zeros(foo.dim), niter=args.round2_niter,
                                     stepsize=args.round2_stepsize, minimizer_kwargs=_minimizer_kwargs,
                                     callback=_callback_global)
        #        res_round2 = op.fmin_powell(obj_near,np.zeros(foo.dim))
        if args.showResult:
            print("result (round 2) with i = ", i)
            print(res_round2)
            print()
        R_star = res_round2.fun
        # change this because I could have used the R_quick.
        X_star = nth_fp_vectorized(res_round2.x, scale(sp, i))
    return X_star, R_star

def run_mcmc_single(args):
    sys.path.insert(0, os.path.join(os.getcwd(), "build/R_square"))
    import foo as foo_square
    sp = np.zeros(foo_square.dim) + args.startPoint
    obj = lambda X: foo_square.R(*X)
    res = op.basinhopping(obj, sp, niter=args.niter, stepsize=args.stepSize, minimizer_kwargs={'method': args.method},
                          callback=_callback_global)
    if args.showResult:
        print("result round 1 with single processor ")
        print(res)
        print()
    if args.round2:
        if args.showTime:
            print("[Xsat] round2 with single processor")
        sys.path.insert(0, os.path.join(os.getcwd(), "build/R_ulp"))
        import foo as foo_ulp
        importlib.reload(foo_ulp)  # necessary because name 'foo' now still points to foo_square
        X_star = [res.x + 0] if res.x.ndim == 0 else res.x
        obj_near = lambda N: foo_ulp.R(*nth_fp_vectorized(N, X_star))
        # print op.fmin_powell(obj_near,np.zeros(foo.dim))
        print("*" * 50)
        print(obj_near(0))
        print("*" * 50)
        res_round2 = op.basinhopping(obj_near, np.zeros(foo_ulp.dim), niter=args.round2_niter, stepsize=100.0,
                                     minimizer_kwargs={'method': args.method}, callback=_callback_global)
        if args.showResult:
            print("result (round 2):")
            print(res_round2)
            print()
            res.fun = res_round2.fun
            res.x = nth_fp_vectorized(res_round2.x, X_star)
    return res

