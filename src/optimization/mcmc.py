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

def R_quick(X,i,f):
    return f(* scale(X,i))

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
    bit_pattern = struct.pack('Q', m | sign_bit)
    return struct.unpack('d', bit_pattern)[0]

@np.vectorize
def nth_fp32_vectorized(n, x):
    if x < 0: return -nth_fp32_vectorized(-n, -x)
    n = int(n)
    x_f32 = np.float32(x)
    s = struct.pack('<f', x_f32)
    i = struct.unpack('<I', s)[0]
    m = i + n
    if m < 0:
        sign_bit = 0x80000000
        m = -m
    else:
        sign_bit = 0
    if m >= 0x7f800000:  # Float32 infinity
        warnings.warn(f"Float32 value out of range, n={n}, x={x}")
        m = 0x7f800000
    bit_pattern = struct.pack('I', m | sign_bit)
    return struct.unpack('f', bit_pattern)[0]

def _init():
    sys.path.insert(0, os.path.join(os.getcwd(), "build/R_ulp"))
    import foo
    importlib.reload(foo)
    with open("build/f32_mask.npy", "rb") as f:
        f32_mask = np.load(f)
    _original_R = foo.R
    f32_indices = np.where(f32_mask)[0]
    f64_indices = np.where(~f32_mask)[0]
    # if np.any(f32_mask):
    #     buffer = np.zeros(foo.dim, dtype=np.float64)
    #     def R_typed(*args):
    #         """Wrapped version that handles type conversion internally"""
    #         # Handle different input formats
    #         if len(args) == 1 and hasattr(args[0], '__len__'):
    #             X = args[0]
    #         else:
    #             X = np.array(args)
    #         if len(f32_indices) > 0:
    #             buffer[f32_indices] = X[f32_indices].astype(np.float32)
    #         if len(f64_indices) > 0:
    #             buffer[f64_indices] = X[f64_indices]
    #         return _original_R(*buffer)
    #     foo.R = R_typed
    return foo, f32_indices, f64_indices

def mcmc(args, i):
    foo, f32_indices, f64_indices = _init()
    nth_fp_dispatchers = [
        nth_fp32_vectorized if j in f32_indices else nth_fp_vectorized
        for j in range(foo.dim)
    ]
    np.random.seed()
    _minimizer_kwargs = dict(method=noop_min) if args.method == 'noop_min' else dict(method=args.method)
    best_X_star = np.zeros(foo.dim)
    best_R_star = float('inf')
    for round_num in range(args.nStartOver):
        if args.showResult:
            print(f"--- Process-{i}, Start-Over Attempt {round_num + 1}/{args.nStartOver} ---")
        sp = np.zeros(foo.dim) + args.startPoint + i + np.random.uniform(-0.5, 0.5, foo.dim)
        res_global = op.basinhopping(
            lambda X: foo.R(*scale(X, i)),
            sp,
            niter=args.niter,
            stepsize=args.stepSize,
            minimizer_kwargs=_minimizer_kwargs,
            callback=_callback_global
        )
        if args.showResult:
            print(f"Result (Global Search, i={i}, attempt={round_num + 1}):")
            print(res_global)
        current_X_star = scale(tr_help(res_global.x), i)
        current_R_star = res_global.fun
        if 0 < current_R_star < args.round2_threshold:
            if args.showTime:
                print(f"[Xsat Process-{i}] Starting local refinement...")
            sp_refine = tr_help(res_global.x)
            X_scaled_cached = scale(sp_refine, i)
            # This objective function now iterates through the N values and applies
            # the appropriate scalar nth_fp function for each element based on the dispatch table.
            def obj_near(N):
                X_moved = np.array([
                    nth_fp_dispatchers[j](n_val, x_base_val)
                    for j, (n_val, x_base_val) in enumerate(zip(N, X_scaled_cached))
                ])
                return foo.R(*X_moved)
            res_refine = op.basinhopping(
                obj_near,
                np.zeros(foo.dim),
                niter=args.round2_niter,
                stepsize=args.round2_stepsize,
                minimizer_kwargs=_minimizer_kwargs,
                callback=_callback_global
            )
            if args.showResult:
                print(f"Result (Local Refinement, i={i}, attempt={round_num + 1}):")
                print(res_refine)
            if res_refine.fun < current_R_star:
                current_R_star = res_refine.fun
                current_X_star = np.array([
                    nth_fp_dispatchers[j](n_val, x_base_val)
                    for j, (n_val, x_base_val) in enumerate(zip(res_refine.x, X_scaled_cached))
                ])
        if current_R_star < best_R_star:
            best_R_star = current_R_star
            best_X_star = current_X_star
            if args.showResult:
                print(f"--- Process-{i} New Best Found: R_star = {best_R_star:.4g} ---")
        if best_R_star == 0:
            if args.showResult:
                print(f"Solution found by Process-{i} in attempt {round_num + 1}.")
            break
    return best_X_star, best_R_star

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
