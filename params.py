from sage.all import (
        next_prime, ceil, log, sqrt, ZZ, QuaternionAlgebra, GF, EllipticCurve,
        matrix, Zmod, floor, next_prime
)

import ec
import constants

# Fixed parameters [Sec 4.2, Appendix B]

# General scheme parameters
p, f, lb, rb, a = [None for _ in range(5)]
D_mix, e_rsp, D_rsp, D_chl, e_chl = [None for _ in range(5)]
HASH_ITERATIONS, FP_ENC_BYTES = None, None

# Quaternion computations and ideal to isogeny
QUAT_primality_num_iter = None
QUAT_equiv_bound_coeff = 64
QUAT_prime_cofactor = None
B, O0 = None, None

# Elliptic curves
Fp, Fp2, Fp2_i = None, None, None
E0, P0, Q0 = None, None, None
iota, frob = None, None
mat_1, mat_i, mat_ij2, mat_1k2 = [None for _ in range(4)]

def set_sqi_params(lvl):
    """
    Return parameters for sqisign at the given level
    """
    global p, f, lb, D_mix, e_rsp, D_rsp, D_chl, e_chl
    global HASH_ITERATIONS, FP_ENC_BYTES

    if lvl == 1:
        p = ZZ(5*2**248 - 1)
        f = ZZ(248)
        lb = ZZ(128)

    elif lvl == 3:
        p = ZZ(65*2**376 - 1)
        f = ZZ(376)
        lb = ZZ(192)

    elif lvl == 5:
        p = ZZ(27*2**500 - 1)
        f = ZZ(500)
        lb = ZZ(256)

    else:
        raise ValueError(f"level {lvl} not recognized")

    D_mix = next_prime(2**(4*lb))
    e_rsp = ceil(log(sqrt(p), 2))
    D_rsp = ZZ(2**e_rsp)
    D_chl = ZZ(2**f)
    e_chl = ZZ(f - e_rsp)

    HASH_ITERATIONS = 2**(32 * ceil(log(p, 2) / 64) - (f - e_rsp))
    FP_ENC_BYTES = 8 * floor((log(p, 2) + 63) / 64)

    # Quaternions
    global B, O0
    global QUAT_prime_cofactor
    B = QuaternionAlgebra(-1, -p)
    _i, _j, _k = B.gens()
    O0 = B.maximal_order(order_basis=(B(1), _i, (_i+_j)/2, (1-_k)/2))
    QUAT_prime_cofactor = next_prime(2**ceil(log(p, 2)))

    # Elliptic Curves
    global Fp, Fp2, Fp2_i
    global E0, P0, Q0, unred_tate
    Fp = GF(p)
    Fp2, Fp2_i = GF(p**2, name="i", modulus=[1, 0, 1]).objgen()
    E0 = EllipticCurve(Fp2, [1, 0])
    E0.set_order((p + 1) ** 2)

    # Hardcoded SQIsign implementation points and endomorphisms actions
    global iota, frob
    global mat_1, mat_i, mat_ij2, mat_1k2
    mat_1 = matrix(Zmod(2**f), 2, [1, 0, 0, 1])

    if lvl == 1:
        P0 = E0(constants.Px1, constants.Py1)
        Q0 = E0(constants.Qx1, constants.Qy1)
        mat_i = matrix(Zmod(2**f), 2, constants.mat_i_1)
        mat_ij2 = matrix(Zmod(2**f), 2, constants.mat_ij2_1)
        mat_1k2 = matrix(Zmod(2**f), 2, constants.mat_1k2_1)

        unred_tate = Fp2(constants.unred_tate1)

    elif lvl == 3:
        P0 = E0(constants.Px3, constants.Py3)
        Q0 = E0(constants.Qx3, constants.Qy3)
        mat_i = matrix(Zmod(2**f), 2, constants.mat_i_3)
        mat_ij2 = matrix(Zmod(2**f), 2, constants.mat_ij2_3)
        mat_1k2 = matrix(Zmod(2**f), 2, constants.mat_1k2_3)

        unred_tate = Fp2(constants.unred_tate3)

    elif lvl == 5:
        P0 = E0(constants.Px5, constants.Py5)
        Q0 = E0(constants.Qx5, constants.Qy5)
        mat_i = matrix(Zmod(2**f), 2, constants.mat_i_5)
        mat_ij2 = matrix(Zmod(2**f), 2, constants.mat_ij2_5)
        mat_1k2 = matrix(Zmod(2**f), 2, constants.mat_1k2_5)

        unred_tate = Fp2(constants.unred_tate3)

def set_prism_params(lvl):
    """
    Return parameters for sqisign at the given level
    """
    global p, f, lb, a, rb, D_mix
    global FP_ENC_BYTES

    if lvl == 1:
        p = ZZ(5*2**248 - 1)
        f = ZZ(248)
        a = ZZ(192)
        rb = ZZ(32)
        lb = ZZ(248)

    elif lvl == 3:
        p = ZZ(65*2**376 - 1)
        f = ZZ(376)
        a = ZZ(256)
        rb = ZZ(20)
        lb = ZZ(376)

    elif lvl == 5:
        p = ZZ(27*2**500 - 1)
        f = ZZ(500)
        a = ZZ(320)
        rb = ZZ(40)
        lb = ZZ(500)

    else:
        raise ValueError(f"level {lvl} not recognized")

    D_mix = next_prime(2**(4*lb))
    FP_ENC_BYTES = 8 * floor((log(p, 2) + 63) / 64)

    # Quaternions
    global B, O0
    global QUAT_prime_cofactor
    B = QuaternionAlgebra(-1, -p)
    _i, _j, _k = B.gens()
    O0 = B.maximal_order(order_basis=(B(1), _i, (_i+_j)/2, (1-_k)/2))
    QUAT_prime_cofactor = next_prime(2**ceil(log(p, 2)))

    # Elliptic Curves
    global Fp, Fp2, Fp2_i
    global E0, P0, Q0
    Fp = GF(p)
    Fp2, Fp2_i = GF(p**2, name="i", modulus=[1, 0, 1]).objgen()
    E0 = EllipticCurve(Fp2, [1, 0])
    E0.set_order((p + 1) ** 2)

    # Hardcoded SQIsign implementation points and endomorphisms actions
    global iota, frob
    global mat_1, mat_i, mat_ij2, mat_1k2
    mat_1 = matrix(Zmod(2**f), 2, [1, 0, 0, 1])

    if lvl == 1:
        P0 = E0(constants.Px1, constants.Py1)
        Q0 = E0(constants.Qx1, constants.Qy1)
        mat_i = matrix(Zmod(2**f), 2, constants.mat_i_1)
        mat_ij2 = matrix(Zmod(2**f), 2, constants.mat_ij2_1)
        mat_1k2 = matrix(Zmod(2**f), 2, constants.mat_1k2_1)

    elif lvl == 3:
        P0 = E0(constants.Px3, constants.Py3)
        Q0 = E0(constants.Qx3, constants.Qy3)
        mat_i = matrix(Zmod(2**f), 2, constants.mat_i_3)
        mat_ij2 = matrix(Zmod(2**f), 2, constants.mat_ij2_3)
        mat_1k2 = matrix(Zmod(2**f), 2, constants.mat_1k2_3)

    elif lvl == 5:
        P0 = E0(constants.Px5, constants.Py5)
        Q0 = E0(constants.Qx5, constants.Qy5)
        mat_i = matrix(Zmod(2**f), 2, constants.mat_i_5)
        mat_ij2 = matrix(Zmod(2**f), 2, constants.mat_ij2_5)
        mat_1k2 = matrix(Zmod(2**f), 2, constants.mat_1k2_5)
