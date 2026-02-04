import hashlib
import logging
import time

logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)
logger_sh = logging.StreamHandler()
formatter = logging.Formatter('%(name)s [%(levelname)s] %(message)s')
logger_sh.setFormatter(formatter)
logger.addHandler(logger_sh)

from sage.all import *
proof.all(False)

try:
    from . import params, ec, hd, misc
    from . import quaternions as qt
    from . import qlapoti as qlpt
except ImportError:
    import params, ec, hd, misc
    import quaternions as qt
    import qlapoti as qlpt

class SQIsign:

    def __init__(self):
        # Key generation
        self.pk, self.sk = self.key_gen()

    def key_gen(self):
        """
        Key generation [Alg 4.1].
        Output:
        - sk: secret key
        - pk: public key
        """
        logger.info('Starting keygen')
        _t0 = time.time()

        I_sk = qt.RandomIdealGivenNorm(params.D_mix, True)
        assert I_sk.left_order() == params.O0 and I_sk.norm() == params.D_mix

        I_sk = qt.RandomEquivalentPrimeIdeal(I_sk)

        _t1 = time.time()
        logger.info(f'- Ideal generation done: {_t1-_t0:.3f}s')

        E_pk, phi_P0, phi_Q0 = qlpt.IdealToIsogeny(I_sk)

        _t2 = time.time()
        logger.info(f'- Qlapoti done: {_t2-_t1:.3f}s')

        P_pk, Q_pk = ec.TorsionBasis(E_pk)
        M_sk = ec.ChangeOfBasis((phi_P0, phi_Q0), (P_pk, Q_pk))
        assert ec.EvalMatrix(M_sk, (phi_P0, phi_Q0)) == (P_pk, Q_pk)

        pk = (E_pk, (P_pk, Q_pk))
        sk = (E_pk, (P_pk, Q_pk), I_sk, M_sk)

        _t3 = time.time()
        logger.info(f'- Change of basis done: {_t3-_t2:.3f}s')
        logger.info(f'Keygen done: {_t3-_t0:.3f}s')

        return pk, sk

    def sign(self, msg):
        """
        Signing [Alg 4.2].
        Input:
        - msg: the message
        Output:
        - sig: a valid signature
        """
        logger.info('Starting signing')
        _t0 = time.time()

        E_pk, (P_pk, Q_pk), I_sk, M_sk = self.sk

        # Commitment
        I_com = qt.RandomIdealGivenNorm(params.D_mix, True)
        I_com = qt.RandomEquivalentPrimeIdeal(I_com)

        E_com, P_com, Q_com = qlpt.IdealToIsogeny(I_com)

        _t1 = time.time()
        logger.info(f'- Commitment done: {_t1-_t0:.3f}s')

        # Challenge
        enc_pk = misc.encode_curve_a(E_pk)
        enc_com = misc.encode_curve_j(E_com)
        if type(msg) == str: msg = msg.encode()

        chl = enc_pk + enc_com + msg
        for _ in range(params.HASH_ITERATIONS - 1):
            chl = hashlib.shake_256(chl).digest(32) # 256 bits
        chl = hashlib.shake_256(chl).digest(ceil(params.e_chl / 8))
        chl = int.from_bytes(chl, 'little') % (2**params.e_chl)

        _t2 = time.time()
        logger.info(f'- Challenge computed: {_t2-_t1:.3f}s')

        # Response
        # - challenge to ideal
        c1, c2 = M_sk.transpose() * vector(Zmod(2**params.f), (1, chl))
        I_chl1 = qt.KernelDecomposedToIdeal(c1, c2)
        I_chl = qt.Pushforward(I_chl1, I_sk)
        assert I_chl.left_order() == I_sk.right_order()

        # - sample response
        I_rsp = I_com.conjugate() * I_sk * I_chl

        # - normalize response
        n_sk, n_com = I_sk.norm(), I_com.norm()
        alpha_rsp = qt.RandomEquivalentQuaternion(I_rsp, n_sk*n_com)

        assert alpha_rsp.reduced_norm() / I_rsp.norm() <= 2**params.e_rsp, "too short"

        alpha_rsp, n_bt = qt.ComputeBacktrackingAndNormalize(alpha_rsp)

        # - find correct degrees for odd and even part
        n_rsp = n_sk * n_com * 2**(params.f - n_bt)
        d_rsp = alpha_rsp.reduced_norm() / n_rsp
        r_rsp = valuation(d_rsp, 2)
        q_rsp = d_rsp / 2**r_rsp

        # - compute commitment + response
        I_com_rsp = params.O0*alpha_rsp.conjugate() + params.O0*(n_com*q_rsp)
        e_rsp1 = params.e_rsp - r_rsp - n_bt # Degree of the 2D isogeny

        _t3 = time.time()
        logger.info(f'- Response ideal computed: {_t3-_t2:.3f}s')

        # - compute and split auxiliary isogeny
        assert e_rsp1 > 0, 'e_rsp1 = 0; this should never happen'

        I_aux = qt.RandomIdealGivenNorm(2**e_rsp1 - q_rsp, False)
        I_cra = I_com_rsp.intersection(I_aux)
        assert I_cra.norm() == n_com * q_rsp * (2**e_rsp1 - q_rsp)

        E_cra, P_cra, Q_cra = qlpt.IdealToIsogeny(I_cra)

        E_chl, (P_chl, Q_chl), E_aux, (P_aux, Q_aux) = \
            hd.SplitAuxiliaryIsogeny(
                E_com, E_cra,
                P_com, Q_com, P_cra, Q_cra,
                q_rsp, e_rsp1, r_rsp
        )

        _t4 = time.time()
        logger.info(f'- Auxiliary isogeny done: {_t4-_t3:.3f}s')

        # - prepare the signature
        if r_rsp > 0:
            E_chl, P_chl, Q_chl = ec.ComputeEvenNonBacktrackingResponse(
                    E_chl, P_chl, Q_chl, alpha_rsp.conjugate(), e_rsp1, r_rsp
            )

        # here the input E_chl must be isomorphic to the actual challenge curve
        E_chl, P_chl, Q_chl = ec.ComputeChallengeIsogeny(
            E_pk, chl, P_pk, Q_pk, E_chl, P_chl, Q_chl, n_bt
        )

        M_chl, B_aux, B_chl = ec.SetChangeOfBasisMatrix(
                E_aux, E_chl, P_aux, Q_aux, P_chl, Q_chl, e_rsp1 + r_rsp
        )

        sigma = (
                E_aux, n_bt, r_rsp, M_chl, chl, B_aux, B_chl
        )

        _t5 = time.time()
        logger.info(f'- Challenge isogeny done: {_t5-_t4:.3f}s')
        logger.info(f'Signature done: {_t5-_t0:.3f}s')

        return sigma

def SQIsign_verify(msg, sigma, pk):
    """
    Verification [Algorithm 4.9]
    Input:
    - msg: message
    - sigma: signature
    - pk: public key
    Output:
    - boolean verification
    """
    logger.info('Starting verification')
    _t0 = time.time()
    E_pk, (P_pk, Q_pk) = pk
    E_aux, n_bt, r_rsp, M_chl, chl, (P_aux, Q_aux), (P_chl, Q_chl) = sigma

    E_pk._order = params.p+1
    E_aux._order = params.p+1

    e_rsp1 = params.e_rsp - r_rsp - n_bt
    if e_rsp1 < 0:
        logger.error('Invalid response size')
        return False

    if M_chl.base_ring().cardinality() != 2**(e_rsp1 + r_rsp):
        logger.error('Invalid matrix M_chl')
        return False

    K = (2**n_bt) * (P_pk + chl*Q_pk)
    K._order = ZZ(2**(params.f - n_bt))

    _t1 = time.time()
    logger.info(f'- Consistency check done: {_t1-_t0:.3f}s')

    phi_chl = E_pk.isogeny(K, model="montgomery")
    E_chl = phi_chl.codomain()
    if not P_chl in E_chl:
        logger.error('Wrong challenge curve')
        return False

    _t2 = time.time()
    logger.info(f'- 1D part done: {_t2-_t1:.3f}s')

    # Scaling the bases
    cff_aux = 2**(params.f - e_rsp1)
    cff_chl = 2**(params.f - e_rsp1 - r_rsp)

    P_aux, Q_aux = cff_aux * P_aux, cff_aux * Q_aux
    P_chl, Q_chl = cff_chl * P_chl, cff_chl * Q_chl

    _a, _b, _c, _d = M_chl.list()
    P_chl, Q_chl = _a*P_chl + _b*Q_chl, _c*P_chl + _d*Q_chl

    # Update the even part
    if r_rsp > 0:
        if (M_chl[0][0] % 2) ==  (M_chl[0][1] % 2) == 0:
            # P_chl is in the kernel of the isogeny -> Q_chl is the dual
            K_rsp = (2**e_rsp1) * Q_chl
        else:
            # Either P_chl is the dual or both are
            K_rsp = (2**e_rsp1) * P_chl

        E_chl._order = params.p + 1
        K_rsp._order = ZZ(2**(r_rsp))
        phi_even = E_chl.isogeny(K_rsp)
        E_chl = phi_even.codomain()
        P_chl, Q_chl = phi_even(P_chl), phi_even(Q_chl)

    _t3 = time.time()
    logger.info(f'- Scaling and even part done: {_t3-_t2:.3f}s')

    # 2D part
    K = ((P_chl, P_aux), (Q_chl, Q_aux))
    Phi = hd.Dim2Iso(K, e_rsp1)

    _t4 = time.time()
    logger.info(f'- 2D part done: {_t4-_t3:.3f}s')

    # Check Fiat-Shamir
    E1, E2 = Phi.codomain()
    # TODO: know in advance which one?
    for E_com in [E1, E2]:
        enc_pk = misc.encode_curve_a(E_pk)
        enc_com = misc.encode_curve_j(E_com)
        if type(msg) == str: msg = msg.encode()

        check_chl = enc_pk + enc_com + msg
        for _ in range(params.HASH_ITERATIONS - 1):
            check_chl = hashlib.shake_256(check_chl).digest(32) # 256 bits
        check_chl = hashlib.shake_256(check_chl).digest(ceil(params.e_chl / 8))
        check_chl = int.from_bytes(check_chl, 'little') % (2**params.e_chl)
        if check_chl == chl:
            _t5 = time.time()
            logger.info(f'- Hashing done: {_t5-_t4:.3f}s')
            logger.info(f'Verification done: {_t5-_t0:.3f}s')
            return True
    logger.error('Hash not matching')
    return False

if __name__ == "__main__":
    # Add logging
    logger.setLevel(logging.DEBUG)

    rr = randint(1, 2**64)

    set_random_seed(rr)
    logger.debug(f'Running with seed {rr}')

    # Setting parameters [Sec 4.2]
    lvl = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    params.set_sqi_params(lvl)

    alice = SQIsign()

    msg = 'Hello world'
    sigma = alice.sign(msg)

    out = SQIsign_verify(msg, sigma, alice.pk)
    print(f'Verification: {out}')


