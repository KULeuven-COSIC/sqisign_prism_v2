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

def hash_to_prime(msg, E_pk, r=None):
    """
    Hash (message || pk || r) into an (a-1)-bit number
    randomly generating r until the output is prime.
    """

    enc_pk = misc.encode_curve_j(E_pk)
    if type(msg) == str: msg = msg.encode()
    if b'&&' in msg:
        # TODO: better domain separation
        raise ValueError('invalid message')

    msg = enc_pk + b'&&' + msg + b'&&'

    if r:
        chl = msg + r
        h = hashlib.sha256(chl).digest()
        q = int.from_bytes(h, 'big') % (2**params.a)
        if is_pseudoprime(q):
            return q, r
        raise ValueError("invalid salt provided")

    while True:
        r = randint(0, 2**params.rb)
        r = int(r).to_bytes((params.rb+7)// 8, 'big')
        # cnt = str(counter).encode()
        chl = msg + r
        h = hashlib.sha256(chl).digest()
        q = int.from_bytes(h, 'big') % (2**params.a)
        if is_pseudoprime(q):
            return q, r

class PRISM:

    def __init__(self):
        # Key generation
        self.pk, self.sk = self.key_gen()

    def key_gen(self):
        """
        Key generation
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
        P_pk *= 2**(params.f - params.a)
        Q_pk *= 2**(params.f - params.a)
        M_sk = ec.ChangeOfBasis((phi_P0, phi_Q0), (P_pk, Q_pk))
        assert ec.EvalMatrix(M_sk, (phi_P0, phi_Q0)) == (P_pk, Q_pk)

        pair_pk = pari.ellweilpairing(E_pk, P_pk, Q_pk, 2**params.a)
        pk = (E_pk, (P_pk, Q_pk), pair_pk)
        sk = (E_pk, I_sk, M_sk)

        _t3 = time.time()
        logger.info(f'- Change of basis done: {_t3-_t2:.3f}s')
        logger.info(f'Keygen done: {_t3-_t0:.3f}s')

        return pk, sk

    def sign(self, msg):
        """
        Signing.
        Input:
        - msg: the message
        Output:
        - sig: a valid signature
        """
        logger.info('Starting signing')
        _t0 = time.time()

        E_pk, I_sk, M_sk = self.sk

        # Hash the message to a prime
        q, r = hash_to_prime(msg, E_pk)
        assert q < 2**params.a

        _t1 = time.time()
        logger.info(f'- Hashing done: {_t1-_t0:.3f}s')

        # Construct the response to the challenge
        n_rsp = q*(2**params.a - q)
        I_rsp = qt.RandomIdealGivenNorm(n_rsp, False)
        I_rsp = qt.Pushforward(I_rsp, I_sk)
        I_cra = I_sk * I_rsp

        _t2 = time.time()
        logger.info(f'- Response ideal done: {_t2-_t1:.3f}s')

        # Compute the corresponding isogeny
        E_rsp, P_cra, Q_cra = qlpt.IdealToIsogeny(I_cra)

        # Scale torsion and multiplication by q^-1
        q_inv = inverse_mod(q, 2**params.a)
        P_rsp, Q_rsp = ec.EvalMatrix(q_inv*M_sk, basis=(P_cra, Q_cra))

        pts_rsp = (P_rsp, Q_rsp)
        sigma = (
            E_rsp, pts_rsp, r
        )

        _t3 = time.time()
        logger.info(f'- Response isogeny done: {_t3-_t2:.3f}s')
        logger.info(f'Signature done: {_t3-_t0:.3f}s')
        return sigma


def PRISM_verify(msg, sigma, pk):
    """
    Verification
    Input:
    - msg: message
    - sigma: signature
    - pk: public key
    Output:
    - boolean verification
    """
    logger.info('Starting verification')
    _t0 = time.time()

    # (P_pk, Q_pk) are 2^a-torsion
    E_pk, (P_pk, Q_pk), pair_pk = pk
    E_rsp, pts_rsp, r = sigma
    q, r1 = hash_to_prime(msg, E_pk, r = r)
    assert r1 == r # Hashing gives a prime

    P_rsp, Q_rsp = pts_rsp

    # Compute a 2D isogeny to check the response
    K = ((P_pk, P_rsp), (Q_pk, Q_rsp))
    Phi = hd.Dim2Iso(K, params.a)

    _t1 = time.time()
    logger.info(f'- 2D isogeny done: {_t1-_t0:.3f}s')

    # Check the degree using pairings (à la SQIsign2D-East)
    P = hd.CouplePoint(P_pk, E_rsp(0))
    Q = hd.CouplePoint(Q_pk, E_rsp(0))

    P1, Q1 = Phi(P)[0], Phi(Q)[0]

    pair = pari.ellweilpairing(P1.curve(), P1, Q1, 2**params.a)
    pair_q = pair_pk ** q
    pair_qinv = pair_q ** (-1)

    if pair in [pair_q, pair_qinv]:
        _t2 = time.time()
        logger.info(f'- Degree checking done: {_t2-_t1:.3f}s')
        logger.info(f'Verification done: {_t2-_t0:.3f}s')
        return True
    logger.error('Degree not matching')
    return False

if __name__ == "__main__":
    # Add logging
    logger.setLevel(logging.DEBUG)

    rr = randint(1, 10000)

    set_random_seed(rr)
    logger.debug(f'Running with seed {rr}')

    # Setting parameters
    lvl = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    params.set_prism_params(lvl)

    alice = PRISM()

    msg = 'Hello world'
    sigma = alice.sign(msg)

    out = PRISM_verify(msg, sigma, alice.pk)
    print(f'Verification: {out}')

