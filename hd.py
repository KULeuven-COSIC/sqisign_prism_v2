from sage.all import pari, inverse_mod

try:
    from . import params
    from .theta.theta_structures.couple_point import CouplePoint
    from .theta.theta_isogenies.product_isogeny import EllipticProductIsogeny
    from .theta.theta_isogenies.product_isogeny_sqrt import EllipticProductIsogenySqrt
except ImportError:
    import params
    from theta.theta_structures.couple_point import CouplePoint
    from theta.theta_isogenies.product_isogeny import EllipticProductIsogeny
    from theta.theta_isogenies.product_isogeny_sqrt import EllipticProductIsogenySqrt

def Dim2Iso(K, n):
    """
    Compute the 2D (2^n,2^n)-isogeny with kernel K. Wrapper for the theta isogeny
    code. Does not handle split first steps.
    Input:
    - K: two pair of points on E1 x E2, the kernel of the isogeny
    - n: log 2 of the (polarized) degree (or the length of the chain)
    Output:
    - phi: the corresponding isogeny
    """
    # TODO: add strategies and zeta

    (P1, P2), (Q1, Q2) = K
    P = CouplePoint(P1, P2)
    Q = CouplePoint(Q1, Q2)
    kernel = (P, Q)

    return EllipticProductIsogenySqrt(kernel, n)

def EvalTwoPoints(phi, X1, X2):
    """
    Given two 2D-points X1 and X2 on E1xE2 and phi: E1xE2 -> E3xE4, evaluate X1
    and X2 consistently fixing the sign.
    Input:
    - phi: a 2D isogeny
    - X1, X2: CouplePoints on the domain of phi
    Output:
    - phi(X1), phi(X2)
    """
    # TODO: there is probably a cleaner way
    X3 = X1 - X2

    Y1 = phi(X1)
    Y2 = phi(X2)
    Y3 = phi(X3)

    # First component
    out = {
        Y1[0] - Y2[0] : (Y1[0], Y2[0]),
        Y1[0] + Y2[0] : (Y1[0], -Y2[0]),
        -Y1[0] - Y2[0] : (-Y1[0], Y2[0]),
        -Y1[0] + Y2[0] : (-Y1[0], -Y2[0])
    }

    P1, Q1 = out[Y3[0]]

    # Second component
    out = {
        Y1[1] - Y2[1] : (Y1[1], Y2[1]),
        Y1[1] + Y2[1] : (Y1[1], -Y2[1]),
        -Y1[1] - Y2[1] : (-Y1[1], Y2[1]),
        -Y1[1] + Y2[1] : (-Y1[1], -Y2[1])
    }

    P2, Q2 = out[Y3[1]]

    return CouplePoint(P1, P2), CouplePoint(Q1, Q2)

def EmbeddedIsogeny(phi, n, d, basis=None):
    """
    Given the (2^n, 2^n)-isogeny phi embedding a 1D isogeny psi of degree d,
    recover the embedded isogeny and evaluate it on a provided basis (P, Q).
    Input:
    - phi: the 2D-isogeny
    - n: (polarized) 2-degree of phi
    - d: degree of psi
    - basis = (P, Q): basis of the 2^n torsion on the starting curve; if not
      provided, the standard basis (P0, Q0) on E0 is used
    Output:
    - E_psi: the codomain of psi
    - P_psi, Q_psi = psi(P), psi(Q)
    """
    if not basis:
        P = params.P0
        Q = params.Q0
    else:
        P, Q = basis

    T = CouplePoint(P, phi.E2(0))
    S = CouplePoint(Q, phi.E2(0))

    # Evaluate 2D iso on the points
    phi_T, phi_S = EvalTwoPoints(phi, T, S)

    # Detect the correct side with pairings
    # TODO: precomputed for E0
    cof = (params.p**2 - 1) // 2**n
    tPQ = pari.elltatepairing(phi.E1, P, Q, 2**n)**(cof * d)
    t1 = pari.elltatepairing(phi.codomain()[0], phi_T[0], phi_S[0], 2**n)**cof

    if tPQ == t1:
        return phi.codomain()[0], phi_T[0], phi_S[0]

    assert pari.elltatepairing(phi.codomain()[1],phi_T[1],phi_S[1],2**n)**cof==tPQ
    return phi.codomain()[1], phi_T[1], phi_S[1]

def SplitAuxiliaryIsogeny(E1, E2, P1, Q1, P2, Q2, q, e1, r):
    """
    There is an isogeny phi : E1 -> E2 of degree q(2^e1 - q); (P2, Q2) =
    (phi(P1), phi(Q1)). Split the Kani diagram, and return the two sides
    together with image points of order 2^(e1 * r), so that after pushing them
    through the 2^r response part they have the correct order [Algorithm 4.5]
    Input:
    - as above
    Output:
    - F1: first curve
    - (S1, R1): (phi_q(P1), phi_q(Q1))
    - F2: fist curve
    - (S2, R2): (phi_(2^e1-q)(P1), phi_(2^e1-q)(Q1))
    """
    # Prepare the kernel
    f = params.f
    q_inv = inverse_mod(q, 2**e1)

    #  Response points to push through the (eventual) 2^r part
    P1, Q1 = 2**(f-e1-r)*P1, 2**(f-e1-r)*Q1
    P2, Q2 = (2**(f-e1-r) * q_inv)*P2, (2**(f-e1-r)*q_inv)*Q2

    # Actual kernel
    P11, Q11 = 2**r * P1, 2**r * Q1
    P22, Q22 = 2**r * P2, 2**r * Q2

    K = ((P11, P22), (Q11, Q22))

    # Phi
    phi = Dim2Iso(K, e1)

    # Compute the two sides of the diagram; similar to EmbeddedIsogeny, except
    # now we need both sides
    T = CouplePoint(P1, E2(0))
    S = CouplePoint(Q1, E2(0))

    # Evaluate 2D iso on the points
    phi_T, phi_S = EvalTwoPoints(phi, T, S)

    # Check sides with pairings
    ePQ = pari.ellweilpairing(E1, P1, Q1, 2**(e1+r)) ** q
    ephi = pari.ellweilpairing(
            phi.codomain()[0], phi_T[0], phi_S[0], 2**(e1+r))

    # if tPQ == t1:
    if ePQ == ephi:
        # Phi[0] is the q-side
        F1, S1, R1 = phi.codomain()[0], phi_T[0], phi_S[0]
        F2, S2, R2 = phi.codomain()[1], phi_T[1], phi_S[1]

    else:
        # Phi[0] is the 2^e - q side
        F1, S1, R1 = phi.codomain()[1], phi_T[1], phi_S[1]
        F2, S2, R2 = phi.codomain()[0], phi_T[0], phi_S[0]

    return F1, (S1, R1), F2, (S2, R2)
