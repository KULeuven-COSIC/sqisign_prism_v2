from sage.all import (
    randint, mod, valuation, kronecker, gcd, ZZ, QQ, is_prime, Matrix, round,
    vector, Zmod, isqrt, ceil, sqrt
)

try:
    from . import params
except ImportError:
    import params

def ModularSQRT(n, m):
    """
    Modular square root [Alg 3.1]
    Input:
    - m: an odd prime
    - n: an integer which is a square mod m
    Output:
    - an integer x s.t. x^2 = n mod m
    """
    if n % m == 0:
        return ZZ(0)
    if m % 4 == 3:
        return ZZ(mod(n, m)**((m+1)/4))

    # Tonelli-Shanks for m = 1 mod 4
    # https://rosettacode.org/wiki/Tonelli-Shanks_algorithm#Python
    s = valuation(m-1, 2)
    q = (m-1)//2**s

    z = 2
    while kronecker(z, m) != -1:
        z += 1

    c = mod(z, m)**q
    n = mod(n, m)
    r = n ** ((q+1)//2)
    t = n ** q

    e = s
    t2 = 0
    while (t-1) != 0:
        t2 = t * t
        for i in range(1, m):
            if t2 - 1 == 0:
                break
            t2 = t2 * t2
        b = c ** (1 << (e - i - 1))
        r = r * b
        c = b * b
        t = t * c
        e = i
    return ZZ(r)

def GeneralizedRepresentInteger(M, om, O):
    """
    Given an order O and a target norm M, compute an element gamma in O with
    norm M [Algorithm 3.12]
    Input:
    - M: odd integer > p
    - om: quaternion such that q = -om^2 is a positive integer and q = 1 mod 4
    - O: a special extremal maximal order containing Z[om] + jZ[om] as
      suborder, with j from the standard basis of B_p,inf
    - isogeny_cond: boolean flag only relevant for -om^2 = 1, enforces a
      modular condition needed for FixedDegreeIsogeny
    Output:
    - gamma in O with nrd(gamma) = M, or False
    """
    q = ZZ(-om**2)
    bound = ceil(4*M / (params.p*sqrt(q)))
    counter = 0

    m = isqrt(4*M/params.p - q)
    while counter < bound:
        counter += 1
        z = randint(1, m)
        m1 = isqrt((4*M - params.p*z**2) / (q*params.p))
        t = randint(-m1, m1)

        M1 = ZZ(4*M - params.p*(z**2 + q*t**2))
        if is_prime(M1):
            res = Cornacchia(q, M1)
            if not res:
                continue
            x, y = res

            j = params.B.gens()[1]
            gamma = x + om*y + j*z + om*j*t
            d = 1
            while True:
                if gamma / (d+1) in O:
                    d += 1
                else:
                    break

            if d != 2:
                continue
            return gamma / d

    return False

def Cornacchia(q, m):
    """
    Solve the equation x^2 + q*y^2 = m [Algorithm 3.11]
    Input:
    - q, m: integers with q prime and 0 <= q <= m
    Output:
    - x, y such that x^2 + qy^2 = m or False if a solution do not exist
    """
    if not kronecker(-q, m):
        return False
    if m == 2:
        if q == 1:
            return 1, 1
        return False

    r = ModularSQRT(-q % m, m)
    s = m

    while r ** 2 > m:
        r, s = s % r, r

    x, y = r, isqrt( (m - r**2) / q)
    if x**2 + q*y**2 == m:
        return x, y
    return False

def RandomIdealGivenNorm(N, prime):
    """
    Compute a random ideal from given norm [Alg 3.10]. It may fail if N is not
    prime.
    Input:
    - N: positive integer coprime with p, norm of the output ideal
    - prime: a boolean indicating whether N is prime
    Output:
    - a random left ideal J of O0 of norm N, or False if not found
    """
    if prime:
        while True:
            g1, g2, g3 = [randint(0, N-1) for _ in range(3)]
            gamma = params.B([0, g1, g2, g3])
            n_gamma = gamma.reduced_norm()
            if kronecker(-n_gamma, N) == 1:
                gamma += ModularSQRT(-n_gamma, N)
                break
    else:
        gamma = GeneralizedRepresentInteger(
                params.QUAT_prime_cofactor*N,
                params.B.gens()[0], params.O0)
        if not gamma:
            return False

    # Scaling by a random beta mod N
    while True:
        x, y, z, w = [randint(0, N-1) for _ in range(4)]
        beta = params.B([x, y, z, w])
        if gcd(beta.reduced_norm(), N) == 1:
            break
    return params.O0 * (gamma * beta) + params.O0 * N

def LLLReducedBasis(I):
    """
    Given an ideal I, computes an LLL-reduced basis of I. Calls pari on the
    Gram matrix of I.
    Input:
    - I: an ideal
    Output:
    - 4 elements of I forming an LLL-reduced basis
    """
    B = I.basis()
    M = []
    for a in B:
        M.append([QQ(2)*(a*b.conjugate()).reduced_trace() for b in B])
    G = Matrix(QQ, M)
    U = G.LLL_gram().transpose()
    return [sum(c*beta for c, beta in zip(row, B)) for row in U]

def RandomEquivalentPrimeIdeal(I):
    """
    Find an equivalent ideal to I with prime and bounded norm [Alg 3.9].
    Input:
    - I: a left O-ideal
    Output:
    - J ~ I of small prime norm, or raise an exception
    """
    B = LLLReducedBasis(I)

    cnt = 0
    b = params.QUAT_equiv_bound_coeff
    N = I.norm()
    while cnt < (2 * b + 1)**4:
        cnt += 1
        c_i = [randint(-b, b) for _ in range(4)]
        beta = sum([c_i[i] * B[i] for i in range(4)])
        if ZZ(beta.reduced_norm()/N).is_prime():
            J = I * beta.conjugate() * (1 / N)
            assert J.norm() == beta.reduced_norm()/N
            assert J.is_left_equivalent(I)
            return J

    raise RuntimeError(f'RandomEquivalentPrimeIdeal failed on {I = }')

def ReducedIdeal(I):
    """
    Given an ideal I compute the ideal equivalent to I with the smallest norm
    Input:
    - I: a left O0-ideal
    Ouptut:
    - J: the ideal equivalent to I with smallest norm
    - beta: the element in I generating J
    """
    B = LLLReducedBasis(I)
    beta = B[0]
    J = I * (beta.conjugate() / I.norm())
    assert J.is_left_equivalent(I)
    return J, beta

def SmallishGenerator(I, N, I_basis):
    """
    Compute a rather small random generator of I
    Input:
    - I: target left O0-ideal
    - N: norm of I
    - I_basis: an LLL-reduced basis of I
    Output:
    - alpha: a generator of I
    """
    while True:
        alpha = sum(randint(1,100000)*gen for gen in I_basis)

        a_alpha = alpha.coefficient_tuple()[0]
        if gcd(2*a_alpha, N) == 1 and gcd(alpha.reduced_norm(), N*N) == N:
            break

    return alpha

def _succ_min(L):
    """
    Compute the first two minima of L
    Input:
    - L: a matrix
    Output:
    - lam1, lam2: the two minima of L
    """
    fourth_root_p = round(params.p**(1/4), 10)
    lam1 = round(L.row(0).norm()/fourth_root_p, 10)
    lam2 = round(L.row(1).norm()/fourth_root_p, 10)
    return lam1, lam2

def KernelDecomposedToIdeal(c1, c2):
    """
    Given c1, c2 defining a point K = c1*P0 + c2*Q0 in E0[2^f] return the ideal
    corresponding to the isogeny with kernel K [Algorithm 3.17]
    Input:
    - c1, c2: integers defining a point K = c1*P0 + c2*Q0 in E0[2^f]
    Output:
    - I: a left O0-ideal corresponding to the isogeny of kernel K
    """
    v = vector(Zmod(2**params.f), (c1, c2))

    # theta = j + (1+k)/2 (and j = -i + 2*((i+j)/2))
    M_theta = -params.mat_i + 2*params.mat_ij2 + params.mat_1k2
    d1, d2 = M_theta.transpose() * v # c1*theta(P0) + c2*theta(Q0) = theta(P)

    M = Matrix(Zmod(2**params.f), 2, [c1, d1, c2, d2])

    a, b = M**-1 * params.mat_i.transpose() * v
    a, b = ZZ(a), ZZ(b)

    # Now aP + b*thetaP == eta*P
    # Sanity:
    # import ec
    # M = a*params.mat_1+(-1-b)*params.mat_i+2*b*params.mat_ij2+b*params.mat_1k2
    # X = ec.EvalMatrix(M)
    # P = c1*X[0] + c2*X[1]
    # assert P == 0

    ker = params.B((a + b/2, -1, b, b/2))
    I = params.O0 * ker + params.O0 * (2**params.f)
    assert I.norm() == 2**params.f
    return I

def Pushforward(I, J):
    """
    Compute the pushforward of I under J, i.e. [J]_* I.
    Input:
    - I, J: two left O-ideals with coprime norm
    Output:
    - the ideal [J]_*I
    """
    assert I.left_order() == J.left_order()
    assert gcd(I.norm(), J.norm()) == 1

    O = J.right_order()
    return J.conjugate()*I + O*I.norm()

def RandomEquivalentQuaternion(I, n_sk_comm):
    """
    Return a uniform element in the (ideal) lattice I.
    Input:
    - I: an ideal
    - n_sk_comm: n(I_sk) * n(I_comm)
    Output:
    - b: uniformly sampled in I such that nrd(b) < D_rsp*n(I_sk)*n(I_comm)*2^f
    """
    B = I.basis()
    M = []
    for a in B:
        M.append([QQ(2)*(a*b.conjugate()).reduced_trace() for b in B])
    G = Matrix(QQ, M)
    U = G.LLL_gram().transpose()
    B = [sum(c*beta for c, beta in zip(row, B)) for row in U]

    # G_inv = G**-1
    # U_star = G_inv.LLL_gram()
    # H = U_star.transpose() * G_inv * U_star

    # I.norm() == 2**params.f * n_sk_comm; params.D_rsp is the target norm,
    # which is always very close to Minkowski
    bound = params.D_rsp*n_sk_comm*2**(params.f)
    # Bi = [isqrt(bound * G_star[i,i]) for i in range(4)]
    Bi = [5 for i in range(4)]
    while True:
        xi = vector(QQ, [randint(-Bi[i], Bi[i]) for i in range(4)])
        v = sum(xi[i] * B[i] for i in range(4))
        if 0 < v.reduced_norm() < bound:
            return v

def ComputeBacktrackingAndNormalize(alpha):
    """
    Remove backtracking part [Algorithm 4.4]
    Input:
    - alpha: a quaternion element
    Output:
    - alpha: the non-backtracking part of the input
    - n: the 2-valuation of the backtracking coefficient
    """
    a0, a1, a2, a3 = alpha
    x0, x1, x2, x3 = a0-a3, a1-a2, a2*2, a3*2
    g = gcd([x0, x1, x2, x3])
    n = valuation(g, 2)
    return alpha/g, n

