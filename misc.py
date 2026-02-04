try:
    from . import params
except ImportError:
    import params

def encode_fp(x):
    """
    Given an element x in Fp, return its bytes representation (little-endian)
    depending on the parameters.
    Input:
    - x: an element of Fp
    Output:
    - x_enc: the encoding of x in bytes
    """
    x = int(x)
    x_enc = x.to_bytes(params.FP_ENC_BYTES, 'little')
    return x_enc

def encode_fp2(x):
    """
    Given an element x in Fp2, return its bytes representation (little-endian)
    depending on the parameters.
    Input:
    - x: an element of Fp2
    Output:
    - x_enc: the encoding of x in bytes
    """
    a, b = x
    x_enc = encode_fp(a) + encode_fp(b)
    return x_enc

def encode_curve_a(E):
    """
    Given a curve E, return the Fp2 encoding of its Montgomery coefficient.
    Input:
    - E: a curve
    Output:
    - enc_a: the Fp2 encoding of the A coefficient of E
    """
    A = E.a_invariants()[1]
    return encode_fp2(A)

def encode_curve_j(E):
    """
    Given a curve E, return the Fp2 encoding of its j-invariant.
    Input:
    - E: a curve
    Output:
    - enc_j: the Fp2 encoding of j(E)
    """
    j = E.j_invariant()
    return encode_fp2(j)
