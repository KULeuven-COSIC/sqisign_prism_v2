# SQIsign / PRISM v2 in Sage

SageMath implementation of SQIsign and PRISM. SQIsign is implemented following
the NIST specification [version
2.0](https://sqisign.org/spec/sqisign-20250205.pdf), with the only exception of
the improved ideal to isogeny subroutine from
[qlapoti](https://eprint.iacr.org/2025/1604). PRISM is implemented using the
same underlying functions.

This implementation is intended as a proof of concept, and hence sometimes
clarity is favored over efficiency. Most notabily, the elliptic curve and
isogeny arithmetic is done using SageMath functions and not x-only arithmetic.
As a result, the 1D isogeny computation in SQIsign takes a significant portion
of the overall computation, unlike it would be in an optimized implementation.

## Usage - SQIsign

To import and use `SQIsign` in `sage` do:
```python
from sqisign import SQIsign, SQIsign_verify
import params

lvl = 1 # 3, 5
params.set_sqi_params(lvl)

alice = SQIsign()

msg = 'Hello world'
sigma = alice.sign(msg)

out = SQIsign_verify(msg, sigma, alice.pk)
print(f'Verification: {out}')
```

The file `sqisign.py` runs the code above and prints some additional timings.
It can be run with `sage --python -O sqisign.py [level]`; if not provided,
`level` is set to `1`.

Similarly a basic benchmarking script can be run with `sage --python -O
bench_sqi.py [level]`. The average timings for key generation, signing and
verification are computed from 10 runs.

## Usage - PRISM

To import and use `PRISM` in `sage` do:
```python
from prism import PRISM, PRISM_verify
import params

lvl = 1 # 3, 5
params.set_prism_params(lvl)

alice = PRISM()

msg = 'Hello world'
sigma = alice.sign(msg)

out = PRISM_verify(msg, sigma, alice.pk)
print(f'Verification: {out}')
```

The file `prism.py` runs the code above and prints some additional timings.  It
can be run with `sage --python -O prism.py [level]`; if not provided, `level`
is set to `1`.

Similarly a basic benchmarking script can be run with `sage --python -O
bench_prism.py [level]`. The average timings for key generation, signing and
verification are computed from 10 runs.

## Project structure

### Main code

- `sqisign.py`: main implementation of SQIsign
- `prism.py`: main implementation of PRISM

### Libraries and functions

- `constants.py`: constants for different security levels
- `ec.py`: elliptic curve helpers
- `hd.py`: wrappers for 2D isogenies
- `misc.py`: encodings
- `params.py`: parameters for different security levels
- `qlapoti.py`: ideal to isogeny algorithm from
  [qlapoti](https://github.com/KULeuven-COSIC/Qlapoti)
- `quaternions.py`: quaternion function, adapted from
  [quaternion_helpers](https://github.com/Jonathke/quaternion_helpers)
- `sos.py`: sum of squares algorithm from
  [qt-pegasis](https://github.com/KULeuven-COSIC/qt-pegasis)
- `theta/`: 2D isogenies in theta coordinates from
  [two-isogenies](https://github.com/ThetaIsogenies/two-isogenies/tree/main/Theta-SageMath)

### Benchmarking

- `bench_sqi.py`: benchmarking for SQIsign
- `bench_prism.py`: benchmarking for PRISM
