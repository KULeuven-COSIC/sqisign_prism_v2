import time
import sys

from sage.all import randint

from sqisign import SQIsign, SQIsign_verify
import params

if __name__ == "__main__":
    assert False, 'run benchmarks with `sage --python -O bench_sqi.py [level]` to skip assertions'

    lvl = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    params.set_sqi_params(lvl)
    print(f'Running SQIsign lvl {lvl}')

    n_runs = 10

    tot_times = [0, 0, 0] # [keygen, sign, verification]

    for i in range(n_runs):
        print(f'\t- run {i+1}/{n_runs}')

        _t0 = time.time()
        sqi = SQIsign()

        _t1 = time.time()
        msg = f'Hello world {randint(1, 100)}'
        sigma = sqi.sign(msg)

        _t2 = time.time()
        out = SQIsign_verify(msg, sigma, sqi.pk)
        assert out
        _t3 = time.time()

        tot_times[0] += _t1 - _t0
        tot_times[1] += _t2 - _t1
        tot_times[2] += _t3 - _t2

    names = ['Keygen', 'Sign', 'Verify']
    print(f'Done')
    for i in range(3):
        tt = tot_times[i] / n_runs
        print(f'{names[i]}:\t{tt:.3f}s')
    tt = (sum(tot_times)) / n_runs
    print(f'Total time:\t{tt:.3f}s')







