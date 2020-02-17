import numpy as np

from numpy.testing import assert_array_equal

from pyhd98 import Hooke


def test_hooke_data():
    lambda_ = 1.2
    mu = 3.4
    hooke = Hooke(lambda_, mu)
    assert hooke.name == b"Hooke"
    assert hooke.num_int_var == 0
    assert hooke.lambda_ == lambda_
    assert hooke.mu == mu


def test_hooke_update():
    lambda_ = 1.2
    mu = 3.4
    hooke = Hooke(lambda_, mu)

    dim = 3
    sym = 6
    I = np.zeros((sym,), dtype=np.float64)
    I[:dim] = 1.0

    eps1 = np.zeros_like(I)
    delta_eps = np.array([1.2, -3.4, 5.6, -7.8, 9.1, -2.3], dtype=np.float64)

    sig2_act = np.zeros_like(I)
    C2_act = np.zeros((sym, sym), dtype=np.float64)
    hooke.update(delta_eps, eps1, np.array(0.0), sig2_act, np.array(0.0), C2_act)

    sig2_exp = np.zeros_like(I)
    tr_eps = delta_eps[:dim].sum()
    sig2_exp = lambda_ * tr_eps * I + 2 * mu * delta_eps

    C2_exp = lambda_ * I[:, None] * I[None, :] + 2 * mu * np.eye(sym)

    assert_array_equal(sig2_act, sig2_exp)
    assert_array_equal(C2_act, C2_exp)
