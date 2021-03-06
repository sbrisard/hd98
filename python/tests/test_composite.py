import numpy as np

from numpy.testing import assert_allclose
from pyhd98 import Hooke, HalmDragon1998, Composite


def create_default_hooke():
    dim = 3
    κ1, μ1 = 76700.0, 41600.0
    return Hooke(κ1 - 2 * μ1 / dim, μ1)


def create_default_halm_dragon_1998():
    dim = 3
    κ, μ = 60700.0, 31300.0
    α, β = 16000.0, 31000.0
    k0, k1 = 0.11, 2.2
    return HalmDragon1998(κ - 2 * μ / dim, μ, α, β, k0, k1, 0)


def test_composite_update(atol=1e-15, rtol=1e-15):
    np.random.seed(20200213)
    dim, sym = 3, 6
    mat0 = create_default_hooke()
    mat1 = create_default_halm_dragon_1998()
    composite = Composite(mat0, mat1)
    n0, n1 = 5, 2
    n = n0 * n1
    phase = np.remainder(np.arange(n, dtype=np.uintp), 2).reshape((n0, n1))

    ω1 = 0.4 * np.random.rand(n0, n1)
    ω2_act = np.empty_like(ω1)

    ε1 = 1e-3 * np.random.rand(n0, n1, sym)
    Δε = 1e-4 * np.random.rand(n0, n1, sym)

    σ2_act = np.empty_like(ε1, dtype=np.float64)
    C2_act = np.empty((n0, n1, sym, sym), dtype=np.float64)

    composite.update(phase, Δε, ε1, ω1, σ2_act, ω2_act, C2_act)

    σ2_exp = np.empty((sym,), dtype=np.float64)
    # We need to create ω2 so as to be modifiable in place
    ω2_exp = np.empty((1,), dtype=np.float64)
    C2_exp = np.empty((sym, sym), dtype=np.float64)

    for i0 in range(n0):
        for i1 in range(n1):
            mat = mat0 if phase[i0, i1] == 0 else mat1
            mat.update(
                Δε[i0, i1], ε1[i0, i1], ω1[i0, i1, None], σ2_exp, ω2_exp, C2_exp
            )
            assert_allclose(σ2_act[i0, i1], σ2_exp, atol=atol, rtol=rtol)
            if mat is mat1:
                # Assert damage value only if material is damageable!
                assert_allclose(ω2_act[i0, i1], ω2_exp, atol=atol, rtol=rtol)
            assert_allclose(C2_act[i0, i1], C2_exp, atol=atol, rtol=rtol)
