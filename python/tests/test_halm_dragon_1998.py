import numpy as np
import pytest

from numpy.testing import assert_allclose

from pyhd98 import HalmDragon1998


def test_halm_dragon_1998_data():
    κ = 60700.0
    μ = 31300.0
    α = 16000.0
    β = 31000.0
    k0 = 0.11
    k1 = 2.2
    λ = κ - 2 * μ / 3
    mat = HalmDragon1998(λ, μ, α, β, k0, k1, 0)
    assert mat.lambda_ == λ
    assert mat.mu == μ
    assert mat.alpha == α
    assert mat.beta == β
    assert mat.k0_sqrt2 == k0 * np.sqrt(2)
    assert mat.k1_sqrt2 == k1 * np.sqrt(2)
    assert mat.stiffness_type == 0


@pytest.mark.parametrize(
    "ε_dot", [[1.0, 1.0, 1.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]]
)
def test_halm_dragon_1998_proportional_strain(ε_dot, atol=1e-15, rtol=1e-15):
    ε_dot = np.asarray(ε_dot)
    dim, sym = 3, 6
    κ, μ = 60700.0, 31300.0
    α, β = 16000.0, 31000.0
    k0, k1 = 0.11, 2.2
    mat = HalmDragon1998(κ - 2 * μ / 3, μ, α, β, k0, k1, 0)

    I2 = np.zeros((sym,), dtype=np.float64)
    I2[:dim] = 1.0

    tr_ε_dot = ε_dot[:dim].sum()
    C_ε_dot = mat.lambda_ * tr_ε_dot * I2 + 2 * mat.mu * ε_dot
    H_ε_dot = 2 * mat.alpha * tr_ε_dot * I2 + 4 * mat.beta * ε_dot
    ε_dot_H_ε_dot = ε_dot.dot(H_ε_dot)

    t_damage = np.sqrt(2 * mat.k0_sqrt2 / ε_dot_H_ε_dot)
    t_max = 2 * t_damage
    num_increments = 10
    num_steps = num_increments * (num_increments + 1)

    Δt = t_max / num_increments
    Δξ = Δt / t_damage

    sign = np.empty((num_steps,), dtype=np.int8)
    ω_exp = np.empty((num_steps,), dtype=np.float64)
    step = 0
    for increment in range(1, num_increments + 1):
        for i in range(1, increment + 1):
            sign[step] = 1
            if i == increment:
                ξ = increment * Δξ
                ω_exp[step] = (
                    0.0 if ξ <= 1.0 else ((ξ ** 2 - 1.0) * mat.k0_sqrt2 / mat.k1_sqrt2)
                )
            else:
                ω_exp[step] = ω_exp[step - 1]
            step += 1
        for i in range(1, increment + 1):
            sign[step] = -1
            ω_exp[step] = ω_exp[step - 1]
            step += 1

    t = 0.0
    ω = np.zeros((1,), dtype=np.float64)
    ω_new = np.zeros_like(ω)
    ε = np.zeros((sym,), dtype=np.float64)
    σ = np.empty_like(ε)
    σ_exp = np.empty_like(ε)
    C_new = np.empty((sym, sym), dtype=np.float64)
    for step in range(num_steps):
        t += sign[step] * Δt
        Δε = sign[step] * Δt * ε_dot
        mat.update(Δε, ε, ω, σ, ω_new, C_new)
        ε += Δε
        ω[0] = ω_new[0]
        σ_exp = t * (C_ε_dot - ω_exp[step] * H_ε_dot)
        assert_allclose(σ, σ_exp, rtol=rtol, atol=atol)
