import numpy as np
import matplotlib.pyplot as plt

from pyhd98 import HalmDragon1998


def print_underlined(text, symbol):
    print("")
    print(text)
    print(symbol * len(text))


def print_h2(text):
    print_underlined(text, "=")


def print_h3(text):
    print_underlined(text, "-")


if __name__ == "__main__":
    dim, sym = 3, 6
    κ, μ = 60700.0, 31300.0
    α, β = 16000.0, 31000.0
    k0, k1 = 0.11, 2.2
    mat = HalmDragon1998(κ - 2 * μ / 3, μ, α, β, k0, k1)

    I2 = np.zeros((sym,), dtype=np.float64)
    I2[:dim] = 1
    I4 = np.eye(sym, dtype=np.float64)
    J4 = 1 / dim * I2[None, :] * I2[:, None]
    K4 = I4 - J4

    C = dim * κ * J4 + 2 * μ * K4
    H = 2 * (mat.alpha * I2[None, :] * I2[:, None] + 2 * mat.beta * I4)

    # Reference data (uniaxial, strain-diven experiment)
    ε_dot = np.zeros((sym,), dtype=np.float64)
    ε_dot[0] = 1.0

    t_d = np.sqrt(2 * mat.k0_sqrt2 / np.dot(ε_dot, H @ ε_dot))

    num_load_steps = 21
    t = np.linspace(0, 2 * t_d, num=num_load_steps)

    ε_ref = t[:, None] * ε_dot
    ω_ref = mat.k0_sqrt2 / mat.k1_sqrt2 * ((t / t_d) ** 2 - 1)
    σ_ref = t[:, None] * ((C - ω_ref[:, None, None] * H) @ ε_dot)

    # Stress-driven experiment
    Δt = np.diff(t)
    Δσ = np.diff(σ_ref, axis=0)
    ε = np.zeros_like(ε_ref)
    ω = np.zeros_like(ω_ref)

    max_iter = 20
    rtol = 1e-15
    atol = 1e-15

    for load_step in range(1, num_load_steps):
        print_h2(f"Load-step #{load_step}")
        print_h3("Target values")
        print(f"- σ_xx = {σ_ref[load_step, 0]}")
        print(f"- σ_yy = σ_zz = {σ_ref[load_step, 1]}")
        print(f"- ε_xx = {ε_ref[load_step, 0]}")
        ε1 = ε[load_step - 1]
        ω1 = np.zeros((mat.num_int_var,), dtype=np.float64)
        ω1[0] = ω[load_step - 1]
        ω2 = np.zeros_like(ω1)
        Δε = np.zeros((sym,), dtype=np.float64)
        σ2 = np.zeros_like(Δε)
        C2 = np.zeros((sym, sym), dtype=np.float64)

        print_h3("Iterations")

        res = 0.
        for i in range(max_iter):
            mat.update(Δε, ε1, ω1, σ2, ω2, C2)
            ρ = σ_ref[load_step] - σ2
            old_res = res
            res = np.linalg.norm(ρ)
            print(f"- Iteration #{i}, res = {res}, rate = {res / old_res**2}")
            tol = rtol * np.linalg.norm(σ_ref[load_step]) + atol
            if res < tol:
                break
            δε = np.linalg.solve(C2, ρ)
            Δε += δε
        else:
            raise RuntimeError("no convergence")

        ε[load_step] = ε[load_step - 1] + Δε
        print_h3("Results")
        print(f"- expected ε_xx = {ε_ref[load_step, 0]}")
        print(f"- actual ε_xx   = {ε[load_step, 0]}")
