"""Basic definitions."""

import ctypes

from ctypes import c_double, c_char_p, c_int, c_size_t, c_void_p

c_double_p = ctypes.POINTER(c_double)
c_size_t_p = ctypes.POINTER(c_size_t)

path = "/home/sbrisard/.local/opt/hd98/lib/libhd98.so"
hd98 = ctypes.cdll.LoadLibrary(path)


class HD98_MaterialType(ctypes.Structure):
    _fields_ = [("name", c_char_p), ("num_int_var", c_size_t), ("free", c_void_p), ("current_state", c_void_p), ("update", c_void_p)]


HD98_MaterialType_p = ctypes.POINTER(HD98_MaterialType)


class HD98_Material(ctypes.Structure):
    _fields_ = [("type", HD98_MaterialType_p), ("data", c_void_p)]


HD98_Material_p = ctypes.POINTER(HD98_Material)

hd98_material_free_t = ctypes.CFUNCTYPE(None, HD98_Material_p)
hd98_material_current_state_t = ctypes.CFUNCTYPE(None, HD98_Material_p, c_double_p, c_double_p, c_double_p)
hd98_material_update_t = ctypes.CFUNCTYPE(
    None,
    HD98_Material_p,
    c_double_p,
    c_double_p,
    c_double_p,
    c_double_p,
    c_double_p,
    c_double_p,
)


class Material:
    def __init__(self, material_p, data_p_t):
        self._material_p = material_p
        self._addr = ctypes.cast(material_p, c_void_p).value
        self._free = ctypes.cast(
            material_p.contents.type.contents.free, hd98_material_free_t
        )
        self._update = ctypes.cast(
            material_p.contents.type.contents.update, hd98_material_update_t
        )
        self._data_p = ctypes.cast(material_p.contents.data, data_p_t)

    def __del__(self):
        self._free(self._material_p)

    def update(self, delta_eps, eps1, omega1, sig2, omega2, C2):
        self._update(
            self._material_p,
            delta_eps.ctypes.data_as(c_double_p),
            eps1.ctypes.data_as(c_double_p),
            omega1.ctypes.data_as(c_double_p),
            sig2.ctypes.data_as(c_double_p),
            omega2.ctypes.data_as(c_double_p),
            C2.ctypes.data_as(c_double_p),
        )

    @property
    def name(self):
        return self._material_p.contents.type.contents.name

    @property
    def num_int_var(self):
        return self._material_p.contents.type.contents.num_int_var


hd98.hd98_global_update.argtypes = [
    c_size_t,
    c_size_t_p,
    c_void_p,
    c_double_p,
    c_double_p,
    c_double_p,
    c_double_p,
    c_double_p,
    c_double_p,
]
hd98.hd98_global_update.restype = None


def global_update(phase, mat, delta_eps, eps1, omega1, sig2, omega2, C2):
    hd98.hd98_global_update(
        phase.size,
        phase.ctypes.data_as(c_size_t_p),
        mat.ctypes.data_as(c_void_p),
        delta_eps.ctypes.data_as(c_double_p),
        eps1.ctypes.data_as(c_double_p),
        omega1.ctypes.data_as(c_double_p),
        sig2.ctypes.data_as(c_double_p),
        omega2.ctypes.data_as(c_double_p),
        C2.ctypes.data_as(c_double_p),
    )


hd98.hd98_solve_polarizations_plus.argtypes = (
    c_size_t,
    c_size_t_p,
    c_void_p,
    c_double,
    c_double,
) + 4 * (c_double_p,)
hd98.hd98_solve_polarizations_plus.restype = c_int


def solve_polarizations_plus(
    phase, mat, lambda0, mu0, delta_tau, eps1, omega1, delta_eps
):
    err = hd98.hd98_solve_polarizations_plus(
        phase.size,
        phase.ctypes.data_as(c_size_t_p),
        mat.ctypes.data_as(c_void_p),
        lambda0,
        mu0,
        delta_tau.ctypes.data_as(c_double_p),
        eps1.ctypes.data_as(c_double_p),
        omega1.ctypes.data_as(c_double_p),
        delta_eps.ctypes.data_as(c_double_p),
    )
    if err != 0:
        raise RuntimeError(err)
