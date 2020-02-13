"""Basic definitions."""

import ctypes

from ctypes import c_double, c_char_p, c_void_p

c_double_p = ctypes.POINTER(c_double)

path = "/home/sbrisard/.local/opt/hd98/lib/libhd98.so"
hd98 = ctypes.cdll.LoadLibrary(path)


class HD98_MaterialType(ctypes.Structure):
    _fields_ = [("name", c_char_p), ("free", c_void_p), ("update", c_void_p)]


HD98_MaterialType_p = ctypes.POINTER(HD98_MaterialType)


class HD98_Material(ctypes.Structure):
    _fields_ = [("type", HD98_MaterialType_p), ("data", c_void_p)]


HD98_Material_p = ctypes.POINTER(HD98_Material)

hd98_material_free_t = ctypes.CFUNCTYPE(None, HD98_Material_p)
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
