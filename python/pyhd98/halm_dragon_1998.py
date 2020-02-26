"""Definition of the material model of Halm and Dragon (1998)."""

import ctypes
import numpy as np
import pyhd98.core

from ctypes import c_double, c_int
from pyhd98.core import c_double_p, HD98_Material_p, Material


class HD98_HalmDragon1998Data(ctypes.Structure):
    _fields_ = [
        ("lambda_", c_double),
        ("mu", c_double),
        ("alpha", c_double),
        ("beta", c_double),
        ("k0_sqrt2", c_double),
        ("k1_sqrt2", c_double),
        ("stiffness_type", c_int),
    ]


HD98_HalmDragon1998Data_p = ctypes.POINTER(HD98_HalmDragon1998Data)

hd98_halm_dragon_1998_new = pyhd98.core.hd98.hd98_halm_dragon_1998_new
hd98_halm_dragon_1998_new.argtypes = 6 * [c_double]+[c_int]
hd98_halm_dragon_1998_new.restype = HD98_Material_p


class HalmDragon1998(Material):
    def __init__(self, lambda_, mu, alpha, beta, k0, k1, stiffness_type=0):
        super().__init__(
            hd98_halm_dragon_1998_new(lambda_, mu, alpha, beta, k0, k1, stiffness_type),
            HD98_HalmDragon1998Data_p,
        )

    @property
    def lambda_(self):
        return self._data_p.contents.lambda_

    @property
    def mu(self):
        return self._data_p.contents.mu

    @property
    def alpha(self):
        return self._data_p.contents.alpha

    @property
    def beta(self):
        return self._data_p.contents.beta

    @property
    def k0_sqrt2(self):
        return self._data_p.contents.k0_sqrt2

    @property
    def k1_sqrt2(self):
        return self._data_p.contents.k1_sqrt2

    @property
    def stiffness_type(self):
        return self._data_p.contents.stiffness_type
