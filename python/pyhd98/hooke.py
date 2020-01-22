"""Definition of the Hooke material."""

import ctypes
import numpy as np
import pyhd98.core

from ctypes import c_double
from pyhd98.core import c_double_p, HD98_Material_p, Material


class HD98_HookeData(ctypes.Structure):
    _fields_ = [("lambda_", c_double), ("mu", c_double), ("C", c_double_p)]


HD98_HookeData_p = ctypes.POINTER(HD98_HookeData)

hd98_hooke_new = pyhd98.core.hd98.hd98_hooke_new
hd98_hooke_new.argtypes = 2 * [c_double]
hd98_hooke_new.restype = HD98_Material_p


class Hooke(Material):
    def __init__(self, lambda_, mu):
        super().__init__(hd98_hooke_new(lambda_, mu), HD98_HookeData_p)

    @property
    def lambda_(self):
        return self._data_p.contents.lambda_

    @property
    def mu(self):
        return self._data_p.contents.mu
