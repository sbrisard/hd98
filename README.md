# hd98

This C-library implements the simple damage model of Halm and Dragon
(1998). It provides functins to compute the stress and damage
increments induced by a prescribed strain increment. It also exposes a
function that performs this update over a whole grid of material
points (for implemetation within an FFT simulation).

Halm, D., & Dragon, A. (1998). An anisotropic model of damage and frictional sliding for brittle materials. European Journal of Mechanics - A/Solids, 17(3), 439–460. <https://doi.org/10.1016/S0997-7538(98)80054-5>

## Installation

Two components need to be installed: the C library itself, and the
python bindings (wich are based on ctypes, since performance is not
critical – only one function is called to update the whole grid).
