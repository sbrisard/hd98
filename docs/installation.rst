************
Installation
************


Installing the C++ library
==========================

This is a CMake_ based project. The dependencies are: Catch2_ for the tests, and
PyBind11_ for the Python bindings. The installation procedure is
standard. First, clone the repository. Then, ``cd`` into the root directory of
the hd98 project. Let ``hd98_INSTALL_PREFIX`` be the path to the directory where
hd98 should be installed::

  $ git clone https://github.com/sbrisard/hd98
  $ cd hd98
  $ mkdir build
  $ cd build
  $ cmake -DCMAKE_INSTALL_PREFIX=hd98_INSTALL_PREFIX ..
  $ cmake --build . --config Release
  $ cmake --install . --config Release

.. note:: The ``--config`` option might not be available, depending on the
   selected generator.

At this point, hd98 should be installed. You can now run the tests::

  $ ctest . -C Release

.. note:: Depending on the system, you might need to add ``hd98_INSTALL_PREFIX``
   to your ``PATH`` environment variable.


Compiling your first hd98 program
=================================

``cd`` into the ``example`` subdirectory. The provided example program should be
compiled and linked against hd98::

  $ mkdir build
  $ cd build
  $ cmake -Dhd98_DIR=hd98_INSTALL_PREFIX/lib/cmake/hd98 ..
  $ cmake --build . --config Release

An executable called ``example_hd98`` should be present in the ``build/Release``
subdirectory.


Building the documentation
==========================

The documentation of hd98 requires Sphinx_. The C++ API docs are built with
Doxygen_ and the Breathe_ extension to Sphinx_.

To build the HTML version of the docs in the ``public`` subdirectory::

  $ cd docs
  $ sphinx-build -b html . ../public

To build the LaTeX version of the docs::

  $ cd docs
  $ make latex


Installing the Python bindings
==============================

To install the hd98 module, ``cd`` into the ``python`` subdirectory and edit the
``setup.cfg`` file. Set the ``include_dir`` to the appropriate paths. These
should be::

  [hd98]
  include_dir = ${CMAKE_INSTALL_PREFIX}/include

Then, issue the following command::

  $ python setup.py install --user

or (if you intend to edit the project)::

  $ python setup.py develop --user

To run the tests with Pytest_::

  $ python -m pytest tests

.. _Breathe: https://breathe.readthedocs.io/
.. _Catch2: https://github.com/catchorg/Catch2
.. _CMake: https://cmake.org/
.. _Doxygen: https://www.doxygen.nl/
.. _PyBind11: https://github.com/pybind/pybind11
.. _Pytest: https://docs.pytest.org/
.. _Sphinx: https://www.sphinx-doc.org/
