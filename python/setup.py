from setuptools import setup, find_packages

setup(
    author="Sébastien Brisard",
    author_email='',
    python_requires='>=3.5',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: ISC License (ISCL)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="Python bindings to the hd98 library",
    license="BSD-3",
    keywords='pyhd98',
    name='pyhd98',
    packages=["pyhd98"],
    test_suite="tests",
    url="https://github.com/sbrisard/hd98",
)
