#!/usr/bin/python3
# -*- encoding: utf-8 -*-

from setuptools import setup

setup(
    name="apes",
    version="dev",
    description="The All-Purpose Equilibrium Solver",
    author="Salvador Blasco",
    author_email="salvador.blasco@gmail.com",
    keywords="Chemistry Equilibria",
    # packages = find_packages(),
    scripts=['supyquad.py', 'apes.py'],
    install_requires=['numpy', 'scipy', 'PyQt5', 'pyqtgraph'],
    setup_requires=['numpy', 'scipy', 'PyQt5', 'pyqtgraph'],
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 3 - Alpha',

        # Indicate who your project is intended for
        'Intended Audience :: Developers',
        'Topic :: Science / Engineering :: Chemistry',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: MIT License',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 3.4',
        ]
)
