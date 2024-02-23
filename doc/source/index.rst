.. APES documentation master file, created by
   sphinx-quickstart on Fri Jul 22 10:51:46 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to APES's documentation!
================================

APES? What is it?
-----------------

**APES** is the **A**\ ll- **P**\ urpose **E**\ quilibrium **S**\ olver. As such, 
**APES** aims to provide the tools necessary to solve general equilibrium problems.

.. image:: monkey.png
    :scale: 25%
    :alt: Monkey with wrench
    :align: center

What is it for?
---------------

**WARNING** : **APES** is experimental software. Bugs, and errors in the
final values are to be expected.

Describing the System
---------------------

APES handled the stoichiometry coefficients in the form of a 2D array where
each column stands for the independent components and each row stands for
the complexes. For instance, consider the following system where L stands
for a ligand, H stands for a proton and L has two protonable groups. We would
write this as

.. math::

    \begin{eqnarray}
     \rm L + H &\leftrightarrow \rm HL \\
     \rm L + 2H &\leftrightarrow \rm H_2L
    \end{eqnarray}

Now we have two :term:`independent components` which are L and H and 
three :term:`complex`\ es which would be HL, H\ :sub:`2`\ L and, provided that 
this equilibrium takes place in an aqueous environment we must not forget
the OH\ :sup:`-`. The resulting :term:`stoichiometry array` would be

.. math::

    \left(
     \begin{array}{rr}
      1 &  1 \\
      1 &  2 \\
      0 & -1 \\
     \end{array}
    \right)

Contents
========

.. toctree::
   :maxdepth: 2
    
   quickstart
   modules

News
====
* Still in development

Glossary
========

.. glossary::

    complex
        A species obtained by combinantion of independent components

    equilibrium constants array
        An 1D array with the values of the equilibrium constants in logarithmic
        units and in the same order than in the rows of the
        :term:`stoichiometry array`.

    free concentrations array
        An Nx(S+E) array that contains the values of the free concentrations

    independent components
        It is each one of the species with which the complexes are made

    model
        A :term:`stoichiometry array` along with an 
        :term:`equilibrium constants array` that describe the system.

    Nernst's equation
        Relation between the potential and the concentration of electroactive
        species. :math:`E=E_0+\frac{RT}{nF}\log[H]`

    standard potential
        The intercept in :term:`Nernst's equation`

    stoichiometry array
        A 2D array of integers which describes the stoichiometry of the system

    total concentrations array
        An 2D array of floats that represents the analytical concentrations
        for the independent components. The numbers can be negative insofar
        the mass balance is mathematically meaningful. Experimental points in
        rows and species in columns.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

