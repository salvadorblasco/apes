'''
This module contains simple miscelaneous routines that are shared by other
functions or classes in other modules.
'''

import collections
import collections.abc
import functools
import re

import numpy as np


__version__ = '0.1'


def assert_array(*args):
    """Check that all arguments provided are arrays.

    Raises:
        :class:`TypeError` if any of the arguments provided is not an array
    """
    assert_type(np.ndarray, *args)


def assert_array_dim(ndim, *args):
    for i in args:
        if not isinstance(i, np.ndarray):
            raise TypeError("Argument is not an array", i)
        if i.ndim != ndim:
            raise ValueError("%s has wrong number of dimmensions" % i)


def assert_shape(shape, *args):
    """Assert that the arrays provided have the desired shape

    Parameters:
        shape (tuple): The shape desired
        *args (arrays): The arrays to be tested
    Raises:
        ValueError: if any of the array is misshaped
    Example:
        # >>> import numpy
        >>> a = numpy.random.rand(5,4)
        >>> assert_shape((5,4), a)
    """
    assert_array(*args)
    for i in args:
        if i.shape != shape:
            raise ValueError("Array misshaped: found " + str(i.shape) +
                             ", expected " + str(shape))


def assert_type(typ, *args):
    """Assert that the arguments provided are instanvce of a given object.

    Parameters:
        type (type): The type to be tested
        *args (objects): The objects to be tested

    Raises:
        TypeError: if any of the objects is not instance of the intended class
    """
    for i in args:
        if not isinstance(i, typ):
            raise TypeError("%s is not %s" % (i, typ))


def assert_same_len(*args):
    """Assert that a number of objects have all the same length.

    Parameters:
        *args: The objects to be tested. They must have __len__

    Raises:
        ValueError: if not all the objects have the same length
    """
    for i in args:
        if '__len__' not in dir(i):
            raise ValueError("Object %s does not have __len__()" % str(i))
    if not all(len(i) == len(args[0]) for i in args):
        raise ValueError("arrays don't have the same length")


def assert_same_shape(*args):
    """Assert that the arrays provided have the same shape.

    Parameters:
        *args (arrays): The arrays to be tested

    Raises:
        ValueError: if not all the arrays have the same shape
    """
    assert_array(*args)
    if len(args) < 2:
        raise ValueError("At least two arguments needed")
    for arg in args:
        if arg.shape != args[0].shape:
            raise ValueError("Arguments must have the same shape")
    return True


def assert_same_shapeN(N, *args):
    assert_array(*args)
    if len(args) < 2:
        raise ValueError("At least two arguments needed")
    for a in args:
        if a.shape[N] != args[0].shape[N]:
            raise ValueError("Args %dth dimmension must be the same: %d ≠ %d"
                             % (N, a.shape[N], args[0].shape[N]))


def assert_BPT_consistency(B, P, T):
    assert_array_dim(2, P)
    assert_array(B, T)

    E, S = P.shape
    if B.ndim == 1:
        assert_shape((E,), B)
        if T.ndim == 1:
            assert_shape((S,), T)
    elif B.ndim == 2:
        N = B.shape[0]
        assert_shape((N, E), B)
        if T.ndim == 1:
            assert_shape((S,), T)
        elif T.ndim == 2:
            assert_shape((N, S), T)
        else:
            raise ValueError("T has wrong number of dimmensions")
    else:
        raise ValueError("B has wrong number of dimmensions")
    return True


def assert_sequence(*args):
    # from collections.abc import Sequence
    for a in args:
        if not isinstance(a, collections.abc.Sequence):
            raise TypeError


# def setkw(kwargs, keyw, default_values):
#     # import warnings
#     # from collections.abc import Iterable
#     # DELETE. Use default keywords
#     # warnings.warn('Use default keywords instead', DeprecationWarning)
#     if isinstance(keyw, collections.abc.Iterable) and \
#        isinstance(default_values, collections.abc.Iterable):
#         retlist = []
#         for k, dv in zip(keyw, default_values):
#             if k in kwargs.keys():
#                 retlist.append(kwargs[k])
#             else:
#                 retlist.append(dv)
# 
#         return retlist
#     else:
#         if keyw in kwargs:
#             return kwargs[keyw]
#         else:
#             return default_values
# 
# 
# def setkwpop(kwargs, keyw, default_values):
#     # from collections.abc import Iterable
#     # DELETE. Use default keywords
#     # import warnings
#     # warnings.warn('Use default keywords instead', DeprecationWarning)
#     if isinstance(keyw, collections.abc.Iterable) and \
#             isinstance(default_values, collections.abc.Iterable):
#         retlist = []
#         for k, dv in zip(keyw, default_values):
#             if k in kwargs.keys():
#                 retlist.append(kwargs.pop(k))
#             else:
#                 retlist.append(dv)
# 
#         return retlist
#     else:
#         if keyw in kwargs:
#             return kwargs.pop(keyw)
#         else:
#             return default_values


def unravel(x, flags):
    """Unravel data according to flags provided.

    This routine takes an array of data and an array of flags of the same
    length and returns another array with only the independet variables.

    Parameters:
        x (iterable): the original data values
        flags (iterable): the flags indicating how to update x.
            Values must int. Accepted values are

            * 0: value is to be kept constant
            * 1: value is to be refined and the corresponding value from x
              will be substituted by the corresponding value from y.
            * >2: value is restrained. All places with the same number are
              refined together and the ratio between them is maintained.

    Returns:
        generator: Values of **x** processed according to **flags**.
    """
    constr_list = []
    for i, f in zip(x, flags):
        if f == 1:
            yield i
        if f > 1:
            if f not in constr_list:
                yield i
                constr_list.append(f)


def ravel(x, y, flags):
    """Update values from one iterable with other iterable according to flags.

    This function does the opposite action than :func:`unravel`.

    Parameters:
        x (iterable): the original array values
        y (iterable): the updated values to be plugged into *x*.
        flags (sequence): flags indicating how to update *x* with *y*. Accepted
            values are:

            * 0: value is to be kept constant
            * 1: value is to be refined and the corresponding value from x
              will be substituted by the corresponding value from y.
            * >2: value is restrained. All places with the same number are
              refined together and the ratio between them is maintained.

    Yields:
        float: Raveled values.
    """
    # indices of the reference parameter for constraining
    ref_index = {i: flags.index(i) for i in range(2, 1+max(flags))}
    ref_val = {}

    ity = iter(y)
    for i, f in enumerate(flags):
        if f == 1:                      # refinable: return new value
            yield next(ity)
        elif f == 0:                    # constant: return old value
            yield x[i]
        else:                           # constrained: return or compute
            if i == ref_index[f]:
                val = next(ity)         # reference value: return new value
                ref_val[f] = val        # and store ref value
                yield val
            else:                       # other: compute proportional value
                yield x[i]*ref_val[f]/x[ref_index[f]]


# TODO rewrite this to do more stuff
# It should return True or False rather than modify in-place the list.
def check_refine_flags(flags):
    """check that there are not isolated constraints which don't make any
    sense. Replace them with code '1'"""
    maxconstr = 9
    for restrn in range(2, maxconstr):
        if flags.count(restrn) == 1:
            flags[flags.index] = 1


def breakdown_counted_flags(counted_flags):
    constant = counted_flags.pop(0, 0)
    refined = counted_flags.pop(1, 0)
    constraints = dict(counted_flags)
    return constant, refined, constraints


def count_flags(*args):
    counters = (collections.Counter(arg) for arg in args)
    sumall = functools.reduce(lambda x, y: x+y, counters)
    return sumall


def total_refine_params(counted_flags):
    """The total number of parameters to refine.

    Parameters:
        counted_flags (:class:`collections.Counter`): the counted parameters
    Returns:
        int: the total number of parameters to be refined.
    """
    return counted_flags.get(1, 0) + sum(x not in {0,1} for x in counted_flags)


def build_multiT_titr3(titration, concat=False):
    """Build a multi T from a titration object.

    Parameters:
        titration (list): Information about the titrations. Each member
            of the list is one titration. Every element in the list is a dict
            for which the accepted keys are:

                - V (sequence): volume of titre in mL)
                - V0 (float): the initial volume in mL
                - T0 (sequence): total amount of initial components in mmol
                - buret (sequence): concentration in the buret
                  for each component in mmol/mL.
                - Tflags (sequence, optional): The refinement flag for T0
                - buret_flags (sequence, optional): The refinement flag for emf0
                - T (sequence): An (N, S)-array with the initial
                  concentrations. If provided, V, V0, T0, buret and Tflags
                  are ignored.

        concat (bool): if True, the list of T calculated will be concatenated
            along axis=1. If False, the list of T will be returned
    Returns:
        list of arrays or array: The T arrays calculated
    .. seealso:: :func:`build_multiT_titr`, :func:`build_multiT_titr2`
    """
    Tret = []
    if concat:
        call = Tret.extend
    else:
        call = Tret.append

    for t in titration:
        if 'T' in t:
            call(t['T'])
        else:
            call(build_T_titr2(t['T0'], t['buret'], t['V0'], t['V']))

    return Tret


def build_multiT_titr(T0, buret, V0, incrV, N):
    """Build multiple titrations.

    Given the starting conditions calculates the total concentrations
    for every observed point. This function calculates that for more than
    one titration series. Every titration is given in every item of the
    arguments which must be lists or tuples of the same length.
    For this particular function the progression
    of the titration is given in the length of the argument *V*.

    Parameters:
        T0 (sequence): initial millimoles
        buret (sequence): buret contents in mmol/mL
        V0 (sequence): initial volume in mL
        incrV (sequence): volume increment in mL
        N (int): number of points for each titration
    Returns:
        tuple: The total concentrations for every titration
    .. seealso:: :func:`build_T_titr2`
    """
    assert_same_len(T0, buret, V0, incrV, N)
    return tuple(build_T_titr(t0, b, v0, iv, n)
                 for t0, b, v0, iv, n in zip(T0, buret, V0, incrV, N))


def build_multiT_titr2(T0, buret, V0, V):
    """Build multiple titrations.

    Given the starting conditions calculates the total concentrations
    for every observed point. This function calculates that for more than
    one titration series. Every titration is given in every item of the
    arguments which must be lists or tuples of the same length.
    For this particular function the progression
    of the titration is given in the length of the argument *V*.

    Parameters:
        T0 (sequence): initial millimoles
        buret (sequence): buret contents in mmol/mL
        V0 (sequence): initial volume in mL
        V (sequence): volume of titre in mL
    Yields:
        The total concentrations for every titration
    .. seealso:: :func:`build_T_titr2`
    """
    # assert_same_len(T0, buret, V0, V)   # they may be iters
    for t0, b, v0, v in zip(T0, buret, V0, V):
        yield build_T_titr2(t0, b, v0, v)


def build_T_titr(T0, buret, V0, incrV, N):
    """Build total concentrations array from titration data.

    Given the starting conditions calculates the total concentrations
    for every observed point.

    Parameters:
        T0 (sequence): initial millimoles
        buret (sequence): buret contents in mmol/mL
        V0 (float): initial volume in mL
        incrV (float): volume increment in mL
        N (int): number of points
    Returns:
        generator: The total concentrations calculated
    """
    return (tuple((t + b*n*incrV)/(n*incrV + V0) for t, b in zip(T0, buret))
            for n in range(N))


def build_T_titr2(T0, buret, V0, V):
    """Build total concentrations array from titration data.

    Given the starting conditions calculates the total concentrations
    for every observed point.

    Arguments:
        T0 (sequence): initial millimoles for every component.
        buret (sequence): buret contents in mmol/mL, one for every component.
        V0 (float): initial volume in mL
        V (sequence): volume of titre in mL.
    Returns:
        generator: The total concentrations calculated
    """
    if V0 <= 0.0:
        raise ValueError("The starting volume cannot be zero or negative")
    return (tuple((t + b*v)/(v+V0) for t, b in zip(T0, buret)) for v in V)


def build_T(T0, pT, N):
    """Build the total concentrations array from the starting values.

    Parameters:
        T0 (tuple): The total initial concentrations in mol/L. The length of
            this tuple must match the number of free components and each
            element of the tuple must be a doublet of floats containing the
            initial and final titre.
        pT (tuple of bools): A flag stating whether each value in **T0** is
            in logarithmic units (True) or in linear unit (False). Length
            must match that of **T0**.
        N (int): The number of points to be generated.
    Returns:
        tuple: The total concentrations array.
    Examples:
        >>> libaux.build_T(((0.1,0.2),(0.1,0.2)),(False,False), 5)
        ((0.1, 0.1), (0.125, 0.125), (0.15, 0.15), (0.175,0.175), (0.2, 0.2))
    """
    assert_same_len(T0, pT)

    aux = []
    for (t0, t1), pt in zip(T0, pT):
        if pt:
            aux.append(tuple(10**(-i) for i in linspace(t0, t1, N)))
        else:
            aux.append(tuple(linspace(t0, t1, N)))

    return tuple(zip(*aux))


def linspace(init, end, N):
    """Provide **N** evenly spaced numbers between **init** and **end**

    Parameters:
        init (float): The initial number
        end (float): The final number
        N (int): The number of values to return
    Yields:
        **N** evenly spaced numbers between **init** and **end**
    """
    for n in range(N):
        yield init + (end-init)*n/(N-1)


def get_lims(data):
    """This function finds the limits of a series of lists in order to join
    and unjoin them as required.
    """
    lims = [0]
    for i in [len(s) for s in data]:
        lims.append(lims[-1]+i)
    return lims


def make_ranges(a: list, offset: int = 0):
    assert isinstance(a, list)
    assert all(isinstance(i, int) for i in a)

    index_toplot = [n+offset for n, i in enumerate(a) if i != 0]
    reduced_index = [y-x for x, y in zip(index_toplot[0:-1], index_toplot[1:])]
    mode = 0
    ret = ""
    for n, i in zip(index_toplot, reduced_index):
        if i != 1 and mode % 2 == 0:
            if mode == 2:
                ret += ","
            ret += str(n) + ","
            continue

        if i == 1 and mode % 2 == 0:
            ret += str(n) + "-"
            mode = 1
            continue

        if i == 1 and mode == 1:
            continue

        if i != 1 and mode == 1:
            ret += str(n) + ","
            mode = 2

    ret += str(index_toplot[-1])
    return ret


def nwe2str(n, e):
    """Format number and error into string.

    Given a list of numbers with corresponding error, return a string
    for each one containing the formatted number.

    Parameters:
        n (iterable): the numbers
        e (iterable): the error associated with **n** It must have the same
            length as **n**.

    Returns:
        generator: returns string representations of value(error)

    Example:
        >>> list(nwe2str([10.1234, 21.4433], [0.03, 0.02]))
        ['10.12(3)', '21.44(2)']

    .. seealso:: :func:`numwerr`
    """

    return (numwerr(n_, e_) for n_, e_ in zip(n, e))


def nwe2arr(its):
    """Extract value/error values from string.

    Given a list of str representing numbers with corresponding error in
    the format "number(error)", return two lists containing the separate
    values and error in float format.

    Parameters:
        its (iterable): strings representing number and error

    Returns:
        list: numbers
        list: errors

    Example:
        >>> nwe2arr(['10.12(3)', '21.44(2)'])
        ([10.12, 21.44], [0.03, 0.02])
    """
    return (inumwerr(s) for s in its)


def numwerr(n, e):
    """Given two floats representing a value and an error to that value this
    routine returns an human-readable string.

    Parameters:
        n (float): a number
        e (float): the error associated with **n**

    Returns:
        str: A string representation of value(error)

    Example:
        >>> numwerr(10.1234, 0.03)
        "10.12(3)"
    """
    import math
    assert_type(float, n, e)
    if e < 0:
        raise ValueError('Negative errors are meaningless')

    if e == 0.0:
        return str(n)

    dp = -math.floor(math.log10(e))
    s = "%%.%df(%%d)" % dp
    return s % (n, int(e*10**dp))


def inumwerr(s):
    """This function does the opposite as :func:`numwerr`. Given a string
    representing a value and an error to that value this
    routine returns both the number and the error as float numbers.

    Parameters:
        s (str): a number with error in the format number(error)

    Returns:
        float: The number
        float: The error

    Example:
        >>> numwerr("10.12(3)")
        (10.1234, 0.03)
    """
    m = re.compile(r'([+-]?\d+(\.\d+)?)(\((\d+)\))?').match(s)
    if not m:
        raise ValueError("Value '%s' could not be processed" % s)

    nt = m.group(1)
    et = m.group(4)
    if et is None:
        e = 0.0
    else:
        q = nt.find('.')
        if q == -1:
            e = float(et)
        else:
            nd = len(nt)-q-1   # number of decimal positions
            e = float(et)*10**(-nd)
    return float(nt), e


def percent_distribution(concentration, stoichiometry, analytc, reference):
    """Transform free concentrations to relative concentrations.

    This function converts absolute concentrations to relative percent of
    concentrations with respect to a given species.

    Parameters:
        concentration (:class:`numpy.ndarray`): The free concentrations array.
            It must be an (*N*, *E* + *S*)-sized array of floats.
        stoichiometry (:class:`numpy.ndarray`): Stoichimetric coefficients
        analytc (:class:`numpy.ndarray`): The total concentration array
        reference (int): the index of the reference species

    Returns:
        :class:`numpy.ndarray`: The concentrations expressed as percent with
            respect to the reference species. The first dimmension is
            unchanged with respect to that of **C** but the second one is
            shorter as all components whose stoichiometric coefficient with
            respect to the reference one is zero have been removed.
    """
    # import numpy as np
    assert_array_dim(2, stoichiometry)
    n_equilibria, n_species = stoichiometry.shape
    n_data = analytc.shape[0]
    assert_shape((n_data, n_species), analytc)
    assert_shape((n_data, n_equilibria+n_species), concentration)

    # r is the reference species
    rP = stoichiometry[:, reference].tolist()   # species related to 'r'

    # find species to plot (all those that contain reference species)
    to_plot = [n+n_species for n, i in enumerate(rP) if i != 0]
    to_plot.insert(0, reference)
    to_plot.sort()
    q_plot = [i for i in rP if i != 0]
    q_plot.insert(0, 1)

    # remove species not to plot
    new_c = concentration[:, to_plot].copy()

    # divide by T
    T_ref = analytc[:, reference]
    T_ref_ext = T_ref.reshape((n_data, 1)) * np.ones((1, len(to_plot)))
    new_c = new_c / T_ref_ext

    # divide by stoichiometric coef.
    new_c = new_c * (np.ones((n_data, 1)) * np.array(q_plot))

    new_c = 100 * new_c  # multiply by 100

    return new_c


def combine_components(stoich, ref):
    """Group components by kind.

    Sometimes it is useful to group the compounds that share a particular
    component. For example, all the compounds that differ in the degree of
    protonation could be considered as one combined 'species'.

    Arguments:
        stoich (sequence): The stoichiometry array
        ref (int): The reference component

    Returns:
        dict: The keys are tuples containing the indices of the components
            that belong to that family. The values are the indices of the
            equilibria that belong to that family.

    >>> P = [[1,0,1], [1,0,2], [1,0,3], [1,0,4], [1,0,5], [1,0,6], [0,1,1],
    ...     [0,1,2], [0,1,3], [1,1,3], [1,1,4], [1,1,5], [1,1,6], [1,1,7],
    ...     [0,0,-1]]
    >>> ref = 2
    >>> combs = combine_components(P, ref)
    {(0,): [0, 1, 2, 3, 4, 5], (1,): [6, 7, 8], (0, 1): [9, 10, 11, 12, 13],
    (): [14]}
    """
    retval = {}
    for equil, row in enumerate(stoich):
        key = tuple(comp_index
                    for comp_index, comp_stoi in enumerate(row)
                    if comp_stoi != 0 and comp_index != ref)
        if key in retval:
            retval[key].append(equil)
        else:
            retval[key] = [equil]
    return retval


def exhaust(generator):
    "Empty a generator."
    for _ in generator:
        pass


def extend_stoich(array):
    """Return extended stoichiometry array."""
    # import numpy as np
    n_comp = array.shape[1]
    return np.vstack((np.eye(n_comp), array))


def extract_groups(array, xref, yref=None):
    """Separate the components in like groups.

    Given a number of components defined by the stoichiometry array, and a
    reference components (usually H), this routine finds subgroups with
    similar compositions and return a list of ints of the same length than
    the total number of species in which each value represents the groups
    assigned to.

    This function is meant to be a help for coloring speciation graphs.

    Parameter:
        array: The *extended* stoichiometry array.
        xref (int): the x reference component
        yref (int, optional): the y reference component
    Yields:
        tuple: the indices of the group each component belongs to
    """
    y = np.delete(array, xref, axis=1)
    if yref is None:
        c = [tuple(_) for _ in y]
    else:
        c = [tuple(a) for n, a in enumerate(y) if array[n, yref]]

    aux_dict = {t: n for n, t in enumerate(set(c))}
    return [aux_dict[q] for  n, q in enumerate(c)]


def related_components(component: int, stoich):
    _stoich = np.array(stoich)
    size = _stoich.shape[1]
    yield component
    yield from (size + i for i in np.flatnonzero(_stoich[:, component]))


def _utf_subscr(string):
    """Convert numbers to unicode subscript."""
    mdict = dict(zip(u"0123456789-", u"₀₁₂₃₄₅₆₇₈₉₋"))
    return ''.join(mdict.get(c, c) for c in string)


def extended_labels(label_seed, stoich, kind='plain'):
    'Name the constants based on the stoichiometry including free species.'
    yield from iter(label_seed)
    yield from constant_labels(label_seed, stoich, kind)


def constant_labels(label_seed, stoich, kind='plain'):
    '''Name the constants based on the stoichiometry.

    This routine takes the names of the components (label_seed) and the matrix
    of stoichiometric coefficients and makes a list of species in equilibrium.
    For instance: if primary species are 'L' and 'H' and 'L' can take up to
    three protons ('H'), the species in equilibrium are (including H-1): 'L',
    'H', 'HL', 'H_2L', 'H_3L' and 'H-1' and that is the list returned.

    Parameters:
        label_seed (list): list of :class:`str` with the bases for making
            the labels.
        stoich (:class:`ndarray`): the stoichiometric coefficients.
        kind (str, optional): the format to output required. Values accepted
            are 'plain', for plain text output (no subindex or formatting);
            'latex' for LaTeX formatting or 'html' for HTML; 'unicode'
            for unicode substring substitution.

    Yields:
        Strings with the constructed labels in the order they appear in *stoich*.
    '''
    templ = {'plain': '{}', 'latex': '{{}}_{{{}}}', 'html': '<sub>{}</sub>'}
    for row in stoich:
        lbl = []
        for l, p in zip(label_seed, row):
            if p == 0:
                a = ""
            elif p == 1:
                a = l
            else:
                if kind == 'unicode':
                    a = l + _utf_subscr(str(p))
                else:
                    a = l + templ[kind].format(p)
            lbl.append(a)
        yield "".join(lbl)
