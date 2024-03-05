#!/usr/bin/python3
# -*- coding: utf-8 -*-

"""Exceptions used by APES.

This module adds some exceptions that are useful for the iterative functions
to inform upstream when the iteration is not working
"""

__version__ = '0.1'


class FailedCalculateConcentrations(Exception):
    """When calculation of concentrations fails."""

    def __init__(self, msg, last_value, **diagn):
        super().__init__()
        self.msg = msg
        self.last_value = last_value
        self.diagnostic = diagn

    def __str__(self):
        return self.msg


class NothingToRefine(Exception):
    "There are no free variables to refine."


class NotConvergenceException(Exception):
    """When convergence is not reached.

    It contains the last value of the iterations in case this information can
    be valuable.
    """

    def __init__(self, msg, last_value):
        super().__init__()
        self.msg = msg
        self.last_value = last_value

    def __str__(self):
        return self.msg


class ModelIncomplete(Exception):
    """The model is incomplete.

    This exception is thrown when a fitting or update is requested but
    the data required to do so has not been given yet or it is incomplete."""


class TooManyIterations(NotConvergenceException):
    """When maximum number of iterations is reached.

    This exception is thrown when the maximum number of iterations
    has been reached without meeting the convergence criteria
    """

class UnstableIteration(NotConvergenceException):
    """When iteration is unstable and must be stopped."""


class DataFileException(Exception):
    "An error ocurred while reading a data file."
