from __future__ import division
from time import sleep

import numpy as np
import pandas as pd
import numpy.ma as ma
import os
import h5py


def dotd(A, B, out=None):
    """Diagonal of :math:`\mathrm A\mathrm B^\intercal`.
    If ``A`` is :math:`n\times p` and ``B`` is :math:`p\times n`, it is done in :math:`O(pn)`.
    Args:
        A (array_like): Left matrix.
        B (array_like): Right matrix.
        out (:class:`numpy.ndarray`, optional): copy result to.
    Returns:
        :class:`numpy.ndarray`: Resulting diagonal.
    """
    A = ma.asarray(A, float)
    B = ma.asarray(B, float)
    if A.ndim == 1 and B.ndim == 1:
        if out is None:
            return ma.dot(A, B)
        return ma.dot(A, B, out)

    if out is None:
        out = ma.empty((A.shape[0], ), float)

    out[:] = ma.sum(A * B.T, axis=1)
    return out

def logdet(X):
    return np.log(np.linalg.det(X))
    # UC = np.linalg.cholesky(X)
    # return 2*sum(np.log(np.diag(UC)))

def sigmoid(X):
    """ Method to compute sigmoid function """
    return np.divide(1.,1.+np.exp(-X))
    # return 1./(1.+np.exp(-X))

def ddot(d, mtx, left=True):
    """Multiply a full matrix by a diagonal matrix.
    This function should always be faster than dot.

    Input:
      d -- 1D (N,) array (contains the diagonal elements)
      mtx -- 2D (N,N) array
      left: is the diagonal matrix on the left or on the right of the product?

    Output:
      ddot(d, mts, left=True) == dot(diag(d), mtx)
      ddot(d, mts, left=False) == dot(mtx, diag(d))
    """
    if left:
        return (d*mtx.T).T
    else:
        return d*mtx

def lambdafn(X):
    return np.tanh(X/2.)/(4.*X)

def nans(shape, dtype=float):
    """ Method to create an array filled with missing values """
    a = np.empty(shape, dtype)
    a.fill(np.nan)
    return a