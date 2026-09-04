"""Time-domain characterisation helpers for Assignment 1.

``summary_stats``, ``seasonal_cycle`` and ``decorrelation_timescale`` characterise a
series in the time domain. Each docstring states exactly what to return.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from amoc_analysis.spectra_filtering.data_io import fill_gaps


def summary_stats(values: np.ndarray) -> dict[str, float]:
    """Basic descriptive statistics of a (possibly gappy) series.

    Parameters
    ----------
    values : numpy.ndarray, shape (N,)
        The series, which may contain ``NaN``.

    Returns
    -------
    stats : dict[str, float]
        Dictionary with keys:

        ``n``
            Total number of samples (including gaps).
        ``n_missing``
            Number of ``NaN`` values.
        ``mean``, ``std``, ``median``, ``min``, ``max``
            Computed over the finite values (NaN-aware).
        ``range``
            ``max - min``.

    Notes
    -----
    Use the NaN-aware reductions (``numpy.nanmean`` etc.) so gaps do not poison the
    statistics. Decide and document whether ``std`` uses ``ddof=0`` or ``1``.

    DONE (student): implement and return the dictionary.
    """

    min = np.nanmin(values)
    max = np.nanmax(values)

    return {
        'n':        len(values),
        'n_missing':int(np.count_nonzero(np.isnan(values))),
        'mean':     np.nanmean(values),
        'std':      np.nanstd(values),
        'median':   np.nanmedian(values),
        'min':      min,
        'max':      max,
        'range':    max - min
    }


def seasonal_cycle(
    time: np.ndarray,
    values: np.ndarray,
    by: str = "month",
) -> pd.DataFrame:
    """Climatological seasonal cycle (mean and median per calendar period).

    Parameters
    ----------
    time : numpy.ndarray, dtype ``datetime64``
        Time coordinate, same length as ``values``.
    values : numpy.ndarray, shape (N,)
        The series.
    by : {"month", "dayofyear"}, optional
        Grouping period. Default ``"month"`` (a 12-point climatology).

    Returns
    -------
    clim : pandas.DataFrame
        Indexed by the period (e.g. months 1–12), with columns ``"mean"`` and
        ``"median"``.

    Notes
    -----
    Build a ``pandas`` object indexed by ``time`` and use a ``groupby`` on the
    period (e.g. ``s.groupby(s.index.month)``) with ``.agg(["mean", "median"])``.
    The **seasonal range** is then ``clim["mean"].max() - clim["mean"].min()``.

    DONE (student): implement with a pandas groupby and return the DataFrame.
    """

    # Create pandas Series with datetime index
    s = pd.Series(values, index=pd.DatetimeIndex(time))

    # Extract grouping key from 'by' parameter
    if by == "month":
        group_key = s.index.month
    elif by == "dayofyear":
        group_key = s.index.dayofyear
    else:
        raise ValueError(f"'by' must be 'month' or 'dayofyear', got '{by}'")

    # Group by period and aggregate
    clim = s.groupby(group_key).agg(["mean", "median"])

    return clim


def decorrelation_timescale(
    values: np.ndarray,
    dt: float,
) -> tuple[float, float]:
    """Integral decorrelation timescale and effective degrees of freedom.

    Parameters
    ----------
    values : numpy.ndarray, shape (N,)
        Evenly sampled series with **no gaps** (interpolate first if needed).
    dt : float
        Sampling interval (same time units you want the timescale in).

    Returns
    -------
    integral_scale : float
        Integral timescale ``tau`` (in units of ``dt``).
    ndof : float
        Effective degrees of freedom, ``N * dt / tau - 1``.

    Notes
    -----
    Integral-timescale method (after E. Frajka-Williams' ``calc_ndof.m``):

    1. Remove the mean.
    2. Form the **normalised autocovariance** ``R`` (autocorrelation, ``R(0) = 1``)
       at lags ``0, dt, 2*dt, ...`` (the non-negative lags).
    3. Integrate ``R`` from lag 0 **until its first zero crossing**, by the
       trapezoidal rule:
       ``tau = sum_i dt * (R[i] + R[i+1]) / 2`` while ``R[i] >= 0``.
       This one-sided integral (not doubled) is the integral timescale.
    4. ``ndof = N * dt / tau - 1``  (``N`` = number of samples).

    A white series gives ``tau`` of order ``dt`` and ``ndof`` close to ``N``; a
    strongly autocorrelated series gives a large ``tau`` and few independent samples.

    TODO (student): implement the autocovariance, the zero-crossing integral, and
    return ``(integral_scale, ndof)``.
    """

    filled_values = fill_gaps(values)
    N = len(filled_values)

    # 1. Remove the mean
    # values_centered = values - np.mean(values)
    values_centered = filled_values - np.mean(filled_values)

    # 2. Form the normalised autocovariance R (autocorrelation) at non-negative lags
    raw_autocov = np.correlate(values_centered, values_centered, mode='full')
    autocov_positive = raw_autocov[N-1:]

    variance = autocov_positive[0]

    if variance == 0:
        return 0.0, 0.0

    R = autocov_positive / variance

    # 3. Integrate R from lag 0 until its first zero crossing
    zero_crossings = np.where(R <= 0)[0]

    if len(zero_crossings) == 0:
        cutoff = len(R) - 1
    else:
        cutoff = zero_crossings[0]

    if cutoff == 0:
        integral_sum = 0.0
    else:
        r_left = R[:cutoff]
        r_right = R[1:cutoff+1]
        integral_sum = np.sum((r_left + r_right) / 2)

    integral_scale = integral_sum * dt

    # 4. Calculate ndof
    if integral_scale == 0:
        ndof = 0.0
    else:
        ndof = (N * dt) / integral_scale - 1.0

    return float(integral_scale), float(ndof)
