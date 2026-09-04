"""Starter package for the spectra-and-filtering lecture demo.

Worked helpers (``data_io``, frequency axes, the synthetic tone, the Parseval
check) are provided so students can focus on the estimators. The core routines
(``raw_periodogram``, ``welch_psd``, ``butterworth_squared_response``) are left as
stubs in the student version — the accompanying ``pytest`` checks encode the
behaviour they must satisfy.
"""

from .analysis import decorrelation_timescale, seasonal_cycle, summary_stats
from .data_io import fill_gaps, load_moc
from .filters import butterworth_squared_response, nyquist_frequency, tukey_lowpass
from .leakage import synthetic_tone
from .spectra import frequency_axis, parseval_ratio, raw_periodogram, welch_psd

__all__ = [
    "load_moc",
    "fill_gaps",
    "nyquist_frequency",
    "butterworth_squared_response",
    "tukey_lowpass",
    "frequency_axis",
    "raw_periodogram",
    "welch_psd",
    "parseval_ratio",
    "synthetic_tone",
    "summary_stats",
    "seasonal_cycle",
    "decorrelation_timescale",
]
