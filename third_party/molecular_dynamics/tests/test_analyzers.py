# third_party/molecular_dynamics/tests/test_analyzers.py
"""Tests for analyzer components."""

import numpy as np
import pytest

from lib.analyzers.base import detect_convergence


def test_detect_convergence_stable():
    """RMSD 一直稳定 → 从 t=0 收敛。"""
    time_ns = np.arange(0, 50, 0.1)  # 50ns, 0.1ns interval
    values = np.ones_like(time_ns) * 2.0 + np.random.normal(0, 0.01, len(time_ns))
    result = detect_convergence(time_ns, values, window_ns=5.0, variance_threshold=0.05)
    assert result >= 0.0
    assert result < 5.0  # should converge early


def test_detect_convergence_drift():
    """RMSD 持续上升 → 不收敛。"""
    time_ns = np.arange(0, 50, 0.1)
    values = time_ns * 0.1  # linearly increasing
    result = detect_convergence(time_ns, values, window_ns=5.0, variance_threshold=0.001)
    assert result == -1.0


def test_detect_convergence_late():
    """RMSD 先大幅波动后稳定 → 晚期收敛。"""
    time_ns = np.arange(0, 50, 0.1)
    # Early phase: large oscillations ensuring high variance in every window
    values = np.where(time_ns < 20, 2.0 + np.sin(time_ns * 2) * 1.5, 2.0)
    values += np.random.normal(0, 0.01, len(values))
    result = detect_convergence(time_ns, values, window_ns=5.0, variance_threshold=0.05)
    assert result >= 19.0  # converges around t=20ns
