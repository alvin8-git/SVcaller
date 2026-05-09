"""Tests for core overlap logic — no file I/O needed."""
import sys
sys.path.insert(0, "/data/alvin/SVcaller/bin")
import pytest

def reciprocal_overlap(a_start, a_end, b_start, b_end) -> float:
    """Fraction of reciprocal overlap between two intervals."""
    overlap = max(0, min(a_end, b_end) - max(a_start, b_start))
    len_a = a_end - a_start
    len_b = b_end - b_start
    if len_a == 0 or len_b == 0:
        return 0.0
    return overlap / min(len_a, len_b)

def test_full_overlap():
    assert reciprocal_overlap(100, 200, 100, 200) == 1.0

def test_no_overlap():
    assert reciprocal_overlap(100, 200, 300, 400) == 0.0

def test_partial_overlap():
    result = reciprocal_overlap(100, 200, 150, 250)
    assert abs(result - 0.5) < 1e-9

def test_contained():
    result = reciprocal_overlap(100, 400, 150, 250)
    assert result == 1.0
