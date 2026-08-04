"""Core algorithm tests for heatwave-impact-assessment."""
import sys
import os

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(__file__))
from conftest import mod


class TestDetectHeatwave:
    def test_explicit_threshold_deterministic(self):
        """显式阈值下，已知热浪被准确检出。"""
        ts = np.zeros((10, 1, 2), dtype=np.float64)
        ts[:, 0, 0] = [30, 30, 35, 36, 37, 38, 30, 30, 30, 30]  # days 2-5 hot (4 days)
        ts[:, 0, 1] = 30.0
        r = mod.detect_heatwave(ts, threshold=33.0, min_duration=3)
        assert r["hw_days"][0, 0] == 4
        assert r["max_duration"][0, 0] == 4
        assert r["n_events"][0, 0] == 1
        assert r["peak_temp"][0, 0] == 38.0
        assert r["hw_days"][0, 1] == 0
        assert r["hw_mask"][0, 0] == 1
        assert r["hw_mask"][0, 1] == 0

    def test_two_events_counted(self):
        """两段分离热浪 → n_events=2。"""
        ts = np.full((12, 1, 1), 25.0)
        ts[1:4, 0, 0] = 40.0   # event 1: days 1-3
        ts[8:11, 0, 0] = 40.0  # event 2: days 8-10
        r = mod.detect_heatwave(ts, threshold=30.0, min_duration=3)
        assert r["n_events"][0, 0] == 2
        assert r["hw_days"][0, 0] == 6

    def test_short_run_not_counted(self):
        """持续 < min_duration 不算热浪。"""
        ts = np.full((10, 1, 1), 25.0)
        ts[3:5, 0, 0] = 40.0  # only 2 days
        r = mod.detect_heatwave(ts, threshold=30.0, min_duration=3)
        assert r["hw_days"][0, 0] == 0
        assert r["n_events"][0, 0] == 0

    def test_percentile_detection(self):
        """百分位阈值路径：注入 3 天热浪应被检出。"""
        rng = np.random.default_rng(0)
        ts = rng.normal(30, 0.3, size=(30, 1, 1))
        ts[20:23, 0, 0] = 40.0  # 3-day heatwave
        r = mod.detect_heatwave(ts, threshold_pct=90.0, min_duration=3)
        assert r["hw_days"][0, 0] >= 3
        assert r["max_duration"][0, 0] >= 3

    def test_wrong_ndim_raises(self):
        with pytest.raises(mod.ValidationError):
            mod.detect_heatwave(np.zeros((10, 10)))

    def test_too_few_dates_raises(self):
        with pytest.raises(mod.ValidationError):
            mod.detect_heatwave(np.zeros((2, 4, 4)), min_duration=3)


class TestArealEvents:
    def test_clean_single_event(self):
        """合成注入热浪 → 单个区域事件，起止日正确。"""
        temp, pop, rh, info = mod.generate_synthetic([116, 39, 117, 40], n_dates=30, seed=42)
        hw = mod.detect_heatwave(temp, threshold_pct=90.0, min_duration=3)
        events = mod.areal_events(hw["hw_day_mask"], min_duration=3, min_fraction=0.05)
        assert len(events) == 1
        ev = events[0]
        inj = info["heatwave"]
        assert ev["start_day"] == inj["start_day"]
        assert ev["end_day"] == inj["end_day"]
        assert ev["duration_days"] == inj["duration_days"]
        assert ev["peak_affected_fraction"] > 0.4  # ~半区受影响

    def test_empty_mask_no_events(self):
        mask = np.zeros((10, 8, 8), dtype=np.uint8)
        assert mod.areal_events(mask, min_duration=3) == []


class TestWetBulb:
    def test_stull_less_than_dry(self):
        """湿球温度应低于干球温度（RH<100%）。"""
        tw = mod.wet_bulb_stull(35.0, 50.0)
        assert tw < 35.0

    def test_stull_increases_with_rh(self):
        """RH 越高湿球温度越高。"""
        tw_dry = mod.wet_bulb_stull(35.0, 30.0)
        tw_wet = mod.wet_bulb_stull(35.0, 90.0)
        assert tw_wet > tw_dry

    def test_stull_known_approx(self):
        """Stull 公式在 40°C/60% 下约 33°C（文献范围）。"""
        tw = float(mod.wet_bulb_stull(40.0, 60.0))
        assert 31.0 < tw < 35.0

    def test_simple_increases_with_rh(self):
        tw_dry = mod.wet_bulb_simple(30.0, 20.0)
        tw_wet = mod.wet_bulb_simple(30.0, 100.0)
        assert tw_wet > tw_dry
        # RH=100% → Tw≈T
        assert abs(float(mod.wet_bulb_simple(30.0, 100.0)) - 30.0) < 1e-6

    def test_dispatch_and_unknown(self):
        a = mod.estimate_wet_bulb(30.0, 50.0, method="stull")
        b = mod.estimate_wet_bulb(30.0, 50.0, method="simple")
        assert np.isfinite(a) and np.isfinite(b)
        with pytest.raises(mod.UsageError):
            mod.estimate_wet_bulb(30.0, 50.0, method="foo")


class TestHeatRisk:
    def test_levels(self):
        tw = np.array([20.0, 27.0, 29.0, 31.0, 33.0])
        levels = mod.heat_risk_level(tw)
        np.testing.assert_array_equal(levels, [0, 1, 2, 3, 4])

    def test_extreme_at_35(self):
        levels = mod.heat_risk_level(np.array([35.0, 40.0]))
        assert np.all(levels == 4)


class TestPopulationExposure:
    def test_exposed_equals_pop_times_mask(self):
        """暴露人口 = 热浪掩膜 × 人口（像元级一致）。"""
        pop = np.arange(16, dtype=np.float64).reshape(4, 4)
        mask = np.zeros((4, 4), dtype=np.uint8)
        mask[0:2, 0:2] = 1
        expo = mod.population_exposure(pop, mask)
        expected = pop * mask
        np.testing.assert_allclose(expo["exposed"], expected, atol=1e-4)
        assert abs(expo["exposed_population"] - float(expected.sum())) < 1e-3
        assert expo["total_population"] == float(pop.sum())

    def test_shape_mismatch_raises(self):
        with pytest.raises(mod.ValidationError):
            mod.population_exposure(np.zeros((4, 4)), np.zeros((3, 3)))


class TestVulnerability:
    def test_range_01(self):
        rng = np.random.default_rng(1)
        hw = rng.integers(0, 6, (8, 8)).astype(np.float64)
        pop = rng.uniform(0, 5000, (8, 8))
        vuln = mod.vulnerability_index(hw, pop)
        assert vuln.min() >= 0.0
        assert vuln.max() <= 1.0

    def test_zero_heatwave_zero_vuln(self):
        hw = np.zeros((4, 4))
        pop = np.ones((4, 4)) * 100
        vuln = mod.vulnerability_index(hw, pop)
        assert np.all(vuln == 0.0)


class TestSynthetic:
    def test_shapes(self):
        temp, pop, rh, info = mod.generate_synthetic([116, 39, 117, 40], n_dates=20)
        assert temp.shape[0] == 20
        assert pop.shape == rh.shape == temp.shape[1:]
        assert "heatwave" in info

    def test_injected_heatwave_detected_east(self):
        """注入热浪的东半部应被大量检出，西半部基本无。"""
        temp, pop, rh, info = mod.generate_synthetic([116, 39, 117, 40], n_dates=30, seed=7)
        hw = mod.detect_heatwave(temp, threshold_pct=90.0, min_duration=3)
        east = int(hw["hw_mask"][:, 48:].sum())
        west = int(hw["hw_mask"][:, :16].sum())
        assert east > 500
        assert east > west * 10


class TestGeoTiffIO:
    def test_write_and_read_roundtrip(self, tmp_path):
        cube = np.random.default_rng(2).uniform(0, 40, (5, 16, 16)).astype(np.float32)
        bbox = [116.0, 39.0, 117.0, 40.0]
        path = str(tmp_path / "test.tif")
        mod.write_geotiff(path, cube, bbox)
        assert os.path.exists(path)
        read_back, read_bbox = mod.read_geotiff(path)
        assert read_back.shape == cube.shape
        np.testing.assert_allclose(read_bbox, bbox, atol=1e-6)
        np.testing.assert_allclose(read_back, cube, atol=1e-5)

    def test_read_missing_file_raises(self):
        with pytest.raises(mod.UsageError):
            mod.read_geotiff("/nonexistent/path/file.tif")
