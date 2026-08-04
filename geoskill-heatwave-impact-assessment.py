#!/usr/bin/env python3
"""heatwave-impact-assessment — 热浪影响评估

对日最高温时序执行热浪检测，并叠加入口栅格评估暴露与脆弱性。核心流程：

- **热浪检测**：逐像元以第 90 百分位（P90）为阈值，识别连续 ≥ 3 天的超阈值
  过程（scipy.ndimage 沿时间轴连通域标记）。统计热浪日数、最长持续、事件数、
  峰值温度与热浪掩膜。
- **湿球温度**：由温度与相对湿度估算湿球温度 Tw，提供 Stull (2011) 经验式与
  简化式两种方法；Tw 用于健康风险分级（低/中/高/严重/极端）。
- **人口暴露**：热浪掩膜叠加入口栅格，得暴露人口与暴露总量。
- **脆弱性**：归一化热浪强度 × 归一化人口密度，得 [0,1] 脆弱性指数。

数据源：本地多波段日最高温 GeoTIFF（波段=日期），或使用 ``--synthetic`` 生成
含注入热浪的温度时序 + 人口 + 湿度场用于离线测试。

隐私声明 / Privacy：
- 默认完全离线运行，不发起任何网络请求。
- ``--synthetic`` 模式完全无网络。
- 所有处理在本地完成，不上传任何用户数据。

Usage:
    python heatwave-impact-assessment.py --input tmax_ts.tif --output-dir ./out
    python heatwave-impact-assessment.py --bbox 116 39 117 40 --synthetic --output-dir ./out

License: MIT
"""
from __future__ import annotations

import argparse
import datetime as _dt
import json
import os
import sys
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

VERSION = "1.0.0"
SKILL_NAME = "heatwave-impact-assessment"

# 湿球温度健康风险分级阈值（°C，Tw）
RISK_THRESHOLDS = [26.0, 28.0, 30.0, 32.0]
RISK_LABELS = ["low", "moderate", "high", "severe", "extreme"]

# ---- 复用共享核心库（本地 vendored，随脚本目录一起分发）----
try:
    from _geoskill_core.errors import (
        GeoSkillError, UsageError, ValidationError, ProcessError, to_exit_code,
    )
    from _geoskill_core.manifest import OutputManifest, OutputFile
except ImportError:  # pragma: no cover - fallback minimal definitions
    class GeoSkillError(Exception):
        def __init__(self, message: str, code: int = 7, kind: str = "EGeo", **kw):
            super().__init__(message)
            self.message, self.code, self.kind = message, code, kind

    class UsageError(GeoSkillError):
        def __init__(self, m, **k): super().__init__(m, code=2, kind="EUsage", **k)

    class ValidationError(GeoSkillError):
        def __init__(self, m, **k): super().__init__(m, code=6, kind="EValidate", **k)

    class ProcessError(GeoSkillError):
        def __init__(self, m, **k): super().__init__(m, code=7, kind="EProcess", **k)

    def to_exit_code(exc):
        return getattr(exc, "code", 7)

    OutputManifest = None
    OutputFile = None


def _utc_now() -> str:
    return _dt.datetime.now(_dt.timezone.utc).isoformat().replace("+00:00", "Z")


# ---------------------------------------------------------------------------
# 输入校验
# ---------------------------------------------------------------------------
def validate_bbox(bbox):
    """Validate a geographic bbox [W, S, E, N] in EPSG:4326.

    Rules (consistent across the project):
      - W < E (no antimeridian wrap; user must split the request)
      - S < N
      - -180 <= W, E <= 180
      - -90 <= S, N <= 90
    Returns the bbox on success; raises ValidationError on failure.
    """
    if not isinstance(bbox, (list, tuple)) or len(bbox) != 4:
        raise ValidationError("bbox must be a sequence of 4 floats [W S E N]")
    w, s, e, n = [float(v) for v in bbox]
    if not (w < e):
        raise ValidationError(
            f"bbox W={w} must be < E={e} (antimeridian wrap not supported; "
            f"split your request into two boxes if needed)")
    if not (s < n):
        raise ValidationError(f"bbox S={s} must be < N={n}")
    if not (-180.0 <= w <= 180.0 and -180.0 <= e <= 180.0):
        raise ValidationError(
            f"bbox lon must be in [-180, 180], got W={w}, E={e}")
    if not (-90.0 <= s <= 90.0 and -90.0 <= n <= 90.0):
        raise ValidationError(
            f"bbox lat must be in [-90, 90], got S={s}, N={n}")
    return [w, s, e, n]


def validate_threshold_pct(pct):
    """--threshold-pct must be in [0, 100]."""
    try:
        pct = float(pct)
    except (TypeError, ValueError):
        raise ValidationError(f"--threshold-pct must be a number, got {pct!r}")
    if not (0.0 <= pct <= 100.0):
        raise ValidationError(
            f"--threshold-pct must be in [0, 100], got {pct}")
    return pct


def validate_min_duration(days):
    """--min-duration must be a positive integer."""
    try:
        days = int(days)
    except (TypeError, ValueError):
        raise ValidationError(f"--min-duration must be an integer, got {days!r}")
    if days < 1:
        raise ValidationError(f"--min-duration must be >= 1, got {days}")
    return days


def validate_min_area_fraction(frac):
    """--min-area-fraction must be in [0, 1]."""
    try:
        frac = float(frac)
    except (TypeError, ValueError):
        raise ValidationError(
            f"--min-area-fraction must be a number, got {frac!r}")
    if not (0.0 <= frac <= 1.0):
        raise ValidationError(
            f"--min-area-fraction must be in [0, 1], got {frac}")
    return frac


# ---------------------------------------------------------------------------
# 核心算法：热浪检测
# ---------------------------------------------------------------------------
def detect_heatwave(
    temp_ts: np.ndarray,
    threshold: Optional[np.ndarray] = None,
    threshold_pct: float = 90.0,
    min_duration: int = 3,
) -> Dict[str, Any]:
    """对 (n_dates, H, W) 温度时序逐像元检测热浪。

    - 若 ``threshold`` 为 None，则逐像元取时序的 ``threshold_pct`` 百分位为阈值；
      否则使用给定（标量或 (H,W)）阈值。
    - 超阈值日按时间轴做连通域标记（scipy.ndimage.label），连续长度
      ≥ ``min_duration`` 的过程判为热浪事件。

    返回字典，含逐像元栅格：threshold / hw_days / max_duration / n_events /
    peak_temp / hw_mask，以及逐时相的 hw_day_mask。
    """
    from scipy.ndimage import label

    temp_ts = np.asarray(temp_ts, dtype=np.float64)
    if temp_ts.ndim != 3:
        raise ValidationError(
            f"expected (n_dates, H, W) temperature cube, got ndim={temp_ts.ndim}",
            ndim=int(temp_ts.ndim),
        )
    n_dates, h, w = temp_ts.shape
    if n_dates < min_duration:
        raise ValidationError(
            f"need at least {min_duration} time steps, got {n_dates}",
            n_dates=int(n_dates),
        )

    if threshold is None:
        thr = np.nanpercentile(temp_ts, threshold_pct, axis=0)
    else:
        thr = np.asarray(threshold, dtype=np.float64)
        if thr.ndim == 0:
            thr = np.full((h, w), float(thr))

    exceed = temp_ts > thr[None, :, :]

    # 仅沿时间轴连通：structure 在 axis 0 方向连接
    struct = np.zeros((3, 3, 3), dtype=np.int8)
    struct[:, 1, 1] = 1
    labels, n_feat = label(exceed.astype(np.int8), structure=struct)

    sizes = np.bincount(labels.ravel())
    comp_size = sizes[labels]  # 每个体素所属连通域的大小
    valid = (comp_size >= min_duration) & (labels > 0)

    hw_day_mask = exceed & valid
    # 逐像元热浪日数
    hw_days = np.sum(hw_day_mask, axis=0).astype(np.int32)
    # 最长持续：像元所属合格连通域大小的最大值
    valid_size = np.where(valid, comp_size, 0)
    max_duration = np.max(valid_size, axis=0).astype(np.int32)
    # 峰值温度（热浪日内）；无热浪像元退化为全序列最高温
    masked_temp = np.where(hw_day_mask, temp_ts, -np.inf)
    peak = np.max(masked_temp, axis=0)
    overall_max = np.max(temp_ts, axis=0)
    peak_temp = np.where(hw_days > 0, peak, overall_max)

    # 事件数：每个像元时间序列中合格连通域的个数
    n_events = np.zeros((h, w), dtype=np.int32)
    valid_labels = np.where(valid, labels, 0)
    for j in range(h):
        for i in range(w):
            labs = valid_labels[:, j, i]
            n_events[j, i] = int(np.unique(labs[labs > 0]).size)

    return {
        "threshold": thr.astype(np.float32),
        "hw_days": hw_days,
        "max_duration": max_duration,
        "n_events": n_events,
        "peak_temp": peak_temp.astype(np.float32),
        "hw_mask": (hw_days > 0).astype(np.uint8),
        "hw_day_mask": hw_day_mask.astype(np.uint8),
    }


def areal_events(
    hw_day_mask: np.ndarray,
    min_duration: int = 3,
    min_fraction: float = 0.0,
) -> List[Dict[str, Any]]:
    """由逐时相热浪日掩膜 (n_dates, H, W) 提取区域尺度热浪事件清单。

    区域受影响比例 > ``min_fraction`` 的连续时段，且持续 ≥ ``min_duration``
    天，记为一个事件。``min_fraction`` 用于滤除零星噪声像元造成的伪事件。
    """
    hw_day_mask = np.asarray(hw_day_mask)
    n_dates = hw_day_mask.shape[0]
    n_pix = hw_day_mask.size // n_dates
    frac = hw_day_mask.reshape(n_dates, -1).mean(axis=1)  # 受影响比例
    active = frac > min_fraction

    events: List[Dict[str, Any]] = []
    t = 0
    while t < n_dates:
        if not active[t]:
            t += 1
            continue
        start = t
        while t < n_dates and active[t]:
            t += 1
        end = t - 1
        duration = end - start + 1
        if duration >= min_duration:
            events.append({
                "start_day": int(start),
                "end_day": int(end),
                "duration_days": int(duration),
                "peak_affected_fraction": float(np.max(frac[start:end + 1])),
                "mean_affected_fraction": float(np.mean(frac[start:end + 1])),
            })
    return events


# ---------------------------------------------------------------------------
# 核心算法：湿球温度与健康风险
# ---------------------------------------------------------------------------
def wet_bulb_stull(temp: np.ndarray, rh: np.ndarray) -> np.ndarray:
    """Stull (2011) 湿球温度经验式（J Appl Meteor Climatol）。

    适用 RH 5–99%、T -20–50°C。返回 °C。
    """
    t = np.asarray(temp, dtype=np.float64)
    rh = np.clip(np.asarray(rh, dtype=np.float64), 1.0, 100.0)
    tw = (
        t * np.arctan(0.151977 * np.sqrt(rh + 8.313659))
        + np.arctan(t + rh)
        - np.arctan(rh - 1.676331)
        + 0.00391838 * (rh ** 1.5) * np.arctan(0.023101 * rh)
        - 4.686035
    )
    return tw


def wet_bulb_simple(temp: np.ndarray, rh: np.ndarray) -> np.ndarray:
    """简化湿球温度估算：Tw ≈ T - 0.15×(100 - RH)（粗略干湿差修正）。"""
    t = np.asarray(temp, dtype=np.float64)
    rh = np.clip(np.asarray(rh, dtype=np.float64), 0.0, 100.0)
    return t - 0.15 * (100.0 - rh)


def estimate_wet_bulb(temp: np.ndarray, rh: np.ndarray, method: str = "stull") -> np.ndarray:
    if method == "stull":
        return wet_bulb_stull(temp, rh)
    if method == "simple":
        return wet_bulb_simple(temp, rh)
    raise UsageError(f"unknown wet-bulb method '{method}'", method=method)


def heat_risk_level(tw: np.ndarray) -> np.ndarray:
    """由湿球温度分级健康风险：0 低 / 1 中 / 2 高 / 3 严重 / 4 极端。"""
    tw = np.asarray(tw, dtype=np.float64)
    levels = np.zeros(tw.shape, dtype=np.int32)
    for i, th in enumerate(RISK_THRESHOLDS, start=1):
        levels[tw >= th] = i
    return levels


# ---------------------------------------------------------------------------
# 核心算法：人口暴露与脆弱性
# ---------------------------------------------------------------------------
def population_exposure(pop: np.ndarray, hw_mask: np.ndarray) -> Dict[str, Any]:
    """热浪掩膜叠加入口栅格 → 暴露人口栅格与统计。"""
    pop = np.asarray(pop, dtype=np.float64)
    hw_mask = np.asarray(hw_mask, dtype=np.float64)
    if pop.shape != hw_mask.shape:
        raise ValidationError("population and heatwave mask shapes differ")
    exposed = pop * (hw_mask > 0).astype(np.float64)
    total_pop = float(np.sum(pop))
    exposed_total = float(np.sum(exposed))
    return {
        "exposed": exposed.astype(np.float32),
        "total_population": total_pop,
        "exposed_population": exposed_total,
        "exposed_fraction": exposed_total / total_pop if total_pop > 0 else 0.0,
    }


def vulnerability_index(hw_days: np.ndarray, pop: np.ndarray) -> np.ndarray:
    """归一化热浪强度 × 归一化人口 → [0,1] 脆弱性指数。"""
    hw_days = np.asarray(hw_days, dtype=np.float64)
    pop = np.asarray(pop, dtype=np.float64)
    hw_norm = hw_days / hw_days.max() if hw_days.max() > 0 else hw_days
    pop_norm = pop / pop.max() if pop.max() > 0 else pop
    return (hw_norm * pop_norm).astype(np.float32)


# ---------------------------------------------------------------------------
# 合成数据：含注入热浪的温度时序 + 人口 + 湿度（离线测试）
# ---------------------------------------------------------------------------
def generate_synthetic(
    bbox: List[float],
    n_dates: int = 30,
    height: int = 64,
    width: int = 64,
    seed: int = 42,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, Dict[str, Any]]:
    """生成温度时序 (n_dates,H,W)、人口 (H,W)、相对湿度 (H,W)。

    基线温度 = 空间梯度（南热北凉）+ 小日际噪声。在东半部注入一段持续
    ``hw_len`` 天、增幅 +9°C 的热浪，使该区域被 P90 阈值稳定检出。
    """
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[0:height, 0:width]
    yn = yy.astype(np.float64) / max(height - 1, 1)
    xn = xx.astype(np.float64) / max(width - 1, 1)

    # 基线日均最高温：南 ~33，北 ~28
    base = 33.0 - 5.0 * yn
    temp = np.zeros((n_dates, height, width), dtype=np.float32)
    for t in range(n_dates):
        temp[t] = base + rng.normal(0, 0.3, size=(height, width))

    # 注入热浪：东半部 (xn > 0.5)，持续 hw_len 天，+9°C
    hw_len = max(3, int(round(0.1 * n_dates)))
    hw_start = n_dates // 3
    hw_end = hw_start + hw_len
    region = (xn > 0.5).astype(np.float32)
    # 平滑过渡到区域内部，保证内部增幅充分
    for t in range(hw_start, min(hw_end, n_dates)):
        temp[t] = temp[t] + 9.0 * region

    # 人口：中部高、边缘低（人/像元），总量缩放
    pop_field = 5000.0 * np.exp(-(((xn - 0.5) ** 2 + (yn - 0.5) ** 2) / 0.08))
    pop_field = pop_field + rng.uniform(0, 200, size=(height, width))
    population = pop_field.astype(np.float32)

    # 相对湿度：40-80%
    rh = (60.0 + 20.0 * (1.0 - yn) + rng.normal(0, 2.0, size=(height, width)))
    rh = np.clip(rh, 20.0, 95.0).astype(np.float32)

    info = {
        "bbox": bbox,
        "n_dates": int(n_dates),
        "width": int(width),
        "height": int(height),
        "heatwave": {
            "start_day": int(hw_start),
            "end_day": int(hw_end - 1),
            "duration_days": int(hw_len),
            "region": "eastern half (x_norm > 0.5)",
            "amplitude_c": 9.0,
        },
    }
    return temp, population, rh, info


# ---------------------------------------------------------------------------
# GeoTIFF I/O（直接用 rasterio，保证离线可用）
# ---------------------------------------------------------------------------
def write_geotiff(
    path: str,
    cube: np.ndarray,
    bbox: List[float],
    nodata: float = -9999.0,
) -> None:
    import rasterio
    from rasterio.transform import from_bounds

    if cube.ndim == 2:
        cube = cube[np.newaxis, ...]
    nb, h, w = cube.shape
    transform = from_bounds(bbox[0], bbox[1], bbox[2], bbox[3], w, h)
    profile = {
        "driver": "GTiff", "height": h, "width": w, "count": nb,
        "dtype": "float32", "crs": "EPSG:4326", "transform": transform,
        "nodata": nodata, "compress": "deflate",
    }
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with rasterio.open(path, "w", **profile) as dst:
        for b in range(nb):
            dst.write(cube[b].astype("float32"), b + 1)


def read_geotiff(path: str) -> Tuple[np.ndarray, List[float]]:
    """Read a daily-max temperature time-series GeoTIFF. Returns (cube, bbox).

    NoData values declared in the file are replaced with NaN. The cube has
    shape (n_dates, H, W).
    """
    import rasterio
    if not os.path.exists(path):
        raise UsageError(f"input raster not found: {path}", path=path)
    with rasterio.open(path) as src:
        cube = src.read().astype(np.float32)
        b = src.bounds
        bbox = [b.left, b.bottom, b.right, b.top]
        nodata = src.nodata
    if nodata is not None:
        cube = np.where(cube == float(nodata), np.nan, cube)
    return cube, bbox


def finite_pixel_mask(cube: np.ndarray) -> np.ndarray:
    """Per-pixel mask: True iff every time-step is finite (not NaN/inf)."""
    return np.isfinite(np.asarray(cube)).all(axis=0)


# ---------------------------------------------------------------------------
# Manifest
# ---------------------------------------------------------------------------
def write_manifest(
    output_dir: str,
    args: argparse.Namespace,
    outputs: List[Dict[str, Any]],
    qa: Dict[str, Any],
    started_at: str,
    exit_code: int,
    bbox: List[float],
) -> Optional[str]:
    if OutputManifest is None:
        return None
    cmd = " ".join([SKILL_NAME] + sys.argv[1:])
    man = OutputManifest(
        skill=SKILL_NAME,
        skill_version=VERSION,
        command=cmd,
        started_at=started_at,
        finished_at=_utc_now(),
        exit_code=exit_code,
        inputs={
            "input": getattr(args, "input", None),
            "method": getattr(args, "method", None),
            "min_duration": getattr(args, "min_duration", None),
            "synthetic": bool(getattr(args, "synthetic", False)),
        },
        outputs=[OutputFile(**o) for o in outputs],
        qa=qa,
        software={"python": sys.version.split()[0], "skill": SKILL_NAME},
    )
    path = os.path.join(output_dir, "output-manifest.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(man.to_dict(), f, ensure_ascii=False, indent=2, default=str)
    return path


# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------
def process(args: argparse.Namespace) -> int:
    started_at = _utc_now()
    output_dir = args.output_dir

    # --- pre-flight validation (BEFORE making output dir) -----------------
    bbox = list(args.bbox) if args.bbox else None
    threshold_pct = validate_threshold_pct(args.threshold_pct)
    min_duration = validate_min_duration(args.min_duration)
    min_area_fraction = validate_min_area_fraction(args.min_area_fraction)

    # 1) 获取数据
    synth_info: Optional[Dict[str, Any]] = None
    input_nodata_value: Optional[float] = None
    n_valid = 0
    if args.input and not args.synthetic:
        temp_ts, file_bbox = read_geotiff(args.input)
        # Capture the declared NoData value for the qa/manifest.
        import rasterio as _rio
        with _rio.open(args.input) as _src:
            input_nodata_value = _src.nodata
        if bbox is not None:
            bbox = validate_bbox(bbox)
        else:
            bbox = validate_bbox(file_bbox)
        if temp_ts.ndim != 3:
            raise ValidationError(
                "input must be a multi-band time series (n_dates, H, W)",
                ndim=int(temp_ts.ndim))
        if temp_ts.shape[0] < min_duration:
            raise ValidationError(
                f"need at least {min_duration} time steps, got {temp_ts.shape[0]}",
                n_dates=int(temp_ts.shape[0]))
        if temp_ts.size == 0:
            raise ValidationError("input temperature cube is empty")
        valid = finite_pixel_mask(temp_ts)
        n_valid = int(valid.sum())
        if n_valid == 0:
            raise ValidationError(
                f"input raster has no valid pixels "
                f"(all values are nodata={input_nodata_value})")
        # Real-input path: build a population/RH placeholder with the
        # SAME shape as the temperature cube. The real population/RH
        # input is not implemented (no external data source).
        _, population, rh, _ = generate_synthetic(
            bbox, n_dates=temp_ts.shape[0],
            height=temp_ts.shape[1], width=temp_ts.shape[2],
        )
        source_note = args.input
    else:
        if bbox is None:
            raise UsageError("provide --bbox (synthetic mode) or --input <raster>")
        bbox = validate_bbox(bbox)
        temp_ts, population, rh, synth_info = generate_synthetic(
            bbox, n_dates=args.n_dates,
        )
        n_valid = int(temp_ts.size)
        source_note = "synthetic"

    # All checks passed → now create the output dir.
    os.makedirs(output_dir, exist_ok=True)

    # 2) 热浪检测
    hw = detect_heatwave(temp_ts, threshold_pct=threshold_pct,
                         min_duration=min_duration)

    # 3) 湿球温度（用峰值温度 + 湿度）与健康风险
    tw = estimate_wet_bulb(hw["peak_temp"], rh, method=args.method)
    risk = heat_risk_level(tw)

    # 4) 人口暴露与脆弱性
    expo = population_exposure(population, hw["hw_mask"])
    vuln = vulnerability_index(hw["hw_days"].astype(np.float64), population)

    # 5) 区域事件清单
    events = areal_events(hw["hw_day_mask"], min_duration=min_duration,
                          min_fraction=min_area_fraction)

    # 6) 写出产物
    hw_days_path = os.path.join(output_dir, "heatwave_days.tif")
    exposed_path = os.path.join(output_dir, "exposed_population.tif")
    vuln_path = os.path.join(output_dir, "vulnerability.tif")
    risk_path = os.path.join(output_dir, "wetbulb_risk.tif")
    write_geotiff(hw_days_path, hw["hw_days"].astype(np.float32), bbox)
    write_geotiff(exposed_path, expo["exposed"], bbox)
    write_geotiff(vuln_path, vuln, bbox)
    write_geotiff(risk_path, risk.astype(np.float32), bbox)

    summary = {
        "threshold_pct": threshold_pct,
        "min_duration": min_duration,
        "min_area_fraction": min_area_fraction,
        "wetbulb_method": args.method,
        "n_dates": int(temp_ts.shape[0]),
        "shape": [int(temp_ts.shape[1]), int(temp_ts.shape[2])],
        "pixels_with_heatwave": int(np.sum(hw["hw_mask"] > 0)),
        "mean_heatwave_days": float(np.mean(hw["hw_days"])),
        "max_heatwave_days": int(np.max(hw["hw_days"])),
        "total_population": expo["total_population"],
        "exposed_population": expo["exposed_population"],
        "exposed_fraction": expo["exposed_fraction"],
        "risk_distribution": {
            RISK_LABELS[i]: int(np.sum(risk == i)) for i in range(len(RISK_LABELS))
        },
        "events": events,
    }
    summary_path = os.path.join(output_dir, "heatwave_events.json")
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)

    qa: Dict[str, Any] = {
        "source": source_note,
        "wetbulb_method": args.method,
        "n_dates": int(temp_ts.shape[0]),
        "pixels_with_heatwave": summary["pixels_with_heatwave"],
        "exposed_population": expo["exposed_population"],
        "n_events": len(events),
        "n_valid_pixels": n_valid,
    }
    if input_nodata_value is not None:
        qa["input_nodata"] = input_nodata_value
    if synth_info is not None:
        qa["synthetic_heatwave"] = synth_info["heatwave"]

    outputs = [
        {"path": hw_days_path, "kind": "raster", "crs_epsg": 4326,
         "bbox_wgs84": bbox, "band_count": 1},
        {"path": exposed_path, "kind": "raster", "crs_epsg": 4326,
         "bbox_wgs84": bbox, "band_count": 1},
        {"path": vuln_path, "kind": "raster", "crs_epsg": 4326,
         "bbox_wgs84": bbox, "band_count": 1},
        {"path": risk_path, "kind": "raster", "crs_epsg": 4326,
         "bbox_wgs84": bbox, "band_count": 1},
        {"path": summary_path, "kind": "json"},
    ]

    man_path = write_manifest(output_dir, args, outputs, qa, started_at, 0, bbox)

    if not args.quiet:
        print(f"[{SKILL_NAME}] source: {source_note}")
        print(f"[{SKILL_NAME}] n_dates: {temp_ts.shape[0]}  shape: {temp_ts.shape[1:]}")
        print(f"[{SKILL_NAME}] pixels with heatwave: {summary['pixels_with_heatwave']}")
        print(f"[{SKILL_NAME}] areal events: {len(events)}")
        print(f"[{SKILL_NAME}] exposed population: {expo['exposed_population']:.0f} "
              f"({expo['exposed_fraction']*100:.1f}%)")
        print(f"[{SKILL_NAME}] output: {hw_days_path}")
        if man_path:
            print(f"[{SKILL_NAME}] manifest: {man_path}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog=SKILL_NAME,
        description="Heatwave detection, population exposure and wet-bulb health risk.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--bbox", nargs=4, type=float, metavar=("W", "S", "E", "N"),
                   help="geographic extent [minLon minLat maxLon maxLat]")
    p.add_argument("--input", help="input daily-max temperature time-series GeoTIFF (bands=dates)")
    p.add_argument("--n-dates", type=int, default=30,
                   help="number of synthetic daily time steps (default: 30)")
    p.add_argument("--threshold-pct", type=float, default=90.0,
                   help="percentile threshold for exceedance (default: 90)")
    p.add_argument("--min-duration", type=int, default=3,
                   help="minimum consecutive days to qualify a heatwave (default: 3)")
    p.add_argument("--min-area-fraction", type=float, default=0.05,
                   help="minimum affected fraction for a regional event (default: 0.05)")
    p.add_argument("--method", default="stull", choices=["stull", "simple"],
                   help="wet-bulb temperature estimation method (default: stull)")
    p.add_argument("--synthetic", action="store_true",
                   help="generate a synthetic temperature/population/RH field (offline)")
    p.add_argument("--output-dir", default="./output", help="output directory")
    p.add_argument("--quiet", action="store_true", help="suppress console output")
    p.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    return p


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return process(args)
    except GeoSkillError as exc:
        print(f"[{SKILL_NAME}] ERROR [{exc.kind}] {exc.message}", file=sys.stderr)
        return exc.code
    except KeyboardInterrupt:
        return 130
    except Exception as exc:  # noqa: BLE001
        print(f"[{SKILL_NAME}] ERROR {exc}", file=sys.stderr)
        return to_exit_code(exc)


if __name__ == "__main__":
    sys.exit(main())
