# Heatwave Impact Assessment (geoskill-heatwave-impact-assessment)

> Heatwave impact assessment: per-pixel P90 quantile threshold with connected-run detection over ≥3 consecutive days, Stull/simplified wet-bulb temperature health risk classification, population exposure and vulnerability mapping. Heatwave impact assessment: per-pixel P90 threshold with connected-run detection, wet-bulb temperature health risk, population exposure and vulnerability mapping. Outputs heatwave-day/exposed-population/vulnerability/wet-bulb-risk GeoTIFF + event list JSON.

---

## 1. Overview

Performs heatwave detection and impact assessment on daily maximum (or daily mean) temperature time-series rasters, for high-temperature disaster risk screening, public health emergency response, and urban thermal environment research. Core algorithms: - **Heatwave detection**: per-pixel threshold at the 90th percentile (P90) of the time series, with connected-component labeling along the time axis (scipy.ndimage.label); runs exceeding the threshold for ≥ min-duration (default 3) consecutive days are classified as heatwave events. Outputs heatwave days, longest duration, event count, peak temperature, and heatwave mask. - **Wet-bulb temperature**: estimates Tw from peak temperature and relative humidity. `stull` uses the Stull (2011) empirical formula (J. Appl. Meteorol. Climatol.), `simple` uses a simplified dry-wet bulb depression approximation. Tw ≥ 35°C approaches the physiological tolerance limit of the human body. - **Health risk classification**: Tw is divided into five levels: low (<26), moderate (26–28), high (28–30), severe (30–32), and extreme (≥32). - **Population exposure**: heatwave mask overlaid with the population raster yields the exposed population raster and the total/ratio of exposed population. - **Vulnerability**: normalized heatwave intensity × normalized population density gives a composite vulnerability index in [0,1]. The `--synthetic` mode generates simulated temperature, population, and humidity fields with injected heatwave events (eastern half sustained +9°C), enabling full offline verification of the pipeline without network access or real data.

## 2. Features

Performs heatwave detection and impact assessment on daily maximum (or daily mean) temperature time-series rasters, for high-temperature disaster risk screening, public health emergency response, and urban thermal environment research. Core algorithms: - **Heatwave detection**: per-pixel threshold at the 90th percentile (P90) of the time series, with connected-component labeling along the time axis (scipy.ndimage.label); runs exceeding the threshold for ≥ min-duration (default 3) consecutive days are classified as heatwave events. Outputs heatwave days, longest duration, event count, peak temperature, and heatwave mask. - **Wet-bulb temperature**: estimates Tw from peak temperature and relative humidity. `stull` uses the Stull (2011) empirical formula (J. Appl. Meteorol. Climatol.), `simple` uses a simplified dry-wet bulb depression approximation. Tw ≥ 35°C approaches the physiological tolerance limit of the human body. - **Health risk classification**: Tw is divided into five levels: low (<26), moderate (26–28), high (28–30), severe (30–32), and extreme (≥32). - **Population exposure**: heatwave mask overlaid with the population raster yields the exposed population raster and the total/ratio of exposed population. - **Vulnerability**: normalized heatwave intensity × normalized population density gives a composite vulnerability index in [0,1]. The `--synthetic` mode generates simulated temperature, population, and humidity fields with injected heatwave events (eastern half sustained +9°C), enabling full offline verification of the pipeline without network access or real data.

## 3. Quick Start

```bash
pip install -r requirements.txt
python geoskill-heatwave-impact-assessment.py --bbox 116 39 117 40 --synthetic --output-dir ./out
```

## 4. CLI Parameters

Run `python <skill>.py --help` for the full list. Common parameters:

| Parameter | Type | Description |
|---|---|---|
| `--bbox` | `float[4]` | WGS84 bounding box `min_lon min_lat max_lon max_lat` |
| `--input` | `path` | Local input file (GeoJSON/GeoTIFF/etc.) |
| `--output-dir` | `path` | Output directory (default `./output`) |
| `--synthetic` | `flag` | Use synthetic data instead of real input |
| `--quiet` | `flag` | Suppress non-essential stdout |

## 5. Input / Output

| File | Format | Description |
|---|---|---|
| `heatwave_days.tif` | GeoTIFF (float32) | Per-pixel heatwave days, EPSG:4326 |
| `exposed_population.tif` | GeoTIFF (float32) | Exposed population (heatwave mask × population) |
| `vulnerability.tif` | GeoTIFF (float32) | Normalized vulnerability index [0,1] |
| `wetbulb_risk.tif` | GeoTIFF (float32) | Wet-bulb temperature health risk level 0-4 |
| `heatwave_events.json` | JSON | Regional event list + exposure/risk statistics |
| `output-manifest.json` | JSON | Run manifest (input/output/QA/software versions) |


## 6. Technical Principle

(see SKILL.md for details)

## 7. Methodology

This skill has been methodologically reviewed. See [`REVIEW.md`](./REVIEW.md) for:

- P0/P1/P2 issue counts and verdicts
- Reproduction commands
- Known limitations and edge cases

## 8. License

MIT License. See [`LICENSE`](./LICENSE) for full text.

---

# 热浪影响评估（中文版）

> 本部分为完整中文文档，与上方英文部分对应。

---
name: geoskill-heatwave-impact-assessment
description: '热浪影响评估：逐像元 P90 分位数阈值 + 连续≥3天连通域热浪检测，Stull/简化湿球温度健康风险分级，人口暴露与脆弱性制图。Heatwave impact assessment: per-pixel P90 threshold with connected-run detection, wet-bulb temperature health risk, population exposure and vulnerability mapping. 输出热浪日数/暴露人口/脆弱性/湿球风险 GeoTIFF + 事件清单 JSON。'
---

# 热浪影响评估 | Heatwave Impact Assessment

对日最高温（或日均温）时序栅格执行热浪检测与影响评估，用于高温灾害风险
筛查、公共卫生应急响应与城市热环境研究。核心算法：

- **热浪检测**：逐像元以时序第 90 百分位（P90）为阈值，沿时间轴做连通域
  标记（scipy.ndimage.label），连续超阈值 ≥ min-duration（默认 3）天的过程判
  为热浪事件，输出热浪日数、最长持续、事件数、峰值温度与热浪掩膜。
- **湿球温度**：由峰值温度与相对湿度估算 Tw。`stull` 采用 Stull (2011) 经验式
  （J. Appl. Meteorol. Climatol.），`simple` 采用简化干湿差近似。Tw ≥ 35°C 接近
  人体生理耐受极限。
- **健康风险分级**：由 Tw 分为低(<26)/中(26-28)/高(28-30)/严重(30-32)/极端(≥32)
  五级。
- **人口暴露**：热浪掩膜叠加入口栅格，得暴露人口栅格与暴露总量/比例。
- **脆弱性**：归一化热浪强度 × 归一化人口密度，得 [0,1] 综合脆弱性指数。

支持 `--synthetic` 模式生成含注入热浪事件（东半部持续 +9°C）的模拟温度场、
人口与湿度场，无需网络和真实数据即可离线验证全流程。

## 依赖

```bash
pip install 'numpy' 'rasterio' 'scipy'
```

## 使用方法

### 基本用法（仅给 bbox，自动合成）

```bash
python geoskill-heatwave-impact-assessment.py --bbox 116.0 39.0 117.0 40.0 --output-dir ./out
```

### 示例 1：合成数据离线评估

```bash
python geoskill-heatwave-impact-assessment.py --bbox 116.0 39.0 117.0 40.0 --synthetic --n-dates 40 --output-dir ./out
```

### 示例 2：简化湿球温度方法

```bash
python geoskill-heatwave-impact-assessment.py --bbox 121.0 31.0 122.0 32.0 --synthetic --method simple --output-dir ./sh
```

### 示例 3：真实温度时序栅格

```bash
python geoskill-heatwave-impact-assessment.py --input tmax_ts.tif --threshold-pct 90 --min-duration 3 --output-dir ./real
```

### 示例 4：更严格的热浪定义（P95 + ≥5 天）

```bash
python geoskill-heatwave-impact-assessment.py --input tmax_ts.tif --threshold-pct 95 --min-duration 5 --output-dir ./strict
```

### 示例 5：Stull vs 简化湿球温度对比

```bash
python geoskill-heatwave-impact-assessment.py --bbox 116 39 117 40 --synthetic --method stull --output-dir ./cmp_stull --quiet
python geoskill-heatwave-impact-assessment.py --bbox 116 39 117 40 --synthetic --method simple --output-dir ./cmp_simple --quiet
```

## 输出

| 文件 | 格式 | 说明 |
|---|---|---|
| `heatwave_days.tif` | GeoTIFF (float32) | 逐像元热浪日数，EPSG:4326 |
| `exposed_population.tif` | GeoTIFF (float32) | 暴露人口（热浪掩膜×人口） |
| `vulnerability.tif` | GeoTIFF (float32) | 归一化脆弱性指数 [0,1] |
| `wetbulb_risk.tif` | GeoTIFF (float32) | 湿球温度健康风险等级 0-4 |
| `heatwave_events.json` | JSON | 区域事件清单 + 暴露/风险统计 |
| `output-manifest.json` | JSON | 运行清单（输入/输出/QA/软件版本） |

## 数据源 / Source

- **输入**：本地多波段日最高温时序 GeoTIFF（每波段一日）；人口与湿度在真实
  模式下用合成占位场（本 skill 不绑定特定人口数据源）
- **合成模式**：本地生成含注入热浪的温度/人口/湿度场，无外部数据源
- **湿球公式**：Stull, R. (2011), J. Appl. Meteorol. Climatol. 50(11)

## 隐私声明 / Privacy

- 默认完全离线运行，不发起任何网络请求
- `--synthetic` 模式不读取任何外部数据
- 所有计算在本地完成，不上传用户数据

## License

MIT
