# contdid

[English](README.md) | [简体中文](README_zh.md)

**面向 Stata 的连续处理双重差分**

[![Stata 16+](https://img.shields.io/badge/Stata-16%2B-blue.svg)](https://www.stata.com/)
[![License: AGPL-3.0](https://img.shields.io/badge/License-AGPL--3.0-blue.svg)](LICENSE)
[![Version: 1.1.0](https://img.shields.io/badge/Version-1.1.0-green.svg)](CHANGELOG.md)

`contdid` 实现 **Callaway、Goodman-Bacon & Sant'Anna (2024/2025) 的连续处理双重差分
（DiD）估计量**，填补 Stata 生态的空白：作者的参考实现（`contdid`）只有 **R** 版。

## 概述

在很多 DiD 设计中，处理并不是简单地"打开"——它有**剂量/强度**（政策力度、补贴金额、
污染暴露、金融发展指数……）。`contdid` 估计**剂量-反应函数** ATT(d) 及其导数
**平均因果响应** ACRT(d)。

在平行趋势下，特定剂量的效应识别为

> **ATT(d) = E[ΔY | D = d] − E[ΔY | D = 0]**，

其中 ΔY 为结果的处理前到处理后变化，D 为连续剂量，D = 0 标记未处理单位。剂量-反应
用 **B 样条**建模（默认线性），ACRT(d) 为拟合曲线的导数。逐点置信区间来自聚类自助法，
可选的**统一置信带**（`cband`）通过 multiplier bootstrap 同时覆盖整条曲线。

配合 **`covariates()`**，平行趋势放松为**条件（强）平行趋势**（CGBS，补充附录 SI.3，
命题 S3）：`ATT_x(d) = E[ΔY | X=x, D=d] − E[ΔY | X=x, D=0]`，聚合为
`ATT(d) = E_X[ATT_x(d) | D>0]`。这使估计量超出 R 参考实现——后者不支持协变量。

> **范围：**两期面板 / 重复截面，或**交错采纳**（`gvar()`）；连续剂量；ATT(d) 与
> ACRT(d)；以及**条件平行趋势**（`covariates()`）。CCK 数据驱动筛估计仍在路线图中
> （见 [CHANGELOG](CHANGELOG.md)）。

## 安装

```stata
net install contdid, from("https://raw.githubusercontent.com/niuniuhaoyu/contdid/main/") replace
```

需要 Stata 16 及以上，无额外依赖。

## 快速上手

```stata
* 载入两期面板：id、t (0/1)、y（结果）、d（连续剂量 ≥ 0）
use contdid_sim.dta, clear

* 线性剂量-反应（默认），并作图
contdid y, unit(id) time(t) dose(d) npoints(20) reps(999) seed(12345) graph

* 三次 B 样条 + 2 个内部节点 + 统一置信带
contdid y, unit(id) time(t) dose(d) degree(3) nknots(2) reps(999) seed(12345) cband graph

* 交错采纳：单位在不同时点受处理（g = 处理期，0 = 从未受处理）
contdid y, unit(id) time(t) dose(d) gvar(g) npoints(20) reps(999) seed(12345)

* 带协变量的条件平行趋势（超出 R 参考实现）
contdid y, unit(id) time(t) dose(d) covariates(x1 x2) cband
```

| 选项 | 默认 | 说明 |
|---|---|---|
| `unit(varname)` | — | 单位 / 面板标识（必填） |
| `time(varname)` | — | 时间变量，两期（必填） |
| `dose(varname)` | — | 连续处理剂量，≥ 0，0 = 未处理（必填） |
| `covariates(varlist)` | — | 条件（强）平行趋势的协变量（仅两期） |
| `dose_est_method(parametric\|dds)` | parametric | `parametric` = 固定 B 样条；`dds` = 数据驱动筛（LOO-CV 选节点；非精确 CCK） |
| `gvar(varname)` | — | 处理时点（0 = 从未受处理）；启用交错采纳 |
| `degree(#)` | 1 | B 样条次数（1 = 线性，2 = 二次，3 = 三次） |
| `nknots(#)` | 0 | 内部节点数（置于剂量分位点） |
| `knots(numlist)` | — | 显式内部节点（覆盖 `nknots`） |
| `npoints(#)` | 20 | 评估剂量的个数 |
| `reps(#)` | 999 | 聚类自助法重复次数（也是 multiplier 抽样次数） |
| `seed(#)` | 12345 | 随机种子 |
| `cluster(varname)` | `unit()` | bootstrap 的聚类变量 |
| `level(#)` | 95 | 置信水平（%） |
| `cband` | 关 | 统一置信带（sup-t，multiplier bootstrap） |
| `graph` | 关 | 绘制 ATT(d) 与 ACRT(d) |

## 方法

### 设定与记号

考虑一个由单位 *i* = 1, …, *n* 在时期 *t* = 1, …, *T* 上观测得到的面板。每个单位有一个
**剂量** *Dᵢ* ≥ 0（不随时间变化；*D* = 0 标记未处理单位）；（交错设计下）还有一个
**处理时点** *Gᵢ*（单位首次受处理的期数；*G* = 0 表示从未受处理）。令 *Yᵢₜ*(*d*) 为剂量
*d* 下的潜在结果。在无预期假设下，观测结果满足：处理前 *Yᵢₜ* = *Yᵢₜ*(0)，处理后
*Yᵢₜ* = *Yᵢₜ*(*Dᵢ*)。记 Δ*Y* = *Y*₂ − *Y*₁ 为处理前到处理后的变化。

### 识别

关键识别假设是**平行趋势**：

> 对所有 *d* > 0，E[*Y*₂(0) − *Y*₁(0) | *D* = *d*] = E[*Y*₂(0) − *Y*₁(0) | *D* = 0]。

在该假设下，**特定剂量的处理组平均处理效应**可识别（论文定理 3.1）：

> **ATT(*d* | *d*) = E[Δ*Y* | *D* = *d*] − E[Δ*Y* | *D* = 0]**。

这里 ATT(*d* | *d*) 是剂量组 *d* 的*局部*效应。*全局*剂量-反应
ATT(*d*) = E[*Y*₂(*d*) − *Y*₂(0) | *D* > 0] 需要更强的**强平行趋势**假设，它排除按剂量
水平的选择。**平均因果响应** ACRT(*d*) 是剂量-反应对 *d* 的导数。

### 估计

`contdid` 用 **B 样条**对剂量-反应做非参建模。在受处理单位（*D* > 0）上拟合

> Δ*Y* = *α* + Σₖ βₖ *B*ₖ(*D*) + *ε*,

其中 *B*ₖ 为 B 样条基函数（次数与内部节点由用户选定）。估计量为

> ATT(*d*) = *α̂* + Σₖ β̂ₖ *B*ₖ(*d*) − *m̂*₀， 其中 *m̂*₀ = mean(Δ*Y* | *D* = 0)，
>
> ACRT(*d*) = Σₖ β̂ₖ *B*′ₖ(*d*)。

取 `degree(1)` 与 `nknots(0)` 时基函数为线性，退化为简单的线性剂量估计
E[Δ*Y* | *D* = *d*] − E[Δ*Y* | *D* = 0]。

### 推断

- **逐点置信区间**来自聚类自助法（默认聚类 = 单位）：对单位有放回重抽样、重新估计，
  报告 bootstrap 估计的百分位区间。
- **统一置信带**（`cband`）同时覆盖整条剂量-反应曲线。它由样条估计量的影响函数配合
  multiplier bootstrap（Rademacher 权重）构造；sup-t 临界值为评估网格上学生化最大值
  的 (1 − *α*) 分位点。

### 交错采纳

使用 `gvar()` 时，单位在可能不同的时点首次受处理。对每个处理组 *g* 与处理后时期 *t*，
`contdid` 构造"2×2"式比较：cohort *g* 对比 not-yet-treated（*G* > *t* 或 *G* = 0），
差分结果 Δ*Y* = *Yₜ* − *Y*_{*g*−1}，用同样的 B 样条估计组-时剂量-反应
ATT(*g*, *t*; *d*)。报告的 ATT(*d*) 与 ACRT(*d*) 按组规模加权，聚合所有 (*g*, *t*) 对。

### 解读注意事项

- ATT(*d*) 是*局部*剂量效应；*全局*解读需要强平行趋势。
- ACRT(*d*) 是导数，估计精度低于 ATT(*d*)。
- 剂量 *D* 必须非负且单位内不变；需要存在 *D* = 0 的组。
- 与任何 DiD 设计一样，可信度取决于平行趋势（处理前趋势平坦）是否可信。

## 使用指南

### 最小两期示例

```stata
* 两期面板：id、t (0/1)、y、d（连续剂量）
use contdid_sim.dta, clear
contdid y, unit(id) time(t) dose(d)
```

`contdid` 报告两张表：**ATT(d)**（剂量-反应）与 **ACRT(d)**（其导数）。剂量-反应是
"接受剂量 *d* 相对剂量 0"的效应；ACRT(d) 是在剂量水平 *d* 处"剂量增加一单位"的效应。
加 `graph` 同时作图；加 `cband` 时置信带同时（而非逐点）覆盖曲线。

### 选择样条

- `degree(1)`（默认）施加线性剂量-反应——最简约，也是好的起点。
- `degree(2)`/`degree(3)` 配 `nknots(1)`–`nknots(3)` 让数据揭示曲率；可目视比较拟合，
  并报告你所依赖的设定。
- `knots(0.3 0.6)` 把节点固定在特定剂量值上，便于跨设定复现。

### 交错采纳

```stata
* 交错：g = 处理期（0 = 从未受处理）
contdid y, unit(id) time(t) dose(d) gvar(g) cband graph
```

### 复现已发表图

```stata
do examples/make_figure.do   // 生成 examples/contdid_figure.png
```

## 可复现性

- `examples/contdid_simdata.do` —— 生成模拟数据（`data/contdid_sim.dta`），其真实
  剂量-反应为 ATT(d) = 0.5·d。
- `examples/contdid_example.do` —— 一键示例。
- `examples/reference/dump_splines2.R` —— 导出 R `splines2` 的 B 样条基/导数，作为 Mata
  实现的黄金参考。
- `examples/reference/run_contdid_check.R` —— 线性估计量的独立 base-R 重实现；与 Stata
  输出匹配到 6 位小数。

## 引用

方法：

```bibtex
@article{callaway2024continuous,
  title   = {Difference-in-Differences with a Continuous Treatment},
  author  = {Callaway, Brantly and Goodman-Bacon, Andrew and Sant'Anna, Pedro H. C.},
  year    = {2024},
  note    = {NBER Working Paper 32117},
  url     = {https://doi.org/10.3386/w32117}
}
```

软件：

```bibtex
@software{niu2026contdid,
  title   = {contdid: Difference-in-Differences with a Continuous Treatment for Stata},
  author  = {Haoyu Niu},
  year    = {2026},
  version = {1.0.0},
  url     = {https://github.com/niuniuhaoyu/contdid}
}
```

## 许可

[AGPL-3.0](LICENSE)

---

[English](README.md) | [简体中文](README_zh.md)
