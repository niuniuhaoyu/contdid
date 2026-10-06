# 研究笔记：Stata↔R 对拍管线验证 + R `contdid` 的 `dvals` 坑

> 日期：2026-10-06
> 目的：确认「Stata `contdid` ↔ R `contdid`」跨语言对拍管线可用，作为 v3a/v3b 验证的前置。
> 环境：Stata 19 MP（`D:\Stata\StataMP-64.exe`）；R 4.6.1（conda，`D:\OpenCode\环境脚本\Rscript.cmd`）；R 包 `contdid` 0.1.1（GitHub 源码安装）。

## 1. 做法

同一份模拟数据 `examples/reference/contdid_sim.csv`（真实 ATT(d)=0.5d），三方估计，评估点统一铺满 treated 剂量范围：

- **Stata**：`contdid y, unit(id) time(t) dose(d) degree(1) nknots(0) npoints(9)`，导出 `r(attd)`（脚本 `export_attd.do`）
- **base-R**：`lm(dy ~ d)` on treated − `mean(dy | d==0)`
- **R `contdid`**：`cont_did(..., degree=1, num_knots=0, aggregation="dose", control_group="nevertreated", dvals=<铺满剂量>)`

## 2. 结果

| 对比 | max\|差\| |
|---|---|
| Stata vs base-R | 1.18e-09 |
| Stata vs R `contdid`（dvals 铺满剂量范围） | 1.75e-09 |

残差来自 `stata_attd.csv` 仅保留约 15 位有效数字，**估计量本身一致到机器精度**（R 与 base-R 直接比时为 5.6e-17）。
脚本：`examples/reference/export_attd.do`、`examples/reference/compare_attd.R`。

## 3. 发现的坑（重要）

R `contdid` 0.1.1 的内部函数 `cont_did_acrt` 里：

```r
bs      <- splines2::bSpline(dose[dose>0], degree, knots)   # 边界节点 = range(dose)
bs_reg  <- lm(dy ~ ., data = bs)
bs_grid <- splines2::bSpline(dvals,        degree, knots)   # 边界节点 = range(dvals) !
att.d   <- predict(bs_reg, newdata = bs_grid) - mean(dy[dose==0])
```

**回归基用 `range(dose)` 作边界节点，预测基却默认用 `range(dvals)`。** 当 `range(dvals) != range(dose[dose>0])` 时，预测基被重新缩放，报告出的 `att.d` 被**系统性重标度**。

触发条件（很容易踩）：
- 传入自定义 `dvals` 未铺满剂量范围（如 `seq(0.1,0.9,0.1)`）→ 斜率被放大 `range(dose)/range(dvals)` 倍；
- **连默认 `dvals` 也会触发**：默认是剂量的分位数（约 10%–99%），其范围仍窄于完整剂量范围。

**规避**：与 R `contdid` 对拍时，**必须让 `dvals` 精确铺满 `range(dose[dose>0])`**（用 `seq(min,max,length.out=...)`）。我们的 Stata 默认网格（dmin→dmax 均匀）本就铺满，不受影响。

## 4. 对项目的意义

- ✅ 对拍管线打通：Stata / R 都能在同一数据上跑，数值可比。
- ✅ v1「对照 R 到 6 位小数」的正确性成立（对齐的是论文那套线性估计量，R 源码复现一致）。
- ⚠️ v3a/v3b 做 R 对拍时，**dvals 必须铺满剂量范围**，否则会把 R 的输出重标度误判为"我们不一致"。
- 备注：R 该行为是否为 bug 未与作者确认；我们的 Stata 实现不受此影响。
