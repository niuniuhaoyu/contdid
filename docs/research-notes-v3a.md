# 研究笔记 v3a：协变量 / 条件（强）平行趋势

> 日期：2026-10-06
> 目的：为 `contdid` v3a 的 `covariates()` 选项定死识别与估计。
> 结论先行：**R `contdid` 0.1.1 不支持协变量**（`cont_did.R:126` 直接 `stop`）；方法在论文**补充附录 SI.3** 有正式结果。故 v3a 是**超出 R 参考实现的贡献**，对拍只能靠模拟。

## 1. R 参考实现的现状（已核实）

- `bcallaway11/contdid` 源码 `R/cont_did.R` 第 126 行：
  `if (xformula != ~1) stop("covariates not currently supported, please use xformula=~1")`
- 即：R 包参数 `xformula` 存在但**未实现**；`est_method` 也被显式拒绝。
- ⇒ 我们无法用 R 对拍协变量；**改为模拟验证**（条件平行趋势 DGP：无条件估计有偏、条件估计一致）。

## 2. 论文的识别（SI.3，Prop S3）

**Assumption SPT-X（条件强平行趋势）**：对所有 d∈𝒟，

```
E[ Y_{t=2}(d) − Y_{t=1}(0) | X=x, D>0 ] = E[ Y_{t=2}(d) − Y_{t=1}(0) | X=x, D=d ]
```

**Proposition S3**：在 Assumption 1/2/3/4(a) + SPT-X 下，

```
ATT_x(d) = E[ ΔY | X=x, D=d ] − E[ ΔY | X=x, D=0 ]
```

**推论（ACRT_x）**：

```
ACRT_x(d) = ∂/∂d E[ ΔY | X=x, D=d ]      （无选择偏误项）
```

**汇总到无条件曲线**（对处理组 X 分布取期望）：

```
ATT(d) = E_X[ ATT_x(d) | D>0 ]
```

## 3. 估计量设计（推荐）

设 B_k(D) 为剂量 D 上的 B 样条基（沿用现有 `bspline_basis`）。用一个**含交互的线性 sieve** 拟合全体样本：

```
ΔY = Σ_k β_k B_k(D) + Σ_j γ_j X_j + Σ_{k,j} δ_{kj} B_k(D)·X_j + ε
```

则

```
ATT_x(d) = Σ_k ( β_k + Σ_j δ_{kj} x_j ) · ( B_k(d) − B_k(0) )
ATT(d)   = Σ_k ( β_k + Σ_j δ_{kj} x̄_j ) · ( B_k(d) − B_k(0) ),   x̄ = E[X | D>0]
ACRT(d)  = Σ_k ( β_k + Σ_j δ_{kj} x̄_j ) · B'_k(d)
```

要点 / 坑：
1. **必须有 D×X 交互**。若只放可加项 X，则 x 项在 `m(x,d)−m(x,0)` 中相消，退化为无条件估计（等于没做）。
2. **基的域**：现有实现基的边界节点取 `range(dose|D>0)`，因此 `B_k(0)` 是**外推**。协变量路径建议基的域取 `[0, dmax]`，使 `B_k(0)` 有定义（需在 U 型验证：全局无 knots/degree=1 时退化为线性）。
3. **基线 `D=0` 组**：`m(x,0)` 由同一模型在 D=0 处求值（非"控制组均值"），这与 Prop S3 的 `E[ΔY|X=x,D=0]` 一致。
4. 交错路径（`gvar()`）暂不做协变量 → 留 v3c。
5. 推断：沿用 v2a 的 IF + multiplier bootstrap；协变量下 IF 需把 X 一起纳入 X'X 与影响函数（Review Focus）。

> 备选估计量（未采用，记录备查）：分别在 treated / untreated 拟合 E[ΔY|X,D=d]、E[ΔY|X,D=0] 再相减（"separate regressions"）；或 IPW/DR（类比 Sant'Anna–Zhao 二元情形）。选"含交互的联合 sieve"因其与现有 Mata 结构最贴合、且天然给出 ACRT。

## 4. 验证方案（无 R 对拍）

1. **条件 DGP**（`examples/v3a_simdata.do`）：X 同时影响处理选择与（未处理）趋势，使**无条件平行趋势失败、条件成立**；真实 `ATT(d)=(0.5+0.2·x̄)·d`。
2. **R 原型**（本笔记配套脚本，用 `lm` 实现上式）：确认 `ATT(d)` 恢复真值，且**无条件估计显著有偏**。
3. 移植进 Mata 后：Stata 结果 == R 原型（逐位），并跑覆盖率检验。

## 5. 参考

- Callaway, Goodman-Bacon & Sant'Anna (2024/2025), *Difference-in-Differences with a Continuous Treatment*，**Supplemental Appendix SI.3**（SPT-X, Prop S3）。
- 现有实现：`contdid.ado` 的 `_contdid_fit`（`bspline.mata`）。
