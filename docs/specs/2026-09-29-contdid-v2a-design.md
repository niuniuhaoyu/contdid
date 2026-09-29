# contdid v2a — B 样条 + ACRT + 统一置信带（Spec）

> 日期：2026-09-29
> 状态：待用户审阅（Draft）
> 作者：Haoyu Niu
> 范围：仍为两期（共同处理时点）；交错处理属 v2b，另立项

---

## 1. 背景与目标

v1 已实现两期、线性-in-dose 的 ATT(d) + 聚类自助逐点 CI。v2a 在不改变设定（两期）的前提下，
把剂量响应建模从「线性」升级为「B 样条」，并补上：

- **ACRT(d)**：平均因果响应（ATT(d) 对剂量的导数，即边际效应）
- **统一置信带**：整条 ATT(d)/ACRT(d) 曲线同时覆盖的置信带（非逐点）

目标：两期设定下追平 R `contdid` 的核心估计功能（B 样条 + ATT/ACRT + 统一带）。

---

## 2. 估计量（以 R `cont_did_acrt` 源码为准，已核实）

给定 treated 单位（D>0），对 ΔY 在 B 样条基上做 OLS：

```
ΔY = α + Σ_{k=1..K} β_k B_k(D)  （treated 子样本，含截距）
```

则对评估剂量 d：

- **ATT(d) = α + Σ β_k B_k(d) − m0**，其中 m0 = mean(ΔY | D=0)
- **ACRT(d) = Σ β_k B'_k(d)**（B'_k 为样条基导数，截距项无贡献）

> 与 R 源码逐行对应：
> `bs_reg <- lm(dy ~ ., data=bs)`；`att.d <- predict(...) - mean(dy[dose==0])`；
> `acrt.d <- dbs(dvals) %*% coef(bs_reg)[-1]`。

## 3. B 样条基

- 在 Mata 中实现 **Cox-de Boor 递归**生成 B 样条基 B_k(d) 及导数 B'_k(d)。
- 参数：`degree`（默认 1 = 线性，退化为 v1），`knots`（内部节点）。
- 节点位置默认按**分位数**（R `choose_knots_quantile`：在 treated 剂量分布的分位数取 num_knots 个内部节点）。
- 边界节点取 treated 剂量的 min/max；基含截距（OLS 另加常数项，与 R `lm(dy ~ .)` 一致）。

## 4. 推断

- **逐点 CI**（保留 v1 的聚类自助）：对 ATT(d)、ACRT(d) 各自出逐点置信区间。
- **统一置信带**（新增 `cband` 选项）：multiplier bootstrap 计算 sup-t 临界值，
  使整条曲线（所有评估点）同时被覆盖。
  - 影响函数结构移植自 R 源码：`Xe = estfun(bs_reg)`（OLS 得分）、`bread = bread(bs_reg)`；
    ACRT 的影响函数 = (导数基·系数 − 总体 ACRT) + Xe·bread·(0, mean(导数基))'。
  - multiplier bootstrap：独立标准正态/指数乘子 × 影响函数，取 ||sup-t|| 的 (1−α) 分位数。

## 5. 命令语法（v2a）

```
contdid depvar [if] [in], unit(varname) time(varname) dose(varname) [options]
```

| 选项 | 默认 | 说明 |
|---|---|---|
| `degree(#)` | 1 | B 样条次数（1=线性，2=二次，3=三次） |
| `knots(numlist)` | 无 | 内部节点（指定时优先于 nknots） |
| `nknots(#)` | 0 | 内部节点个数（按分位数放置；0=全局多项式） |
| `npoints(#)` | 20 | 评估点数量 |
| `reps(#)` / `seed(#)` / `level(#)` / `cluster()` | 同 v1 | 自助法 |
| `cband` | off | 输出统一置信带（sup-t） |
| `graph` | off | 画 ATT(d) 与 ACRT(d)（含带） |

- 向后兼容：不指定 `degree`/`knots`/`nknots` 时，等价于 v1 线性估计。
- `nknots(0)` = 全局多项式（R 默认），`degree(1)` + `nknots(0)` = 线性（v1）。

## 6. 返回值

- `r(attd)`：npoints × 5（d、ATT、se、lb、ub）
- `r(acrt)`：npoints × 5（d、ACRT、se、lb、ub）
- `r(b)` / `r(beta)`：样条系数；`r(m0)`：未处理基线
- `cband` 时额外返回 `r(crit_att)`、`r(crit_acrt)`（sup-t 临界值）及带上下界

## 7. 包结构（增量）

- 修改：`ado/contdid.ado`（新增 degree/knots/cband 分支）
- 新增：`ado/mata_bspline.mata`（或内联 Mata：B 样条基 + 导数）
- 修改：`sthlp/contdid.sthlp`、`README.md`、`CHANGELOG.md`
- 新增测试：`examples/_test_bspline.do`、`examples/_test_acrt.do`、`examples/_test_cband.do`
- 新增参考：`examples/reference/run_splines2_check.R`（用 R `splines2` 对拍 B 样条基/导数）

## 8. 验证方案

1. **B 样条基正确性**：与 R `splines2::bSpline` / `dbs` 在同一节点/次数下逐项对拍（容差 1e-6）。
2. **ATT(d)/ACRT(d) 正确性**：
   - 二次剂量效应 DGP（真实 ATT(d) 非线性）下，`degree(2)`/`degree(3)` 恢复真实曲线（误差随样本增大收敛）。
   - 与 R 独立实现（`splines2` + `lm`）对拍数值。
3. **ACRT 正确性**：线性 DGP 下 ACRT(d) ≈ 常数（= 线性斜率）；二次 DGP 下 ACRT(d) ≈ 真实导数。
4. **统一带覆盖**：Monte Carlo 检验 sup-t 带的覆盖 ≥ 名义水平（如 95%）。
5. **向后兼容**：`degree(1) nknots(0)` 的结果与 v1 完全一致。

## 9. 边界与不做（v2a）

- ❌ 交错处理、多期（v2b）
- ❌ CCK 非参数筛、协变量
- ❌ 事件研究聚合
- ✅ 保留 v1 的线性逐点路径作为默认

## 10. 实施顺序（预估任务）

1. 研究：精读 R 源码 `cont_did_acrt` + `splines2` 的 B 样条/导数定义；装 `splines2` 作对拍
2. Mata B 样条基 + 导数（TDD：与 `splines2` 对拍）
3. B 样条估计 ATT(d)（替换线性，向后兼容）
4. ACRT(d)
5. 统一置信带（影响函数 + multiplier bootstrap）
6. 文档 + 测试 + push（v0.2.0）
