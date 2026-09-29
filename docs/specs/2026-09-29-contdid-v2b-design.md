# contdid v2b — 交错处理（Staggered Adoption）Spec

> 日期：2026-09-29
> 状态：Draft（执行中）
> 作者：Haoyu Niu
> 范围：把 contdid 从「两期」扩展到「多期交错处理」，保留连续剂量 + B 样条 + ACRT

---

## 1. 目标

在现有两期估计（v2a）基础上，支持**交错处理**：不同单位在不同时间开始接受（连续剂量的）处理。
产出：

- **组-时剂量反应 ATT(g,t; d)**：cohort g 在时期 t 的剂量反应曲线。
- **剂量反应聚合 ATT(d)**：跨 (g,t) 加权平均的剂量反应。
- **事件研究聚合**（可选，后续）。
- 推断：聚类自助（沿用 v2a 的 bootstrap，扩展到交错）。

## 2. 设定与记号

- 面板：unit i，时间 t ∈ {1,...,T}。
- 处理时点 G_i（首次被处理的时期；G=0 表示从未处理）。
- 剂量 D_i（连续，≥0，单位内时间不变；D=0 为未处理）。
- 结局 Y_{i,t}。

## 3. 识别与估计（对齐 R `cont_two_by_two_subset` + `cont_did_acrt`）

对每个 (g, t)，g 为处理 cohort、t ≥ g（处理后时期）：

1. **子样本**：cohort g（G=g）∪ 控制组（`notyettreated`：G>t 或 G=0；或 `nevertreated`：G=0）。
2. **基期**：`g-1`（处理前一期；默认 varying）。
3. **ΔY = Y_t − Y_{g-1}**。
4. **有效剂量**：D_active = D · 1[G=g]（仅 cohort g 的剂量生效；控制组 D_active=0）。
5. **估计**：在子样本上，用 B 样条拟合 ΔY ~ 剂量（treated = G=g 且 D>0），
   ATT(g,t; d) = fitted(d) − mean(ΔY | D_active=0)；ACRT(g,t; d) = 导数。

## 4. 聚合

- **剂量反应**：ATT(d) = Σ_{g,t} w_{g,t} · ATT(g,t; d)，权重按 cohort 规模 / 剂量分布（对齐论文的 dose-density 权重）。
- **事件研究**（后续）：对每个 event time e = t−g，平均 ATT(g,g+e; d)。

## 5. 命令语法（扩展）

```
contdid depvar, unit(id) time(t) dose(D) gvar(G) [degree() nknots() knots() npoints() reps() seed() cluster() cband graph]
```

- 新增 `gvar(varname)`：处理时点（0 = 从未处理）。
- 有 `gvar` 时进入交错路径；无 `gvar` 时保持两期路径（向后兼容）。

## 6. 实现策略（分阶段）

1. **组-时子样本 + 估计**：`_contdid_gt` Mata 函数，对给定 (g,t) 子样本估计 ATT(g,t;d)、ACRT(g,t;d)。
2. **聚合**：跨 (g,t) 加权平均 → 剂量反应 ATT(d)。
3. **事件研究**：按 event time 聚合。
4. **推断**：聚类自助扩展到交错。

## 7. 验证

- 交错 DGP（多 cohort、连续剂量、已知真实 ATT(d)），恢复真实剂量反应。
- 与 R `contdid`（`aggregation="dose"`, `control_group="notyettreated"`）对拍（若可装）。
- 向后兼容：无 `gvar` 时两期结果不变。

## 8. 边界

- v2b 首版做「剂量反应聚合」；事件研究、交错统一置信带留待后续迭代。
