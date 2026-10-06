# contdid v3 — 协变量（条件平行趋势）+ CCK 数据驱动 sieve（Spec）

> 日期：2026-10-06
> 状态：待用户审阅（Draft）
> 作者：Haoyu Niu
> 范围：在 v1/v2a/v2b 之上补齐 README roadmap 的两项——**协变量**与 **CCK 数据驱动 sieve**
> 依据：Callaway, Goodman-Bacon & Sant'Anna (2024/2025) *Difference-in-Differences with a Continuous Treatment*；R 参考实现 `bcallaway11/contdid`

---

## 1. 背景与目标

contdid 已到 v1.0.0：

- v1：两期、线性-in-dose、聚类自助逐点 CI
- v2a：B 样条 ATT(d) + ACRT(d) + 统一置信带（cband）
- v2b：交错处理（`gvar()`）组-时剂量反应聚合

README 的 roadmap 还剩两项未做，本 spec 覆盖它们。二者难度差一个量级，**拆成 v3a / v3b 两档**，避免被 CCK 拖死 v3a。

- **v3a｜协变量（conditional parallel trends）**：识别从"无条件平行趋势"放宽到"条件平行趋势"。这是绝大多数真实应用的刚需（单位特征不同、趋势本就不同）。
- **v3b｜CCK 数据驱动 sieve**：把"节点个数/次数由用户拍定"升级为"由数据自适应选择"（Chen, Christensen & Kankanala, ReStud 2025；R 里经 `npiv` 包实现，`dose_est_method="cck"`）。**非参数 IV sieve，难度最高**，故独立立项。

---

## 2. 方法概述（高层）

### 2.1 v3a 协变量

- 平行趋势改为**条件**版本：给定协变量 X，处理组与对照组在无处理反事实下趋势相同。
- 命令提供 `covariates(varlist)`；估计量对齐 R `cont_did(..., xformla=~X)` 的默认实现（先精读源码定死：是回归调整、IPW，还是双重稳健）。
- 推断沿用 v2a 的影响函数 + multiplier bootstrap 框架；**带协变量后影响函数的构成必须与 R 源码逐行对齐**（这是最容易错的地方）。

### 2.2 v3b CCK sieve

- 剂量反应不再预设 B 样条次数/节点，而是**数据驱动**地选择 sieve 维度，得到 dose-specific effects 的非参数估计。
- 对齐 R `dose_est_method="cck"`（其底层为 `npiv` + Chen–Christensen–Kankanala 方法）。
- 交付形态：新增 `dose_est_method(parametric|cck)`，默认 `parametric`（保持向后兼容）。

> 精确公式、收敛率、节点选择准则**以论文原文 + R `npiv` 源码为准**；实现前先精读，产出 `research-notes-v3a.md` / `research-notes-v3b.md`，禁止凭印象硬编。

---

## 3. 范围

### 做（v3a）

- `covariates(varlist)` 选项（两期路径；交错路径的协变量见第 9 节）
- 条件平行趋势下的 ATT(d) / ACRT(d)（不预设具体 DR 配方，以 R 源码为准）
- 统一置信带（cband）在含协变量下仍可用
- 输入校验：协变量含缺失 / 常数列 / 与处理完全共线时明确报错

### 做（v3b，独立）

- `dose_est_method(cck)` 选项
- 数据驱动 sieve 维度的选择与估计
- 与 `parametric` 路径并存、可切换

### 不做（v3）

- 协变量 × 交错（`gvar()`）的组合（留待 v3c，见第 9 节）
- 事件研究聚合（v2b 已列后续）
- 面板固定效应（FE）版本的协变量调整
- 协变量下的 ACRT **统一带**（先做点估计与逐点带）

---

## 4. 命令语法（v3）

```
contdid depvar [if] [in], unit(varname) time(varname) dose(varname) ///
        [covariates(varlist) dose_est_method(parametric|cck) ///
         degree(#) knots(numlist) nknots(#) npoints(#) ///
         reps(#) seed(#) cluster(varname) level(#) cband graph]
```

| 新增选项 | 默认 | 说明 |
|---|---|---|
| `covariates(varlist)` | 无 | 条件平行趋势的协变量 X（v3a） |
| `dose_est_method(...)` | `parametric` | 剂量反应估计法：`parametric`=B 样条（现状）；`cck`=数据驱动 sieve（v3b） |

- 向后兼容：不带新选项时，结果与 v1.0.0 完全一致（逐位）。

---

## 5. 估计量与推断

### 5.1 v3a（协变量）

1. 精读 R `cont_did` 在 `xformla` 非空时的估计路径（回归调整 / IPW / DR），定死 ATT(d) 公式。
2. 精读其影响函数（`estfun` / `bread` 的协变量扩展），定死 sup-t 临界值构成。
3. 拟步骤：在条件平行趋势下，用 X 调整 ΔY 与剂量关系（具体形式以 notes 为准），ATT(d) = 调整后 fitted(d) − 调整后基线。
4. 推断：影响函数 + multiplier bootstrap；点对点 CI + `cband` 统一带。

### 5.2 v3b（CCK）

1. 精读论文数据驱动 sieve 章节 + R `npiv` 的 CCK 实现。
2. 用 Mata 实现（或包一层）sieve 维度自适应选择与 dose-specific 估计。
3. 推断：按 notes 对齐。
4. 验收：在已知 DGP 下恢复真实剂量反应；若可装 R，则与 `dose_est_method="cck"` 对拍。

---

## 6. 包结构（增量）

- 修改：`contdid.ado`（新增 `covariates()`、`dose_est_method()` 分支）
- 新增：`contdid_cck.mata`（或内联 Mata：CCK sieve 估计）
- 修改：`contdid.sthlp`、`README.md`、`CHANGELOG.md`
- 新增测试：`examples/_test_v3a_cov.do`、`examples/_test_v3a_cband_cov.do`、`examples/_test_v3b_cck.do`
- 新增参考：`examples/reference/run_contdid_cov_check.R`、`examples/reference/run_cck_check.R`（需 R 环境）
- 新增研究笔记：`docs/research-notes-v3a.md`、`docs/research-notes-v3b.md`

> 注意：包布局保持**扁平**（.ado/.mata/.sthlp 在根目录），以支持 `net install`（见提交 `47854d4`）。

---

## 7. 验证方案

1. **协变量正确性（v3a）**：造一个"平行趋势只有条件于 X 才成立"的 DGP；无条件估计有偏、有协变量估计无偏（误差随 N 收敛）。
2. **R 对拍（v3a）**：与 R `cont_did(xformla=~X)` 点估计逐位对齐（容差 ≤ 1e-6）——**需要 R**（当前本机未装，见"还要做"）。
3. **统一带覆盖（v3a）**：Monte Carlo 检验 sup-t 带覆盖 ≥ 名义水平。
4. **CCK 正确性（v3b）**：与 R `dose_est_method="cck"` 对拍；无 R 时至少用模拟 DGP 的自洽性 + 与 `parametric` 在大样本下收敛到同一曲线。
5. **向后兼容**：不带新选项时全部旧测试通过、结果逐位一致。

---

## 8. 交付物与验收标准

- [ ] `covariates()` 可用且有文档（v3a）
- [ ] 条件平行趋势 DGP 下估计无偏（模拟证据）
- [ ] 与 R `contdid` 协变量路径对齐（需 R）
- [ ] `dose_est_method(cck)` 可用且有文档（v3b）
- [ ] 全部旧测试回归通过、向后兼容
- [ ] README / sthlp / CHANGELOG 更新；版本升至 v0.4.0（v3a）与 v0.5.0（v3b）
- [ ] push 到 `niuniuhaoyu/contdid`

---

## 9. 边界与未来

- ❌ v3 不做"协变量 × 交错处理"组合 → 列为 **v3c**。
- ❌ v3 不做协变量下的 ACRT 统一带。
- ✅ v3b 的 CCK 若工作量过大，可先交付"接口 + 研究笔记 + 最小可用实现"，逐步补全。

---

## 10. 实施顺序

1. 研究：精读 CGBS 论文协变量章节 + R `cont_did` 的 `xformla` 路径 → `research-notes-v3a.md`
2. 实现 v3a：`covariates()` 估计 + 校验 + 测试（TDD）
3. 实现 v3a：含协变量的影响函数与 `cband`
4. 文档 + 回归 + push（v0.4.0）
5. 研究：论文 sieve 章节 + R `npiv` → `research-notes-v3b.md`
6. 实现 v3b：`dose_est_method(cck)` + 测试
7. 文档 + 回归 + push（v0.5.0）
