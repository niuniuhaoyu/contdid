# contdid v3 — 协变量 + CCK sieve Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 给 contdid 补上 README roadmap 剩余两项——(v3a) 条件平行趋势的协变量调整、(v3b) CCK 数据驱动 sieve——并保持对 v1.0.0 的完全向后兼容。

**Architecture:** 在扁平包布局下扩展 `contdid.ado`：新增 `covariates()` 与 `dose_est_method()` 分支；协变量调整与影响函数对齐 R `cont_did(xformla=~X)`；CCK 用 Mata 实现（或包 `npiv`）数据驱动 sieve，独立文件 `contdid_cck.mata`。正确性靠"条件平行趋势 DGP 模拟 + R 对拍"双重验证。

**Tech Stack:** Stata 19 MP（`D:\Stata\StataMP-64.exe`）+ Mata；R（仅对拍，需补装）；Git/GitHub。

**Spec:** `contdid/docs/specs/2026-10-06-contdid-v3-design.md`

## Global Constraints

- Stata 16+；批处理跑法见 `SOP/SOP_环境与工具链.md`
- 向后兼容：不带新选项时，结果与 v1.0.0 **逐位一致**
- 估计量/影响函数以论文 + R `contdid` 源码为准，**禁止凭印象硬编**
- 包布局**扁平**（.ado/.mata/.sthlp 在根目录），支持 `net install`
- 许可证 AGPL-3.0；文档英文；v3a 升 v0.4.0、v3b 升 v0.5.0
- 遵守 `SOP/SOP_AI辅助编程核验.md`：对照基准 → 查依赖 → 实跑留痕 → 手工复核关键数字

## Review Focus

最易踩坑、必须逐一配测试：

1. **协变量缺失/常数列/完全共线**——应报错而非静默丢样本或给错结果。（v3a Task 2）
2. **影响函数漏掉协变量部分**——sup-t 临界值会偏小、带过窄、覆盖不足。（v3a Task 3）
3. **无条件 vs 条件 DGP 造错**——DGP 必须保证"无条件平行趋势不成立、条件成立"，否则测不出协变量的作用。（v3a Task 1）
4. **向后兼容**——`covariates()` 缺省、`dose_est_method(parametric)` 缺省时必须与 v1.0.0 逐位一致。（v3a Task 2、v3b Task 1）
5. **CCK 的收敛/维度选择**——sieve 维度准则实现错会导致不收敛或过拟合；需有收敛失败时的明确报错。（v3b Task 2）
6. **R 环境缺失**——对拍步骤依赖 R；本机当前未装 R，须先补装或降级为"模拟自洽 + 大样本收敛"验证并标注。（v3a Task 4、v3b Task 4）

---

## v3a — 协变量（条件平行趋势）

### Task 1: 研究 + 条件 DGP

**Files:**
- Create: `contdid/docs/research-notes-v3a.md`
- Create: `contdid/examples/v3a_simdata.do`

**Interfaces:**
- Consumes: 无
- Produces: `research-notes-v3a.md`（协变量下 ATT(d)/影响函数精确定义）；可生成"仅条件平行趋势成立"的模拟数据

- [ ] **Step 1: 精读 CGBS 论文协变量章节 + R `cont_did` 的 `xformla` 路径**，定死：协变量进入的估计形式（回归调整/IPW/DR）、ATT(d) 公式、影响函数构成。写入 notes。
- [ ] **Step 2: 写条件 DGP**（`v3a_simdata.do`）：X 同时影响剂量与结局趋势，使无条件平行趋势失败、条件成立；已知真实 ATT(d)。
- [ ] **Step 3: 验证 DGP**：无条件估计显著有偏、`covariates(X)` 估计接近真值（小样本粗查即可）。
- [ ] **Step 4: Commit** `docs: v3a 研究笔记 + 条件平行趋势 DGP`

### Task 2: `covariates()` 估计 + 输入校验

**Files:**
- Modify: `contdid.ado`
- Test: `contdid/examples/_test_v3a_cov.do`

**Interfaces:**
- Consumes: `research-notes-v3a.md`
- Produces: 选项 `covariates(varlist)`；`r(attd)`/`r(acrt)` 语义不变

- [ ] **Step 1: 写失败测试**：条件 DGP 下 `contdid y, ..., covariates(X)` 恢复真实 ATT(d)（误差 < 阈值）；无条件版本超阈值。
- [ ] **Step 2: 运行，预期失败**（选项未实现）。
- [ ] **Step 3: 实现输入校验**（Review Focus #1）：协变量缺失、常数列、与剂量完全共线 → `display as error` + `exit`。
- [ ] **Step 4: 实现协变量估计**（按 notes 公式）。
- [ ] **Step 5: 运行，预期通过**；并跑 v1.0.0 旧测试确认无 `covariates()` 时逐位一致（Review Focus #4）。
- [ ] **Step 6: Commit** `feat: covariates() 条件平行趋势估计`

### Task 3: 含协变量的影响函数 + cband

**Files:**
- Modify: `contdid.ado`
- Test: `contdid/examples/_test_v3a_cband_cov.do`

- [ ] **Step 1: 写覆盖检验**：重复小样本，统计含协变量下 sup-t 统一带对整条 ATT(d) 的覆盖率 ≥ 名义。
- [ ] **Step 2: 实现影响函数扩展**（Review Focus #2）。
- [ ] **Step 3: 运行覆盖检验**，预期覆盖 ≥ 0.95（或落入 [0.90, 1.0]）。
- [ ] **Step 4: Commit** `feat: 含协变量的统一置信带`

### Task 4: R 对拍 + 文档 + 发布（v0.4.0）

- [ ] **Step 1: 补装 R**（见下方"环境前置"），写 `examples/reference/run_contdid_cov_check.R`，与 R `cont_did(xformla=~X)` 逐位对拍（≤1e-6）。若无法补装 R，则记录降级验证方式（Review Focus #6）。
- [ ] **Step 2: 更新 README / sthlp / CHANGELOG**，版本升 v0.4.0。
- [ ] **Step 3: 全量回归**：v1/v2a/v2b 所有旧测试 + 新测试全通过。
- [ ] **Step 4: Commit + push**（`git push origin main`）。

---

## v3b — CCK 数据驱动 sieve

### Task 1: 研究（论文 sieve 章节 + R `npiv`）

**Files:**
- Create: `contdid/docs/research-notes-v3b.md`

- [ ] **Step 1: 精读 CGBS 论文数据驱动 sieve 章节**，以及 Chen–Christensen–Kankanala (ReStud 2025) 方法要点。
- [ ] **Step 2: 精读 R `contdid` 的 `dose_est_method="cck"` 与 `npiv` 源码**，定死 sieve 维度选择准则、估计步骤、推断。
- [ ] **Step 3: 写 notes**；评估"纯 Mata 实现 vs 依赖 `npiv`（R/其他）"的取舍，给出建议。
- [ ] **Step 4: Commit** `docs: v3b CCK 研究笔记`

### Task 2: `dose_est_method(cck)` 最小可用实现

**Files:**
- Create: `contdid/contdid_cck.mata`
- Modify: `contdid.ado`
- Test: `contdid/examples/_test_v3b_cck.do`

- [ ] **Step 1: 写失败测试**（选项未实现）。
- [ ] **Step 2: 实现 sieve 维度自适应选择 + dose-specific 估计**（按 notes）。
- [ ] **Step 3: 收敛失败等异常路径给出明确报错**（Review Focus #5）。
- [ ] **Step 4: 运行测试**：大样本下 `cck` 与 `parametric`（足够灵活）收敛到同一曲线。
- [ ] **Step 5: Commit** `feat: CCK 数据驱动 sieve（dose_est_method）`

### Task 3: 验证 + 文档 + 发布（v0.5.0）

- [ ] **Step 1: 与 R `dose_est_method="cck"` 对拍**（需 R）；无 R 则记录降级验证。
- [ ] **Step 2: 更新 README / sthlp / CHANGELOG**，版本升 v0.5.0。
- [ ] **Step 3: 全量回归 + push**。

---

## 环境前置（不在本 plan 之外，但阻塞对拍）

- [ ] **安装 R**（本机当前缺失）——用于 v3a/v3b 与 R `contdid` 对拍。装好后写入 `SOP/SOP_环境与工具链.md`。

---

## Self-Review 结论

- **Spec 覆盖**：§3 做/不做 → v3a/v3b 分档；§4 语法 → v3a Task 2 / v3b Task 2；§5 估计与推断 → v3a Task 2/3、v3b Task 2；§7 验证 → 各 Task；§8 交付 → Task 4/3。
- **占位符**：无 TBD；精确公式通过 Task 1 research-notes 落定（前置研究交付物）。
- **类型一致性**：`r(attd)`/`r(acrt)` 保持 npoints×5 不变；新选项缺省时走旧路径。
- **Review Focus**：6 条分别落到 v3a Task 1/2/3、v3b Task 1/2、以及"环境前置"，均有对应测试或明确降级说明。
