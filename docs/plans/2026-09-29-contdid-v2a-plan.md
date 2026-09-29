# contdid v2a — B 样条 + ACRT + 统一置信带 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax.

**Goal:** 在 contdid（两期）中把剂量响应升级为 B 样条，新增 ACRT(d) 边际效应与统一置信带，向后兼容 v1 线性版本。

**Architecture:** 在 `contdid.ado` 中新增 `degree()`/`knots()`/`nknots()`/`cband` 分支；B 样条基与导数在 Mata 中实现（Cox-de Boor 递归）；估计量仍为「treated 上 OLS ΔY~样条基，ATT(d)=fitted(d)−m0」，ACRT(d)=导数基·系数；统一带用 multiplier bootstrap sup-t。

**Tech Stack:** Stata（ado + Mata）、R（仅 `splines2` 包作对拍）、Git/GitHub。

**Spec:** `contdid/docs/specs/2026-09-29-contdid-v2a-design.md`（本计划从该 spec 展开；执行者两份都读）

## Global Constraints

- Stata 16+
- 向后兼容：`degree(1)` + `nknots(0)` 必须与 v1 线性结果完全一致（逐位）
- 估计量、样条基约定、影响函数结构以 R `contdid` 源码（`cont_did_acrt`）与 R `splines2` 为准
- 许可证 AGPL-3.0；文档英文；版本升至 0.2.0
- v2a 只做两期；交错处理不做

## Review Focus

以下输入/情形最易踩坑，必须逐一配测试（括号标注归属任务）：

1. **样条基约定不一致**——R `splines2::bSpline` 默认 `intercept=FALSE`（K = nknots+degree 个基），而 `lm` 再补截距；Mata 实现必须与 `splines2` 在同一 degree/knots 下逐位一致，否则 ATT(d) 整体错。（Task 2）
2. **边界与内部节点**——内部节点个数 `nknots` 与总节点（含边界）的关系、边界节点取 treated 剂量的 min/max，必须与 R `choose_knots_quantile` 对齐。（Task 1、2）
3. **ACRT 不含截距**——ACRT(d)=导数基·系数，**不含截距项**（`coef[-1]`）。（Task 4）
4. **`degree(1) nknots(0)` 必须退化为 v1**——回归系数、ATT(d) 与 v1 完全一致；任何浮点差异都要查。（Task 3）
5. **统一带的影响函数**——m0（未处理均值）的影响函数贡献不能漏，否则 sup-t 临界值偏小、带过窄、覆盖不足。（Task 5）
6. **节点重复/在边界**——若 `knots(numlist)` 提供的节点落在 treated 剂量范围外或重复，应报错而非静默错算。（Task 2）

---

### Task 1: 研究 + 对拍准备

**Files:**
- Create: `contdid/docs/research-notes-v2a.md`
- Create: `contdid/examples/reference/run_splines2_check.R`

**Interfaces:**
- Consumes: 无
- Produces: `research-notes-v2a.md`（样条基/导数精确约定、节点约定、影响函数公式）；R 对拍脚本

- [ ] **Step 1: 精读 R `splines2` 的 bSpline/dbs 定义**

读 `splines2` 文档（https://cran.r-project.org/package=splines2）：`bSpline(x, degree, knots, intercept=FALSE)` 的基个数公式、边界节点默认值、`dbs` 导数定义。把**基个数 = nknots + degree**、边界节点 = range(x)（或默认）、`intercept=FALSE` 的确切行为写进 notes。

- [ ] **Step 2: 复读 R `cont_did_acrt` 源码**

确认估计量三行（`bs <- bSpline(dose[dose>0], degree, knots)`；`lm(dy ~ .)`；`att.d <- predict - mean(dy[dose==0])`）与 ACRT 行（`acrt.d <- dbs(dvals) %*% coef[-1]`），以及影响函数结构（`Xe=estfun`、`bread=bread`、`inffunc` 组合）。写进 notes。

- [ ] **Step 3: 安装 R `splines2`**

Run: `Rscript -e 'install.packages("splines2", repos="https://cloud.r-project.org")'`
Expected: `requireNamespace("splines2")` → TRUE。

- [ ] **Step 4: 写对拍脚本骨架**

`run_splines2_check.R`：生成一组 x、degree=3、nknots=2 的分位数节点，输出 `bSpline` 基矩阵与 `dbs` 导数矩阵到 csv，供 Task 2 对拍。

- [ ] **Step 5: 写 research-notes-v2a.md**

写入：①B 样条 Cox-de Boor 递归与导数公式；②splines2 的基个数/边界/截距约定；③节点选择（分位数）；④ATT(d)/ACRT(d)/影响函数公式；⑤手算小例子。

- [ ] **Step 6: Commit**

```bash
cd "/Users/niuhaoyu/Documents/open code/contdid"
git add docs/research-notes-v2a.md examples/reference/run_splines2_check.R
git commit -m "docs: v2a 研究笔记 + splines2 对拍脚本"
```

---

### Task 2: Mata B 样条基 + 导数

**Files:**
- Create: `contdid/ado/bspline.mata`（Mata 函数：`bspline_basis(x, degree, knots)`、`bspline_deriv(x, degree, knots)`）
- Test: `contdid/examples/_test_bspline.do`

**Interfaces:**
- Consumes: `research-notes-v2a.md`（基约定）
- Produces: Mata 函数 `bspline_basis(X, degree, knots)` → n×K 基矩阵；`bspline_deriv(X, degree, knots)` → n×K 导数矩阵（K = nknots + degree）

- [ ] **Step 1: 写失败测试**

`_test_bspline.do`：在 Stata 里对一组固定 x（如 0.1..0.9）和 degree=3、nknots=2（分位数节点），调用 `bspline_basis`，把结果矩阵写 csv；然后 R `splines2::bSpline` 同一设定写 csv；两文件逐元素 |差|≤1e-6 才算过。先用占位（Mata 函数未定义）让脚本报错。

- [ ] **Step 2: 运行，预期失败**

Run: stata_do `_test_bspline.do`
Expected: `bspline_basis` not found（Mata 未定义）。

- [ ] **Step 3: 实现 Cox-de Boor 递归**

在 `bspline.mata` 实现：
```
real matrix bspline_basis(real colvector x, real scalar degree, real colvector knots)
  // knots = 内部节点；构造扩展节点向量（两端各补 degree+1 个边界点，取 min/max 重复）
  // B_{i,0}(x) = 1[t_i <= x < t_{i+1}]；递归 B_{i,p} = (x-t_i)/(t_{i+p}-t_i) B_{i,p-1} + (t_{i+p+1}-x)/(t_{i+p+1}-t_{i+1}) B_{i+1,p-1}
  // 返回 n×K，K = length(knots)+degree
real matrix bspline_deriv(...)
  // B'_{i,p} = p/(t_{i+p}-t_i) B_{i,p-1} - p/(t_{i+p+1}-t_{i+1}) B_{i+1,p-1}
```

- [ ] **Step 4: 运行，与 splines2 对拍通过**

Run: stata_do `_test_bspline.do` + R 脚本
Expected: 逐元素 |差| ≤ 1e-6。

- [ ] **Step 5: 节点校验（Review Focus #6）**

在 Mata 入口校验：内部节点须落在 (min(x), max(x)) 内且不重复，否则 `error()`。

- [ ] **Step 6: Commit**

```bash
git add ado/bspline.mata examples/_test_bspline.do
git commit -m "feat: Mata B 样条基 + 导数（对拍 splines2）"
```

---

### Task 3: B 样条估计 ATT(d)（向后兼容）

**Files:**
- Modify: `contdid/ado/contdid.ado`（加 `degree`/`knots`/`nknots` 选项，样条分支）
- Test: `contdid/examples/_test_bspline_att.do`

**Interfaces:**
- Consumes: `bspline_basis`（Task 2）
- Produces: 命令 `contdid` 新选项 `degree(#)` `knots(numlist)` `nknots(#)`；`r(attd)` 仍为 npoints×5

- [ ] **Step 1: 写向后兼容测试**

`_test_bspline_att.do`：`degree(1) nknots(0)` 跑出的 ATT(d) 与 v1（旧 `r(attd)` 黄金值）逐位一致。

- [ ] **Step 2: 写非线性恢复测试**

二次剂量 DGP（真实 ATT(d)=0.3d+0.4d²），`degree(2) nknots(0)`（全局二次）恢复真实曲线，npoints 处最大误差 < 阈值。

- [ ] **Step 3: 实现样条分支**

在 ado 中：解析 `degree`/`knots`/`nknots`；默认 degree=1、nknots=0（走原线性路径）；否则调用 `bspline_basis`，`mata` 里 OLS ΔY~基（treated），ATT(d)=φ(d)'β−m0。

- [ ] **Step 4: 运行测试**

Run: stata_do 两个测试
Expected: 向后兼容逐位一致；非线性恢复误差在阈值内。

- [ ] **Step 5: Commit**

```bash
git add ado/contdid.ado examples/_test_bspline_att.do
git commit -m "feat: B 样条估计 ATT(d)（degree/knots 选项，向后兼容）"
```

---

### Task 4: ACRT(d)

**Files:**
- Modify: `contdid/ado/contdid.ado`（加 `r(acrt)`）
- Test: `contdid/examples/_test_acrt.do`

**Interfaces:**
- Consumes: `bspline_deriv`（Task 2）
- Produces: `r(acrt)`（npoints×5）；ACRT(d)=导数基·系数（**不含截距**）

- [ ] **Step 1: 写失败测试**

`_test_acrt.do`：线性 DGP（真实 ATT(d)=0.5d）下 ACRT(d) 应≈常数 0.5（各评估点 |ACRT−0.5|<阈值）。

- [ ] **Step 2: 运行，预期失败**

Expected: `r(acrt)` 不存在。

- [ ] **Step 3: 实现 ACRT**

ACRT(d)=Σ β_k B'_k(d)，系数为样条回归不含截距的部分（对应 R `coef(bs_reg)[-1]`）。

- [ ] **Step 4: 运行，预期通过**

Expected: 线性 DGP 下 ACRT≈0.5；二次 DGP 下 ACRT(d)≈0.3+0.8d（导数）。

- [ ] **Step 5: Commit**

```bash
git add ado/contdid.ado examples/_test_acrt.do
git commit -m "feat: ACRT(d) 边际因果响应"
```

---

### Task 5: 统一置信带（cband）

**Files:**
- Modify: `contdid/ado/contdid.ado`（加 `cband`）
- Test: `contdid/examples/_test_cband.do`

**Interfaces:**
- Consumes: Task 3/4 的估计量 + 影响函数（research-notes-v2a）
- Produces: `cband` 选项；`r(crit_att)`、`r(crit_acrt)` 及带上下界

- [ ] **Step 1: 写覆盖检验**

`_test_cband.do`：重复 100 次小样本，统计统一带（sup-t）对整条 ATT(d) 曲线的覆盖率 ≥ 名义（如 0.95）。

- [ ] **Step 2: 实现影响函数 + multiplier bootstrap**

按 notes 移植 R 结构：`Xe`（OLS 得分）、`bread`；ATT(d) 与 ACRT(d) 的 IF；multiplier（标准正态）× IF → sup-t 统计量 → 取 (1−α) 分位数为临界值。

- [ ] **Step 3: 运行覆盖检验**

Run: stata_do `_test_cband.do`
Expected: 覆盖 ≥ 0.95（或落入 [0.90, 1.0]）。

- [ ] **Step 4: Commit**

```bash
git add ado/contdid.ado examples/_test_cband.do
git commit -m "feat: 统一置信带 (cband, multiplier bootstrap)"
```

---

### Task 6: 文档 + 版本 + push

**Files:**
- Modify: `contdid/README.md`、`contdid/sthlp/contdid.sthlp`、`contdid/CHANGELOG.md`、`contdid/examples/contdid_example.do`

- [ ] **Step 1: 更新 README（加 degree/knots/ACRT/cband 说明与示例）**
- [ ] **Step 2: 更新 sthlp（新选项表）**
- [ ] **Step 3: CHANGELOG 记 v0.2.0 + ado 顶部版本号**
- [ ] **Step 4: 更新示例（用 degree(2)/degree(3) + cband + graph 画 ATT 与 ACRT）**
- [ ] **Step 5: 全量回归（跑 v1 全部旧测试 + 新测试，确认向后兼容）**
- [ ] **Step 6: Commit + push**

```bash
git add -A && git commit -m "docs: v0.2.0 B样条+ACRT+统一带 文档"
git push origin main
```

---

## Self-Review 结论

- **Spec 覆盖**：§3 样条基 → Task 2；§2 ATT/ACRT → Task 3/4；§4 统一带 → Task 5；§5 语法 → Task 3；§6 返回 → Task 3/4/5；§8 验证 → Task 1-5；§10 顺序 → 6 任务对应。
- **占位符**：无 TBD；精确公式通过 Task 1 research-notes 落定（前置研究交付物）。
- **类型一致性**：`bspline_basis(X, degree, knots)` 返回 n×K（K=nknots+degree）在 Task 2 定义、Task 3/4 复用；`r(attd)`/`r(acrt)` 均 npoints×5，一致。
- **Review Focus**：6 条分别落在 Task 2（#1/#6）、Task 1（#2）、Task 4（#3）、Task 3（#4）、Task 5（#5），均已配测试。
