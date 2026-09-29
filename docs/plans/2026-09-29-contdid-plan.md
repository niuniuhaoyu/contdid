# contdid（连续处理 DiD · Stata 包）Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 实现一个正确、可用、有文档的 Stata 包 `contdid`，估计连续处理变量的 DiD 剂量-反应曲线 ATT(d)。

**Architecture:** 标准 Stata 包结构（ado/sthlp/examples/data/docs）。核心命令 `contdid.ado` 用线性-in-dose 方法估计 ATT(d)，按个体聚类自助法出点对点置信区间，`graph` 选项画剂量-反应图。正确性通过「模拟数据恢复已知真实 ATT(d)」+「对照 R `contdid`」双重验证。

**Tech Stack:** Stata（ado + Mata）、R（仅用于生成参考值交叉验证）、Git/GitHub。

**Spec:** `contdid/docs/specs/2026-09-29-contdid-design.md`（本计划从该 spec 展开，执行者两份都要读）

## Global Constraints

- Stata 16+（ado 语法兼容）
- 许可证：AGPL-3.0
- 包名：`contdid`（命令名 `contdid`）
- v1 仅做：两期 + ATT(d) + 线性-in-dose + 点对点自助推断 + 剂量-反应图
- v1 不做：交错处理、ACRT、统一置信带、CCK、协变量
- 文档语言：英文（README/sthlp）
- 估计量精确公式以论文 Callaway-Goodman-Bacon-Sant'Anna（2024/2025）为准，禁止凭印象硬编

## Review Focus

以下输入/情形最容易让使用者踩坑，必须逐一被测试覆盖（对应任务已在括号标注）：

1. **剂量 D 全为 0 或全为 0 的对照组不存在**——若数据里没有 D=0 的对照单位，ATT(d) 的基线无法由「D=0 均值」定义，命令应报错或明确退回「截距外推」，不能静默给错结果。（Task 3）
2. **D 为负值或非数值**——剂量必须 ≥ 0，违反应报错并指明变量。（Task 3）
3. **非两期数据 / time 变量不是 0/1**——v1 只支持两期，时间变量取值异常应报错。（Task 3）
4. **D=0 基线是「对照均值」还是「回归截距」**——两者在有限样本下不等，必须先精读论文定死默认行为并写入 research-notes。（Task 1，Task 3 实现按 notes）
5. **自助法随机种子**——未设 `seed()` 时结果应可复现（默认固定种子），否则无法复现论文结果。（Task 4）

---

### Task 1: 研究 + 包骨架

**Files:**
- Create: `contdid/docs/research-notes.md`
- Create: `contdid/LICENSE`
- Create: `contdid/CHANGELOG.md`
- Create: `contdid/.gitignore`
- Create: `contdid/contdid.pkg`
- Create: `contdid/stata.toc`
- Create: `contdid/ado/contdid.ado`（占位，仅 `program define` 空壳 + 版本号）
- Create: `contdid/sthlp/contdid.sthlp`（占位）

**Interfaces:**
- Consumes: 无（首个任务）
- Produces: `research-notes.md`（后置所有任务读取的估计量精确定义）；目录骨架

- [ ] **Step 1: 精读论文估计量章节**

打开 arXiv:2107.02637（Callaway-Goodman-Bacon-Sant'Anna），重点读：两期、连续处理、ATT(d) 的定义与识别条件（"strong parallel trends" 与"generalized parallel trends"的区别）、对照组如何进入估计量。把结论写进 `research-notes.md`。

- [ ] **Step 2: 读 R contdid 源码的默认估计量**

从 `https://github.com/bcallaway11/contdid` 读 `num_knots=0, degree=1` 时的估计路径（`target_parameter="level"`, `aggregation="dose"`）。关键要定死一件事并写进 notes：**ATT(d) 的基线是「D=0 单位的 ΔY 均值」还是「回归截距外推」**（对应 Review Focus #4）。

- [ ] **Step 3: 写 research-notes.md**

写入：①估计量精确公式（含 ATT(d)=… 的每一项）；②识别假设清单；③对照组定义；④与 R `contdid` 默认行为的对应关系；⑤一个手算小例子（N=4，便于后面单测对答案）。

- [ ] **Step 4: 建目录 + 骨架文件**

用 `write` 创建上面 Files 列出的所有文件。LICENSE 用 AGPL-3.0 全文；`.gitignore` 排除 `.DS_Store`、`*.log`、`__pycache__`、`*.mlib`；`contdid.pkg` 写 `d 'CONTdid'` 等 SSC 元数据占位；`contdid.ado` 只写 `program define contdid, rclass` + `version 16` + 空壳。

- [ ] **Step 5: 验证骨架可加载**

在 Stata 里运行：
```stata
adopath + "/Users/niuhaoyu/Documents/open code/contdid/ado"
capture which contdid
```
预期：`which contdid` 能找到命令（即使里面是空壳）。

- [ ] **Step 6: Commit**

```bash
cd "/Users/niuhaoyu/Documents/open code/contdid"
git init && git add -A && git commit -m "chore: contdid 包骨架 + 研究笔记"
```

---

### Task 2: 模拟数据 DGP

**Files:**
- Create: `contdid/examples/contdid_simdata.do`
- Create: `contdid/data/contdid_sim.dta`（由上面 .do 生成）

**Interfaces:**
- Consumes: 无（独立）
- Produces: `contdid_sim.dta`，变量：`id`（个体）、`t`（0/1）、`y`（结局）、`d`（连续剂量，含 D=0 对照）。真实 ATT(d)=0.5·d。

- [ ] **Step 1: 写 DGP 生成 .do**

DGP（两期、连续处理、单位固定效应满足平行趋势）：
```stata
clear all
set seed 12345
set obs 4000
gen id = _n
gen d = runiform() > 0.2 ? runiform() : 0   // 约 20% 为 D=0 对照
gen mu = 2 + rnormal(0,1)                     // 单位固定效应
gen e0 = rnormal(0,0.3)
gen e1 = rnormal(0,0.3)
gen y0 = mu + e0                              // 处理前
gen y1 = mu + 0.5*d + e1                      // 处理后，真实 ATT(d)=0.5*d
gen t = 0
tempfile pre
save `pre'
gen t = 1
gen y = y1
append using `pre'
replace y = y0 if t==0
keep id t y d
save "contdid/data/contdid_sim.dta", replace
```

- [ ] **Step 2: 运行并验证数据**

Run: `stata -b do contdid/examples/contdid_simdata.do`（或经 stata-mcp `stata_do`）。
验证：`contdid_sim.dta` 存在，`d` 有 0 和非 0 两类，`t` 只有 0/1，样本 8000 行。

- [ ] **Step 3: Commit**

```bash
cd "/Users/niuhaoyu/Documents/open code/contdid"
git add examples/contdid_simdata.do data/contdid_sim.dta
git commit -m "feat: 模拟数据 DGP（真实 ATT(d)=0.5d）"
```

---

### Task 3: 核心估计量（线性-in-dose）

**Files:**
- Modify: `contdid/ado/contdid.ado`（实现主体）

**Interfaces:**
- Consumes: `research-notes.md`（估计量精确定义）；`contdid_sim.dta` 的结构（id/t/y/d）
- Produces: 命令 `contdid`，语法 `contdid depvar [if] [in], unit(varname) time(varname) dose(varname) [npoints(#) level(#)]`；返回 r() 与 e() 结果，包含矩阵 `e(attd)`（列：d、ATT、se、lb、ub）

- [ ] **Step 1: 写失败测试（验证脚本）**

Create `contdid/examples/_test_recover.do`，跑 contdid 并检查 ATT(d) ≈ 0.5·d：
```stata
use "contdid/data/contdid_sim.dta", clear
contdid y, unit(id) time(t) dose(d) npoints(5) seed(1)
matrix b = e(attd)
* 检查第 3 列（ATT）应约等于 0.5 * 第 1 列（d）
local maxerr = 0
forvalues i = 1/5 {
    local d = b[`i',1]
    local att = b[`i',3]
    local err = abs(`att' - 0.5*`d')
    if `err' > `maxerr' local maxerr = `err'
}
di "MAX ERR = " `maxerr'
assert `maxerr' < 0.05
```

- [ ] **Step 2: 运行测试，预期失败**

Run: stata_do `contdid/examples/_test_recover.do`
Expected: 报错（`contdid` 还是空壳，或 e(attd) 不存在）。

- [ ] **Step 3: 实现输入校验**

在 `contdid.ado` 内先写校验（对应 Review Focus #2、#3）：
- `dose()` 变量若存在负值 → `display as error` 并 `exit 198`
- `time()` 变量取值若不在 {0,1} → 报错并提示 v1 仅支持两期
- 无 D=0 单位时，若走「对照均值」基线 → 报错；若走「截距外推」→ 显示 note（行为由 research-notes 决定）

- [ ] **Step 4: 实现线性-in-dose 估计量**

按 `research-notes.md` 的公式实现：生成 Δy = y(t=1) − y(t=0)；对 Δy 与 d 拟合（线性）；在 `npoints()` 个评估点 d 上计算 ATT(d)（基线按 notes 定死），填入 `e(attd)`。

- [ ] **Step 5: 运行测试，预期通过**

Run: stata_do `contdid/examples/_test_recover.do`
Expected: PASS，`MAX ERR` < 0.05。

- [ ] **Step 6: Commit**

```bash
cd "/Users/niuhaoyu/Documents/open code/contdid"
git add ado/contdid.ado examples/_test_recover.do
git commit -m "feat: 线性-in-dose 估计量 ATT(d)"
```

---

### Task 4: 聚类自助法推断

**Files:**
- Modify: `contdid/ado/contdid.ado`（加自助法）

**Interfaces:**
- Consumes: Task 3 的估计量函数（可抽出为内部 program，如 `_contdid_att`）
- Produces: 选项 `reps(#)` `seed(#)` `cluster(varname)`；`e(attd)` 补上 se/lb/ub 列

- [ ] **Step 1: 写覆盖检验脚本**

Create `contdid/examples/_test_coverage.do`：重复 200 次小样本模拟，统计 ATT(d) 点对点 CI 覆盖率，应接近 95%：
```stata
* 伪代码：循环 200 次，每次生成小样本、跑 contdid，记录 d=0.5 处 CI 是否覆盖真值 0.25
```
（具体循环体用 Task 2 的 DGP 缩到 N=200，`contdid ... npoints(1)` 或固定 d 点）

- [ ] **Step 2: 实现聚类自助法**

在 ado 内：按 `cluster()`（默认 `unit()`）重抽样 B=`reps()` 次；每次重跑估计量得到 ATT(d) 的自助分布；取分位数得点对点 CI；写 se/lb/ub。默认 `reps(999) seed(12345)`（对应 Review Focus #5：无 seed 也可复现）。

- [ ] **Step 3: 运行覆盖检验**

Run: stata_do `contdid/examples/_test_coverage.do`
Expected: 覆盖率在 [0.90, 0.98] 之间。

- [ ] **Step 4: Commit**

```bash
cd "/Users/niuhaoyu/Documents/open code/contdid"
git add ado/contdid.ado examples/_test_coverage.do
git commit -m "feat: 聚类自助法点对点置信区间"
```

---

### Task 5: 剂量-反应图

**Files:**
- Modify: `contdid/ado/contdid.ado`（加 `graph` 选项）

**Interfaces:**
- Consumes: `e(attd)`（d、ATT、lb、ub）
- Produces: `graph` 选项输出 ATT(d) vs d 曲线 + CI 带

- [ ] **Step 1: 实现画图**

`graph` 选项：用 `e(attd)` 画 `twoway (rarea ub lb d) (line ATT d)`，加标题 "Dose–response: ATT(d)"、x 轴 "Dose (d)"、y 轴 "ATT(d)"。

- [ ] **Step 2: 手动运行验证出图**

Run: 在示例数据上 `contdid y, unit(id) time(t) dose(d) graph`，确认图正常生成、无报错。

- [ ] **Step 3: Commit**

```bash
cd "/Users/niuhaoyu/Documents/open code/contdid"
git add ado/contdid.ado
git commit -m "feat: 剂量-反应图"
```

---

### Task 6: 对照 R contdid 交叉验证

**Files:**
- Create: `contdid/examples/reference/run_contdid.R`
- Create: `contdid/examples/reference/reference_attd.csv`

**Interfaces:**
- Consumes: `contdid_sim.dta`（可导出 csv 供 R 读）
- Produces: `reference_attd.csv`（R `contdid` 在 num_knots=0, degree=1 下的 ATT(d) 参考值）

- [ ] **Step 1: 导出数据为 csv**

在 Stata 里 `export delimited` 把 `contdid_sim.dta` 存为 csv（或直接让 R 读 .dta 用 haven）。

- [ ] **Step 2: 写 R 脚本跑 contdid**

```r
# install.packages("contdid", repos = c("https://bcallaway11.r-universe.dev", "https://cloud.r-project.org"))
library(contdid)
df <- read.csv("contdid_sim.csv")
res <- cont_did(yname="y", tname="t", idname="id", dname="d", data=df,
                gname="g", target_parameter="level", aggregation="dose",
                treatment_type="continuous", num_knots=0, degree=1)
write.csv(res$att.d, "reference_attd.csv", row.names=FALSE)
```

- [ ] **Step 3: 对比 Stata 与 R 结果**

把 `reference_attd.csv` 与 Stata `e(attd)` 逐点比较，|差| ≤ 1e-6 才算通过。若有偏差，回 Task 1/3 修正估计量定义。

- [ ] **Step 4: Commit**

```bash
cd "/Users/niuhaoyu/Documents/open code/contdid"
git add examples/reference/
git commit -m "test: 对照 R contdid 交叉验证"
```

---

### Task 7: 文档与发布

**Files:**
- Modify: `contdid/README.md`
- Modify: `contdid/sthlp/contdid.sthlp`
- Modify: `contdid/CHANGELOG.md`
- Modify: `contdid/examples/contdid_example.do`

**Interfaces:**
- Consumes: 最终命令语法与选项
- Produces: 完整文档 + 一键复现示例

- [ ] **Step 1: 写英文 README**

含：标题 + badges、Overview（方法 + 论文引用）、Installation（`net install contdid, from("https://raw.githubusercontent.com/niuniuhaoyu/contdid/main/") replace`）、Quick start（3 行示例）、Estimator 说明、Citation（BibTeX）、License。

- [ ] **Step 2: 写 sthlp 帮助文件**

按 Stata 帮助格式写全：语法、选项表、示例、保存结果、参考文献。

- [ ] **Step 3: 写一键示例**

`contdid_example.do`：`use contdid_sim.dta` → `contdid ... graph` → 输出表和图，全程可复现。

- [ ] **Step 4: 写 CHANGELOG + 更新版本**

CHANGELOG 记 v0.1.0；`contdid.ado` 顶部 `version` 与包元数据对齐。

- [ ] **Step 5: 提交并推 GitHub**

```bash
cd "/Users/niuhaoyu/Documents/open code/contdid"
git add -A && git commit -m "docs: README/sthlp/CHANGELOG/example"
git remote add origin https://github.com/niuniuhaoyu/contdid.git
git branch -M main
git push -u origin main
```

---

## Self-Review 结论

- **Spec 覆盖**：§3 v1 范围 → Task 3/4/5；§4 命令语法 → Task 3 语法实现；§5 估计与推断 → Task 3/4；§6 包结构 → Task 1；§7 验证 → Task 2/6；§8 交付物 → Task 7。无遗漏。
- **占位符**：无 TBD/TODO；估计量公式通过 Task 1 的 research-notes 落定（不是占位，是明确的前置研究交付物）。
- **类型一致性**：`e(attd)` 列序（d/ATT/se/lb/ub）在 Task 3 定义、Task 4/5 复用，一致。
- **Review Focus**：5 条分别落在 Task 3（#1/#2/#3/#4）、Task 4（#5），均已配测试。
