# contdid — Stata 包设计（Spec）

> 日期：2026-09-29
> 状态：待用户审阅（Draft）
> 作者：Haoyu Niu
> 定位：连续处理变量（Continuous Treatment）的 Difference-in-Differences，Stata 实现

---

## 1. 背景与目标

**方法**：Callaway, Goodman-Bacon & Sant'Anna（2024/2025）*"Difference-in-Differences with a Continuous Treatment"*（NBER WP 32117 / arXiv:2107.02637）。

**缺口**：官方仅提供 R 包 `contdid`（alpha 版），**无 Stata 实现**。大量使用 Stata 的实证研究者面对连续处理（政策剂量、金融发展指数、补贴金额等）时没有现成工具。

**目标**：实现一个**正确、可用、有文档、可复现**的 Stata 包 `contdid`，填补该缺口，作为作者 GitHub 作品集的核心项目。

**个人价值对齐**：作者的研究方向（如「金融发展 → 收入不平等」）中，金融发展本身是连续变量，本包可直接服务于其论文。

---

## 2. 方法概述（高层）

- 连续处理变量 D ≥ 0（D = 0 表示未处理）。
- 核心参数：**ATT(d)** —— 剂量为 d 的处理组，其平均处理效应（treatment effect on the treated, 作为剂量的函数）。
- 识别依赖**广义平行趋势**（generalized parallel trends）。
- 剂量-反应关系 ATT(d) 描述「处理强度每增加一单位，效应如何变化」。

> 精确的估计量公式、识别条件、假设以论文原文为准；实现阶段先精读论文，再逐项对照 R `contdid` 验证数值正确性。

---

## 3. v1 范围

### 做（v1）

- 两期（或共同处理时点）面板 / 重复截面，连续处理变量 D ≥ 0
- 估计剂量-反应曲线 **ATT(d)**
- 估计方法：**线性-in-dose**（对应 R 包 `num_knots=0, degree=1` 的默认），代码接口预留 B 样条扩展
- 推断：**点对点置信区间**（按个体聚类自助法 + 稳健标准误）
- 输出：ATT(d) 数值表 + 剂量-反应图
- 配套：模拟示例数据、英文 README、命令帮助文件（sthlp）、LICENSE、CHANGELOG

### 不做（v2+，明确列为未来）

- 交错处理（staggered adoption）
- ACRT（平均因果响应，ATT 的导数 / slope 参数）
- 统一置信带（uniform confidence bands）
- CCK 非参数筛（data-driven sieve）估计
- 协变量条件平行趋势

---

## 4. 命令语法（设计）

```
contdid depvar [if] [in] , unit(varname) time(varname) dose(varname) [options]
```

| 项 | 说明 |
|---|---|
| `depvar` | 结局变量 |
| `unit(varname)` | 面板个体 id（必填） |
| `time(varname)` | 时间变量，两期（0=处理前 / 1=处理后） |
| `dose(varname)` | 连续处理剂量 D（≥0，0=未处理） |

### v1 选项

| 选项 | 默认 | 说明 |
|---|---|---|
| `npoints(#)` | 20 | 评估点数量（ATT(d) 曲线上的 d 取值） |
| `reps(#)` | 999 | 自助法重复次数 |
| `seed(#)` | 无 | 随机种子 |
| `cluster(varname)` | `unit()` | 聚类变量 |
| `graph` | 否 | 输出剂量-反应图 |
| `level(#)` | 95 | 置信水平 |
| `method(linear)` | linear | v1 仅 linear，预留 spline |

---

## 5. 估计量与推断（v1）

1. 构造每个个体两期的结局变化 ΔY = Y_post − Y_pre。
2. 用 ΔY 对剂量 D 做（线性）拟合，得到剂量-反应函数 f(d)。
3. 对每个评估剂量 d：**ATT(d) = f(d) − f(0)**（处理组剂量=d 相对对照基线的效应）。
4. 点对点标准误 / 置信区间：按 `unit()` 聚类自助法（默认 999 次）；同时提供稳健（HC）标准误作为解析选项。
5. 输出：d、ATT(d)、SE、pointwise CI 下界/上界；`graph` 选项画出 ATT(d) vs d 曲线与置信区间。

> 注意：上述 1–3 是**高层直觉**；精确估计量（含对照组如何进入、是否含权重）以论文 §3 为准，实现时逐项核对 R `contdid` 的默认行为（`num_knots=0, degree=1`）。

---

## 6. 包结构

```
contdid/
├── README.md              # 英文，含安装/快速上手/方法/引用
├── LICENSE                # AGPL-3.0（与作者既有包一致）
├── CHANGELOG.md
├── contdid.pkg            # SSC 包元数据
├── stata.toc
├── ado/
│   └── contdid.ado        # 主命令
├── sthlp/
│   └── contdid.sthlp      # 帮助文件
├── examples/
│   ├── contdid_example.do # 可运行的复现示例
│   └── contdid_simdata.do # 模拟数据生成
├── data/
│   └── contdid_sim.dta    # 模拟示例数据
└── docs/
    └── specs/
        └── 2026-09-29-contdid-design.md
```

---

## 7. 验证方案（数值正确性）

1. 用固定种子的 DGP 生成模拟面板数据（含已知真实 ATT(d) 曲线）。
2. 在 R 中用 `contdid`（`num_knots=0, degree=1`）对同一数据估计，保存参考值。
3. 在 Stata 中用 `contdid` 对同一数据估计。
4. 逐评估点比较 ATT(d)，容差 ≤ 1e-6。
5. 附上 Monte Carlo 覆盖检验：重复模拟，检查点对点 CI 的覆盖率接近名义水平。

---

## 8. 交付物与验收标准

- [ ] `contdid.ado` + `contdid.sthlp` 可安装、可运行
- [ ] 模拟示例数据可生成，`examples/contdid_example.do` 一键复现
- [ ] ATT(d) 与 R `contdid` 对齐（容差 ≤ 1e-6）
- [ ] 英文 README 完整（安装 / 快速上手 / 方法 / 引用 / 许可证）
- [ ] LICENSE（AGPL-3.0）+ CHANGELOG
- [ ] 上传 GitHub 仓库 `niuniuhaoyu/contdid`

---

## 9. 边界与未来（YAGNI）

- ❌ v1 不做交错处理、ACRT、统一置信带、CCK、协变量（见 §3）
- ✅ v2 候选：B 样条、交错处理、ACRT、统一置信带
- ✅ 与作者既有「现代 DiD」定位一致，可作为作品集第二个核心包

---

## 10. 实施顺序

1. 精读论文 + 读 R `contdid` 源码（搞清估计量精确公式）
2. 搭包骨架（目录 + 模拟数据 + 文档占位）
3. 实现 `contdid.ado`（线性-in-dose + 聚类自助法）
4. 对照 R `contdid` 验证数值
5. 写 README / sthlp / CHANGELOG
6. 上传 GitHub
