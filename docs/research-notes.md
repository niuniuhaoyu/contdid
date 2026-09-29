# contdid 研究笔记（估计量精确定义）

> 来源：Callaway, Goodman-Bacon & Sant'Anna (2024/2025) arXiv:2107.02637（v8, 2025-12-31）
> 及官方 R 包 `bcallaway11/contdid` 源码（`cont_did.R`、`cont_did_acrt`）
> 用途：Task 3 实现严格照本文实现，禁止凭印象。

---

## 1. 设定（v1：两期，无交错）

- 两期面板：t=1（处理前，无人被处理）、t=2（处理后，部分单位接受剂量 D）。
- 剂量 D：支持 `{0} ∪ D_+`，其中 D=0 表示未处理，D_+ ⊆ (0,∞)。要求 P(D=0) > 0（存在未处理对照）。
- 观测：`{Y_{i,2}, Y_{i,1}, D_i}` iid。
- ΔY = Y_{t=2} − Y_{t=1}。

## 2. 关键识别结果（Theorem 3.1）

在 **Assumption PT（平行趋势）** 下：

> **ATT(d|d) = E[ΔY | D=d] − E[ΔY | D=0]**

其中 ATT(d|d) 是「剂量组 d 的平均水平处理效应」（local ATT，效果相对于 0 剂量）。

**注意区分**：
- **ATT(d|d)**（local，剂量组 d 自身的效应）——标准 PT 下**可识别**。
- **ATT(d)**（global，对所有 treated 单位）——标准 PT 下**不可识别**，需要「强平行趋势」(strong PT)。
- **ACRT(d)**（因果响应/slope）——v1 不做。

## 3. 精确估计量（对应 R `cont_did_acrt`，degree=1, knots=0 = 线性）

源码核心三行（已核实）：
```r
bs <- splines2::bSpline(dose[dose > 0], degree = degree, knots = knots)  # 仅 treated
bs$dy <- dy[dose > 0]
bs_reg <- lm(dy ~ ., data = bs)                    # 仅 treated 上回归 ΔY ~ 剂量
att.d <- predict(bs_reg, newdata = bs_grid) - mean(dy[dose == 0])  # 减 untreated 均值
```

翻译为线性（degree=1, knots=0）时，B 样条基退化为线性项，即：

1. 取 **treated 单位（D > 0）** 的子样本；
2. 对 ΔY 做 **OLS：ΔY = α + β·D**（含截距），得 `(α̂, β̂)`；
3. 计算 **m0 = mean(ΔY | D = 0)**（未处理单位 ΔY 的样本均值）；
4. 对每个评估剂量 d：
   **ATT(d) = (α̂ + β̂·d) − m0**

**关键裁决（Review Focus #4）**：基线是 **mean(ΔY | D=0)**，**不是**回归截距 α̂。
- 这对应「D=0 对照均值」方案，已由源码 `- mean(dy[dose == 0])` 定死。
- 若数据无 D=0 单位（违反 Assumption 2），命令必须报错（见 Task 3 校验）。

## 4. 识别假设清单

- Assumption 1：iid 随机抽样。
- Assumption 2：t=1 无人被处理；t=2 剂量 D ∈ {0}∪D_+；P(D=0)>0。
- Assumption 3：无预期效应，Y_{i,1}=Y_{i,1}(0)，Y_{i,2}=Y_{i,2}(D_i)。
- Assumption PT：对所有 d>0，E[Y_2(0)−Y_1(0) | D=d] = E[Y_2(0)−Y_1(0) | D=0]。

## 5. 手算小例子（供 Task 3 单测对答案）

假设 4 个单位（两期），ΔY = Y2 − Y1：

| 单位 | D | ΔY |
|---|---|---|
| 1 | 0  | 2.0 |
| 2 | 0  | 1.0 |
| 3 | 1  | 4.0 |
| 4 | 2  | 6.0 |

- treated（D>0）子样本：单位 3、4 → ΔY = α + β·D
  - β̂ = cov(D,ΔY)/var(D) = [ (1−1.5)(4−5) + (2−1.5)(6−5) ] / [ (1−1.5)²+(2−1.5)² ]
    = [ (−0.5)(−1) + (0.5)(1) ] / [ 0.25+0.25 ] = [ 0.5 + 0.5 ] / 0.5 = 1.0 / 0.5 = **2.0**
  - α̂ = mean(ΔY) − β̂·mean(D) = 5 − 2.0·1.5 = 5 − 3 = **2.0**
- m0 = mean(ΔY | D=0) = (2.0+1.0)/2 = **1.5**
- 因此：ATT(d) = (2.0 + 2.0·d) − 1.5 = **0.5 + 2.0·d**
- 验证点：ATT(0)=0.5、ATT(1)=2.5、ATT(2)=4.5

> 注：这个例子里 treated 只有 2 个点，线性拟合完全确定（R²=1）。Task 3 单测可构造 N 更大的 DGP（真实 ATT(d)=0.5d，见 Task 2），断言 ATT(d)≈0.5d。

## 6. v1 设计选择（与 R 的差异，需记录）

| 项 | R `contdid` | 本包 v1 |
|---|---|---|
| 估计量（点估计） | 同上（linear, degree=1, knots=0） | **完全一致** |
| 评估点 dvals | 默认 treated 剂量分位数 | `npoints()` 个在 [min D_+, max D_+] 等距 |
| 推断 | multiplier bootstrap + influence function + uniform band | **聚类自助法（按 unit）点对点 CI**（简化） |
| 交错处理 | 支持 | 不支持（v2） |
| ACRT | 支持 | 不支持（v2） |

> 点估计与 R 逐点一致是验证硬指标；推断方式 v1 简化为聚类自助，属已记录的偏离。

## 7. 实现伪代码（contdid.ado 核心）

```
1. 校验：dose 无负值、time ∈ {0,1}、存在 D=0 单位
2. 生成 dy = y(t=1) - y(t=0)   （注意 R 的 get_first_difference 方向：post - pre）
3. treated 子样本（D>0）：regress dy on D → b[D], b[_cons]
4. m0 = mean(dy | D==0)
5. for d in npoints 评估点: att(d) = b[_cons] + b[D]*d - m0
6. 聚类自助：按 unit 重抽样 reps 次，重算 3-5，取分位数得 SE/CI
```

> ΔY 方向：R 源码 `get_first_difference` 取的是 post − pre（t=2 减 t=1）。本包 `dy = y(post) - y(pre)`，与之一致。
