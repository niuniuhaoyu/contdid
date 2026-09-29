# contdid v2a 研究笔记（B 样条基 / 导数 / 影响函数）

> 来源：R `splines2` 0.5.4（`bSpline`/`dbs`）实际输出 + `contdid` R 源码 `cont_did_acrt`
> 用途：Task 2/3/4/5 严格照本文实现，禁止凭印象。

---

## 1. B 样条基约定（已用 R 实测核实）

- 输入：评估点 x、次数 degree、内部节点 knots（nk 个）。
- **边界节点** = range(x)（min、max），在两端各重复 **degree+1** 次。
- 扩展节点向量：`t = [min× (degree+1), knots, max× (degree+1)]`，长度 = 2(degree+1) + nk。
- **完整基函数个数** = length(t) − (degree+1) = **degree + 1 + nk**。
- **`intercept=FALSE`**（R `bSpline` 默认，也是 R contdid 用的）：**丢弃第 1 个基函数**（左端那个「近似常数」的），
  故返回 **K = degree + nk** 列。
- 实测：degree=3, knots=c(0.4,0.7) → K = 5（= 2+3）✓。

### Cox-de Boor 递归

```
B_{i,0}(x) = 1 若 t_i <= x < t_{i+1}，否则 0
B_{i,p}(x) = (x − t_i)/(t_{i+p} − t_i) · B_{i,p−1}(x)
           + (t_{i+p+1} − x)/(t_{i+p+1} − t_{i+1}) · B_{i+1,p−1}(x)
```

### 导数递归

```
B'_{i,p}(x) = p/(t_{i+p} − t_i) · B_{i,p−1}(x) − p/(t_{i+p+1} − t_{i+1}) · B_{i+1,p−1}(x)
```

> 约定：分母为 0 时该项取 0（节点重复处的极限）。

## 2. 估计量（v2a，与 v1 同构）

treated（D>0）上 OLS：

```
ΔY = α + Σ_{k=1..K} β_k B_k(D)      （B_k 为「intercept=FALSE」后的 K 列基，OLS 另加常数项 α）
```

- **ATT(d) = α + Σ β_k B_k(d) − m0**，m0 = mean(ΔY | D=0)
- **ACRT(d) = Σ β_k B'_k(d)**（导数基，**不含截距**，对应 R `coef(bs_reg)[-1]`）

> 退化：degree=1、nk=0 → K=1，B_1(d)=d（线性），ATT(d)=α+βd−m0 = v1。

## 3. 节点选择

- `nknots(#)`：在 treated 剂量分布的分位数处取 nk 个内部节点
  （R `choose_knots_quantile`：`quantile(x, probs=seq(0,1,length=nk+2))[-c(1,nk+2)]`）。
- `knots(numlist)`：用户显式指定，优先于 nknots。

## 4. 影响函数与统一置信带（cband）

以 R `cont_did_acrt` 为蓝本：

- 记设计阵 X = [1, B_1(D), …, B_K(D)]（treated，含截距列），`Xe` = OLS 得分（n×(K+1)），
  `bread` = (X'X/n)^{-1}。
- ATT(d) 的 IF：`φ(d)' · (Xe · bread)' − 1[D=0 的 IF 项]`，其中 φ(d) = [1, B_1(d), …, B_K(d)]，
  未处理项 IF = 1[D=0]·(ΔY−m0)/p0（p0 = P(D=0)）。
- ACRT(d) 的 IF：`φ'(d)' · (Xe · bread)'`，φ'(d) = [0, B'_1(d), …, B'_K(d)]。
- **multiplier bootstrap**：抽 B 次独立标准正态乘子 ξ_i，对每点 d 计算
  `max_d |Σ_i ξ_i·IF_i(d)| / se(d)`，取 (1−α) 分位数为 sup-t 临界值；
  统一带 = 点估计 ± 临界值 × se(d)。

> 精确的 m0 未处理项 IF、以及 se(d) 的逐点估计，实现时以 Task 5 测试（覆盖 Monte Carlo）为准绳校准。

## 5. 手算小例子（供单测）

degree=3, knots=c(0.4,0.7), x 网格 0.1..0.9（9 点）：基/导数黄金值见
`examples/reference/splines2_basis.csv`、`splines2_deriv.csv`（由 `dump_splines2.R` 生成）。

- K = 5（列名 1..5）
- x=0.1 行：基 = [0,0,0,0,0]，导数 = [10, 0,0,0,0]
- x=0.2 行：基 = [0.564815, 0.131944, 0.006944, 0, 0]
