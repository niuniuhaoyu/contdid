{smcl}
{* *! version 1.1.0 07oct2026}{...}
{hline}
{p 4 8 2}{bf:contdid} —— 连续处理双重差分{right:版本 1.1.0}
{hline}

{p 4 4 2}{it:英文帮助：} {help contdid}

{title:标题}

{p 4 4 2}
{cmd:contdid} —— 估计连续处理双重差分设计下的剂量-反应函数 ATT(d) 与平均因果响应 ACRT(d)。

{title:语法}

{p 8 12 2}
{cmd:contdid} {it:depvar} {ifin}, {cmdab:unit:(}{it:varname}{cmd:)}
{cmdab:time:(}{it:varname}{cmd:)} {cmdab:dose:(}{it:varname}{cmd:)}
{cmd:[}{cmdab:covariates:(}{it:varlist}{cmd:)} {cmdab:gvar:(}{it:varname}{cmd:)}
{cmdab:degree:(}{it:#}{cmd:)}
{cmdab:nknots:(}{it:#}{cmd:)}
{cmdab:dose_est_method:(}{it:parametric|c dds}{cmd:)} {cmdab:maxknots:(}{it:#}{cmd:)}
{cmdab:knots:(}{it:numlist}{cmd:)} {cmdab:npoints:(}{it:#}{cmd:)}
{cmdab:level:(}{it:#}{cmd:)} {cmdab:reps:(}{it:#}{cmd:)}
{cmdab:seed:(}{it:#}{cmd:)} {cmdab:cluster:(}{it:varname}{cmd:)}
{cmd:cband} {cmd:graph}{cmd:]}

{title:描述}

{p 4 4 2}
{cmd:contdid} 实现 Callaway、Goodman-Bacon 与 Sant'Anna (2024, NBER WP 32117) 的连续处理
双重差分估计量。在平行趋势下，特定剂量的处理组平均处理效应识别为

{p 12 12 2}
ATT(d) = E[ΔY | D = d] − E[ΔY | D = 0]

{p 4 4 2}
其中 ΔY 为 {it:depvar} 的处理前到处理后变化，D 为连续剂量（{cmd:dose()}），D = 0 标记
未处理单位。{cmd:contdid} 用 B 样条对剂量-反应建模（在处理组内对 ΔY 关于 D 的 B 样条基
做 OLS，再减去未处理组 ΔY 均值），估计 ATT(d) 及其导数 ACRT(d)，并报告聚类自助法的逐点
置信区间；{cmd:cband} 通过 multiplier bootstrap 追加统一置信带。

{p 4 4 2}
配合 {cmd:covariates()}，识别改用条件（强）平行趋势（CGBS 2024/2025，补充附录 SI.3）：
ATT(d) = E_X[ ATT_x(d) | D>0 ]，其中 ATT_x(d) = E[ΔY | X=x, D=d] − E[ΔY | X=x, D=0]。
这使估计量超出 R 参考实现——后者不支持协变量。

{title:选项}

{p 4 4 2}{cmdab:unit:(}{it:varname}{cmd:)} 面板单位标识。
{p 4 4 2}{cmdab:time:(}{it:varname}{cmd:)} 时间变量，必须恰好取两个不同值。
{p 4 4 2}{cmdab:dose:(}{it:varname}{cmd:)} 连续处理剂量。必须非负，0 表示未处理；需要存在
{cmd:dose} = 0 的组。
{p 4 4 2}{cmdab:covariates:(}{it:varlist}{cmd:)} 以协变量为条件，把平行趋势放松为{bf:条件
（强）平行趋势}（CGBS 2024/2025，补充附录 SI.3，命题 S3）。它估计含剂量×协变量交互的联合
筛，其中 ATT(d) 为 E[ΔY | X, D=d] − E[ΔY | X, D=0] 对处理组协变量分布的期望。协变量须为
单位层面（时间不变）。
{p 4 4 2}{cmdab:gvar:(}{it:varname}{cmd:)} 处理时点（0 = 从未受处理）；指定后启用多处理组
的交错采纳。
{p 4 4 2}{cmdab:degree:(}{it:#}{cmd:)} B 样条次数（默认 1 = 线性；2 = 二次；3 = 三次）。
{p 4 4 2}{cmdab:nknots:(}{it:#}{cmd:)} 内部节点数，置于剂量分位点（默认 0 = 全局多项式）。
{p 4 4 2}{cmdab:dose_est_method:(}{it:parametric|dds}{cmd:)} 剂量-反应估计方法：
{cmd:parametric}（默认；固定 B 样条次数/节点）或 {cmd:dds}（数据驱动筛：内部节点由
留一交叉验证选择，上限 {cmd:maxknots()}）。{cmd:dds} 是可落地的数据驱动筛，不是精确的
Chen-Christensen-Kankanala/npiv 估计量，也不以对齐 R {cmd:contdid} 的
{cmd:dose_est_method="cck"} 为目标。
{p 4 4 2}{cmdab:maxknots:(}{it:#}{cmd:)} {cmd:dose_est_method(dds)} 尝试的最大内部节点数
（默认 5）。
{p 4 4 2}{cmdab:knots:(}{it:numlist}{cmd:)} 显式指定内部节点（覆盖 {cmd:nknots()}）。
{p 4 4 2}{cmdab:npoints:(}{it:#}{cmd:)} 评估剂量的个数（默认 20，最小 2）。
{p 4 4 2}{cmdab:level:(}{it:#}{cmd:)} 置信水平（%），默认 95。
{p 4 4 2}{cmdab:reps:(}{it:#}{cmd:)} bootstrap 重复次数（默认 999）；也用于 {cmd:cband}
的 multiplier 抽样。
{p 4 4 2}{cmdab:seed:(}{it:#}{cmd:)} 随机种子（默认 12345）。
{p 4 4 2}{cmdab:cluster:(}{it:varname}{cmd:)} 按聚类（而非单位）重抽样（默认用单位标识）。
{p 4 4 2}{cmd:cband} 报告统一置信带（sup-t），经 multiplier bootstrap。
{p 4 4 2}{cmd:graph} 绘制 ATT(d) 与 ACRT(d)。

{title:存储结果}

{p 4 4 2}
{cmd:r(attd)} 与 {cmd:r(acrt)} 为 {it:npoints} × 5 矩阵，列为
{cmd:d}、{cmd:ATT}/{cmd:ACRT}、{cmd:se}、{cmd:lb}、{cmd:ub}。
{cmd:r(degree)}、{cmd:r(dmin)}、{cmd:r(dmax)} 返回样条次数与剂量范围。
{p 4 4 2}
使用 {cmd:cband} 时：{cmd:r(crit_att)} 与 {cmd:r(crit_acrt)} 为 sup-t 临界值；
{cmd:r(cb_att)} 与 {cmd:r(cb_acrt)} 为 {it:npoints} × 4 矩阵，列为
{cmd:d}、估计值、{cmd:cb_lb}、{cmd:cb_ub}。

{title:示例}

{p 4 4 2}
{cmd:y} 关于连续剂量 {cmd:d} 的线性剂量-反应，并作图：
{p 8 8 2}{cmd:. contdid y, unit(id) time(t) dose(d) npoints(20) reps(999) seed(1) graph}
{p 4 4 2}
三次 B 样条 + 2 个内部节点 + 统一置信带：
{p 8 8 2}{cmd:. contdid y, unit(id) time(t) dose(d) degree(3) nknots(2) reps(999) seed(1) cband}
{p 4 4 2}
带协变量 x1、x2 的条件（强）平行趋势 + 统一带：
{p 8 8 2}{cmd:. contdid y, unit(id) time(t) dose(d) covariates(x1 x2) cband}

{title:参考文献}

{p 4 4 2}
Callaway, B., A. Goodman-Bacon, and P. H. C. Sant'Anna. 2024.
Difference-in-Differences with a Continuous Treatment. NBER Working Paper 32117.

{p 4 4 2}
Niu, H. 2026. contdid: Difference-in-Differences with a Continuous Treatment
for Stata. Version 1.1.0.

{title:另见}

{p 4 4 2}
英文帮助：{help contdid}；{help csdid}（二值交错 DiD）、{help drdid}（双重稳健 DiD）。
