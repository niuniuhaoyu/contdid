*=====================================================================
* contdid 演示：从数据到曲线
*=====================================================================
clear all
set more off
adopath + "D:/OpenCode/contdid"

*---------------------------------------------------------------------
* 第 1 步：数据长什么样
*---------------------------------------------------------------------
use "D:/OpenCode/contdid/data/contdid_sim.dta", clear
display as txt "{hline 78}"
display as txt "1. 数据结构"
display as txt "{hline 78}"
describe, short
list id t y d in 1/8, sepby(id) noobs
quietly count if d == 0
local n0 = r(N)/2
quietly count if d > 0
local n1 = r(N)/2
display "没被处理的单位 (d=0): `n0' 个；  有剂量的单位 (d>0): `n1' 个"
quietly summarize d if d > 0
display "剂量 d 的范围: " %5.3f r(min) " 到 " %5.3f r(max)

*---------------------------------------------------------------------
* 第 2 步：真相 + 把数据压成「每个单位一个 ΔY」
*---------------------------------------------------------------------
display as txt "{hline 78}"
display as txt "2. 真相与原始 DiD 分组均值"
display as txt "{hline 78}"
display "这个模拟数据的真实 ATT(d) = 0.5 * d  （论文里手算的那个对照）"
display ""

preserve
    bysort id (t): gen double dy = y[2] - y[1]
    by id: keep if _n == 1
    quietly summarize dy if d == 0
    local base = r(mean)
    display "未处理组的 ΔY 均值（基线）= " %8.4f `base'
    display ""
    display as txt "  剂量区间        n      E[ΔY|d]     E[ΔY|d]-基线   真实ATT(d)"
    display as txt "  ------------------------------------------------------------------"
    local cut1 = 0
    local cut2 = 0.25
    local cut3 = 0.5
    local cut4 = 0.75
    foreach hi in 0.25 0.50 0.75 1.00 {
        quietly summarize dy if d > `cut1' & d <= `hi'
        local m = r(mean)
        quietly count if d > `cut1' & d <= `hi'
        local nnn = r(N)
        local mid = (`cut1' + `hi')/2
        display "  (`cut1', `hi']  " %6.0f `nnn' "   " %9.4f `m' "   " ///
            %10.4f (`m'-`base') "   " %10.4f (0.5*`mid')
        local cut1 = `hi'
    }
restore

*---------------------------------------------------------------------
* 第 3 步：跑 contdid（线性版，对应手算）
*---------------------------------------------------------------------
display as txt "{hline 78}"
display as txt "3. contdid（默认线性）"
display as txt "{hline 78}"
contdid y, unit(id) time(t) dose(d) npoints(7) seed(1)

matrix A = r(attd)
matrix C = r(acrt)
display ""
display as txt "  评估剂量 d     ATT(d)        标准误      95%下限      95%上限      真实0.5d     ACRT(d)"
display as txt "  -----------------------------------------------------------------------------------------"
forvalues k = 1/7 {
    display "  " %10.4f A[`k',1] "  " %11.4f A[`k',2] "  " %10.4f A[`k',3] ///
        "  " %11.4f A[`k',4] "  " %11.4f A[`k',5] "  " %10.4f (0.5*A[`k',1]) ///
        "  " %10.4f C[`k',2]
}
display ""
display "注意：ATT(0) 这一行 ~ 0，这是归一化——效应是「相对未处理组」算的。"

*---------------------------------------------------------------------
* 第 4 步：灵活曲线 + 一致置信带 + 出图
*---------------------------------------------------------------------
display as txt "{hline 78}"
display as txt "4. B 样条曲线 + 一致置信带 + 图"
display as txt "{hline 78}"
contdid y, unit(id) time(t) dose(d) npoints(30) degree(3) nknots(2) cband reps(199) seed(1) graph
matrix B = r(attd)
matrix CB = r(cb_att)
display "B 样条曲线点数: " rowsof(B)
display "一致置信带 (cb_att) 首行:  d=" %6.3f CB[1,1] "  lb=" %7.4f CB[1,2] "  ub=" %7.4f CB[1,3]
graph export "D:/OpenCode/contdid/examples/dose_response_demo.png", replace width(1600)

display as txt "{hline 78}"
display "演示结束"
