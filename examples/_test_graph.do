*! _test_graph.do — graph option runs without error and exports a file
version 16
clear all
cd "/Users/niuhaoyu/Documents/open code/contdid"
adopath + "/Users/niuhaoyu/Documents/open code/contdid"
use "data/contdid_sim.dta", clear

contdid y, unit(id) time(t) dose(d) npoints(20) reps(99) seed(1) graph

* confirm the graph was created
graph export "examples/contdid_dose_response.png", width(1200) replace
di as result "PASS: graph option ran and exported"
