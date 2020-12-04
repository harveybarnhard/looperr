* This program produces timing results of the regressby function to
* benchmark against the fastols_by function in R (benchmarking.R)

/* timer helpers */
* https://github.com/matthieugomez/benchmark-stata-r/blob/master/code/3-benchmark-stata.do
cap program drop Tic
program define Tic
	syntax, n(integer)
	timer on `n'
end

cap program drop Toc
program define Toc
	syntax, n(integer)
	timer off `n'
end

set processors 2

local i = 0
foreach n in 10 100 1000 10000 100000 1000000 {
	foreach grps in 10 100 1000 10000 100000 1000000 {
		if(`n'*`grps'>10000000){
			continue
		}
		import delimited using "${looperr}/data-raw/benchmark`n'_`grps'", clear
		Tic, n(`++i')
		regressby y x, by(g)
		Toc, n(`i')
	}
}

drop _all
gen seconds = .
gen nObs = .
gen nGroups = .
set obs `i'
timer list
forval j = 1/`i'{
	replace seconds = r(t`j') if _n == `j'
}
local j = 1
foreach n in 10 100 1000 10000 100000 1000000 {
	foreach grps in 10 100 1000 10000 100000 1000000 {
		if(`n'*`grps'>10000000){
			continue
		}
		replace nObs = `n' if _n == `j'
		replace nGroups = `grps' if _n == `j'
		local j = `j' + 1
	}
}
gen expr = "regressby"
outsheet using "${looperr}/data-raw/results_stata.csv", replace
