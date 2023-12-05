
* looping through cancer diagnosis and cancer death analysis files and analysing
foreach outcome in "cancer_diag" "cancer_death"{
	foreach fup in "med_term" "long_term"{
		
		local temp "`outcome'_`fup'"
		display "`temp'"
		import delimited "./data/final_data/stata/`temp'.csv", varnames(1) clear
	
		stset follow, failure(outcome_num = 1)
		
		* looping and extracting results
		foreach male in 0 1 {
			cap frame drop results
			frame create results
	
			stcrreg i.case c.age_start i.diab i.ht i.smok if male == `male', compete(outcome_num == 2)
			matrix list r(table)
	* save transposed (') table of results
			matrix result_matrix = r(table)'
			matrix list result_matrix
	
			frame results {
				svmat result_matrix, names(col)
				list
				local rownames : rownames result_matrix
				display "`rownames'"
		
				gen var = ""
				local i 1
				foreach rowname in `rownames'{
					replace var = "`rowname'" in `i'
					local ++i
				}
				order var
				list
		
				export delimited using ./results/stata/`temp'_male`male'.csv, replace
			}
		}
	}
}

/*
gen byte outcome_new = mod(outcome_num, 3)
tab outcome_new outcome_num
label define outcome_new 0 "Censor" 1 "Event" 2 "Competing", replace
label values outcome_new outcome_new

gen byte male = sex_lb == "Men"

gen random = runiform()
keep if random < 0.01
stset follow, failure(outcome_num = 2)
stcrreg i.case c.age_start i.diab i.ht i.smok if male, compete(outcome_num == 3)
desc

export delimited "test.csv", replace
*/


*tab case outcome_new, row

/*
stset follow, failure(outcome_num = 1)


* looping and extracting results
foreach male in 0 1 {
	cap frame drop results
	frame create results
	
	stcrreg i.case c.age_start i.diab i.ht i.smok if male == `male', compete(outcome_num == 2)
	matrix list r(table)
	* save transposed (') table of results
	matrix result_matrix = r(table)'
	matrix list result_matrix
	
	frame results {
		svmat result_matrix, names(col)
		list
		local rownames : rownames result_matrix
		display "`rownames'"
		
		gen var = ""
		local i 1
		foreach rowname in `rownames'{
			replace var = "`rowname'" in `i'
			local ++i
		}
		order var
		list
		
		export delimited using ./results/stata/cancer_diag_med_male`male'.csv
	}
}
*/
