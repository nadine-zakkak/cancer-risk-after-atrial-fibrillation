
* looping through filenames for cancer diagnosis and cancer death analysis files
foreach outcome in "cancer_region" "cancer_organsys"{
* foreach outcome in "cancer_organsys"{
	foreach fup in "med_term" "long_term"{
		
		local filename "`outcome'_`fup'"
		display "`filename'"
		
		
		import delimited "./data/final_data/stata/`filename'.csv", varnames(1) clear
		
		/*
		encode alc_desc, generate(alc_num)
		
		
		sort subtype male
		by subtype male: gen random = runiform()
		by subtype male: keep if random < 0.01
		
		stset follow, failure(outcome_num = 1)
	*/
		levelsof subtype, local(subtype_value)

		* looping and extracting results
		foreach subtype in `subtype_value'{
			foreach male in 0 1 {
				import delimited "./data/final_data/stata/`filename'.csv", varnames(1) clear
				
				
				keep if subtype == `subtype' & male == `male'
				stset follow, failure(outcome_num = 1)

				
				/*
				gen random = runiform()
				keep if random < 0.01
				*/
				
				count if subtype == `subtype' & male == `male'
				local combination_exists = r(N)
				
				* only run analysis if subtype-sex combination exsits in data
				if `combination_exists' > 0 {
					display "`filename'_subtype`subtype'_male`male'"
					
					encode alc_desc, generate(alc_num)
					stset follow, failure(outcome_num = 1)
					
					cap frame drop results
					frame create results
			
					* add alcohol as covariate if analysing digestive system cancers
					if "`outcome'" == "cancer_organsys" & `subtype' == 7 {
						stcrreg i.case c.age_start i.diab i.ht i.smok i.alc_num ///
						if male == `male' & subtype == `subtype', ///
						compete(outcome_num == 2)	
					} 
					else {
						stcrreg i.case c.age_start i.diab i.ht i.smok if ///
						male == `male' & subtype == `subtype', ///
						compete(outcome_num == 2)	
					}
				
					matrix list r(table)
					*save transposed (') table of results
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
						* add subtype column to table of results
						gen subtype = `subtype'
						order var
						list
				
						export delimited using ///
						./results/stata/`filename'_male`male'_subtype`subtype'.csv, ///
						replace
					}
				}
			}
		}
	}
}


