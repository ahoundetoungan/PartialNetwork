
clear
import sasxport "~/tmpah/SchoolData/Inschool.xpt"
keep sqid aid sschlcde s1 s2 s3 s4 s6a s6b s6c s6e s9 s10* s11 s12 s14 s15 s17 s18 s20 s21 s34a-s43e s44* s47 s48 s59a
destring sqid, replace
destring aid, replace
destring sschlcde, replace
sort sqid
drop if sqid == 999999
save ~/tmpah/SchoolData/cleandta, replace

clear
import sasxport "~/tmpah/SchoolData/sfriend.xpt"
destring sqid, replace
sort sqid
drop if sqid == 999999
merge 1:1 sqid using ~/tmpah/SchoolData/cleandta
drop _merge
save ~/tmpah/SchoolData/cleandta, replace

drop if aid == .
// Age 
replace s1 = . if s1 == 99
gen age = s1
gen age2 = age^2
// Sex, 1 is male and 2 is female
replace s2 = . if s2 == 9
gen male = (s2 == 1)
gen female = (s2 == 2)
replace male = . if s2 == .
replace female = . if s2 == .
// Grade
replace s3 = . if s3 == 13 | s3 == 99
gen grade = s3
// Hispanic
replace s4 = 0 if s4 == 8 | s4 == 9
gen hispanic = s4
// Race
gen racewhite = (s6a == 1)
gen raceblack = (s6b == 1)
gen raceasian = (s6c == 1)
gen raceother = (s6e == 1)
// Year in school
replace s9 = . if s9 == 99
gen yearinschl = s9
// GPA
gen nsubject = (s10a < 5) + (s10b < 5) + (s10c < 5) + (s10d < 5) 
replace s10a = 0 if s10a > 4
replace s10b = 0 if s10b > 4
replace s10c = 0 if s10c > 4
replace s10d = 0 if s10d > 4
egen gpa = rowtotal(s10a-s10d)
replace gpa = gpa/nsubject
replace gpa =. if nsubject == 0
drop nsubject
// Live with mother
replace s11 = . if s11 == 9
gen withmom = s11
// Mother education
gen melhigh = (s12 == 1 | s12 == 2 | s12 == 9 | s12 == 10)
gen mehigh = (s12 == 3 | s12 == 4)
gen memhigh = (s12 >= 5 & s12 <=8)
gen memiss = (s12 >= 11)
// Mother job
gen mjhome = (s14 == 1 | (s14 >= 16 & s14 <= 18))
gen mjprof = (s14 == 2 | s14 == 3)
gen mjother = (s14 >= 4 & s14 <= 15) | s14 == 19
gen mjmiss = (s14 >= 20)
// Mother works for pay
replace s15 = 0 if s15 > 1
gen mworkpay = s15
// Live with father
replace s17 = . if s17 == 9
gen withdad = s17
//Live with both
gen withbothpar = (withmom == 1 & withdad == 1)
replace withbothpar = . if (withmom == . | withdad == .)
// father education
gen felhigh = (s18 == 1 | s18 == 2 | s18 == 9 | s18 == 10)
gen fehigh = (s18 == 3 | s18 == 4)
gen femhigh = (s18 >= 5 & s18 <=8)
gen femiss = (s18 >= 11)
// Father job
gen fjhome = (s20 == 1 | (s20 >= 16 & s20 <= 18))
gen fjprof = (s20 == 2 | s20 == 3)
gen fjother = (s20 >= 4 & s20 <= 15) | s20 == 19
gen fjmiss = (s20 >= 20)
// father works for pay
replace s21 = 0 if s21 > 1
gen fworkpay = s21
// Number of cluds
egen nclubs = rowtotal(s44*)
// TV hours
replace s47 = . if s47==9
gen TVhours = s47
// Do my work
replace s48 = . if s48==9
gen DomyworkS48 = s48
// cigarettes
replace s59a = . if s59a == 99
gen cigarettes = s59a

drop s1* s2* s3* s4* s5* s6* s9*


//keep if inlist(sschlcde,4, 175, 90, 45, 24, 79, 145, 162, 34, 62, 158, 75, 65, 267, 93, 59, 369, 73, 47, 259, 269, 29, 147, 83, 87, 186, 27, 106, 58, 86, 12, 187, 57, 18, 71, 72, 31, 121, 194, 42, 43, 115, 60, 193, 19, 268, 183, 13, 16, 105, 53, 23, 41, 6, 271, 14, 9, 2, 76, 22, 171, 74, 20, 168, 191, 78, 141, 50, 114, 82, 170, 160, 66, 35, 89, 91, 371, 7, 33, 52, 77, 67, 54, 21, 119, 68, 64, 120, 15, 176) // keep only schools where >=90% of the students have grades

save ~/tmpah/SchoolData/cleandta, replace


