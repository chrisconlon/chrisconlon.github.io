clear
set more off
log using tables.txt, text replace
use daily-wide
drop if elapsed ==. | dex_read_code == .
recode year .=2004

regress servicevisit totalinv
regress servicevisit totalso
regress servicevisit totalinv totalso 
regress servicevisit totalinv elapsed
regress servicevisit totalso elapsed
regress servicevisit totalinv totalso elapsed
xi i.machine_code
regress servicevisit totalinv elapsed _I*
regress servicevisit totalso elapsed _I*
regress servicevisit totalinv totalso elapsed _I*

  


probit servicevisit totalinv
est store prob1
probit servicevisit totalso
est store prob2
probit servicevisit totalinv totalso 
est store prob3
probit servicevisit totalinv elapsed
est store prob4
probit servicevisit totalso elapsed
est store prob5
probit servicevisit totalinv totalso elapsed
est store prob6
probit servicevisit totalinv elapsed _I*
est store prob7
probit servicevisit totalso elapsed _I*
est store prob8
probit servicevisit totalinv totalso elapsed _I*
est store prob9
outreg2  [prob1 prob2 prob3 prob4 prob5 prob6 prob7 prob8 prob9] using elapsed.xls, excel replace label rdec(4) nonotes

probit servicevisit totalinv year
est store prob1
probit servicevisit totalso year
est store prob2
probit servicevisit totalinv totalso year
est store prob3
probit servicevisit totalinv elapsed year
est store prob4
probit servicevisit totalso elapsed year
est store prob5
probit servicevisit totalinv totalso elapsed year
est store prob6
probit servicevisit totalinv elapsed _I* year
est store prob7
probit servicevisit totalso elapsed _I* year
est store prob8
probit servicevisit totalinv totalso elapsed _I* year
est store prob9
outreg2  [prob1 prob2 prob3 prob4 prob5 prob6 prob7 prob8 prob9] using elapsedy.xls, excel replace label rdec(4) nonotes




gen totalso03 = (year==2003) * totalso
gen totalso04 = (year==2004) * totalso

gen totalinv03 = (year==2003) * totalinv
gen totalinv04 = (year==2004) * totalinv

gen elapsed03 = (year==2003) * elapsed
gen elapsed04 = (year==2004) * elapsed

drop elapsed totalso totalinv

probit servicevisit totalinv* year
est store prob1
probit servicevisit totalso* year
est store prob2
probit servicevisit totalinv* totalso* year
est store prob3
probit servicevisit totalinv* elapsed* year
est store prob4
probit servicevisit totalso* elapsed* year
est store prob5
probit servicevisit totalinv* totalso* elapsed* year
est store prob6
probit servicevisit totalinv* elapsed* _I* year
est store prob7
probit servicevisit totalso* elapsed* _I* year
est store prob8
probit servicevisit totalinv* totalso* elapsed* _I* year
est store prob9
outreg2  [prob1 prob2 prob3 prob4 prob5 prob6 prob7 prob8 prob9] using elapsedyear.xls, excel replace label rdec(4) nonotes



stset visitday , id(revenueid  ) failure(servicevisit )
*duplicates tag revenueid endday, generate(
streg totalso* , dist(exp)
est store haz1
test _b[totalso03]=_b[totalso04]

streg totalinv*, dist(exp)
est store haz2
test _b[totalinv03]=_b[totalinv04]

streg totalso* totalinv*, dist(exp)
est store haz3
test _b[totalinv03]=_b[totalinv04]
test _b[totalso03]=_b[totalso04]

streg year totalso* , dist(exp) nocons
est store haz4
test _b[totalso03]=_b[totalso04]

streg year totalinv*, dist(exp) nocons
est store haz5
test _b[totalinv03]=_b[totalinv04]

streg year totalso* totalinv*, dist(exp) nocons
est store haz6
test _b[totalinv03]=_b[totalinv04]
test _b[totalso03]=_b[totalso04]

streg year totalso* totalinv* elapsed* , dist(exp) nocons
est store haz7
test _b[totalinv03]=_b[totalinv04]
test _b[totalso03]=_b[totalso04]
test _b[elapsed03]=_b[elapsed04]

outreg2  [haz1 haz2 haz3 haz4 haz5 haz6 haz7] using hazard.xls, excel replace label rdec(4) nonotes  eform 

xi i.machine_code

streg totalso* _I* , dist(exp)
est store haz1
test _b[totalso03]=_b[totalso04]

streg totalinv* _I*, dist(exp)
est store haz2
test _b[totalinv03]=_b[totalinv04]

streg totalso* totalinv* _I*, dist(exp)
est store haz3
test _b[totalinv03]=_b[totalinv04]
test _b[totalso03]=_b[totalso04]

streg year totalso* _I*, dist(exp) nocons
est store haz4
test _b[totalso03]=_b[totalso04]

streg year totalinv* _I*, dist(exp) nocons
est store haz5
test _b[totalinv03]=_b[totalinv04]

streg year totalso* totalinv* _I*, dist(exp) nocons
est store haz6
test _b[totalinv03]=_b[totalinv04]
test _b[totalso03]=_b[totalso04]

streg year totalso* totalinv* elapsed* _I*, dist(exp) nocons
est store haz7
test _b[totalinv03]=_b[totalinv04]
test _b[totalso03]=_b[totalso04]
test _b[elapsed03]=_b[elapsed04]

outreg2  [haz1 haz2 haz3 haz4 haz5 haz6 haz7] using hazard2.xls, excel replace label rdec(4) nonotes  eform 
log close

