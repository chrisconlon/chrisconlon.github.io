clear
set more off
insheet using visits.csv, comma
gen starttime2 = Clock(starttime,"YMDhms")
gen endtime2 = Clock(endtime,"YMDhms")
drop starttime endtime
rename starttime2 starttime
rename endtime2 endtime
format starttime endtime %tC
encode product_name, generate(product)
*encode machinerank, generate(machrank)
drop product_name
compress
gen netrev = (price-cost) * vends
gen grossrev = (price) * vends
gen year = year(dofC(endtime))
order revenueid revenueidlag interval new_product_code machine_code starttime endtime beginv endinv servinv vends vendout stockout price cost newpar maxpar  product netrev grossrev year  servicevisit
keep revenueid revenueidlag new_product_code machine_code interval  beginv endinv servinv vends vendout price cost newpar maxpar starttime endtime product netrev grossrev year servicevisit stockout
save visits.dta, replace

collapse (sum) vends beginv endinv netrev grossrev (min) stockout vendout, by(revenueid revenueidlag starttime endtime interval machine_code product year servicevisit)
keep revenueid-interval machine_code starttime endtime product  vends netrev endinv grossrev year servicevisit stockout vendout
reshape wide vends endinv netrev grossrev stockout vendout, i(revenueid-interval machine_code starttime endtime year servicevisit) j(product)
collapse (sum) vends* endinv* netrev* grossrev* (min) interval stockout* vendout*, by(revenueid machine_code endtime year servicevisit)
compress
egen totalinv = rsum(endinv*)
egen totalvends = rsum(vends*)
egen totalgrossrev = rsum(grossrev*)
egen totalrev = rsum(netrev*)
drop grossrev*
saveold  visits-wide, replace

drop vends* netrev* endinv*
save temp, replace
log using table.txt, text replace
table year machine, c(n totalgrossrev p25 totalgrossrev median totalgrossrev p75 totalgrossrev )
table year machine, c(n totalgrossrev p25 interval median interval p75 interval )
log close

clear
use visits
collapse (median) price cost , by(revenueid new_product_code )
save  pricelist, replace

clear
use visits
contract new_product_code product
save prodlist, replace

clear
use visits
contract revenueid starttime endtime
drop _freq
gen year = year(dofc(endtime))
rename starttime visitstart
rename endtime visitend
save times, replace


clear
insheet using daily.csv, names
gen endtime = Clock(dexdatetime,"YMDhms")
gen visit_day = dofc(clock(visitday,"YMDhms"))
format %tC endtime
format %td visit_day
*encode product_name, generate(product)
drop lagdextime dexdatetime visitday
rename visit_day visitday
format %td visitday

sort revenueid new_product_code
merge revenueid new_product_code using pricelist, nokeep
tab _merge
drop _merge


gen grossrev = price*vends
gen netrev = (price-cost)*vends
save  daily, replace

clear
use daily
keep revenueid machine_code dex_read_code new_product_code vends currentinv netrev grossrev endtime servicevisit stockout vendout visitday
egen duplicates = count(new_product_code) , by(revenueid machine_code new_product_code visitday servicevisit)
tab duplicates
collapse (sum) vends currentinv grossrev netrev (min) stockout vendout, by(revenueid machine_code dex_read_code new_product_code servicevisit visitday)
reshape wide vends currentinv grossrev netrev stockout vendout ,i(revenueid machine_code dex_read_code servicevisit visitday) j(new_product_code )
egen totalinv = rsum(currentinv*)
egen totalvends = rsum(vends*)
egen totalgrossrev = rsum(grossrev*)
egen totalrev = rsum(netrev*)
egen totalso = rsum(stockout*)
compress
sort revenueid
merge revenueid using times, nokeep
tab _merge
drop _merge
gen elapsed = day(visitday - dofc(visitstart))
save  daily-wide.dta, replace



clear
insheet using fulldataset.csv, names
gen starttime = Clock(lagdextime,"YMDhms")
gen endtime = Clock(dexdatetime,"YMDhms")
format %tC starttime endtime
*encode product_name, generate(product)
drop lagdextime dexdatetime

sort revenueid new_product_code
merge m:1 revenueid new_product_code using pricelist, keep(match master)
tab _merge
drop _merge


gen grossrev = price*vends
gen netrev = (price-cost)*vends
gen elapsed = hours(endtime-starttime)
save  fulldataset, replace

keep revenueid machine_code dex_read_code new_product_code vends currentinv netrev grossrev vendout stockout
egen duplicates = count(new_product_code) , by(revenueid machine_code dex_read_code new_product_code)
tab duplicates
collapse (sum) vends currentinv grossrev netrev (min) vendout stockout, by(revenueid machine_code dex_read_code new_product_code)
reshape wide vends currentinv grossrev netrev vendout stockout,i(revenueid machine_code dex_read_code ) j(new_product_code )
egen totalinv = rsum(currentinv*)
egen totalvends = rsum(vends*)
egen totalgrossrev = rsum(grossrev*)
egen totalrev = rsum(netrev*)
egen totalso = rsum(stockout*)
compress
sort dex_read_code
save  dataset-wide.dta, replace

use fulldataset
contract dex_read_code starttime endtime elapsed
merge dex_read_code using dataset-wide
drop _merge
save dataset-wide, replace

use fulldataset
drop if elapsed < 0
replace stockout = stockout*100
replace vendout = vendout*100
gen stockoutlb = stockout * (1-marginal)
gen vendoutlb = vendout * (1-marginal)
collapse (mean) stockout* vendout* [aweight=elapsed ], by(new_product_code )
merge new_product_code using prodlist
list product stockoutlb stockout vendoutlb vendout , clean noobs


