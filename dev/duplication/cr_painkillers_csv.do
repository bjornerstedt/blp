use painkillers, clear
drop year?*
drop productname
ren product productname
gen product = productname
order year month product
export delimited using "painkillers.csv", quote replace
