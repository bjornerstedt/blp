cd "~/GitHub/blp/dev/logit"
import delimited "~/GitHub/blp/dev/logit/data.csv", clear

egen obs = group(marketid personid)
clogit sh p x , group(obs)

nlogitgen tt = type(fast:1, slow:2)

nlogit sh p x || tt: || productid: , case(obs)
