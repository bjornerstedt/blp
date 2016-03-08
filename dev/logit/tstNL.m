%% Test 1: NLIndividualDemand
display '**********************  Test 1  *************************'
test = false;

m = SimMarket();
m.model.markets = 50;
if test
    m.demand = NLDemand
    m.demand.alpha = .3;
    sresults = m.create()
    return
else
    m.demand = NLIndividualDemand
    m.demand.settings.nind = 1000;
    m.demand.alpha = .3;
    sresults = m.create();
end
display(m.model)
data = m.data;
data.q = [];
data = repmat(data, m.demand.settings.nind, 1);
data.personid = reshape(repmat(1:m.demand.settings.nind, size(m.data, 1), 1), [], 1);
data.sh = reshape(m.data.q, [], 1);
writetable(data, 'data.csv');