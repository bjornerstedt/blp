% Import painkillers data 

pk = readtable('painkillers.csv','TreatAsEmpty','NA');
catvar = {'firm','brand','substance','form'};

for i = 1:length(catvar) 
    pk.(catvar{i}) = categorical(pk.(catvar{i}));
end
pk.date = datetime( pk.year,pk.month,1);
pk.date.Format = 'yyyy-MM';

save 'painkillers' pk 