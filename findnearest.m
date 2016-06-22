function [n1, n2, ni1, ni2, distvector] = findnearest(p,data)
maxindex = size(data,2);
ppp = repmat(p,1,maxindex);
dada = data - ppp;
dada = dada.*dada;
distvector = sum(dada).^(1/2);
[~,ni1] = min(distvector);
n1 = data(:,ni1);
pushdist = distvector(ni1);
distvector(ni1) = NaN; 
[~,ni2] = min(distvector);
n2 = data(:,ni2);
distvector(ni1) = pushdist;
end
