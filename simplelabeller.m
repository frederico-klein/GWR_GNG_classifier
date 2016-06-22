function labels = simplelabeller(nodes, data, y)
% Usage:
%
% labels = simplelabeller(nodes, data, y)
try
    [~,ni] = pdist2(data',nodes', 'euclidean', 'Smallest',1); %will fail if stats package is not installed
    labels = y(:,ni);
catch
    [labels, ~ ]= labelling(nodes, data, y);
end
end
function [labels, ni1 ]= labelling(nodes, data, y)
maxmax = size(nodes,2);
labels = zeros(1,maxmax);
for i = 1:maxmax
    [~, ~, ni1 , ~ , ~] = findnearest(nodes(:,i), data); 
    labels(i) = y(ni1);
end
end