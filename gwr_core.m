function [gas]= gwr_core(eta,gas)
%eta = data(:,k); % this the k-th data sample
[gas.gwr.ws, ~, gas.gwr.s, gas.gwr.t, ~] = findnearest(eta, gas.A); %step 2 and 3
if gas.C(gas.gwr.s,gas.gwr.t)==0 %step 4
    gas.C = spdi_bind(gas.C,gas.gwr.s,gas.gwr.t);
else
    gas.C_age = spdi_del(gas.C_age,gas.gwr.s,gas.gwr.t);
end
gas.a = exp(-norm((eta-gas.gwr.ws).*gas.awk)); %step 5

[gas.gwr.neighbours] = findneighbours(gas.gwr.s, gas.C);
gas.gwr.num_of_neighbours = size(gas.gwr.neighbours,2);

if (gas.a < gas.params.at) && (gas.r <= gas.params.nodes) %step 6
    gas.gwr.wr = 0.5*(gas.gwr.ws+eta); %too low activity, needs to create new node r
    gas.A(:,gas.r) = gas.gwr.wr;
    gas.C = spdi_bind(gas.C,gas.gwr.t,gas.r);
    gas.C = spdi_bind(gas.C,gas.gwr.s,gas.r);
    gas.C = spdi_del(gas.C,gas.gwr.s,gas.gwr.t);
    gas.r = gas.r+1;
else %step 7
    gas.gwr.wi = gas.A(:,gas.gwr.neighbours);
    gas.A(:,gas.gwr.neighbours) = gas.gwr.wi + gas.params.en*(repmat(eta,1,gas.gwr.num_of_neighbours)-gas.gwr.wi).*repmat(gas.h(gas.gwr.neighbours),size(eta,1),1);
    gas.A(:,gas.gwr.s) = gas.gwr.ws + gas.params.eb*gas.h(gas.gwr.s)*(eta-gas.gwr.ws); 
end
%step 8 : age edges with end at s
%first we need to find if the edges connect to s
gas.C_age = spdi_add(gas.C_age,gas.gwr.s,gas.gwr.neighbours);
%step 9: again we do it inverted, for loop first
%%%% this strange check is a speedup for the case when the algorithm is static
if gas.params.STATIC % skips this if algorithm is static
    gas.h = gas.hizero;
    gas.h(gas.gwr.s) = gas.hszero;
else
    for i = 1:gas.r
        gas.h(i) = gas.hi(gas.time,gas.params);
    end
    gas.h(gas.gwr.s) = gas.hs(gas.time,gas.params);
    gas.time = (cputime - gas.t0)*1;
end
%step 10: check if a node has no edges and delete them
if gas.r > gas.params.nodes
    R = gas.params.nodes;
else
    R = gas.r;
end
if gas.r>2 % don't remove everything 
    [gas.C(1:R,1:R), gas.C_age(1:R,1:R) ] = removeedge(gas.C(1:R,1:R), gas.C_age(1:R,1:R), gas.params.amax);
    [gas.C(1:R,1:R), gas.A(:,1:R), gas.C_age(1:R,1:R), gas.h, gas.r ] = removenode(gas.C(1:R,1:R), gas.A(:,1:R), gas.C_age(1:R,1:R), gas.h, gas.r);  %inverted order as it says on the algorithm to remove points faster
end
end
function sparsemat = spdi_add(sparsemat, a, b)
sparsemat(a,b) = sparsemat(a,b) + 1;
sparsemat(b,a) = sparsemat(a,b) + 1;
end
function sparsemat = spdi_bind(sparsemat, a, b) % adds a 2 way connection
sparsemat(a,b) = 1;
sparsemat(b,a) = 1;
end
function sparsemat = spdi_del(sparsemat, a, b) % removes a 2 way connection.
sparsemat(a,b) = 0;
sparsemat(b,a) = 0;
end
function [C, C_age ] = removeedge(C, C_age, amax) 
[row, col] = find(C_age > amax);
a = size(row,2);
if ~isempty(row)
    for i = 1:a
        C_age(row(i),col(i)) = 0;
        C_age(col(i),row(i)) = 0;
        C(row(i),col(i)) = 0;
        C(col(i),row(i)) = 0;
    end
end
end
function [C, A, C_age, h,r ] = removenode(C, A, C_age, h,r) %depends only on C operates on everything
[row,~] = find(C);
maxa = max(row);
pointstoremove = find(any(bsxfun(@eq, row, 1:maxa))==0);
if ~isempty(pointstoremove)
    numnum = length(pointstoremove); 
    C = clipsimmat(C,pointstoremove,numnum);
    A = clipA(A,pointstoremove,numnum);
    C_age = clipsimmat(C_age,pointstoremove,numnum);
    h = clipvect(h,pointstoremove,numnum);
    r = r-numnum;    
end
end
function C = clipsimmat(C,i,n)
C(i,:) = [];
C(:,i) = [];
C = padarray(C,[n n],'post');
end
function V = clipvect(V, i, n)
V(i) = [];
V = [V zeros(1,n)];
end
function A = clipA(A, i,n)
A(:,i) = [];
ZERO = zeros(size(A,1),n);
A = [A ZERO];
end
function neighbours = findneighbours(s,C)
ne = find(C(s,:));
neighbours = ne(ne~=0);
end
