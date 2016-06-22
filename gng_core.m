function gasgas = gng_core(eta, gasgas)

PARAMS = gasgas.params;
% Unsupervised Self Organizing Map. Growing Neural Gas (GNG) Algorithm.

% Main paper used for development of this neural network was:
% Fritzke B. "A Growing Neural Gas Network Learns Topologies", in
%                         Advances in Neural Information Processing Systems, MIT Press, Cambridge MA, 1995.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NumOfEpochs   = 5;
%NumOfSamples = fix(size(Data,2)/NumOfEpochs);
%age_inc               = PARAMS.age_inc;
%max_age             = PARAMS.amax;
max_nodes               = PARAMS.nodes;
%eb                         = PARAMS.eb;
%en                         = PARAMS.en;
lamda                   = PARAMS.lambda;%3;
alpha                    = PARAMS.alpha;%.5;     % q and f units error reduction constant.
d                           = PARAMS.d;%.99;   % Error reduction factor.

%RMSE                  = zeros(1,NumOfEpochs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the params vector where the GNG algorithm parameters are stored:
params = [ PARAMS.age_inc;
    PARAMS.amax;
    PARAMS.nodes;
    PARAMS.eb;
    PARAMS.en;
    PARAMS.alpha;
    PARAMS.lambda;
    PARAMS.d;
    0;   % Here, we insert the sample counter.
    1;   % This is reserved for s1 (bmu);
    2;]; % This is reserved for s2 (secbmu);

n = gasgas.n;
ages = gasgas.C_age;
nodes = gasgas.A;
edges = gasgas.C;
errorvector = gasgas.errorvector;


NumOfNodes = size(nodes,2);
params(9) = gasgas.n;

% Step 2. Find the two nearest units s1 and s2 to the new data sample.
[s1, s2, distances] = findTwoNearest(eta,nodes);
params(10) = s1;
params(11) = s2;

% Steps 3-6. Increment the age of all edges emanating from s1 .
[nodes, edges, ages, errorvector] = edgeManagement(eta,nodes,edges,ages,errorvector, distances, params);

% Step 7. Dead Node Removal Procedure.
[nodes, edges, ages, errorvector] = removeUnconnected(nodes, edges,ages,errorvector);

% Step 8. Node Insertion Procedure.
if mod(n,lamda)==0 && size(nodes,2)<max_nodes
    [nodes,edges,ages,errorvector] = addNewNeuron(nodes,edges,ages,errorvector,alpha);
end

% Step 9. Finally, decrease the error of all units.
gasgas.errorvector = d*errorvector;

%update n
n = n+1;
%%%% writing stuff to gas:

gasgas.C_age = ages;
gasgas.A = nodes;
gasgas.C = edges;
gasgas.n = n;
gasgas.a = norm(errorvector)/sqrt(NumOfNodes);
gasgas.r = NumOfNodes+1;

end

function [s1 s2 distances] = findTwoNearest(Input,nodes) %#codegen

NumOfNodes = size(nodes,2);

distances = zeros(1,NumOfNodes);
for i=1:NumOfNodes
    distances(i) = norm(Input - nodes(:,i));
end

[sdistances indices] = sort(distances);

s1 = indices(1);
s2 = indices(2);

% MEX Code Generation:

% mexcfg = coder.config('mex');
% mexcfg.DynamicMemoryAllocation = 'AllVariableSizeArrays'; 
% codegen -config mexcfg findTwoNearest -args {coder.typeof(In(:,n),[Inf 1]),coder.typeof(nodes,[Inf Inf])}

end

function [nodes, edges, ages,error] = edgeManagement(Input,nodes,edges,ages,error,distances,params)  %#codegen

age_inc   = params(1);
max_age = params(2);
            eb = params(4);
            en = params(5);
            s1 = params(10);
            s2 = params(11);

NumOfNodes = size(nodes,2);

% Step 3. Increment the age of all edges emanating from s1. 
s1_Neighbors = find(edges(:,s1)==1);
SizeOfNeighborhood = length(s1_Neighbors);

ages(s1_Neighbors,s1) = ages(s1_Neighbors,s1) + age_inc;
ages(s1,s1_Neighbors) = ages(s1_Neighbors,s1);

% Step 4. Add the squared distance to a local error counter variable:
error(s1) = error(s1) + distances(s1)^2;

% Step 5. Move s1 and its topological neighbors towards î.
nodes(:,s1) = nodes(:,s1) + eb*(Input-nodes(:,s1));
nodes(:,s1_Neighbors) = nodes(:,s1_Neighbors) + en*(repmat(Input,[1 SizeOfNeighborhood])-nodes(:,s1_Neighbors));

% Step 6.
% If s1 and s2 are connected by an edge, set the age of this edge to zero.
% If such an edge does not exist, create it.
edges(s1,s2) = 1;
edges(s2,s1) = 1;
ages(s1,s2) = 0;
ages(s2,s1) = 0;

% Step 7. Remove edges with an age>max_age.
[DelRow DelCol] = find(ages>max_age);
SizeDeletion = length(DelRow);
for i=1:SizeDeletion
    edges(DelRow(i),DelCol(i)) = 0;
      ages(DelRow(i),DelCol(i)) = NaN;
end

% MEX-code generation:

% mexcfg = coder.config('mex');
% mexcfg.DynamicMemoryAllocation = 'AllVariableSizeArrays'; 
% codegen -config mexcfg edgeManagement.m -args {coder.typeof(In(:,n),[Inf 1]),coder.typeof(nodes,[Inf Inf]),coder.typeof(edges,[Inf Inf]),coder.typeof(ages,[Inf Inf]),coder.typeof(error,[1 Inf]),coder.typeof(distances,[1 Inf]),coder.typeof(params,[11,1])}
end

function  [nodes,edges,ages,error] = addNewNeuron(nodes,edges,ages,error,alpha)  %#codegen 
                                                                                                                                    % Checks whether this function is
                                                                                                                                    % suitable for automatic .mex code generation.
NumOfNodes = size(nodes,2);
    
[max_error q] = max(error);
       
% Find q-Neighborhood
  q_Neighbors = find(edges(:,q)==1);
  
% Find the neighbor f with the largest accumulated error. 
  [value index] = max(error(q_Neighbors));
  f = q_Neighbors(index); 
    
% Add the new node half-way between nodes q and f: 
   nodes = [nodes .5*(nodes(:,q) + nodes(:,f))];
   
% Remove the original edge between q and f.
   edges(q,f) = 0;
   edges(f,q) = 0;
   ages(q,f) = NaN;
   ages(f,q) = NaN;
   
   NumOfNodes = NumOfNodes + 1;
   r = NumOfNodes;
   
   % Insert edges connecting the new unit r with units q anf f. 
   edges = [edges  zeros(NumOfNodes-1,1)];
   edges = [edges; zeros(1,NumOfNodes)];
     
   edges(q,r) = 1;
   edges(r,q) = 1;
   edges(f,r) = 1;
   edges(r,f) = 1;
  
   ages = [ages  NaN*ones(NumOfNodes-1,1)];
   ages = [ages; NaN*ones(1,NumOfNodes)];
   
   ages(q,r) = 0;
   ages(r,q) = 0;
   ages(f,r) = 0;
   ages(r,f) = 0;
   
   error(q) = alpha*error(q);
   error(f) = alpha*error(f);
   
   error = [error error(q)];

% MEX-code Generation:  
   
% mexcfg = coder.config('mex');
% mexcfg.DynamicMemoryAllocation = 'AllVariableSizeArrays'; 
% codegen -config mexcfg addNewNeuron.m -args {coder.typeof(nodes,[Inf Inf]), coder.typeof(edges,[Inf Inf]), coder.typeof(ages,[Inf Inf]), coder.typeof(error,[1 Inf]), double(0)}
end
function [nodes, edges, ages, error] = removeUnconnected(nodes, edges,ages,error) %#codegen

NumOfNodes = size(nodes,2);

i = 1;
while NumOfNodes >= i
    if any(edges(i,:)) == 0
        
        edges = [edges(1:i-1,:); edges(i+1:NumOfNodes,:);];
        edges = [edges(:,1:i-1)  edges(:,i+1:NumOfNodes);];
       
        ages = [ages(1:i-1,:); ages(i+1:NumOfNodes,:);];
        ages = [ages(:,1:i-1)  ages(:,i+1:NumOfNodes);];

        nodes = [nodes(:,1:i-1) nodes(:,i+1:NumOfNodes);];
        error = [error(1,1:i-1) error(1,i+1:NumOfNodes);];
        
        NumOfNodes = NumOfNodes - 1;
        
        i = i -1; 
    end
    i = i+1;
end

% MEX-code generation

% mexcfg = coder.config('mex');
% mexcfg.DynamicMemoryAllocation = 'AllVariableSizeArrays'; 
% codegen -config mexcfg removeUnconnected -args {coder.typeof(nodes,[Inf Inf]),coder.typeof(edges,[Inf Inf]),coder.typeof(ages,[Inf Inf]),coder.typeof(error,[1 Inf])}
end
