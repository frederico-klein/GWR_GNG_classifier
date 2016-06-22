function [A, C ,outparams] = gas_wrapper(data,varargin)
global VERBOSE LOGIT
VERBOSE = true;
LOGIT = true;
if nargin == 3
    params = varargin{1};
    gastype = varargin{2};
    arq_connect = struct();
else
    arq_connect = varargin{1};
    gastype = [];
end
if isfield(arq_connect, 'params')&&~isempty(arq_connect.params)
    params = arq_connect.params;
end
if isfield(arq_connect, 'method')&&~isempty(arq_connect.method)
    gastype = arq_connect.method;
elseif 0
    gastype = 'gwr';
end
if ~isfield(params, 'savegas')
    params.savegas.name = 'gas';
    params.savegas.resume = false;
    params.savegas.path = '~/';
end
if isempty(params)    
    NODES = 100;    
    params = struct();    
    params.use_gpu = false;
    params.PLOTIT = true; %
    params.RANDOMSTART = false; % if true it overrides the .startingpoint variable
    params.RANDOMSET = false;
    params.savegas.name = 'gas';
    params.savegas.resume = false;
    params.savegas.path = '~/';
    params.savegas.gasindex = 1;
    n = randperm(size(data,2),2);
    params.startingpoint = [n(1) n(2)];
    
    params.amax = 500; %greatest allowed age
    params.nodes = NODES; %maximum number of nodes/neurons in the gas
    params.en = 0.006; %epsilon subscript n
    params.eb = 0.2; %epsilon subscript b
    
    %Exclusive for gwr
    params.STATIC = true;
    params.MAX_EPOCHS = 1; % this means data will be run over twice
    params.at = 0.80; %activity threshold
    params.h0 = 1;
    params.ab = 0.95;
    params.an = 0.95;
    params.tb = 3.33;
    params.tn = 3.33;
    
    %Exclusive for gng
    params.age_inc                  = 1;
    params.lambda                   = 3;
    params.alpha                    = .5;     % q and f units error reduction constant.
    params.d                           = .99;   % Error reduction factor.
else
   %    
end
if  isfield(arq_connect, 'name')&&params.savegas.resume %%%resuming is not working.
     params.savegas.name = strcat(arq_connect.name,'-n', num2str(params.nodes), '-s',num2str(size(data,1)),'-q',num2str(params.q),'-i',num2str(params.savegas.gasindex));
elseif params.savegas.resume
    error('Strange arq_connect definition. ''.name'' field is needed.')
end

MAX_EPOCHS = params.MAX_EPOCHS;
PLOTIT = params.PLOTIT;

%%% things that are specific for skeletons:
if isfield(params,'skelldef')
    skelldef = params.skelldef;
else
    skelldef = [];
end
if isfield(params,'layertype')
    layertype = params.layertype;
else
    layertype = [];
end

if isfield(params,'plottingstep')        
    if params.plottingstep == 0
        plottingstep = size(data,2);
    else
        plottingstep = params.plottingstep;
    end    
else
    plottingstep = fix(size(data,2)/20);
end
if ~isfield(params,'use_gpu')|| gpuDeviceCount==0
    params.use_gpu = false;
end
if PLOTIT
    figure
    plotgas() % clears plot variables
end
datasetsize = size(data,2);
errorvect = nan(1,MAX_EPOCHS*datasetsize);
epochvect = nan(1,MAX_EPOCHS*datasetsize);
nodesvect = nan(1,MAX_EPOCHS*datasetsize);
if  params.use_gpu %do not use!
    data = gpuArray(data);
    errorvect = gpuArray(errorvect);
    epochvect = gpuArray(epochvect);
    nodesvect = gpuArray(nodesvect);
end
switch gastype
    case 'gwr'
        gasfun = @gwr_core;
        gasgas = gas;
        gasgas = gasgas.gwr_create(params,data);
    case 'gng'
        gasfun = @gng_core;
        gasgas = gas;
        gasgas = gasgas.gng_create(params,data);
    otherwise
        error('Unknown method.')
end
if isfield(params, 'savegas')&&params.savegas.resume
    %disabled.
end
therealk = 0; 
for num_of_epochs = 1:MAX_EPOCHS
    if params.RANDOMSET
        kset = randperm(datasetsize);
    else
        kset = 1:datasetsize;
    end
    for k = kset %step 1
        therealk = therealk +1;        
        gasgas = gasfun(data(:,k), gasgas);        
        errorvect(therealk) = gasgas.a;
        epochvect(therealk) = therealk;
        nodesvect(therealk) = gasgas.r-1;
        if PLOTIT&&mod(k,plottingstep)==0&&numlabs==1 
            plotgas(gasgas.A,gasgas.C,errorvect,epochvect,nodesvect, skelldef, layertype)
            drawnow
        end        
    end
    gasgas = gasgas.update_epochs(num_of_epochs);
    if isfield(params, 'savegas')&&isfield(params.savegas,'save') &&params.savegas.save
        save(strcat(savegas,'-e',num2str(num_of_epochs)), 'gasgas')
    end
end
outparams.graph.errorvect = errorvect;
outparams.graph.epochvect = epochvect;
outparams.graph.nodesvect = nodesvect;
outparams.accumulatedepochs = gasgas.params.accumulatedepochs;
outparams.initialnodes = [gasgas.ni1,gasgas.ni2];
A = gasgas.A;
C = gasgas.C;
end
