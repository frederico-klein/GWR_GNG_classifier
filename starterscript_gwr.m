load('Shapes.mat')
%Data = simplecluster_dataset;
%pkg load image %statistics
NODES = 100;

params = struct();

params.use_gpu = false;
params.PLOTIT = true; 
params.RANDOMSTART = false; % if true it overrides the .startingpoint variable
params.RANDOMSET = true;
n = randperm(size(Data,2),2);
params.startingpoint = [n(1) n(2)];

params.amax = 50; %greatest allowed age
params.nodes = NODES; %maximum number of nodes/neurons in the gas
params.en = 0.006; %epsilon subscript n
params.eb = 0.2; %epsilon subscript b
params.MAX_EPOCHS = 10; % this means data will be run over MAX_EPOCHS times

%%%Exclusive for gwr
params.STATIC = true;
params.at = 0.80; %activity threshold
params.h0 = 1;
params.ab = 0.95;
params.an = 0.95;
params.tb = 3.33;
params.tn = 3.33;

%%%Exclusive for gng
params.age_inc                  = 1;
params.lambda                   = 3;
params.alpha                    = .5;     % q and f units error reduction constant.
params.d                           = .99;   % Error reduction factor.

startergasfun(Data, params, 'gwr', T, ColorLabels);