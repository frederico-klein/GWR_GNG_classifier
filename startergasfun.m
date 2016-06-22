function startergasfun(Data, params, gastype, varargin)
if nargin >3
    T = varargin{1};
    ColorLabels = varargin{2};
end
tic
A = gas_wrapper(Data,params,gastype);
if params.PLOTIT
    subplot(1,2,1)
    hold on
   if exist('ColorLabels', 'var')
       scatter(Data(1,:), Data(2,:),[],ColorLabels,'filled')
       scatter( A(1,:)', A(2,:)', 'r','filled')
   else
       plot(Data(1,:),Data(2,:), '.g', A(1,:)', A(2,:)', '.r')
   end
end

if exist('T', 'var') % if there is a target try to classify it
%%% simple labeller:
nodes_Y = simplelabeller(A, Data, T);
Y = simplelabeller(Data, A, nodes_Y);

figure
plotconfusion(T, Y);
end
toc

