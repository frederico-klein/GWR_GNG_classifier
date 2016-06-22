function plotgas(varargin)
if nargin == 0
    theactualplot('','','', '', 'clear')
    return
else
    [A,C, error_vect, epoch_vect, nodes, skelldef, layertype] = varargin{:};
end
[row,col] = find(C);
ax = A(1,row);
ay = A(2,row);
bx = A(1,col);
by = A(2,col);
X = reshape([ax;bx;NaN*ones(size(ax))],size(ax,2)*3,1)';
Y = reshape([ay;by;NaN*ones(size(ax))],size(ax,2)*3,1)'; 
if size(A,1)>=3
    az = A(3,row);
    bz = A(3,col);
    Z = reshape([az;bz;NaN*ones(size(ax))],size(ax,2)*3,1)';
    theactualplot(A, error_vect, epoch_vect, nodes, X, Y, Z)
else
    theactualplot(A, error_vect, epoch_vect, nodes, X, Y) 
end
end
function theactualplot(A, error_vect, epoch_vect, nodes, varargin)
persistent hdl_main hdl_main_p hdl_error hdl_nodes axmain %axmaincopy
if strcmp(varargin{1},'clear')
    clear hdl_main hdl_main_p hdl_error hdl_nodes axmain
    return
end
if isempty(hdl_main)||isempty(hdl_error)||isempty(hdl_nodes) % initialize the plot window
    axmain = subplot(1,2,1);
end
if length(varargin)==3
    if isempty(hdl_main)
        hdl_main = plot3(axmain,varargin{1},varargin{2},varargin{3});
        hold on
        hdl_main_p = plot3(axmain,A(1,:),A(2,:),A(3,:), '.r');
        hold off
    else
        set(hdl_main, 'XData',varargin{1},'YData',varargin{2},'ZData',varargin{3});
        set(hdl_main_p, 'XData',A(1,:),'YData',A(2,:),'ZData',A(3,:));
    end
else
    if isempty(hdl_main)
        hdl_main = plot(axmain,varargin{1},varargin{2});
        hold on
        hdl_main_p = plot(axmain,A(1,:),A(2,:), '.r');
        hold off
    else
        set(hdl_main, 'XData',varargin{1},'YData',varargin{2});
        set(hdl_main_p, 'XData',A(1,:),'YData',A(2,:));
    end    
end
set(gca,'box','off')
%%% will set the same scale for axmain and axmaincopy
%%% Now I will plot the error
axerror = subplot(2,2,2);
if ~isempty(hdl_error)
    set(hdl_error, 'XData',epoch_vect,'YData',error_vect);
else
    title('Activity or RMS error')
    hdl_error = plot(axerror, epoch_vect, error_vect);
end
axnodes = subplot(2,2,4);
if ~isempty(hdl_nodes)
    set(hdl_nodes, 'XData',epoch_vect,'YData',nodes);
else
    hdl_nodes = plot(axnodes, epoch_vect, nodes);
    title('Number of Nodes')
end
end
function B = threedeeA(A)
C = makefatskel(A);
B = [];
for i = 1:size(C,3)
    B = cat(1, C(:,:,i), B);
end
end
