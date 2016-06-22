function plotgas(varargin)
if nargin == 0
    theactualplot('','','', '', 'clear')
    return
else
    [A,C, error_vect, epoch_vect, nodes, skelldef, layertype] = varargin{:};
end
if size(A,1) >50 
    A = rebuild(A,skelldef,layertype);
end
[row,col] = find(C);
ax = A(1,row);
ay = A(2,row);
bx = A(1,col);
by = A(2,col);
X = reshape([ax;bx;NaN*ones(size(ax))],size(ax,2)*3,1)';
Y = reshape([ay;by;NaN*ones(size(ax))],size(ax,2)*3,1)'; 
if size(A,1)>=3&&size(A,1)<75&&size(A,1)~=72
    az = A(3,row);
    bz = A(3,col);
    Z = reshape([az;bz;NaN*ones(size(ax))],size(ax,2)*3,1)';
    theactualplot(A, error_vect, epoch_vect, nodes, X, Y, Z)
elseif size(A,1) == 75||size(A,1) == 72
    if size(A,1) == 72
        tdskel = zeros(24,3,size(A,2));
        for k = 1:size(A,2)
            for i=1:3
                for j=1:24
                    tdskel(j,i,k) = A(j+24*(i-1),k);
                end
            end
        end
        tdskel = cat(1,zeros(1,3,size(A,2)), tdskel);     
    else
        tdskel = zeros(25,3,size(A,2));
        for k = 1:size(A,2)
            for i=1:3
                for j=1:25
                    tdskel(j,i,k) = A(j+25*(i-1),k);
                end
            end
        end
    end
    if all(size(tdskel) ~= [25 3 size(A,2)])
        error('wrong skeleton building procedure!')
    end
    moresticks = [];
    for i=1:size(tdskel,1)
        for j=1:size(row)
            moresticks = cat(2,moresticks,[tdskel(i,:,row(j));tdskel(i,:,col(j)); [NaN NaN NaN]]');
        end
    end    
    SK = skeldraw(A(:,1),0);
    for i = 2:size(A,2)
       SK = [SK skeldraw(A(:,i),0)];      
    end
    T = [SK moresticks];    
    % make A into a sequence of 3d points
    A = threedeeA(A);
    theactualplot(A', error_vect, epoch_vect, nodes, T(1,:),T(2,:),T(3,:))
       
else
    theactualplot(A, error_vect, epoch_vect, nodes, X, Y) 
end
end
function AA = rebuild(A, skelldef, layertype) 
a = size(A,1)/3;
c = size(A,2);
switch layertype
    case 'pos'
        assumed_q = a/(size(skelldef.pos,2)/3);
    case 'vel'
        assumed_q = a/(size(skelldef.vel,2)/3);
    case 'all'
         if a > 25&&((size(skelldef.elementorder,2))<=a) % maybe this is useless
             assumed_q = a/size(skelldef.elementorder,2)*3;
         else
             assumed_q = 1;
             %disp('Warning. Unexpected condition while reshaping skeletons.')
         end
end
B = reshape(A,ceil(a/assumed_q),3,[]); 
if a == 49
    A = B(1:24,:,:);
    AA = reshape(A,72,c);
    return
else    
    AA = zeros(skelldef.length/6,3,c*assumed_q);
    switch layertype
        case 'pos'
            AA(skelldef.elementorder(1:a/assumed_q),:,:) = B;
            AA = reshape(AA,skelldef.length/2,c*assumed_q);
            return
        case 'vel'
             ini = size(skelldef.pos,2)/3+1;
             een = size(skelldef.vel,2)/3+ini-1;
             AA(skelldef.elementorder(ini:een)-skelldef.length/6,:,:) = B;
             AA = reshape(AA,skelldef.length/2,c*assumed_q);
         case 'all'
            a = size(skelldef.pos,2);
            AA(skelldef.elementorder(1:a/3),:,:) = B(1:a/3,:,:);
            AA = reshape(AA,[],c*assumed_q);
        otherwise
            error('unknown layer type.')
    end
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
