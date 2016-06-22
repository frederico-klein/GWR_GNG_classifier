classdef gas
    properties
        params
        A
        C
        C_age
        errorvector
        hizero
        hszero
        time
        t0
        n
        r
        h
        a
        awk
        ni1
        ni2
        n1n2
        gwr
    end
    methods
        function X = hs(gasgas, t)
            h0 = gasgas.params.h0;%= 1;
            ab = gasgas.params.ab;%= 0.95;
            tb = gasgas.params.tb;%= 3.33;
            X = h0 - gasgas.S(t)/ab*(1-exp(-ab*t/tb));
        end
        function X = hi(gasgas, t)
            h0 = gasgas.params.h0;%= 1;
            an = gasgas.params.an;%= 0.95;
            tn = gasgas.params.tn;%= 3.33;
            X = h0 - gasgas.S(t)/an*(1-exp(-an*t/tn));
        end
        function X = S(~,t)
            X = 1;
        end
        function [n1n2, ni1,ni2] = initialnodes(gasgas, data)            
            RANDOMSTART = gasgas.params.RANDOMSTART;            
            % (1)
            % pick n1 and n2 from data
            ni1 = 1;
            ni2 = 2;
            if RANDOMSTART
                nn = randperm(size(data,2),2);
                ni1 = nn(1);
                ni2 = nn(2);
            elseif isfield(gasgas.params,'startingpoint')&&~isempty('params.startingpoint')
                ni1 = gasgas.params.startingpoint(1);
                ni2 = gasgas.params.startingpoint(2);
            end
            if ni1> size(data,2)
                dbgmsg('n1 overflow. using last data point',num2str(ni1),1)
                n1 = data(:,end);
            else
                n1 = data(:,ni1);
            end
            if ni2> size(data,2)
                dbgmsg('n2 overflow. using last data point',num2str(ni2),1)
                n2 = data(:,end);
            else
                n2 = data(:,ni2);
            end
            n1n2 = [n1, n2];
        end
        function gasgas = gas_create(gasgas, params,data)
            gasgas.params = params;
            gasgas.params.accumulatedepochs = 0;
            gasgas = gasgas.gas_finalize;
            [gasgas.n1n2, gasgas.ni1,gasgas.ni2] = initialnodes(gasgas, data);            
            %%%% 
            if gasgas.params.use_gpu  
                gasgas.n1n2 = gather(gasgas.n1n2);
            end  
            if isfield(gasgas.params, 'skelldef')&&isfield(gasgas.params.skelldef, 'awk')&&isfield(params, 'layertype')
                switch gasgas.params.layertype
                    case 'pos'
                        gasgas.awk = repmat(params.skelldef.awk.pos,1,params.q(1));
                    case 'vel'
                        gasgas.awk = repmat(params.skelldef.awk.vel,1,params.q(1));
                end
                try
                    n1.*gasgas.awk;
                catch
                    dbgmsg('Tried to use your awk definition, but I failed.',1)
                    gasgas.awk = ones(size(gasgas.n1n2,1),1);
                end
            else
                gasgas.awk = ones(size(gasgas.n1n2,1),1);
            end
        end
        function gasgas = gng_create(gasgas,params,data)
            
            gasgas = gasgas.gas_create(params,data);            
            %%% gng specific init                        
            gasgas.A = zeros(size(gasgas.n1n2,1),2);
            gasgas.A(:,[1 2]) = gasgas.n1n2;                        
            % Initial connections (edges) matrix.
            edges = [0  1;
                1  0;];
            gasgas.C = edges;
            % Initial ages matrix.
            ages = [ NaN  0;
                0  NaN;];
            gasgas.C_age = ages;
            % Initial Error Vector.
            gasgas.errorvector = [0 0];            
            gasgas.n = 0; %samplercounter
            gasgas = gasgas.gas_finalize;
        end
        function gasgas = gwr_create(gasgas, params, data)            
            gasgas = gasgas.gas_create(params,data);
            %%% creating gwr structure
            gasgas.gwr = struct('s',[],'t',[],'wi',[],'wr',[],'ws',[],'num_of_neighbours',[],'neighbours',[]);
            gasgas.A = zeros(size(gasgas.n1n2,1),params.nodes);
            gasgas.A(:,[1 2]) = gasgas.n1n2;            
            %%%gwr specific initialization
            %test some algorithm conditions:
            if ~(0 < params.en || params.en < params.eb || params.eb < 1)
                error('en and/or eb definitions are wrong. They are as: 0<en<eb<1.')
            end            
            % (2)
            % initialize empty set C            
            gasgas.C = sparse(params.nodes,params.nodes); % this is the connection matrix.
            gasgas.C_age = gasgas.C;            
            gasgas.r = 3; %the first point to be added is the point 3 because we already have n1 and n2
            gasgas.h = zeros(1,params.nodes);%firing counter matrix            
            %%% SPEEDUP CHANGE
            if params.STATIC
                gasgas.hizero = gasgas.hi(0)*ones(1,params.nodes);
                gasgas.hszero = gasgas.hs(0);
            else
                gasgas.time = 0;
            end
            gasgas.t0 = cputime; % my algorithm is not necessarily static!
            gasgas = gasgas.gas_finalize;
        end
        function gasgas = gas_finalize(gasgas)
            if isfield(gasgas.params, 'use_gpu')&&gasgas.params.use_gpu %%%not working/very slow
                gasgas.C = gpuArray(full(gasgas.C));
                gasgas.C_age = gpuArray(full(gasgas.C_age));
                gasgas.A = gpuArray(gasgas.A);
                gasgas.errorvector = gpuArray(gasgas.errorvector);
                gasgas.hizero = gpuArray(gasgas.hizero);
                gasgas.hszero= gpuArray(gasgas.hszero);
                gasgas.time= gpuArray(gasgas.time);
                gasgas.t0= gpuArray(gasgas.t0);
                gasgas.n= gpuArray(gasgas.n);
                gasgas.r= gpuArray(gasgas.r);
                gasgas.h= gpuArray(gasgas.h);
                gasgas.a= gpuArray(gasgas.a);
                gasgas.awk= gpuArray(gasgas.awk);
                gasgas.n1n2 = gpuArray(gasgas.n1n2);
                gasgas.params.en = gpuArray(gasgas.params.en);
                gasgas.params.eb = gpuArray(gasgas.params.eb);
                gasgas.params.at = gpuArray(gasgas.params.at);
                gasgas.params.ab = gpuArray(gasgas.params.ab);
                gasgas.params.an = gpuArray(gasgas.params.an);
                gasgas.params.tb = gpuArray(gasgas.params.tb);
                gasgas.params.tn = gpuArray(gasgas.params.tn);
                gasgas.params.h0 = gpuArray(gasgas.params.h0);
                gasgas.params.amax = gpuArray(gasgas.params.amax);
                gasgas.params.STATIC = gpuArray(gasgas.params.STATIC);
                gasgas.params.eb = gpuArray(gasgas.params.eb);
                gasgas.params.nodes = gpuArray(gasgas.params.nodes);                
                gasgas.gwr.s = gpuArray();
                gasgas.gwr.t = gpuArray();
                gasgas.gwr.wi = gpuArray();
                gasgas.gwr.wr = gpuArray();
                gasgas.gwr.ws = gpuArray();
                gasgas.gwr.num_of_neighbours = gpuArray();
                gasgas.gwr.neighbours = gpuArray();
            end
        end
        function gasgas = update_epochs(gasgas,epochs)
            if isfield(gasgas.params, 'accumulatedepochs')
                gasgas.params.accumulatedepochs = gasgas.params.accumulatedepochs + epochs;
            else
                gasgas.params.accumulatedepochs = epochs;
            end
        end
    end
end
