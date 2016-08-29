classdef DuctNetwork < handle
    % DuctNetwork implements all algorithms for balancing duct network
    
    properties
        n; % number of nodes without ground
        n_NodeDescription; % cell array of size n by 1, text description on each node
        P; % array of size n by 1, pressure vector on each node, Pa
        b; % number of branches in the network
        b_BranchDescription; % cell array of size b by 1, text description on each branch
        Q; % array of size b by 1, flow rate vector on each branch, m^3/s
        A; % Associate matrix of size n by b
        % A(i,j)==1 means branch j leaves node i
        % A(i,j)==0 means branch j is not connected with node i
        % A(i,j)==-1 means branch j enters node i
        t; % dimension of null space
        U; % null space basis of matrix A, of size b by t, AU=0, U'U=1
        b_Pdrop; %cell array of size b by 1 for the pressure drop in terms of q and s
        %if there exist multiple Pdrop functions, use a cell array to store each of them
        %only record the original Pdrop function for the fitting, handle the fitting direction in pressure drop calculation
        b_dPdQ; %cell array of size b by 1 for the partial pressure drop over q in terms of q and s,
        %only record the original dPdQ function for the fitting, handle the fitting direction in pressure drop calculation
        b_dPdS; %cell array of size b by 1 for the partial pressure drop over s in terms of q and s
        %only record the original dPdS function for the fitting, handle the fitting direction in pressure drop calculation
        b_Qidx; %cell array of size b by 1 for the index of dependent Q of Pdrop, dPdQ and dPdS functions, so use Q(abs(b_Qidx{b})).*sign(b_Qidx{b}) in these functions.
        % Qidx{ii}[jj] can be negative, the sign represents the fitting direction w.r.t. branch direction
        b_Sidx; %cell array of size b by 1 for the index of dependent S of Pdrop, dPdQ and dPdS functions, so use S(b_Sidx{b}) in these functions.
        s; % number of parameters in the model
        s_ParamDescription; % cell array of size b by 1, text description on each parameter
        S; % array of size s by 1, parameter vector for whole system
        s_m; % array of size s by 1, s_m(i)=0 means no need to identify S(i), s_m(i)>0 is the multiplicity of this parameter that need to be identified
        s_MultiS; % cell array of size s by 1, containing values of multiplicities of each parameter in row vectors, Use [s_MultiS{:}]' to obtain all in one column vector
        
        n_trail;
    end
    
    methods
        function obj = DuctNetwork()
            obj.n = 0;
            obj.n_NodeDescription = cell(0,1);
            obj.P = zeros(0,1);
            obj.b = 0;
            obj.b_BranchDescription = cell(0,1);
            obj.Q = zeros(0,1);
            obj.A = zeros(0,0);
            obj.t = 0;
            obj.U = zeros(0,0);
            obj.b_Pdrop = cell(0,1);
            obj.b_dPdQ = cell(0,1);
            obj.b_dPdS = cell(0,1);
            obj.b_Qidx = cell(0,1);
            obj.b_Sidx = cell(0,1);
            obj.s = 0;
            obj.s_ParamDescription = cell(0,1);
            obj.S = zeros(0,1);
            obj.s_m = zeros(0,1);
            obj.s_MultiS = cell(0,1);
            
            obj.n_trail = 0;
        end

        function Branch_Idx = AddBranch(obj,FromNode,ToNode,varargin)
            [~, b_Curr]=size(obj.A);
            Branch_Idx = b_Curr+1;
            SetNode(FromNode,1);
            SetNode(ToNode,-1);
            [obj.n,obj.b]=size(obj.A);
            obj.U = null(obj.A);
            obj.t = size(obj.U,2);
            obj.b_Pdrop{Branch_Idx,1}=cell(0,1);
            obj.b_dPdQ{Branch_Idx,1}=cell(0,1);
            obj.b_dPdS{Branch_Idx,1}=cell(0,1);
            obj.b_Qidx{Branch_Idx,1}=cell(0,1);
            obj.b_Sidx{Branch_Idx,1}=cell(0,1);
            obj.b_BranchDescription{Branch_Idx,1}=cell(0,1);
            if nargin>3
                obj.AddFitting(Branch_Idx,varargin{:});
            end
            function SetNode(Node,a)
                if ischar(Node), Node = find(cellfun(@(S)strcmp(S,Node), obj.n_NodeDescription),1); end;
                if isempty(Node), Node=0; end;
                if Node>0, obj.A(Node,Branch_Idx) = a; end;
            end
        end
        
        function Node=AddNode(obj,varargin)
            Node = obj.n+1;
            obj.n = Node;
            if nargin>1
                obj.n_NodeDescription{Node,1} = varargin{1};
            end
        end
        
        function AddFitting(obj,Branch,Model,Param) %Branch has direction, positve means fitting direction is same as the branch direction, negative means opposite
            if isa(Model,'function_handle')
                if Model('Is_Junction')
                    % separate the junction into three/four junction branches, and add each of them separately, parameters are shared by all of them.
                    CenterNode=find(arrayfun(@(i)isempty(setdiff(find(obj.A(i,:)),abs(Branch))),1:obj.n),1);
                    if isempty(CenterNode)
                        disp('no such junction exist, fail to create junction'); return
                    end    
                    Branch = -obj.A(CenterNode,abs(Branch)).*abs(Branch);
                    branch_function_handle = Model('Get_Branches');
                    branch_assignment = Model('Branch_Assignment');
                    param_assignment = Model('Parameter_Assignment');
                    ParamDescription = Model('Parameter_Description');
                    bShared = Model('Is_Shared_Parameter');
                    Param_Idx = zeros(size(ParamDescription));
                    for ii = 1:length(ParamDescription)
                        if bShared(ii) % the parameter is shared by all models
                            s_Idx = find(strcmp(obj.s_ParamDescription,ParamDescription{ii}),1);
                            if isempty(s_Idx) % shared parameter is first included in the model
                                obj.s = obj.s+1; %add the parameter at the end
                                obj.s_ParamDescription{obj.s,1} = ParamDescription{ii}; % assign parameter description
                                Param_Idx(ii) = obj.s;
                            else
                                Param_Idx(ii) = s_Idx;
                            end
                        else % the parameter is unique
                            obj.s = obj.s+1; %add the parameter at the end
                            obj.s_ParamDescription{obj.s} = ParamDescription{ii}; % assign parameter description
                            Param_Idx(ii) = obj.s;
                        end
                    end
                    for ii = 1:length(branch_function_handle)
                        BranchList = Branch(branch_assignment{ii});
                        PrimaryBranch = abs(BranchList(1));
                        obj.b_BranchDescription{PrimaryBranch} = [obj.b_BranchDescription{PrimaryBranch};{branch_function_handle{ii}('Model_Description')}];
                        obj.b_Pdrop{PrimaryBranch} = [obj.b_Pdrop{PrimaryBranch};{@(q,s)branch_function_handle{ii}({'Pdrop','dPdQ','dPdS'},q,s)}];
                        obj.b_dPdQ{PrimaryBranch} = [obj.b_dPdQ{PrimaryBranch};{@(q,s)branch_function_handle{ii}('dPdQ',q,s)}];
                        obj.b_dPdS{PrimaryBranch} = [obj.b_dPdS{PrimaryBranch};{@(q,s)branch_function_handle{ii}('dPdS',q,s)}];
                        obj.b_Qidx{PrimaryBranch} = [obj.b_Qidx{PrimaryBranch};{BranchList}];
                        obj.b_Sidx{PrimaryBranch} = [obj.b_Sidx{PrimaryBranch};{Param_Idx(param_assignment{ii})}];
                    end
                    obj.S(Param_Idx(1:length(Param))) = Param;
                    obj.s_m(Param_Idx) = Model('Is_Identified_Parameter');
                    obj.s_MultiS(Param_Idx(1:length(Param))) = arrayfun(@(a,idx)a*ones(1,obj.s_m(idx)),Param,Param_Idx(1:length(Param)),'UniformOutput',false);
                else
                    ParamDescription = Model('Parameter_Description');
                    bShared = Model('Is_Shared_Parameter');
                    Param_Idx = zeros(size(ParamDescription));
                    for ii = 1:length(bShared)
                        if bShared(ii) % the parameter is shared by all models
                            s_Idx = find(cellfun(@(S)strcmp(S,ParamDescription{ii}),obj.s_ParamDescription),1);
                            if isempty(s_Idx)
                                obj.s = obj.s+1; %add the parameter at the end
                                obj.s_ParamDescription{obj.s,1} = ParamDescription{ii}; % assign parameter description
                                Param_Idx(ii) = obj.s;
                            else
                                Param_Idx(ii) = s_Idx;
                            end
                        else % the parameter is unique
                            obj.s = obj.s+1; %add the parameter at the end
                            obj.s_ParamDescription{obj.s,1} = ParamDescription{ii}; % assign parameter description
                            Param_Idx(ii) = obj.s;
                        end
                    end
                    
                    BranchList = Branch;
                    PrimaryBranch = abs(Branch(1));
                    obj.b_BranchDescription{PrimaryBranch} = [obj.b_BranchDescription{PrimaryBranch};{Model('Model_Description')}];
                    obj.b_Pdrop{PrimaryBranch} = [obj.b_Pdrop{PrimaryBranch};{@(q,s)Model({'Pdrop','dPdQ','dPdS'},q,s)}];
                    obj.b_dPdQ{PrimaryBranch} = [obj.b_dPdQ{PrimaryBranch};{@(q,s)Model('dPdQ',q,s)}];
                    obj.b_dPdS{PrimaryBranch} = [obj.b_dPdS{PrimaryBranch};{@(q,s)Model('dPdS',q,s)}];
                    obj.b_Qidx{PrimaryBranch} = [obj.b_Qidx{PrimaryBranch};{BranchList}];
                    obj.b_Sidx{PrimaryBranch} = [obj.b_Sidx{PrimaryBranch};{Param_Idx}];
                    
                    obj.S(Param_Idx(1:length(Param))) = Param;
                    obj.s_m(Param_Idx) = Model('Is_Identified_Parameter');
                    obj.s_MultiS(Param_Idx(1:length(Param))) = arrayfun(@(a,idx)a*ones(1,obj.s_m(idx)),Param, Param_Idx(1:length(Param)),'UniformOutput',false);
                end
            elseif ischar(Model)
                Model = DuctNetwork.FittingDatabase(Model);
                obj.AddFitting(Branch,Model,Param);
            elseif iscell(Model)
                cellfun(@(a,b)obj.AddFitting(Branch,a,b),Model,Param);
            end
        end

        function Param_Idx = AddParameter(obj, Description, ParamValue, varargin)
            if iscell(Description)
                VAR=cellfun(@(a)num2cell(a),varargin,'UniformOutput',false);
                Param_Idx = cellfun(@(a,b,varargin)obj.AddParameter(a,b,varargin{:}),Description, num2cell(ParamValue),VAR{:});
            elseif isa(Description,'function_handle')
                Model = Description;
                Description = Model('Parameter_Description');
                Is_Identified_Parameter = Model('Is_Identified_Parameter');
                Param_Idx = AddParameter(obj, Description, ParamValue, Is_Identified_Parameter);
            elseif ischar(Description)
                obj.s = obj.s + 1;
                Param_Idx = obj.s;
                obj.s_ParamDescription{Param_Idx,1} = Description;
                obj.S(Param_Idx,1) = ParamValue;
                if nargin>=4
                    obj.s_m(Param_Idx,1) = varargin{1}; 
                    obj.s_MultiS{Param_Idx,1} = ParamValue*ones(1,varargin{1});
                else
                    obj.s_m(Param_Idx,1) = 0;
                    obj.s_MultiS{Param_Idx,1} = [];
                end
            end
        end
            
        function varargout = BranchPressureDrop(obj,Branch_Idx,Q,S)  %[dP, dPdX, dPdS]
            N = length(obj.b_BranchDescription{Branch_Idx});
            dP = zeros(N,1);dPdQ = zeros(N,obj.b);dPdX = zeros(N,obj.t);dPdS = zeros(N,obj.s);
            if nargout>=3
                for k = 1:N
                    dir = sign(obj.b_Qidx{Branch_Idx}{k});
                    idx_Q = abs(obj.b_Qidx{Branch_Idx}{k});
                    idx_S = obj.b_Sidx{Branch_Idx}{k};
                    [dP(k), dPdQ(k,idx_Q), dPdS(k,idx_S)] = obj.b_Pdrop{Branch_Idx}{k}(Q(idx_Q),S(idx_S));
                    dP(k) = dir(1)*dP(k);
                    dPdX(k,:) = (dir(1)*dPdQ(k,idx_Q).*dir)*obj.U(idx_Q,:);
                    dPdS(k,:) = dir(1)*dPdS(k,:);
                end
                varargout{1} = sum(dP,1);
                varargout{2} = sum(dPdX,1);
                varargout{3} = sum(dPdS,1);
            elseif nargout==2
                for k = 1:N
                    dir = sign(obj.b_Qidx{Branch_Idx}{k});
                    idx_Q = abs(obj.b_Qidx{Branch_Idx}{k});
                    idx_S = obj.b_Sidx{Branch_Idx}{k};
                    [dP(k), dPdQ(k,idx_Q)] = obj.b_Pdrop{Branch_Idx}{k}(Q(idx_Q),S(idx_S));
%                 Jac2 = [dPdQ(k,idx_Q), dPdS(k,idx_S)];
%                 Jac1 = jacobianest(@(x)obj.b_Pdrop{Branch_Idx}{k}(x(1:length(idx_Q)),x(length(idx_Q)+1:end)),[Q(idx_Q);S(idx_S)]);
%                 if any(abs(Jac2-Jac1)./Jac1>1e-8)
%                     disp(['Jacobian Estimation Error in ',obj.b_BranchDescription{Branch_Idx}{k}]);
%                     disp([find(abs(Jac2-Jac1)./Jac1>1e-6), max(abs(Jac2-Jac1)./Jac1), dir]);
%                 end
                    dP(k) = dir(1)*dP(k);
                    dPdX(k,:) = (dir(1)*dPdQ(k,idx_Q).*dir)*obj.U(idx_Q,:);
                end
                varargout{1} = sum(dP,1);
                varargout{2} = sum(dPdX,1);
            elseif nargout==1
                for k = 1:N
                    dir = sign(obj.b_Qidx{Branch_Idx}{k});
                    idx_Q = abs(obj.b_Qidx{Branch_Idx}{k});
                    idx_S = obj.b_Sidx{Branch_Idx}{k};
                    dP(k) = obj.b_Pdrop{Branch_Idx}{k}(Q(idx_Q),S(idx_S));
                    dP(k) = dir(1)*dP(k);
                end
                varargout{1} = sum(dP,1);
            end
        end
            
        function varargout = res_StateEquation(obj,X,S)
            % Input:
            % X is obj.t dimensional vector of internal state
            % S is obj.s dimensional vector of system parameter used in simulation
            % Output:
            % e is the residual of State Equations of size obj.t by 1
            % dedX is the Jacobian of e wrt internal state X of size obj.t by obj.t
            % dedS is the Jacobian of e wrt system parameter S of size obj.t by obj.s
            X = real(X);
            if nargout>=3
                [dP, dPdX, dPdS]=arrayfun(@(Branch_idx)obj.BranchPressureDrop(Branch_idx,obj.U*X,S),(1:obj.b)','UniformOutput',false);
                dP = cell2mat(dP); dPdX = cell2mat(dPdX); dPdS = cell2mat(dPdS);
                e = obj.U'*dP;
                dedX = obj.U'*dPdX;
                dedS = obj.U'*dPdS;
                varargout = cell(1,nargout);
                varargout{1}=e; varargout{2}=dedX; varargout{3}=dedS;
            elseif nargout==2
                [dP, dPdX]=arrayfun(@(Branch_idx)obj.BranchPressureDrop(Branch_idx,obj.U*X,S),(1:obj.b)','UniformOutput',false);
                dP = cell2mat(dP); dPdX = cell2mat(dPdX);
                e = obj.U'*dP; dedX = obj.U'*dPdX; varargout={e,dedX};
            elseif nargout==1
                dP=arrayfun(@(Branch_idx)obj.BranchPressureDrop(Branch_idx,obj.U*X,S),(1:obj.b)','UniformOutput',false);
                dP = cell2mat(dP);
                e = obj.U'*dP; varargout={e};
            end
        end
        
        function [X,Q,P] = Sim(obj,Param_Value, Param_Idx)
            if nargin==1
                S_value = obj.S;
            elseif nargin==2
                S_value = Param_Value;
            elseif nargin==3
                S_value = obj.S;
                if iscell(Param_Idx), Param_Idx = cellfun(@(Str) find(strcmp(Str,obj.s_ParamDescription)),Param_Idx); end
                S_value(Param_Idx) = Param_Value;
            end
            
            options = optimoptions(@lsqnonlin,'Display','none',...
                'Algorithm','trust-region-reflective',...
                'FunctionTolerance',1e-6,'StepTolerance',1e-6,...
                'MaxIterations',obj.t*5,...
                'SpecifyObjectiveGradient',true,'CheckGradients',false,...
                'FiniteDifferenceType','forward','FiniteDifferenceStepSize',1e-10);
            exitflag = -1;
            while exitflag<=0
                obj.n_trail = obj.n_trail+1;
                X0 = rand(obj.t,1);
                [X,~,~,exitflag] = lsqnonlin(@(x) obj.res_StateEquation(x,S_value),X0,[],[],options);
                %options = optimoptions(options,'SpecifyObjectiveGradient',false);
            end
            dP=arrayfun(@(Branch_idx)obj.BranchPressureDrop(Branch_idx,obj.U*X,S_value),(1:obj.b)','UniformOutput',false);
            dP=cell2mat(dP);
            Q = obj.U*X;
            P = (obj.A*obj.A')\obj.A*dP;
            obj.Q = Q;
            obj.P = P;
        end
        
        function varargout = Measurements(obj,ctrl_ParamIdx, proc_ctrl_ParamValue, ob_SensorType, ob_Uncertainty, proc_ob_SensorPosition)
            % Input:
            % ctrl_ParamIdx is parameterID vector of length ctrl. 
            % proc_in_ParamValue is parameter value matrix of size proc by ctrl.  
            % ob_SensorType is the type of each sensor, currently either pressure (p) or flow (q), cell vector of length ob 
            % ob_Uncertainty is the uncertainty of each sensor, vector of size ob, can be cell array containing function handles or numbers for Gaussian  
            % proc_ob_SensorPosition is sensor position in each step of experiment in matrix of size proc by ob 
            % Output
            % nargout should be the same as number of sensors, ob
            % varargout{ii} is columne vector of length proc recording the readings of each sensor. 
            
            [proc, ctrl] = size(proc_ctrl_ParamValue);
            ob = length(ob_SensorType);
            if isreal(ob_Uncertainty)
                ob_UncertaintyArray = ob_Uncertainty;
                ob_Uncertainty = cell(size(ob_Uncertainty));
                for ii=1:ob
                    ob_Uncertainty{ii} = @(x) x+randn()*ob_UncertaintyArray(ii);
                end
            end
            varargout  = cell(1,ob);
            for ii = 1:proc
                obj.Sim(proc_ctrl_ParamValue(ii,:),ctrl_ParamIdx);
                for jj = 1:ob
                    switch ob_SensorType
                        case {'P','Pressure'}
                            varargout{jj}(ii,1) = ob_Uncertainty{ii}(obj.P(proc_ob_SensorPosition(ii,jj)));
                        case {'q','Flowrate'}
                            varargout{jj}(ii,1) = ob_Uncertainty{ii}(obj.Q(proc_ob_SensorPosition(ii,jj)));
                    end
                end
            end
        end
    end
    
    methods (Static = true)
        function ModelFcn = FittingDatabase(ModelStr)
            ModelFcn=str2func(['DuctNetwork.',ModelStr]);
        end
        
        function varargout = FanQuadratic(query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {}; return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Model_Description'
                        varargout{ii}='Fan Using Quadratic Q-dP Relationship';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.FanQuadratic};
                    case 'Branch_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:2};
                    case 'Parameter_Description'
                        varargout{ii}={'Max Pressure(Pa)','Max Flow(m^3/s)'};
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[1,1];
                    case 'Pdrop'
                        q = varargin{1};s = varargin{2};
                        if q<0
                            varargout{ii} = -s(1);
                        else
                            varargout{ii} = s(1)*(q^2/s(2)^2-1);
                        end
                    case 'dPdQ'
                        q = varargin{1};s = varargin{2};
                        if q<0
                            varargout{ii} = 0;
                        else
                            varargout{ii} = s(1)*2*q/s(2)^2;
                        end
                    case 'dPdS'
                        q = varargin{1};s = varargin{2};
                        if q<0
                            varargout{ii} = [1,0];
                        else
                            varargout{ii} = [q^2/s(2)^2-1, -s(1)*2*q^2/s(2)^3];
                        end
                    otherwise
                        varargout{ii}=[];
                end
            end
        end
        
        function varargout = DuctQuadratic (query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {};return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Model_Description'
                        varargout{ii}='Duct Using Quadratic Q-dP Relationship';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.DuctQuadratic};
                    case 'Branch_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Description'
                        varargout{ii}={'Resistance(Pa/(m^3/s)^2)'};
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[1];
                    case 'Pdrop'
                        q = varargin{1};s = varargin{2};
                        varargout{ii} = s(1)*q*abs(q);
                    case 'dPdQ'
                        q = varargin{1};s = varargin{2};
                        varargout{ii} = s(1)*2*abs(q);
                    case 'dPdS'
                        q = varargin{1};s = varargin{2};
                        varargout{ii} = q*abs(q);
                    otherwise
                        varargout{ii}=[];
                end
            end
        end
        
        function varargout = CircularDarcyWeisbach (query, varargin)
            if ischar(query)
                query = {query};
            end
            varargout = cell(1,nargout);
            [varargout{:}] = DuctNetwork.CircularDarcyWeisbachHaaland (query, varargin{:});
        end
        
        function varargout = CircularDarcyWeisbachHaaland (query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {};return
            end
            varargout = cell(1,n);
            for ii = 1:nargout
                switch query{ii}
                    case 'Pdrop'
%                         if (~exist('dP','var')), DataInitiation(); end;
                        q = varargin{1};
                        s = varargin{2};
                        L = s(1);
                        D = s(2);
                        rho = s(3);
                        e = s(4);
                        nu = s(5);

                        Area = pi*(D/2)^2;
                        V = q/Area;
                        Re =abs(V)*D/nu;

                        lambda = 1/(1+exp(-(Re-3750)/250));
                        Cf_lam = 64/Re;
                        A = (e/3.7/D)^3.33;
                        B = (6.9/Re)^3;
                        T = log10(A+B);
                        Cf_turb = (-0.6*T)^(-2);
                        Cf = Cf_lam*(1-lambda)+lambda*Cf_turb;
                        %dP = Cf*L/D*rho/2*V^2;
                        T5 = L/D*rho*abs(V);
                        dP = Cf*T5*V/2;
                        varargout{ii}=dP;
                    case 'dPdQ'
%                         if (~exist('dP','var')), DataInitiation(); end;
                        dCf_lamdq = -16*nu*pi*D/q^2;
%                         if (~exist('dCf_turbdAB','var')) dCf_turbdAB = -1/0.18/T^3/log(10)/(A+B); end;
                        dCf_turbdAB = -1/0.18/T^3/log(10)/(A+B); 
                        dCf_turbdq = -dCf_turbdAB*3*B/abs(q);
                        dlambdadq = lambda*(1-lambda)*4/nu/pi/D;
                        dCfdq = dCf_lamdq*(1-lambda) + lambda*dCf_turbdq + (Cf_turb-Cf_lam)*dlambdadq;
                        dPdq = T5*(Cf/Area+abs(V)/2*dCfdq);
                        varargout{ii}=dPdq;
                    case 'dPdS'
%                         if (~exist('dP','var')), DataInitiation(); end;
%                         if (~exist('dCf_turbdAB','var')) dCf_turbdAB = -1/0.18/T^3/log(10)/(A+B); end;
                        dPdL = dP/L;
                        dCf_lamdD = Cf_lam/D;
                        dCf_turbdD = dCf_turbdAB*(-3.33*A+3*B)/D;
                        dlambdadD = -lambda*(1-lambda)*Re/D;
                        dCfdD = dCf_lamdD*(1-lambda) + lambda*dCf_turbdD + (Cf_turb-Cf_lam)*dlambdadD;
                        dPdD = -4*dP/D+dP/Cf*dCfdD;
                        dPdrho = dP/rho;
                        dCf_turbde = dCf_turbdAB*3.33*A/e;
                        dPde = dP/Cf*lambda*dCf_turbde;
                        dCfdnu = lambda*(1-lambda)/nu*(Cf_lam/lambda+dCf_turbdAB*3*B/(1-lambda)-(Cf_turb-Cf_lam)*Re);
                        dPdnu = dP/Cf*dCfdnu;
                        varargout{ii}=[dPdL, dPdD, dPdrho, dPde, dPdnu];
                    case 'Model_Description'
                        varargout{ii}='Circular Straight Duct Using Darcy Weisbach Equation by Haaland Approximation';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.CircularDarcyWeisbach};
                    case 'Branch_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:5};
                    case 'Parameter_Description'
                        varargout{ii}={'Length(m)','Diameter(m)','Density(kg/m^3)','Roughness(mm)','Dynamic Viscosity(m^2/s)'};
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,true,true,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[1,0,0,0,0];
                    otherwise
                        varargout{ii}=[];
                end
            end
        end
        
        function varargout = CircularDarcyWeisbachChurchill (query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {};return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = varargin{1};
                        s = varargin{2};
                        L = s(1);
                        D = s(2);
                        rho = s(3);
                        e = s(4);
                        nu = s(5);

                        Area = pi*(D/2)^2;
                        V = q/Area;
                        Re =abs(V)*D/nu;
                        T21 = power((7/Re),0.9);
                        T2 = T21 +(0.27*e/D);
                        T1 = -2.457*log(T2);
                        A = T1^16;
                        B = power((37530/Re),16);
                        T3 = power((8/Re),12);
                        T4 = 1/power(A+B,1.5);
                        Cf = 8*power(T3+T4,1/12);
                        T5 = L/D*rho*abs(V);
                        %dP = Cf*L/D*rho/2*V^2;
                        dP = Cf*T5/2*V;
                        varargout{ii}=dP;
                    case 'dPdQ'
                        dAdq = 35.3808*T1^15*T21/T2/abs(q);
                        dBdq = -16*B/abs(q);
                        dCfdq = Cf/(T3+T4)*(-T3/abs(q)-(dAdq+dBdq)/8/(A+B)^2.5);
                        dPdq = T5*(Cf/Area+abs(V)/2*dCfdq);
                        varargout{ii}=dPdq;
                    case 'dPdS'
                        dAdT2 = -2.457*16*A/T1/T2;
                        dAdD = dAdT2*(0.9*T21-0.27*e/D)/D;
                        dBdD = 16*B/D;
                        G = Cf/12/(T3+T4);
                        H1 = 12*T3;
                        H2 = 1.5/(A+B)^2.5;
                        dPdrho = dP/rho;
                        dPdL = dP/L;
                        dCfdD = G*(H1/D-H2*(dAdD+dBdD));
                        dPdD = dP*(dCfdD/Cf-5/D);
                        dAdnu = dAdT2*(0.9*T21/nu);
                        dBdnu = 16*B/nu;
                        dCfdnu = G*(H1/nu-H2*(dAdnu+dBdnu));
                        dPdnu = dP/Cf*dCfdnu;
                        dPde = 1.32678*dP*T1^15/T2/(T3+T4)/(A+B)^2.5/D;
                        dPds = [dPdL, dPdD, dPdrho, dPde, dPdnu];
                        varargout{ii}=dPds;
                    case 'Model_Description'
                        varargout{ii}='Circular Straight Duct Using Darcy Weisbach Equation by Churchill Approximation';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.CircularDarcyWeisbach};
                    case 'Branch_Assignment'
                        varargout{ii}={1};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:5};
                    case 'Parameter_Description'
                        varargout{ii}={'Length(m)','Diameter(m)','Density(kg/m^3)','Roughness(mm)','Dynamic Viscosity(m^2/s)'};
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,true,true,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[1,0,0,0,0];
                    otherwise
                        varargout{ii}=[];
                end
            end
            
        end
        
        function varargout = CircularTJunction( query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {}; return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        varargout{ii}=0;
                    case 'dPdQ'
                        varargout{ii}=zeros(3,1);
                    case 'dPdS'
                        varargout{ii}=zeros(4,1);
                    case 'Model_Description'
                        varargout{ii}='Circular T-Junction using ED5-3,ED5-4,SD5-18,SD5-9';
                    case 'Is_Junction'
                        varargout{ii}=true;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.CircularTJunction_Horizontal,@DuctNetwork.CircularTJunction_Horizontal,@DuctNetwork.CircularTJunction_Vertical};
                    case 'Branch_Assignment'
                        varargout{ii}={[1,2,3];[2,1,3];[3,1,2]};
                    case 'Parameter_Assignment'
                        varargout{ii}={[1,2,3,4];[2,1,3,4];[3,1,2,4]};
                    case 'Parameter_Description'
                        varargout{ii}={'Main 1 Diameter(m)','Main 2 Diameter(m)','Branch Diameter(m)','Density(kg/m^3)'};
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,false,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[0,0,0,0];
                end
            end
        end
        
        function varargout = CircularTJunction_Horizontal( query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {}; return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = reshape(varargin{1},1,[]);
                        s = reshape(varargin{2},1,[]);
                        [dP,dPdQ,dPdS] = DataInitiation(q,s);
                        varargout{ii} = dP;
                    case 'dPdQ'
                        varargout{ii} = dPdQ;
                    case 'dPdS'
                        varargout{ii} = dPdS;
                    case 'Model_Description'
                        varargout{ii}='Main of Circular T-Junction using ED5-3,ED5-4,SD5-18,SD5-9';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.CircularTJunction_Horizontal};
                    case 'Branch_Assignment'
                        varargout{ii}={[1,2,3]};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:4};
                    case 'Parameter_Description'
                        varargout{ii}={'Main 1 Diameter(m)','Main 2 Diameter(m)','Branch Diameter(m)','Density(kg/m^3)'};
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,false,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[0,0,0,0,0];
                end
            end
            function [dP,dPdQ,dPdS]=DataInitiation(q,s)
                dir = sign(q);
                switch (q>0)*[4;2;1]
                    case 1 %[0,0,1]*[4;2;1], - - +, bullhead diverge SD5-18, use Cb
                        [dP, dPdQ([1,2,3]), dPdS([1,2,3,4])] = DuctNetwork.Calc_SD5_18(abs(q([1,2,3])), s([1,2,3,4]),'b');
                    case 2 %[0,1,0]*[4;2;1], - + -, T diverge SD5-9 at downstream side, use Cs
                        [dP, dPdQ([1,3,2]), dPdS([1,3,2,4])] = DuctNetwork.Calc_SD5_9(abs(q([1,3,2])),s([1,3,2,4]),'s');
                    case 3 %[0,1,1]*[4;2;1], - + +, flow inward from the branch, T converge ED5-3 at downsteram side, no pressure drop
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                    case 4 %[1,0,0]*[4;2;1], + - -, flow outward from the branch, T diverge SD5-9 at upstream side, no pressure drop
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                    case 5 %[1,0,1]*[4;2;1], + - +, flow inward from the branch, T converge ED5-3 at upstream side, use Cs
                        [dP, dPdQ([1,3,2]), dPdS([1,3,2,4])] = DuctNetwork.Calc_ED5_3(abs(q([1,3,2])),s([1,3,2,4]),'s');
                    case 6 %[1,1,0]*[4;2;1], + + -, flow inward from the opposite main, bullhead converge ED5-4, use Cb
                        if s(1)>=s(2) % D1>D2, use Cb1
                            [dP, dPdQ([1,2,3]), dPdS([1,2,3,4])] = DuctNetwork.Calc_ED5_4(abs(q([1,2,3])), s([1,2,3,4]), 'b1');
                        else % D1<D2, use Cb2
                            [dP, dPdQ([2,1,3]), dPdS([2,1,3,4])] = DuctNetwork.Calc_ED5_4(abs(q([2,1,3])), s([2,1,3,4]), 'b2');
                        end
                    case 7 %[1,1,1]*[4;2;1], + + +, impossible
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                    case 0 %[0,0,0]*[4;2;1], - - -, impossible                        
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                end
                dP = dir(1)*dP;dPdQ = dir(1)*dPdQ.*dir;dPdS = dir(1)*dPdS;
            end
        end
        
        function varargout = CircularTJunction_Vertical( query, varargin)
            if ischar(query)
                n=1; query = {query};
            elseif iscell(query)
                n = length(query);
            else
                varargout = {}; return
            end
            varargout = cell(1,n);
            for ii = 1:n
                switch query{ii}
                    case 'Pdrop'
                        q = reshape(varargin{1},1,[]);
                        s = reshape(varargin{2},1,[]);
                        [dP,dPdQ,dPdS] = DataInitiation(q,s);
                        varargout{ii} = dP;
                    case 'dPdQ'
                        varargout{ii} = dPdQ;
                    case 'dPdS'
                        varargout{ii} = dPdS;
                    case 'Model_Description'
                        varargout{ii}='Branch of Circular T-Junction using ED5-3,ED5-4,SD5-18,SD5-9';
                    case 'Is_Junction'
                        varargout{ii}=false;
                    case 'Get_Branches'
                        varargout{ii}={@DuctNetwork.CircularTJunction_Vertical};
                    case 'Branch_Assignment'
                        varargout{ii}={[1,2,3]};
                    case 'Parameter_Assignment'
                        varargout{ii}={1:4};
                    case 'Parameter_Description'
                        varargout{ii}={'Branch Diameter(m)','Main 1 Diameter(m)','Main 2 Diameter(m)','Density(kg/m^3)'};
                    case 'Is_Shared_Parameter'
                        varargout{ii}=[false,false,false,true];
                    case 'Is_Identified_Parameter'
                        varargout{ii}=[0,0,0,0,0];
                end
            end
            function [dP,dPdQ,dPdS]=DataInitiation(q,s)
                dir = sign(q);
                switch (q>0)*[4;2;1]
                    case 1 %[0,0,1]*[4;2;1], - - +, T diverge SD5-9 at downstream side, use Cb
                        [dP, dPdQ([2,1,3]), dPdS([2,1,3,4])] = DuctNetwork.Calc_SD5_9(abs(q([2,1,3])),s([2,1,3,4]),'b');
                    case 2 %[0,1,0]*[4;2;1], - + -, T diverge SD5-9 at downstream side, use Cb
                        [dP, dPdQ([3,1,2]), dPdS([3,1,2,4])] = DuctNetwork.Calc_SD5_9(abs(q([3,1,2])),s([3,1,2,4]),'b');
                    case 3 %[0,1,1]*[4;2;1], - + +, bullhead converge ED5-4 at downsteram side, no pressure drop
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                    case 4 %[1,0,0]*[4;2;1], + - -, bullhead diverge SD5-18 at upstream side, no pressure drop
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                    case 5 %[1,0,1]*[4;2;1], + - +, T converge ED5-3 at branch side, use Cb
                        [dP, dPdQ([3,1,2]), dPdS([3,1,2,4])] = DuctNetwork.Calc_ED5_3(abs(q([3,1,2])),s([3,1,2,4]),'b');
                    case 6 %[1,1,0]*[4;2;1], + + -, flow inward from the main 1, T converge ED5-3 at branch side, use Cb
                        [dP, dPdQ([2,1,3]), dPdS([2,1,3,4])] = DuctNetwork.Calc_ED5_3(abs(q([2,1,3])),s([2,1,3,4]),'b');
                    case 7 %[1,1,1]*[4;2;1], + + +, impossible
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                    case 0 %[0,0,0]*[4;2;1], - - -, impossible                        
                        dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
                end
                dP = dir(1)*dP;dPdQ = dir(1)*dPdQ.*dir;dPdS = dir(1)*dPdS;
            end
        end
    end
    
    properties (Constant)
        ED5_3 = load('FittingData/ED5_3.mat');
        SD5_9 = load('FittingData/SD5_9.mat');
        ED5_4 = load('FittingData/ED5_4.mat');
        SD5_18 = load('FittingData/SD5_18.mat');
    end
    
    methods (Static = true)
        function [dP, dPdQ, dPdS]=Calc_ED5_3(q,s,Selection)
            gExp = @(q,s) 0.5*s(4)*(q(3)/(pi*s(3)^2/4))^2;
            dgdq = @(q,s) [0,0,s(4)*q(3)/(pi*s(3)^2/4)^2];
            dgds = @(q,s) [0,0,-2/s(3),1/s(4)]*gExp(q,s);
            switch Selection
                case 's'
                    if s(3)<=0.25 %Dc<= 0.25m
                        Cs_Table = DuctNetwork.ED5_3.Cs_part1;
                    else
                        Cs_Table = DuctNetwork.ED5_3.Cs_part2;
                    end
                    GridVec = {DuctNetwork.ED5_3.QsQc,DuctNetwork.ED5_3.AbAc,DuctNetwork.ED5_3.AsAc};
                    ZExp = @(q,s)[q(1)/q(3);(s(2)/s(3))^2;(s(1)/s(3))^2];
                    dZdq = @(q,s)[1/q(3),0,-q(1)/q(3)^2;0,0,0;0,0,0];
                    dZds = @(q,s)[0,0,0,0;0,2*s(2)/s(3)^2,-2*s(2)^2/s(3)^3,0;2*s(1)/s(3)^2,0,-2*s(1)^2/s(3)^3,0];
                    [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, Cs_Table, ZExp, dZdq, dZds, gExp, dgdq, dgds, q, s);
                case 'b'
                    if s(3)<=0.25 %Dc<= 0.25m
                        Cb_Table = DuctNetwork.ED5_3.Cb_part1;
                    else
                        Cb_Table = DuctNetwork.ED5_3.Cb_part2;
                    end
                    GridVec = {DuctNetwork.ED5_3.QbQc,DuctNetwork.ED5_3.AbAc,DuctNetwork.ED5_3.AsAc};
                    ZExp = @(q,s)[q(2)/q(3);(s(2)/s(3))^2;(s(1)/s(3))^2];
                    dZdq = @(q,s)[0,1/q(3),-q(2)/q(3)^2;0,0,0;0,0,0];
                    dZds = @(q,s)[0,0,0,0;0,2*s(2)/s(3)^2,-2*s(2)^2/s(3)^3,0;2*s(1)/s(3)^2,0,-2*s(1)^2/s(3)^3,0];
                    [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec, Cb_Table, ZExp, dZdq, dZds, gExp, dgdq, dgds, q, s);
                otherwise
                    dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
            end
        end
        
        function [dP, dPdQ, dPdS]=Calc_SD5_9(q,s,Selection)
            gExp = @(q,s) 0.5*s(4)*(q(3)/(pi*s(3)^2/4))^2;
            dgdq = @(q,s) [0,0,s(4)*q(3)/(pi*s(3)^2/4)^2];
            dgds = @(q,s) [0,0,-2/s(3),1/s(4)]*gExp(q,s);
            switch Selection
                case 'b'
                    GridVec = {DuctNetwork.SD5_9.QbQc,DuctNetwork.SD5_9.AbAc};
                    ZExp = @(q,s)[q(2)/q(3);(s(2)/s(3))^2];
                    dZdq = @(q,s)[0,1/q(3),-q(2)/q(3)^2;0,0,0];
                    dZds = @(q,s)[0,0,0,0;0,2*s(2)/s(3)^2,-2*s(2)^2/s(3)^3,0];
                    [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec,DuctNetwork.SD5_9.Cb, ZExp, dZdq, dZds, gExp, dgdq, dgds, q, s);
                case 's'
                    GridVec = {DuctNetwork.SD5_9.QsQc,DuctNetwork.SD5_9.AsAc};
                    ZExp = @(q,s)[q(1)/q(3);(s(1)/s(3))^2];
                    dZdq = @(q,s)[1/q(3),0,-q(1)/q(3)^2;0,0,0];
                    dZds = @(q,s)[0,0,0,0;2*s(1)/s(3)^2,0,-2*s(1)^2/s(3)^3,0];
                    [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec,DuctNetwork.SD5_9.Cs, ZExp, dZdq, dZds, gExp, dgdq, dgds, q, s);
                otherwise
                    dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
            end
        end
        
        function [dP, dPdQ, dPdS]=Calc_ED5_4(q,s,Selection)
            gExp = @(q,s) 0.5*s(4)*(q(3)/(pi*s(3)^2/4))^2;
            dgdq = @(q,s) [0,0,s(4)*q(3)/(pi*s(3)^2/4)^2];
            dgds = @(q,s) [0,0,-2/s(3),1/s(4)]*gExp(q,s);
            switch Selection
                case 'b1'
                    GridVec = {DuctNetwork.ED5_4.QbQc,DuctNetwork.ED5_4.AbAc,DuctNetwork.ED5_4.AbAc};
                    ZExp = @(q,s)[q(1)/q(3);(s(1)/s(3))^2;(s(2)/s(3))^2];
                    dZdq = @(q,s)[1/q(3),0,-q(1)/q(3)^2;0,0,0;0,0,0];
                    dZds = @(q,s)[0,0,0,0;2*s(1)/s(3)^2,0,-2*s(1)^2/s(3)^3,0;0,2*s(2)/s(3)^2,-2*s(2)^2/s(3)^3,0];
                    [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec,DuctNetwork.ED5_4.Cb1, ZExp, dZdq, dZds, gExp, dgdq, dgds, q, s);
                case 'b2'
                    GridVec = {DuctNetwork.ED5_4.QbQc,DuctNetwork.ED5_4.AbAc,DuctNetwork.ED5_4.AbAc};
                    ZExp = @(q,s)[q(2)/q(3);(s(1)/s(3))^2;(s(2)/s(3))^2];
                    dZdq = @(q,s)[0,1/q(3),-q(2)/q(3)^2;0,0,0;0,0,0];
                    dZds = @(q,s)[0,0,0,0;2*s(1)/s(3)^2,0,-2*s(1)^2/s(3)^3,0;0,2*s(2)/s(3)^2,-2*s(2)^2/s(3)^3,0];
                    [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec,DuctNetwork.ED5_4.Cb2, ZExp, dZdq, dZds, gExp, dgdq, dgds, q, s);
                otherwise
                    dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
            end
        end
        
        function [dP, dPdQ, dPdS]=Calc_SD5_18(q,s,Selection)
            gExp = @(q,s) 0.5*s(4)*(q(3)/(pi*s(3)^2/4))^2;
            dgdq = @(q,s) [0,0,s(4)*q(3)/(pi*s(3)^2/4)^2];
            dgds = @(q,s) [0,0,-2/s(3),1/s(4)]*gExp(q,s);
            switch Selection
                case 'b'
                    GridVec = {DuctNetwork.SD5_18.QbQc,DuctNetwork.SD5_18.AbAc};
                    ZExp = @(q,s)[q(1)/q(3);(s(1)/s(3))^2];
                    dZdq = @(q,s)[1/q(3),0,-q(1)/q(3)^2;0,0,0];
                    dZds = @(q,s)[0,0,0,0;2*s(1)/s(3)^2,0,-2*s(1)^2/s(3)^3,0];
                    [dP, dPdQ, dPdS] = DuctNetwork.Interp_Gradient(GridVec,DuctNetwork.SD5_18.Cb, ZExp, dZdq, dZds, gExp, dgdq, dgds, q, s);
                otherwise
                    dP = 0;dPdQ = [0,0,0];dPdS = [0,0,0,0];
            end
        end
        
        function [f,dfdq,dfds] = Interp_Gradient(GridVector, InterpTable, ZExp, dZdq, dZds, gExp, dgdq, dgds, q,s)
            % Z = ZExp(q,s)
            % C = Interpolation(GridVector,InterpTable,Z)
            % g = gExp(q,s)
            % f = g*C;
            % dfdq = g*dCdZ*dZdq + C*dgdq;
            % dfds = g*dCdZ*dZds + C*dgds;
            Coeff_Interpolant = griddedInterpolant(GridVector,InterpTable,'linear','nearest');
            Z = reshape(ZExp(q,s),1,[]);
            NZ = length(Z);
            Cf = Coeff_Interpolant(Z);
            StepSize = cell(1,NZ);IndexInTable = cell(1,NZ);
            for ii = 1:NZ
                if Z(ii)>GridVector{ii}(end)
                    IndexInTable{ii}=length(GridVector{ii})*[1;1];
                    StepSize{ii}=1;
                elseif Z(ii)<=GridVector{ii}(1)
                    IndexInTable{ii}=[1;1];
                    StepSize{ii}=1;
                else
                    tmp=find(Z(ii)>GridVector{ii},1,'last');
                    IndexInTable{ii}=[tmp;tmp+1];
                    StepSize{ii}=GridVector{ii}(tmp+1)-GridVector{ii}(tmp);
                end
            end
            cell_Grad = cell(1,NZ);
            [cell_Grad{:}] = gradient(InterpTable(IndexInTable{:}),StepSize{:});
            dCfdZ = cellfun(@(M)M(1),cell_Grad,'UniformOutput',true);
            g = gExp(q,s);
            dZdq = dZdq(q,s);
            dZds = dZds(q,s);
            dgdq = dgdq(q,s);
            dgds = dgds(q,s);
            f = Cf*g;
            dfdq = g*dCfdZ*dZdq + Cf*dgdq;
            dfds = g*dCfdZ*dZds + Cf*dgds;
        end
        
        function [f,dfdq,dfds] = Interp_Gradient2(GridVector, InterpTable, ZExp, dZdq, dZds, gExp, dgdq, dgds, q,s)
            % Z = ZExp(q,s)
            % C = Interpolation(GridVector,InterpTable,Z)
            % g = gExp(q,s)
            % f = g*C;
            % dfdq = g*dCdZ*dZdq + C*dgdq;
            % dfds = g*dCdZ*dZds + C*dgds;
            Z = reshape(ZExp(q,s),1,[]);
            NZ = length(Z);
            StepSize = cell(1,NZ);IndexInTable = cell(1,NZ);ZMesh = cell(1,NZ);IndexInMesh = cell(1,NZ);
            for ii = 1:NZ
                if Z(ii)>GridVector{ii}(end)
                    IndexInTable{ii}=length(GridVector{ii})*[1;1];
                    StepSize{ii} = 2*(Z(ii)-GridVector{ii}(end));
                    ZMesh{ii}=[GridVector{ii}(end);Z(ii);2*Z(ii)-GridVector{ii}(end)];
                elseif Z(ii)<=GridVector{ii}(1)
                    IndexInTable{ii}=[1;1];
                    StepSize{ii} = -2*(Z(ii)-GridVector{ii}(1));
                    ZMesh{ii}=[2*Z(ii)-GridVector{ii}(1);Z(ii);GridVector{ii}(1)];
                else
                    IndexInTable{ii}=find(Z(ii)<=GridVector{ii},1)*[1;1]+[-1;0];
                    StepSize{ii} = range(GridVector{ii}(IndexInTable{ii}));
                    ZMesh{ii}=[GridVector{ii}(IndexInTable{ii}(1));GridVector{ii}(IndexInTable{ii}(end))];
                end
                IndexInMesh{ii}=[1;2];
            end
            CfMesh = InterpTable(IndexInTable{:});
            lambda = cellfun(@(v,x)(x-v(1))/(v(2)-v(1)),ZMesh,num2cell(Z));
            for ii=1:NZ
                Index_1 = IndexInMesh; Index_1{ii}=1;
                Index_2 = IndexInMesh; Index_2{ii}=2;
                Index_3 = IndexInMesh; Index_3{ii}=3;
                CfMesh(Index_3{:}) = CfMesh(Index_2{:});
                CfMesh(Index_2{:}) = CfMesh(Index_1{:})*lambda(ii)+CfMesh(Index_3{:})*(1-lambda(ii));
                IndexInMesh{ii}=[1;2;3];
            end
            Index_mid = num2cell(2*ones(1,NZ));
            Cf=CfMesh(Index_mid{ii});
            cell_Grad = cell(1,NZ);
            [cell_Grad{:}] = gradient(CfMesh,StepSize{:});
            dCfdZ = cellfun(@(M)M(Index_mid{:}),cell_Grad);
            g = gExp(q,s);
            dZdq = dZdq(q,s);
            dZds = dZds(q,s);
            dgdq = dgdq(q,s);
            dgds = dgds(q,s);
            f = Cf*g;
            dfdq = g*dCfdZ*dZdq + Cf*dgdq;
            dfds = g*dCfdZ*dZds + Cf*dgds;
        end
    end
end

