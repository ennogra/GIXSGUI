% ***********************************************
% Copyright (c) 2020 UChicago Argonne, LLC
% See LICENSE file.
% ***********************************************
%
% LINEFIT Tool to fit 1D data.
%
%   OBJ = LINEFIT(DATA) loads 2-column data and create a LINEFIT object
%   with default settings.

%   Zhang Jiang
%   $Revision: 0.1$  $Date: 2016/10/06 $
%   $Revision: 0.2$  $Date: 2016/10/25 $ Correct analytical properties.

classdef linefit
    properties  % public
        % (Public) Input data. Must be a 2-column data [XData,YData]
        Data = [];
        % (Public) Index of curve models. -2/-1/[1]/2/3/4: Custom/Mult-Model/[Lorentzian]/...
        CurveModelIndex = 1;
        % (Public) Index of background models. -2/-1/[0]/1/2/3/...: Custom/Power/[No]/Constant/Linear/Quadratic/...
        BkgdModelIndex = 0;
        % (Public) Custom curve model
        CustomCurveModel    = struct;
        % (Public) Custom background model
        CustomBkgdModel     = struct;
        % (Public) Complete model for fitting
        Model = struct;
        % (Public) Flag to force a non-negative background. 0/[1]: no/[yes]
        FlagNonNegBkgd = 1;
        % (Public) Flag to fit in log scale for YData. [0]/1: [no]/yes
        FlagLogScale = 0;        
        % (Public) Optimzation opitions for slover lsqcurvefit
        FitOptions = optimoptions('lsqcurvefit',...
            'Display','iter-detailed',...
            'TolX',1e-9,...
            'TolFun',1e-9,...
            'MaxFunEvals',10000,...
            'MaxIter',10000);
        % (Public) Store fitting result
        FitOutput = struct;        
    end % End of public properties
    properties (SetAccess = private)
        % (Private) List of built-in modles
        ModelList = struct;
        % (Private) Peak value and positions ([value, position]) for the initialization for peak-shaped models
        PeaksFound = [];
    end % End of properties with SetAccess = private
    properties (Dependent = true, SetAccess = private)
        % (Private) Base models to be replicated for multiple curves
        ModelBase = struct;                
        % (Private) Number of multiple curves.
        NOfModelCurves = 1;
    end % End of properties with Dependent = true, SetAccess = private
    properties (Hidden = true)
        % (Hidden) Set to pass parameters and fitting flags (internal use only)
        FitParamsSet = struct;        
    end        
    methods % Dynamic methods
        % --- interface constructor
        function obj = linefit(varargin)
            if nargin == 1
                obj.Data = varargin{1};
            end
            % ==============================
            % --- initialize model list
            % ==============================            
            % --- symmetric without support                        
            obj.ModelList    = obj.createmodel(@linefit.lorentzian,3,...
                [1,0,0.01],[0,-Inf,0],[Inf,Inf,Inf],[1,1,1],...
                'Cauchy/Lorentzian','lorentzian(x,[amp,x0,gamma])');
            obj.ModelList(2)    = obj.createmodel(@linefit.ilorentzian,3,...
                [1,0,0.01],[0,-Inf,0],[Inf,Inf,Inf],[1,1,1],...
                'Intermeidate Lorentzian','ilorentzian(x,[amp,x0,gamma])');
            obj.ModelList(3)    = obj.createmodel(@linefit.mlorentzian,3,...
                [1,0,0.01],[0,-Inf,0],[Inf,Inf,Inf],[1,1,1],...
                'Modified Lorentzian','mlorentzian(x,[amp,x0,gamma])');            
            obj.ModelList(4)    = obj.createmodel(@linefit.gaussian,3,...
                [1,0,0.01],[0,-Inf,0],[Inf,Inf,Inf],[1,1,1],...
                'Gaussian','gaussian(x,[amp,x0,sigma])');
            obj.ModelList(5)    = obj.createmodel(@linefit.generalizednorm1,4,...
                [1,0,0.01,2],[0,-Inf,0,0],[Inf,Inf,Inf,Inf],[1,1,1,1],...
                'Generalized Normal 1','generalizednorm1(x,[amp,a,b,c])');            
            obj.ModelList(6)    = obj.createmodel(@linefit.laplace,3,...
                [1,0,0.01],[0,-Inf,0],[Inf,Inf,Inf],[1,1,1],...
                'Laplace','laplace(x,[amp,a,b])');
            obj.ModelList(7)    = obj.createmodel(@linefit.logistic,3,...
                [1,0,0.01],[0,-Inf,0],[Inf,Inf,Inf],[1,1,1],...
                'Logistic','logistic(x,[amp,a,b])');            
            obj.ModelList(8)    = obj.createmodel(@linefit.pearsonVII,4,...
                [1,0,0.01,2],[0,-Inf,0,0],[Inf,Inf,Inf,Inf],[1,1,1,1],...
                'PearsonVII','pearsonVII(x,[amp,a,b,c])');
            obj.ModelList(9)    = obj.createmodel(@linefit.pseudovoigtTCH,4,...
                [1,0,0.01,0.01],[0,-Inf,0,0],[Inf,Inf,Inf,Inf],[1,1,1,1],...
                'Pseudo-Voigt (TCH)','pseudovoigtTCH(x,[amp,x0,sigma,gamma])');
            obj.ModelList(10)    = obj.createmodel(@linefit.pseudovoigtIAT,4,...
                [1,0,0.01,0.01],[0,-Inf,0,0],[Inf,Inf,Inf,Inf],[1,1,1,1],...
                'Pseudo-Voigt (IAT)','pseudovoigtIAT(x,[amp,x0,sigma,gamma])');
            obj.ModelList(11)    = obj.createmodel(@linefit.pseudovoigtLLHGD,4,...
                [1,0,0.01,0.01],[0,-Inf,0,0],[Inf,Inf,Inf,Inf],[1,1,1,1],...
                'Pseudo-Voigt (LLHGD)','pseudovoigtLLHGD(x,[amp,x0,sigma,gamma])');            
            obj.ModelList(12)    = obj.createmodel(@linefit.uniform,3,...
                [1,0,0.01],[0,-Inf,-Inf],[Inf,Inf,Inf],[1,1,1],...
                'Uniform','uniform(x,[amp,a,b])');            
            obj.ModelList(13)    = obj.createmodel(@linefit.voigt,4,...
                [1,0,0.01,0.01],[0,-Inf,0,0],[Inf,Inf,Inf,Inf],[1,1,1,1],...
                'Voigt','voigt(x,[amp,x0,sigma,gamma]');
            % --- skewe without support            
            obj.ModelList(14)    = obj.createmodel(@linefit.skewlaplace,4,...
                [1,0,0.01,2],[0,-Inf,0,0],[Inf,Inf,Inf,Inf],[1,1,1,1],...
                'Skew-Laplace','skewlaplace(x,[amp,a,b,c])');
            obj.ModelList(15)    = obj.createmodel(@linefit.skewlogistic,4,...
                [1,0,1,2],[0,-Inf,0,0],[Inf,Inf,Inf,Inf],[1,1,1,1],...
                'Skew-Logistic','skewlogistic(x,[amp,a,b,c])');
            obj.ModelList(16)    = obj.createmodel(@linefit.skewnorm,4,...
                [1,0,0.01,0],[0,-Inf,0,-Inf],[Inf,Inf,Inf,Inf],[1,1,1,1],...
                'Skew-Normal','skewnorm(x,[amp,a,b,c])');            
            obj.ModelList(17)    = obj.createmodel(@linefit.skewpvSB,5,...
                [1,0,0.01,0.5,0],[0,-Inf,0,0,-Inf],[Inf,Inf,Inf,1,Inf],[1,1,1,1,1],...
                'Skew-Pseudo-Voigt (SB)','skewpvSB(x,[amp,x0,fwhm0,eta,a])');
            obj.ModelList(18)    = obj.createmodel(@linefit.skewpvSSG,6,...
                [1,0,0.01,0.5,0,0],[0,-Inf,0,0,-Inf,-Inf],[Inf,Inf,Inf,1,Inf,Inf],[1,1,1,1,1,1],...
                'Skew-Pseudo-Voigt (SSG)','skewpvSSG(x,[amp,x0,fwhm0,eta,a,b])');
            % --- skew with support
            obj.ModelList(19)    = obj.createmodel(@linefit.burrXII,5,...
                [1,0,1,2,1],[0,-Inf,0,0,0],[Inf,Inf,Inf,Inf,Inf],[1,0,1,1,1],...
                'Burr XII','burrXII(x,[amp,x0,a,b,c])');                
            obj.ModelList(20)    = obj.createmodel(@linefit.expdist,3,...
                [1,0,1],[0,-Inf,0],[Inf,Inf,Inf],[1,1,1],...
                'Exponential','expdist(x,[amp,a,b])');
            obj.ModelList(21)    = obj.createmodel(@linefit.gammadist,4,...
                [1,0,1,1],[0,-Inf,0,0],[Inf,Inf,Inf,Inf],[1,1,1,1],...
                'Gamma','gammadist(x,[amp,a,b,c])');            
            obj.ModelList(22)    = obj.createmodel(@linefit.inversegamma,4,...
                [1,0,1,1],[0,-Inf,0,0],[Inf,Inf,Inf,Inf],[1,0,1,1],...
                'Inverse Gamma','inversegamma(x,[amp,x0,a,b])');                        
            obj.ModelList(23)    = obj.createmodel(@linefit.inversenorm,4,...
                [1,0,1,1],[0,-Inf,0,0],[Inf,Inf,Inf,Inf],[1,0,1,1],...
                'Inverse Normal','inversenorm(x,[amp,x0,a,b])');
            obj.ModelList(24)    = obj.createmodel(@linefit.levy,3,...
                [1,0,1],[0,-Inf,0],[Inf,Inf,Inf],[1,1,1],...
                'Levy','levy(x,[amp,a,b])');
            obj.ModelList(25)    = obj.createmodel(@linefit.logcauchy,4,...
                [1,0,0,1],[0,-Inf,-Inf,0],[Inf,Inf,Inf,Inf],[1,0,1,1],...
                'Log-Cauchy','logcauchy(x,[amp,x0,a,b])');            
            obj.ModelList(26)    = obj.createmodel(@linefit.lognorm,4,...
                [1,0,0,1],[0,-Inf,-Inf,0],[Inf,Inf,Inf,Inf],[1,0,1,1],...
                'Log-Normal','lognorm(x,[amp,x0,a,b])');
            obj.ModelList(27)    = obj.createmodel(@linefit.paretoI,3,...
                [1,1,1],[0,0,0],[Inf,Inf,Inf],[1,1,1],...
                'Pareto Type I','paretoI(x,[amp,a,b])');
            obj.ModelList(28)    = obj.createmodel(@linefit.weibull,4,...
                [1,1,1,1],[0,-Inf,0,0],[Inf,Inf,Inf,Inf],[1,1,1,1],...
                'Weibull','weibull(x,[amp,a,b,c])');  
            obj.ModelList(29)    = obj.createmodel(@linefit.atan,3,...
                [1,0,1],[0,-Inf,-Inf],[Inf,Inf,Inf],[1,1,1],...
                'Atan','atan(x,[amp,a,b])');  
            obj.ModelList(30)    = obj.createmodel(@linefit.erf,3,...
                [1,0,1],[0,-Inf,-Inf],[Inf,Inf,Inf],[1,1,1],...
                'Error','erf(x,[amp,a,b])');  
            obj.ModelList(31)    = obj.createmodel(@linefit.powerlaw,3,...
                [1,0,1],[0,-Inf,-Inf],[Inf,Inf,Inf],[1,0,1],...
                'Power Law','powerlaw(x,[amp,a,b])');              
            obj.Model = obj.ModelBase;  % pass Model base to Model
        end
        %         % --- Copy object (No need for non-handle class)
        %         function new = copy(this)
        %             new = feval(class(this));
        %             m = metaclass(this);
        %             %mp = findobj([m.Properties{:}],'SetAccess','public');
        %             mp = findobj([m.Properties{:}],'Dependent',0);
        %             p = {mp.Name};
        %             for ii=1:length(p)
        %                 new.(p{ii}) = this.(p{ii});
        %             end
        %         end
        function obj = set.Data(obj,data)
            dim_data = size(data);
            if ~isnumeric(data) || ~ismatrix(data) || (dim_data(2)~=2 && dim_data(2)~=3)
                obj.Data = [];
                return;
            end
            obj.Data = data;
            obj.PeaksFound = []; %#ok<MCSUP>
            obj.FitOutput = struct;               %#ok<MCSUP>
        end
        function obj = set.CurveModelIndex(obj,modelindex)
            if ~isscalar(modelindex) || ~isnumeric(modelindex) || ~nnz(modelindex == [-2,-1,1:length(obj.ModelList)]) %#ok<MCSUP>
                error('Invalid CurveModelIndex.')
            end
            n = obj.NOfModelCurves;     %#ok<MCSUP> % get number of curves
%             old_index = obj.CurveModelIndex;
%             if old_index ~= modelindex
%                 obj.CurveModelIndex = modelindex;
%                 obj.Model = obj.ModelBase; %#ok<MCSUP>
%             end
            obj.CurveModelIndex = modelindex;
            obj.Model = obj.ModelBase; %#ok<MCSUP>
            obj = addcurves(obj,n-1);
        end
        function obj = set.BkgdModelIndex(obj,modelindex)
            if ~isscalar(modelindex) || ~isnumeric(modelindex) || floor(modelindex) ~= modelindex || modelindex<-2 || obj.BkgdModelIndex == modelindex
                return;
            end
            obj.BkgdModelIndex = modelindex;
            tmp_obj = obj.ModelBase; %#ok<MCSUP>
            obj.Model.BkgdModel = tmp_obj.BkgdModel; %#ok<MCSUP>
            obj.FitOutput = struct; %#ok<MCSUP>
            %obj.Model = obj.ModelBase; %#ok<MCSUP>
        end
        function obj = set.CustomCurveModel(obj,model)
            if ~isstruct(model) || ~isfield(model,'ModelFcnHandle')
                error('Invalid CustomCurveModel. Please use createmodel method to create a curve model.');
            end
            obj.CustomCurveModel = model;
            if ismember(obj.CurveModelIndex,[-2,-1]) %#ok<MCSUP>
                n = obj.NOfModelCurves;  %#ok<MCSUP> % get model number of curves                
                obj.Model = obj.ModelBase; %#ok<MCSUP>
                obj = addcurves(obj,n-1);
            end
        end
        function obj = set.CustomBkgdModel(obj,model)
            if ~isstruct(model) || ~isfield(model,'ModelFcnHandle')
                error('Invalid CustomBkgdModel. Please use createmodel method to create a background model.');
            end
            obj.CustomBkgdModel = model;
            if obj.BkgdModelIndex == -2 %#ok<MCSUP>
                obj.Model = obj.ModelBase; %#ok<MCSUP>
            end
        end
        function obj = set.Model(obj,model)
            obj.Model = model;
            obj.FitOutput = struct; %#ok<MCSUP>
        end
        function s = get.NOfModelCurves(obj)
            s = obj.Model.CurveModel.NOfParams(2);
        end
        function s = get.ModelBase(obj)
            if ismember(obj.CurveModelIndex,[-2,-1])
                s1 = obj.CustomCurveModel;
                if ~isstruct(s1) || ~isfield(s1,'ModelFcnHandle') || ~isa(s1.ModelFcnHandle,'function_handle')
                    error('Invalid CustomCurveModel.')
                end
            else
                s1 = obj.ModelList(obj.CurveModelIndex);
            end
            s1.NOfParams = [s1.NOfParams,1];
            switch obj.BkgdModelIndex
                case -2
                    s2 = obj.CustomBkgdModel;
                    if ~isstruct(s2) || ~isfield(s2,'ModelFcnHandle')
                        error('Invalid CustomBkgdModel.')
                    end
                case -1
                    s2 = obj.createmodel(...
                        @(x,p)(p(1)+p(2)*x.^p(3)),3,...
                        [0,0,-1],[-Inf,-Inf,-Inf],[Inf,Inf,Inf],[1,1,1],...
                        'Power');
                case 0
                    s2 = obj.createmodel(...
                        @(x,p)(zeros(size(x))),1,0,0,0,0,'None');
                otherwise
                    s2 = obj.createmodel(...
                        @(x,p)(polyval(p(:),x)),obj.BkgdModelIndex,...
                        zeros(obj.BkgdModelIndex,1),-Inf*ones(obj.BkgdModelIndex,1),...
                        Inf(obj.BkgdModelIndex,1),ones(obj.BkgdModelIndex,1),...
                        ['Polynomial of degree ',num2str(obj.BkgdModelIndex-1),'']);
            end
            s.CurveModel = s1;
            s.BkgdModel  = s2;
        end
        function s = get.FitParamsSet(obj)
            fitparam = [...
                [obj.Model.CurveModel.StartParams(:), ...
                obj.Model.CurveModel.LowerBounds(:), ...
                obj.Model.CurveModel.UpperBounds(:),...
                obj.Model.CurveModel.FitFlags(:)]; ...
                [obj.Model.BkgdModel.StartParams(:), ...
                obj.Model.BkgdModel.LowerBounds(:), ...
                obj.Model.BkgdModel.UpperBounds(:),...
                obj.Model.BkgdModel.FitFlags(:)]
                ];
            x       = fitparam(:,1);
            lb      = fitparam(:,2);
            ub      = fitparam(:,3);
            fitflag = fitparam(:,4);
            x1 = x(fitflag==1);
            x2 = x(fitflag==0);
            lb1 = lb(fitflag==1);
            ub1 = ub(fitflag==1);
            s.fitparam  = fitparam;
            s.x         = x;
            s.x1        = x1;
            s.x2        = x2;
            s.lb1       = lb1;
            s.ub1       = ub1;
            s.fitflag   = fitflag;
        end
        function obj = addcurves(obj,varargin)
        % ADDCURVES Add curves of identical ModelBase.CurveModel to Model.CurveModel
        %   OBJ2=ADDCURVES(OBJ1) or OBJ2=OBJ1.ADDCURVES Add one curve to
        %   the end. 
        %
        %   OBJ2=ADDCURVES(OBJ1,N) or OBJ2=OBJ1.ADDCURVES(N) Add N curves
        %   to the end. 
        %
        %   OBJ2=ADDCURVES(OBJ1,N,POS) or OBJ2=OBJ1.ADDCURVES(N,POS) Add N
        %   curves to position POS. POS must be a non-zero integer. The
        %   value of POS (i.e. |POS|) is the position; the sign of POS is
        %   the direction (i.e. + means insertion after |POS| and - means
        %   insertion before |POS|). 
        %
        %   To add more curves to current object OBJ1 without creating a
        %   new object, set OBJ2 to OBJ1.
            old_n = obj.NOfModelCurves;
            switch nargin
                case 1
                    n   = 1;
                    pos = old_n;
                case 2
                    n   = varargin{1};
                    pos = old_n;
                case 3
                    n   = varargin{1};
                    pos = varargin{2};
                otherwise
                    error('Invalid input argument.');
            end
            if ~isscalar(n) || ~isnumeric(n) || n<0 || floor(n)~=n
                error('Invalid additional number of curves.');
            end
            if ~isscalar(pos) || ~isnumeric(pos) || pos == 0 || floor(pos)~=pos || abs(pos)>old_n
                error('Invalid insert position.');
            end
            obj.Model.CurveModel.NOfParams(2) = old_n+n;
            if pos < 0  % covert "insertion before" to "insertion after"
                pos = abs(pos)-1;
            end
            obj.Model.CurveModel.StartParams    = [...
                obj.Model.CurveModel.StartParams(:,1:pos),...
                repmat(obj.ModelBase.CurveModel.StartParams,1,n),...
                obj.Model.CurveModel.StartParams(:,pos+1:end)];
            obj.Model.CurveModel.LowerBounds    = [...
                obj.Model.CurveModel.LowerBounds(:,1:pos),...
                repmat(obj.ModelBase.CurveModel.LowerBounds,1,n),...
                obj.Model.CurveModel.LowerBounds(:,pos+1:end)];
            obj.Model.CurveModel.UpperBounds    = [...
                obj.Model.CurveModel.UpperBounds(:,1:pos),...
                repmat(obj.ModelBase.CurveModel.UpperBounds,1,n),...
                obj.Model.CurveModel.UpperBounds(:,pos+1:end)];
            obj.Model.CurveModel.FitFlags       = [...
                obj.Model.CurveModel.FitFlags(:,1:pos),...
                repmat(obj.ModelBase.CurveModel.FitFlags,1,n),...
                obj.Model.CurveModel.FitFlags(:,pos+1:end)];
            obj.Model.CurveModel.FitStdErrors       = [...
                obj.Model.CurveModel.FitStdErrors(:,1:pos),...
                repmat(obj.ModelBase.CurveModel.FitStdErrors,1,n),...
                obj.Model.CurveModel.FitStdErrors(:,pos+1:end)];
        end
        function obj = reducecurves(obj,varargin)
        % REDUCECURVES Reduce curves from Model.CurveModel
        %   OBJ2=REDUCECURVES(OBJ1) or OBJ2 = OBJ1.REDUCECURVES Reduce the
        %   last one curve at the end. 
        %
        %   OBJ2=REDUCECURVES(OBJ1,INDEX) or OBJ2=OBJ1.REDUCECURVES(INDEX)
        %   Reduce curves with postion INDEX. INDEX must be a scalar or a
        %   vector, e.g., [2,5:7] removes curve #2, #5, #6, and #7.
            old_n = obj.NOfModelCurves;
            switch nargin
                case 1
                    ind = old_n;
                case 2
                    ind = varargin{1};
                otherwise
                    error('Invalid input argument.');
            end
            if ~isvector(ind) || ~isequal(floor(ind),ind) || numel(ind)>=old_n || min(ind)<1 || max(ind)>old_n
                error('Invalid curve indices to be reduced.');
            end
            ind = ind(:);
            new_n = old_n-length(ind);
            obj.Model.CurveModel.NOfParams(2) = new_n;
            obj.Model.CurveModel.StartParams(:,ind)  = [];
            obj.Model.CurveModel.LowerBounds(:,ind)  = [];
            obj.Model.CurveModel.UpperBounds(:,ind)  = [];
            obj.Model.CurveModel.FitFlags(:,ind)     = [];
            obj.Model.CurveModel.FitStdErrors(:,ind) = [];
        end
        function obj = sortcurves(obj,row,varargin)
        % SORTCURVES Sort the order of curves in Model.CurveModel
        %   OBJ2=SORTCURVES(OBJ1,ROW) or OBJ2=OBJ1.SORTCURVES(ROW) Sort the
        %   curves in the default MODE ('ascend') manner according to the 
        %   values in the row number ROW.
        %
        %   OBJ2=SORTCURVES(OBJ1,ROW,MODE) or OBJ2=OBJ1.SORTCURVES(ROW,MODE)
        %   MODE = 'ascend' or 'descent'. 
            if nargin == 1
                error('Invalid input argument.');
            end
            if ~isscalar(row) || ~isnumeric(row) || floor(row) ~= row || row<1 || row>obj.Model.CurveModel.NOfParams(1)
                error('Invalid row number.');
            end
            if nargin == 2
                mode = 'ascend';
            elseif nargin == 3 && ischar(varargin{1}) && ...
                    (strcmpi(varargin{1},'ascend') || strcmpi(varargin{1},'descend'))
                mode = lower(varargin{1});
            else
                error('Invalid sorting mode.');
            end
            startparams = (obj.Model.CurveModel.StartParams)';
            [~,ind] = sortrows(startparams,row);
            if strcmpi(mode,'descend')
                ind = flipud(ind(:));
            end
            obj.Model.CurveModel.StartParams     = obj.Model.CurveModel.StartParams(:,ind);
            obj.Model.CurveModel.LowerBounds     = obj.Model.CurveModel.LowerBounds(:,ind);
            obj.Model.CurveModel.UpperBounds     = obj.Model.CurveModel.UpperBounds(:,ind);
            obj.Model.CurveModel.FitFlags        = obj.Model.CurveModel.FitFlags(:,ind);
            obj.Model.CurveModel.FitStdErrors    = obj.Model.CurveModel.FitStdErrors(:,ind);
        end
        function obj = applypeaks(obj,varargin)
        % APPLYPEAKS Take peak positions in PeaksFound and pass to the Model.CurveModel.StartParams, and create or reduce number of curves to the number of found peaks.
        %   OBJ2=APPLYPEAKS(OBJ1) or OBJ2=OBJ1.ACCEPTPEAKS Pass to the 1st and 2nd
        %   input parameter (peak positions) for curve models #1-18
        %   (peak-shaped models).
        %   
        %   OBJ2=APPLYPEAKS(OBJ1,ROWS) or OBJ2=OBJ1.ACCEPTPEAKS(ROWS) Pass
        %   to ROWS of Model.CurveModel.StartParams. For example,
        %   ROWS=[1,2], where the 1 is the amplitude row, 2 is
        %   the peak location row.
            if isempty(obj.PeaksFound)
                return;
            end
            n = size(obj.PeaksFound,2);
            if nargin == 1
                switch obj.CurveModelIndex
                    case num2cell(1:18), row = [1;2];
                    otherwise
                        error('Row of Model.CurveModel.StartParamas has to be specified.')
                end
            elseif nargin == 2
                row = varargin{1};
                if ~isvector(row) || ~isnumeric(row) || numel(row)~=2 || ~isequal(floor(row),row) || ...
                        max(row)>obj.Model.CurveModel.NOfParams(1) || min(row)<1 || row(1) == row(2)
                    error('Invalid row numbers of Model.CurveModel.StartParams');
                end
            else
                error('Invalid number of input arguments.');
            end
            old_n = obj.NOfModelCurves;
            if old_n<n
                obj = addcurves(obj,n-old_n);
            elseif old_n>n
                obj = reducecurves(obj,n+1:old_n);
            end
            obj.Model.CurveModel.StartParams(row(1),:) = obj.PeaksFound(1,:)*0.1;
            obj.Model.CurveModel.StartParams(row(2),:) = obj.PeaksFound(2,:);
        end
        function varargout = plot(obj,varargin)
        % PLOT Plot model calculations or fitted result.
        %   PLOT(OBJ) or OBJ.PLOT Preview the model calcuations.
        %
        %   PLOT(OBJ,MODE) or OBJ.PLOT(MODE) Plot model calculations
        %   (MODE='preview') or fitted result (MODE='fit'). Default is
        %   MODE='preview'. In the 'fit' mode, fitting residual is plotted
        %   and chi2 is labeled. 
        %
        %   HLINES=OBJ.PLOT(...) or HLINES=PLOT(OBJ,...) Return handles of
        %   the plotted lines. 
            if nargin == 1
                mode = 'preview';
            elseif nargin == 2
                mode = varargin{1};
                if ~ischar(mode) || (~strcmpi(mode,'preview') && ~strcmpi(mode,'fit'))
                    error('Invalid plot mode.');
                end
            end
            if strcmpi(mode,'preview')
                obj2 = obj;
                legend_str = {'Data';'Preview';'Background'};
            else
                obj2 = acceptfit(obj);
                legend_str = {'Data';'Fit';'Background'};
            end
            [ydata_cal,iydata_cal,ydata_bkgd] = evalmodel(obj2);            
            % --- start plotting
            hfig = figure;
            %clf(hfig);
            if strcmpi(mode,'preview')
                hax = axes('Parent',hfig);
            else
                hax = subplot(5,1,(1:4));
            end
            hline = plot(obj2.Data(:,1),obj2.Data(:,2),'o','Parent',hax);
            hold on; box on;
            plot(obj2.Data(:,1),ydata_cal,'-','Parent',hax,'linewidth',1.5);
            plot(obj2.Data(:,1),ydata_bkgd,'-','Parent',hax,'linewidth',1);
            if obj2.CurveModelIndex == -1 
                obj3_list = splitmultimodel(obj2);  % split mutlimodels to a list array of single models                
            end
            if obj2.NOfModelCurves>1 && obj2.CurveModelIndex ~= -1 
                plot(obj2.Data(:,1),iydata_cal,'--');
                legend_str = [legend_str;cellfun(@(x)['Curve # ',x],cellstr(num2str((1:obj2.NOfModelCurves)')),'UniformOutput',false)];
            elseif obj2.NOfModelCurves==1 && obj2.CurveModelIndex == -1 
                new_legend_str = cell(length(obj3_list),1);
                for jj=1:length(obj3_list)
                    [my_cal,~,~] = evalmodel(obj3_list(jj));
                    plot(obj3_list(jj).Data(:,1),my_cal,'--');
                    new_legend_str{jj} = obj3_list(jj).ModelBase.CurveModel.ModelName;
                end
                legend_str = [legend_str; new_legend_str];
            elseif obj2.NOfModelCurves>1 && obj2.CurveModelIndex == -1 
                new_legend_str = cell(length(obj3_list)*obj2.NOfModelCurves,1);
                for jj=1:length(obj3_list)
                     [~,miy_cal,~] = evalmodel(obj3_list(jj));
                     plot(obj3_list(jj).Data(:,1),miy_cal,'--');
                     ind_beg = (jj-1)*obj2.NOfModelCurves+1;     % for legend string
                     ind_end = ind_beg+obj2.NOfModelCurves-1;                     
                     new_legend_str(ind_beg:ind_end) = cellfun(@(x)['Curve # ',x,' - ',obj3_list(jj).ModelBase.CurveModel.ModelName],cellstr(num2str((1:obj3_list(jj).NOfModelCurves)')),'UniformOutput',false);
                end   
                legend_str = [legend_str; new_legend_str];                
            end
            hold off;
            if obj2.FlagLogScale==1     % check for log scale
                set(hax,'yscale','log');
            end
            legend(legend_str,'Parent',hfig,'Location','best');
            set(hax,'xminortick','on','yminortick','on');
            title(hax,['Curve: ',obj2.Model.CurveModel.ModelName,' | Background: ' obj2.Model.BkgdModel.ModelName]);
            % plot residue
            figure(hfig);            
            if strcmpi(mode,'fit')
                set(hax,'xticklabel',[]);
                hax2 = subplot(5,1,5);
                plot(obj2.Data(:,1),obj2.FitOutput.residual,'o','parent',hax2);
                set(hax2,'xminortick','on','yminortick','on');
                ylabel(hax2,'Residual');
                linkaxes([hax,hax2],'x');
                chi2 = sum((obj2.FitOutput.residual).^2./obj2.FitOutput.ydata);
                title(hax2,['\chi^2 = ',num2str(chi2),...
                    ' | Resnorm = ',num2str(obj2.FitOutput.resnorm),...
                    ' | Exitflag = ',num2str(obj2.FitOutput.exitflag)]);
            end
            if nargout == 1
                varargout{1} = hline;
            end
        end
        function obj = acceptfit(obj)
        % ACCEPTFIT Accept and pass fit result to Model.(CurveModel/BkgdModel).(StartParams/FitStdErrors)
        %   OBJ2=ACCEPT(OBJ1) or OBJ2=OBJ1.ACCEPT Accept and create a new
        %   object. Set OBJ2 to OBJ1 to update the current object.
            if isempty(fieldnames(obj.FitOutput))
                error('No fitting result exits.')
            end
            % assign fitted result and std errors
            fitflag = obj.FitParamsSet.fitflag;
            x=nan(1,length(fitflag));
            x(fitflag==1) = obj.FitOutput.fitted(:,1);
            x(fitflag==0) = obj.FitParamsSet.x2;
            x_curve = reshape(x(1:prod(obj.Model.CurveModel.NOfParams)),obj.Model.CurveModel.NOfParams(1),[]);
            x_bkgd = x(end+1-obj.Model.BkgdModel.NOfParams(1):end);
            err=-1*ones(1,length(fitflag));
            err(fitflag==1) = obj.FitOutput.fitted(:,2);
            err_curve = reshape(err(1:prod(obj.Model.CurveModel.NOfParams)),obj.Model.CurveModel.NOfParams(1),[]);
            err_bkgd = err(end+1-obj.Model.BkgdModel.NOfParams(1):end);
            fresult = obj.FitOutput;    % get fit result
            obj.Model.CurveModel.StartParams = x_curve;
            obj.Model.BkgdModel.StartParams  = x_bkgd(:);
            fitflag_curves  = logical(obj.Model.CurveModel.FitFlags(:));
            fitflag_bkdg    = logical(obj.Model.BkgdModel.FitFlags(:));
            obj.Model.CurveModel.FitStdErrors(fitflag_curves) = err_curve(fitflag_curves);
            obj.Model.BkgdModel.FitStdErrors(fitflag_bkdg)  = err_bkgd(fitflag_bkdg);
            obj.FitOutput = fresult;    % keep fit result after accepting fit
        end
        function [ydata_cal,iydata_cal,ydata_bkgd] = evalmodel(obj,varargin)
        % EVALMODEL Evaluate Model
        %   EVALMODEL(OBJ) or Y=OBJ.EVALMODEL Evaluate using XDATA in
        %   OBJ.DATA.
        %
        %   EVALMODEL(OBJ,XV) or OBJ.EVALMODEL(XV) Evaluate using a vector
        %   XV as the input value.
        %
        %   Y=EVALMODEL(OBJ,...) or Y=OBJ.EVALMODEL(...) Return the total
        %   evaluated values to Y.
        %
        %   [Y,IY]=EVALMODEL(OBJ,...) or [Y,IY]=OBJ.EVALMODEL(...) Also
        %   return the value for each curve.
        %
        %   [Y,IY,BKGD]=EVALMODEL(OBJ,...) or [Y,IY,BKGD]=OBJ.EVALMODEL(...)
        %   Returns the background values.
            if nargin == 1
                xdata = obj.Data(:,1);
            elseif nargin == 2 && isnumeric(varargin{1}) && isvector(varargin{1})
                xdata = varargin{1}(:);
            else
                error('Invalid input arguments');
            end
            [ydata_cal,iydata_cal,ydata_bkgd] = fcn_fit_core(...
                obj.FitParamsSet.x1,xdata,...
                obj.FitParamsSet,...
                obj.Model,...
                0,...
                obj.FlagNonNegBkgd);   % calcualte in linear scale
        end
        function splitobj= splitmultimodel(obj)
        % SPLITMULTIMODEL Split multi-models into LINEFIT objects each having a single model
        %   OBJS=SPLITMULTIMODEL(OBJ) or OBJS=OBJ.SPLITMULTIMODEL Split
        %   into an array of new objects.
            if obj.CurveModelIndex ~= -1
                error('The current model is not multi-model.');
            end
            modellist = obj.getmodellist;
            modellist(1:(-modellist.CurveModelIndex(1)+1),:) = [];  % remove [-2:0] models
            modelindex_str = obj.Model.CurveModel.ModelName;
            modelindex = str2num(modelindex_str(strfind(modelindex_str,':')+1:end)); %#ok<ST2NM>
            nofp = modellist.NOfParams(modelindex);
            StartParams = obj.Model.CurveModel.StartParams;
            LowerBounds = obj.Model.CurveModel.LowerBounds;
            UpperBounds = obj.Model.CurveModel.UpperBounds;
            FitFlags    = obj.Model.CurveModel.FitFlags;
            FitStdErrors = obj.Model.CurveModel.FitStdErrors;
            splitobj = repmat(obj,length(modelindex),1);
            for jj=1:length(modelindex)
                new_obj = obj;     % create a new object for each model
                new_obj.CurveModelIndex = modelindex(jj);
                new_obj.Model.CurveModel.StartParams    = StartParams(1:nofp(jj),:);
                new_obj.Model.CurveModel.LowerBounds    = LowerBounds(1:nofp(jj),:);
                new_obj.Model.CurveModel.UpperBounds    = UpperBounds(1:nofp(jj),:);
                new_obj.Model.CurveModel.FitFlags       = FitFlags(1:nofp(jj),:);
                new_obj.Model.CurveModel.FitStdErrors   = FitStdErrors(1:nofp(jj),:);
                splitobj(jj) = new_obj;
                StartParams(1:nofp(jj),:)   = [];
                LowerBounds(1:nofp(jj),:)   = [];
                UpperBounds(1:nofp(jj),:)   = [];
                FitFlags(1:nofp(jj),:)      = [];
                FitStdErrors(1:nofp(jj),:)  = [];
            end
        end
        function s = getmodelprops(obj,varargin)
        % GETMODELPROPS Get model properties
        %   S=GETMODELPROPS(OBJ) or S=OBJ.GETMODELPROPS(OBJ) Get the
        %   properties of the model using analytical method and return to a
        %   structure S. 
        %
        %   S=GETMODELPROPS(OBJ,XV) or S=OBJ.GETMODELPROPS(XV) Get the
        %   properties numerically using input values XV. 
            if nargin == 1             % --- check for method
                if obj.CurveModelIndex == -2
                    error('Analytical method is not defined for customized models. Please use numerical method.');
                else
                    getprops_method = 1;    % analytical
                end
            elseif nargin == 2
                getprops_method = 2;    % numerical
                xdata = varargin{1}(:);
                if ~isnumeric(xdata) || ~isvector(xdata)
                    error('Invalid x values for model evaulatoin.');
                end
                xdata = xdata(:);
            else
                error('Invalid number of input arguments.');
            end
            % --- check for obj
            if isempty(fieldnames(obj.FitOutput))
                obj2 = obj;
            else
                obj2 = acceptfit(obj);
            end
            if obj2.CurveModelIndex == -1   % for multi-models
                obj3_list = splitmultimodel(obj2);  % split mutlimodels to a list array of single models
                s = cell(length(obj3_list),1);  % initialize output cell
            end
            if getprops_method == 1     % analytical method
                if obj2.CurveModelIndex == -1   % for multi-models
                    for ii=1:length(s)
                        s{ii} = getmodelprops(obj3_list(ii));
                    end
                else % for other defined single model
                    props = [];
                    for ii=1:obj2.NOfModelCurves
                        fh_model = obj2.Model.CurveModel.ModelFcnHandle;
                        [~,tmp_props] = fh_model(obj2.Data(1,1),obj2.Model.CurveModel.StartParams(:,ii));
                        props = [props;tmp_props(:)']; %#ok<AGROW>
                    end
                    fname = fieldnames(props);
                    for ii=1:length(fname)
                        s.(fname{ii}) = cell2mat({props.(fname{ii})});
                    end
                end
            elseif  getprops_method == 2    % numerical method
                if obj2.CurveModelIndex == -1   % for multi-models
                    for ii=1:length(s)
                        s{ii} = getmodelprops(obj3_list(ii),xdata);
                    end
                else
                    [~,iydata_cal,~] = evalmodel(obj2,xdata);
                    % initialize
                    s.peak      = nan(3,obj2.NOfModelCurves);
                    s.area      = nan(1,obj2.NOfModelCurves);
                    s.mean      = s.area;
                    s.variance  = s.area;
                    s.skewness  = s.area;
                    s.kurtosis  = s.area;
                    s.fwhm      = s.area;
                    s.fwhm_x    = s.area;
                    for ii = 1:obj2.NOfModelCurves
                        tmp = obj2.lineprops([xdata,iydata_cal(:,ii)]);
                        s.peak(:,ii)    = tmp.peak(:);
                        s.area(ii)      = tmp.area;
                        s.mean(ii)      = tmp.mean;
                        s.variance(ii)  = tmp.variance;
                        s.skewness(ii)  = tmp.skewness;
                        s.kurtosis(ii)  = tmp.kurtosis;                        
                        s.fwhm(ii)      = tmp.fwhm;
                        s.fwhm_x(ii)    = tmp.fwhm_x;
                    end
                end
            end
        end
        function  obj = startfit(obj)
        % STARTFIT Start fitting
        %   OBJ2=STARTFIT(OBJ1) or OBJ2=OBJ1.STARTFIT Start fitgint and
        %   store fit result in FitOutput of OBJ2.
        %
        %   Set OBJ2 to OBJ1 to store fit result back to OBJ1.
            yfit = obj.Data(:,2);
            if obj.FlagLogScale == 1
                yfit = log10(yfit);
            end
            [fittedX1,resnorm,residual,exitflag,output,lambda,jacobian] = ...
                lsqcurvefit(@fcn_fit_core,obj.FitParamsSet.x1,obj.Data(:,1),yfit,...
                obj.FitParamsSet.lb1,obj.FitParamsSet.ub1,obj.FitOptions,...
                obj.FitParamsSet,obj.Model,obj.FlagLogScale,obj.FlagNonNegBkgd);
            % --- error estimation
            s2=resnorm/(length(residual) - 3);
            [Q,R]=qr(jacobian,0); %#ok<ASGLU>
            Rinv=inv(R);
            sigmaest=(Rinv*Rinv')*s2;
            stderrors=sqrt(diag(sigmaest));
            % --- assign result
            fresult.fitted      = [fittedX1,full(stderrors)];
            fresult.resnorm     = resnorm;
            fresult.residual    = residual;
            fresult.exitflag    = exitflag;
            fresult.output      = output;
            fresult.lambda      = lambda;
            fresult.jacobian    = jacobian;
            fresult.ydata       = residual+yfit;
            if obj.FlagLogScale == 1
                fresult.ydata   = 10.^fresult.ydata;
            end
            obj.FitOutput = fresult;
        end
        function obj = findpeakauto(obj,varargin)
        % FINDPEAKAUTO Automatically search for peaks
        %   OBJ2=FINDPEAKAUTO(OBJ1) or OBJ2=OBJ1.FINDPEAKAUTO Use
        %   PEAKFINDER (external function) to find peaks and store them to
        %   OBJ2.PEAKSFOUND.
        %
        %   OBJ2=FINDPEAKAUTO(OBJ1,ARG1,ARG2,ARG3,ARG4,ARG5) or
        %   OBJ2=OBJ1.FINDPEAKAUTO(,ARG1,ARG2,ARG3,ARG4,ARG5) Use arguement
        %   defined in PEAKFINDER
        %
        %   PEAKFINDER.m (File ID: #25500) can be downloaded from Matlab
        %   Central File Exchange.
            a = cell(1,5);
            a(1:nargin-1) = varargin;
            if obj.FlagLogScale == 1
                ind = peakfinder(log10(obj.Data(:,2)),a{1},a{2},a{3},a{4},a{5});
            else
                ind = peakfinder(obj.Data(:,2),a{1},a{2},a{3},a{4},a{5});
            end
            if isempty(ind)
                return;
            end
            x = obj.Data(ind,1)';
            y = obj.Data(ind,2)';
            figure
            hold on; box on;
            plot(obj.Data(:,1),obj.Data(:,2),'-');
            if obj.FlagLogScale == 1
                set(gca,'yscale','log');
            end
            set(gca,'xminortick','on','yminortick','on');            
            ymax = max(obj.Data(:,2));
            ymin = min(obj.Data(:,2));
            for ii=1:length(ind)
                plot([x(ii),x(ii)],[ymin,ymin+1.1*(ymax-ymin)],'r:');
                text(x(ii),ymin+1.1*(ymax-ymin),num2str(ii),'color','red');
            end
            hold off;
            obj.PeaksFound = [y;x];
        end
        function obj = findpeakmanual(obj)
        % FINDPEAKMANUAL Manually and interactively search for peaks
        %   OBJ2=FINDPEAKMANUAL(OBJ1) or OBJ2=OBJ1.FINDPEAKMANUAL Use
        %   interactive tool GETPTS (requires Matlab Image Processing Toolbox)
        %   to label peaks.
            if nargin ~= 1 || isempty(obj.Data)
                error('Invalid number of input arguments or no data exists.');
            end
            figure; box on;
            plot(obj.Data(:,1),obj.Data(:,2),'-');
            set(gca,'xminortick','on','yminortick','on');                        
            if obj.FlagLogScale == 1
                set(gca,'yscale','log');
            end
            if exist('getpts') == 2 %#ok<EXIST>
                [x,~]=getpts(gcf);
            elseif exist('gu_getpts') == 2 %#ok<EXIST>
                [x,~]=gu_getpts(gcf);
            else
                error('Need Matlab Image Processing Toolbox');
            end
            if isempty(x)
                return;
            end
            x = x';
            y = interp1(obj.Data(:,1),obj.Data(:,2),x,'linear','extrap');
            ymax = max(obj.Data(:,2));
            ymin = min(obj.Data(:,2));
            hold on; box on;
            for ii=1:length(x)
                plot([x(ii),x(ii)],[ymin,ymin+1.1*(ymax-ymin)],'r:');
                text(x(ii),ymin+1.1*(ymax-ymin),num2str(ii),'color','red');
            end
            hold off;            
            obj.PeaksFound = [y;x];
        end
    end % End of dynamic methods
    methods(Static) % static method
        function s = createmodel(varargin)
        % CREATEMODEL Create model stucture for LINEFIT
        %   S=CREATEMODEL(FCN_HANDLE,NOFP,STARTPARAMS,LOWERBOUNDS,UPPERBOUNDS,FITFLAGS)
        %   FCN_HANDLE is a function handle for the model. NOFP is the
        %   number of total parameters in the model. STARTPARAMS, LOWERBOUNDS,
        %   and UPPERBOUNDS are initial parameters, lower and upper bounds
        %   for the fitting. FITFLAGS are logical flags (false/true or 0/1)
        %   to indicate whether the parameter is to be fixed or fitted
        %   during the fitting. STARTPARAMS,LOWERBOUNDS,UPPERBOUNDS and
        %   FITFLAGS must be either vector of the dimension NOFPx1, or
        %   empty to use default values (random values for STARTPARAMS,
        %   -Inf, Inf and 1 for LOWERBOUNDS, UPPERBOUNDS, and FITFLAGS).
        %
        %   Note on FCN_HANDLE: it must have the format @(x,p)UserFcn(x,p)
        %   where x is the input x value, p is a vector of parameters; or
        %   it can be an anonymous function, e.g.,
        %   f=@(x,p)p(1)+p(2)*x+p(3)*sin(x), which gives 3.7296 when
        %   evaluate x=4 for p=[2,1,3], i.e. f(4,[2,1,3] = 3.7296
        %
        %   For Multi-Model (CurveModelIndex=-1) where the Model.CurveModel is
        %   constructed as a sum of built-in modles, use
        %   S=CREATEMODEL(MODELS,[],STARTPARAMS,LOWERBOUNDS,UPPERBOUNDS,FITFLAGS)
        %   where MODELS is a list of index of built-in models, e.g.,
        %   MODELS=[1,4,13] constructs a multi-model using the sum of
        %   Lorentzian (Model #1), Gaussian (Model #4) and Voigt (Model
        %   #13). The default of NOPF will be set automatically to the sum
        %   of input parameters of the built-in models.
        %
        %   S=CREATEMODEL(...,FITFLAGS,ModelName) Gives name (string) for
        %   the model.
        %
        %   S=CREATEMODEL(...,FITFLAGS,ModelName,MODELUSAGE) Gives useage
        %   hint (string) for the model.
            if nargin ~= 6 && nargin ~= 7 && nargin ~= 8
                error('Invalid number of input arguments.');
            end
            if isa(varargin{1},'function_handle')
                flag_model = 1;     % single customized or predefined model
                s.ModelFcnHandle    = varargin{1};
            elseif isnumeric(varargin{1})
                modellist = linefit.getmodellist;
                modellist(1:(-modellist.CurveModelIndex(1)+1),:) = [];  % remove [-2:0] models
                modelindex = varargin{1}(:);
                if nnz(ismember(modelindex,1:(size(modellist,1))))~=length(modelindex)
                    error('Invalid model index.')
                end
                s.ModelFcnHandle = @(x,p)multimodels(x,p);
                flag_model = 2;     % stacked predefined models
            else
                error('Invalid input arguments.')
            end
            s.NOfParams         = varargin{2};
            s.StartParams       = varargin{3}(:);
            s.LowerBounds       = varargin{4}(:);
            s.UpperBounds       = varargin{5}(:);
            s.FitFlags          = varargin{6}(:);
            if flag_model == 2
                nofp = modellist.NOfParams(modelindex);
                s.NOfParams = sum(nofp);
            elseif isempty(s.NOfParams) || ~isscalar(s.NOfParams) || ~isnumeric(s.NOfParams) || s.NOfParams<0 || floor(s.NOfParams)~=s.NOfParams
                error('Invalid total number of parameters.');
            end
            if isempty(s.StartParams)
                if flag_model == 2
                    s.StartParams = cell2mat(modellist.StartParams(modelindex)')';
                else
                    s.StartParams = rand(s.NOfParams,1);
                end
            elseif numel(s.StartParams) ~= s.NOfParams
                error('Size of start parameters does not match total number of parameters.');
            end
            if isempty(s.LowerBounds)
                if flag_model == 2
                    s.LowerBounds = cell2mat(modellist.LowerBounds(modelindex)')';
                else                
                    s.LowerBounds = -Inf*ones(s.NOfParams,1);
                end
            elseif numel(s.LowerBounds) ~= s.NOfParams
                error('Size of lower bounds does not match total number of parameters.');
            end
            if isempty(s.UpperBounds)
                if flag_model == 2
                    s.UpperBounds = cell2mat(modellist.UpperBounds(modelindex)')';
                else                                
                    s.UpperBounds = Inf*ones(s.NOfParams,1);
                end
            elseif numel(s.UpperBounds) ~= s.NOfParams
                error('Size of upper bounds does not match total number of parameters.');
            end
            if isempty(s.FitFlags)
                if flag_model == 2
                    s.FitFlags = cell2mat(modellist.FitFlags(modelindex)')';
                else                 
                    s.FitFlags = true(s.NOfParams,1);
                end
            elseif numel(s.FitFlags) ~= s.NOfParams
                error('Size of fit flags does not match total number of parameters.');
            end
            if nargin == 6
                s.ModelName     = 'Custom';
                s.ModelUsage    = 'Not specified.';
            elseif nargin == 7
                s.ModelName     = varargin{7};
                s.ModelUsage    = 'Not specified.';
            else
                s.ModelName     = varargin{7};
                s.ModelUsage    = varargin{8};
            end
            if flag_model == 2
                s.ModelName = ['Multi-Model:',sprintf('% d',(modelindex'))];  % use model name to store stacked model index
            end
            s.FitStdErrors  = -1*ones(s.NOfParams,1);
            function varargout = multimodels(x,p)
                y = zeros(size(x));
                for ii=1:length(modelindex)
                    y = y+modellist.ModelFcnHandle{modelindex(ii)}(x,p(1:nofp(ii)));
                    p(1:nofp(ii)) = [];
                end
                varargout{1} = y;
            end
        end
        % --- method to numerically find line statistics (peak, com, fwhm)
        function s = lineprops(data)
        % LINEPROPS Return line properties numerically
        %   S=OBJ.LINEPROPS(DATA) or S=LINEFIT.LINEPROPS(DATA) Calculate
        %   area, peak, mean, fwhm and the position of fwhm of a 2-column
        %   DATA, and returns a structure S with fields: area, peak, mean,
        %   fwhm, and fwhm_x. 
        %
        %   Area: numerically integrated area (with Matlab trapz function)
        %   Peak: the [X,Y] values, and position index of the peak
        %   Mean: mean value
        %   FWHM: full-width-half-maximum
        %   FWHM_X: the center X position of FWHM
            data = data(:,1:2);
            % --- peak
            s.area = trapz(data(:,1),data(:,2));
            [peak,peakIndex] = max(data(:,2));
            s.peak = [data(peakIndex,:),peakIndex];
            % --- Mean
            s.mean = trapz(data(:,1),data(:,1).*data(:,2))/trapz(data(:,1),data(:,2));
            %            s.com = trapz(data(:,1),data(:,2:end).*repmat(data(:,1),1,size(data,2)-1)) ./  trapz(data(:,1),data(:,2:3));
            % --- variance
            s.variance =  trapz(data(:,1),(data(:,1)-s.mean).^2.*data(:,2))/trapz(data(:,1),data(:,2));
            % --- skewness
            s.skewness =  trapz(data(:,1),(data(:,1)-s.mean).^3.*data(:,2))/trapz(data(:,1),data(:,2))/s.variance^1.5;            
            % --- kurtosis
            s.kurtosis =  trapz(data(:,1),(data(:,1)-s.mean).^4.*data(:,2))/trapz(data(:,1),data(:,2))/s.variance^2-3;            
            % --- FWHM
            try
                leftError = 0;
                halfLeftIndex = find(data(1:peakIndex,2)<peak/2);
                left_data = data(halfLeftIndex(end):peakIndex,:);
                [~,m,~] = unique(left_data(:,2));
                left_data = left_data(m,:);
                left_data = sortrows(left_data,1);
                left = interp1(...
                    left_data(:,2),left_data(:,1),peak/2);
            catch
                leftError = 1;
            end
            try
                rightError = 0;
                halfRightIndex = find(data(peakIndex:end,2)<peak/2)+peakIndex-1;
                right_data = data(peakIndex:halfRightIndex(1),:);
                [~,m,~] = unique(right_data(:,2));
                right_data = right_data(m,:);
                right_data = sortrows(right_data,1);
                right = interp1(...
                    right_data(:,2),right_data(:,1),peak/2);
            catch
                rightError = 1;
            end
            if leftError == 0 && rightError == 0
                s.fwhm = right-left;
                s.fwhm_x = (left+right)/2;
            elseif leftError == 0 && rightError == 1
                s.fwhm = data(end,1)-left;
                s.fwhm_x = left;
            elseif leftError == 1 && rightError == 0
                s.fwhm = right-data(1,1);
                s.fwhm_x = right;
            else
                s.fwhm = data(end,1)-data(1,1);
                s.fwhm_x = (data(end,1)+data(1,1))/2;
            end
        end
        % --- collect built-in models
        function s = getmodellist
        % GETMODELLIST Get the list of built-in models
        %   S=OBJ.GETMODELLIST or S=LINEFIT.GETMODELLIST return the list to
        %   table S.
            if nargin>0
                error('No input argument is needed');
            end
            tmp_obj = linefit;
            ModelName   = ['User-Customized';'Predefined Multi-Model';'Reserved';{tmp_obj.ModelList.ModelName}'];
            ModelUsage  = ['User-Customized';'Predefined Multi-Model';'Reserved';{tmp_obj.ModelList.ModelUsage}'];
            ModelFcnHandle = [{@()[];@()[];@()[]};{tmp_obj.ModelList.ModelFcnHandle}'];                        
            NOfParams =   [NaN;NaN;NaN;cell2mat({tmp_obj.ModelList.NOfParams}')];            
            StartParams = [NaN;NaN;NaN;cellfun(@(x)(x'),{tmp_obj.ModelList.StartParams}','UniformOutput',false)];
            LowerBounds = [NaN;NaN;NaN;cellfun(@(x)(x'),{tmp_obj.ModelList.LowerBounds}','UniformOutput',false)];
            UpperBounds = [NaN;NaN;NaN;cellfun(@(x)(x'),{tmp_obj.ModelList.UpperBounds}','UniformOutput',false)];
            FitFlags    = [NaN;NaN;NaN;cellfun(@(x)(x'),{tmp_obj.ModelList.FitFlags}','UniformOutput',false)];            
            CurveModelIndex = (1:length(ModelName))'-3;
            s = table(CurveModelIndex,ModelName,ModelUsage,ModelFcnHandle,NOfParams,StartParams,LowerBounds,UpperBounds,FitFlags);
        end
        % =========================================
        % --- built-in model functions
        % =========================================
        function varargout = gaussian(x,p) 
        % GAUSSIAN (Built-in Model) Gaussian or normal distribution
        %   Y=(OBJ/LINEFIT).GAUSSIAN(x,[AMP,X0,SIGMA])
        %       AMP: Amplitude, X0: Location, SIGMA: Scale
            amp    = p(1); 
            a       = p(2);    
            b       = p(3);         
            varargout{1} = amp/sqrt(2*pi*b^2)*exp(-(x-a).^2/(2*b^2));
            if nargout == 2
                props.amp      = amp;
                props.x0         = a;
                props.sigma      = b;
                props.mean      = a;
                props.variance  = b^2;
                props.skewness  = 0;
                props.kurtosis  = 0;
                props.mode      = a;
                props.median    = a;      
                props.fwhm      = 2*sqrt(2*log(2))*b;
                varargout{2} = props;
            end
        end
        function varargout = generalizednorm1(x,p) 
        % GENERALIZEDNORM1 (Built-in Model) Genearalized normal version 1 (or exponential power) distribution
        %   Y=(OBJ/LINEFIT).GENERALIZEDNORM1(X,[AMP,A,B,C])     
        %       AMP: Amplitude, A: Location, B: Scale, C: Shape
        %
        %   Reference:
        %   (1) https://en.wikipedia.org/wiki/Generalized_normal_distribution
            amp     = p(1); 
            a       = p(2);    
            b       = p(3);         
            c       = p(4);   
            varargout{1} = amp*c/(2*b*gamma(1/c))*exp(-(abs(x-a)/b).^c);
            if nargout == 2
                props.amp      = amp;
                props.a         = a;
                props.b         = b;
                props.c         = c;
                props.mean      = a;
                props.variance  = b^2*gamma(3/c)/gamma(1/c);
                props.skewness  = 0;
                props.kurtosis  = gamma(5/c)*gamma(1/c)/gamma(3/c)^2-3;
                props.mode      = a;
                props.median    = a;
                varargout{2} = props;
            end
        end
        function varargout = lorentzian(x,p) 
        % LORENTZIAN (Built-in Model) Lorentzian or Cauchy distribution
        %   Y=(OBJ/LINEFIT).LORENTZIAN(X,[AMP,X0,GAMMA])
        %       AMP: Amplitude, X0: Location, GAMMA: Scale or HWHM
            amp    = p(1); 
            a       = p(2);    
            b       = p(3);             
            varargout{1} = amp/(pi*b)./(1+((x-a)/b).^2);
            if nargout == 2
                props.amp      = amp;
                props.x0         = a;
                props.gamma      = b;
                props.mode      = a;
                props.median    = a;
                props.fwhm      = 2*b;
                varargout{2} = props;
            end            
        end
        function varargout = ilorentzian(x,p) 
        % ILORENTZIAN (Built-in Model) Intermediate Lorentzian
        %   Y=(OBJ/LINEFIT).ILORENTZIAN(X,[AMP,X0,GAMMA])
        %       AMP: Amplitude, X0: Location, GAMMA: Scale
        % 
        %   Reference: 
        %   (1) C.P. Khattak and D.E. Cox, J. App. Cryst. 10, 405 (1977)
            amp    = p(1); 
            a       = p(2);    
            b       = p(3);             
            varargout{1} = amp*sqrt(2^(2/3)-1)/(2*b)./(1+(2^(2/3)-1)/b^2*(x-a).^2).^1.5;
            if nargout == 2
                props.amp      = amp;
                props.x0        = a;
                props.gamma     = b;
                props.fwhm      = 2*b;
                varargout{2} = props;                
            end
        end
        function varargout = mlorentzian(x,p)
        % MLORENTZIAN (Built-in Model) Modified Lorentzian
        %   Y=(OBJ/LINEFIT).MLORENTZIAN(X,[AMP,X0,GAMMA])
        %       AMP: Amplitude, X0: Location, GAMMA: Scale
        %
        %   Reference:
        %   (1) G. Malmros and J.O. Thomas, J. Appl. Cryst. 10, 7 (1977)
            amp    = p(1); 
            a       = p(2);    
            b       = p(3);             
            varargout{1} = amp*2*sqrt(sqrt(2)-1)/pi/b./(1+(sqrt(2)-1)/b^2*(x-a).^2).^2;
            if nargout == 2
                props.amp      = amp;
                props.x0        = a;
                props.gamma     = b;
                props.fwhm      = 2*b;
                varargout{2} = props;                
            end            
        end
        function varargout = pearsonVII(x,p) 
        % PEARSONVII (Built-in Model) Pearson type VII or student't distribution
        %   Y=(OBJ/LINEFIT).PEARSONVII(X,[AMP,A,B,C])
        %       AMP: Amplitude, A: Location, B: Scale, C: Shape
        %
        %   Reference:
        %   (1) M.P. McLaughlin, Compendium of Common Probability
        %   Distributions, 2nd Edition,
        %   http://www.causascientia.org/math_stat/Dists/Compendium.pdf 
        %   (2) https://en.wikipedia.org/wiki/Student's_t-distribution
            amp    = p(1); 
            a       = p(2);    
            b       = p(3);  
            c       = p(4);
            varargout{1} = amp/(b*sqrt(c)*beta(c/2,0.5))*(1+(x-a).^2/(b^2*c)).^(-(c+1)/2);
            if nargout == 2
                props.amp      = amp;
                props.a         = a;
                props.b         = b;
                props.c         = c;
                if c>1
                    props.mean      = a;
                else
                    props.mean      = NaN;
                end
                if c>2
                    props.variance  = c*b^2/(c-2);
                elseif c<=2 && c>1
                    props.variance  = Inf;
                else
                    props.variance   = NaN;
                end
                if c>3
                    props.skewness  = 0;
                else
                    props.skewness  = NaN;
                end
                if c>4
                    props.kurtosis  = 6/(c-4);
                elseif c<=4 && c>2
                    props.kurtosis  = Inf;
                else
                    props.kurtosis  = NaN;
                end
                props.mode      = a;
                props.median    = a;
                varargout{2} = props;                
            end 
        end
        function varargout = voigt(x,p)
        % VOIGT (Built-in Model) Voigt distribution - convoluation of Gaussian and Lorentzian
        %   Y=(OBJ/LINEFIT).VOIGT(X,[AMP,X0,SIGMA,GAMMA])
        %       AMP: Amplitude, X0: Location, SIGMA: Scale of Gaussian,
        %       GAMMA: Scale of Lorentzian
        %
        %   Note: requires Steven G. Johnson's Faddeeva Package (File ID:
        %   #38787) on Matlab Central File Exchange
        %
        %   Reference:
        %   (1) Steven G. Johnson, Faddeeva Package: complex error functions
        %   http://ab-initio.mit.edu/Faddeeva
            amp    = p(1); 
            x0      = p(2);
            sigma   = p(3);
            gamma   = p(4);
            z = (x-x0+1i*gamma)/sigma/sqrt(2);
            w = Faddeeva_w(z);
            varargout{1} = amp*real(w)/sigma/sqrt(2*pi);
            if nargout == 2
                props.amp      = p(1);
                props.x0        = p(2);
                props.sigma     = p(3);
                props.gamma     = p(4);
                props.fwhm      = 0.5346*(2*gamma)+sqrt(0.2166*(2*gamma).^2+(2*sigma*sqrt(2*log(2))).^2);
                varargout{2} = props;
            end
        end
        function varargout = pseudovoigtTCH(x,p) 
        % PSEUDOVOIGTTCH (Built-in Model) Pseudo-voigt type TCH    
        %   Y=(OBJ/LINEFIT).PSEUDOVOIGTTCH(X,[AMP,X0,SIGMA,GAMMA])
        %       AMP: Amplitude, X0: Location, SIGMA: Scale of Gaussian,
        %       GAMMA: Scale of Lorentzian        
        %
        %   Reference:
        %   (1) P. Thompson, et al., J. Appl. Cryst. 20, 79 (1987)
            amp    = p(1); 
            x0      = p(2);
            sigma   = p(3);
            gamma   = p(4);
            % FWHM of original G and L
            fG = 2*sigma*sqrt(2*log(2));
            fL = 2*gamma;
            % FWHM and eta
            f = (fG.^5 + 2.69269*fG.^4.*fL + 2.42843*fG.^3.*fL.^2 ...
                + 4.47163*fG.^2.*fL.^3 + 0.07842*fG.*fL.^4 + fL.^5).^0.2;
            eta = 1.36603*(fL./f) - 0.47719*(fL./f).^2 + 0.11116*(fL./f).^3;
            % new sigma and gamma for new Gaussian and Lorentzian
            sigma = f/2/sqrt(2*log(2));
            gamma = f/2;
            G = 1/sqrt(2*pi)./sigma.*exp(-(x-x0).^2/2./sigma.^2);
            L = 1/pi*gamma./((x-x0).^2+gamma.^2);
            varargout{1} = amp*(eta.*L+(1-eta).*G);
            if nargout == 2
                props.amp      = p(1);
                props.x0        = p(2);
                props.sigma     = p(3);
                props.gamma     = p(4);
                props.fwhm      = f;
                varargout{2} = props;
            end
        end
        function varargout = pseudovoigtIAT(x,p) 
        % PSEUDOVOIGTIAT (Built-in Model) Pseudo-voigt type IAT    
        %   Y=(OBJ/LINEFIT).PSEUDOVOIGTIAT(X,[AMP,X0,SIGMA,GAMMA])
        %       AMP: Amplitude, X0: Location, SIGMA: Scale of Gaussian,
        %       GAMMA: Scale of Lorentzian        
        %
        %   Reference:
        %   (1) T. Ida, et al., J. Appl. Cryst. 33, 1311 (2000)
            amp    = p(1); 
            x0      = p(2);
            sigma   = p(3);
            gamma   = p(4);            
            pcoeff = [...
                0.66        -0.42179	1.19913     1.10186     -0.30165	0.25437     1.01579
                0.15021     -1.25693	1.43021     -0.47745	-1.38927	-0.14107	1.50429
                -1.24984	10.30003	-15.36331	-0.68688	9.3155      3.23653     -9.21815
                4.74052     -23.45651	47.06071	2.76622     -24.10743	-11.09215	23.59717
                -9.48291	29.14158	-73.61822	-4.55466	34.96491	22.10544	-39.71134
                8.48252     -16.50453	57.92559	4.05475     -21.18862	-24.12407	32.83023
                -2.95553	3.19974     -17.80614	-1.26571	3.7029      9.76947     -10.02142
                ];
            fG  = 2*sigma*sqrt(2*log(2));
            fL  = 2*gamma;            
            rho = fL/(fG+fL);
            n = (0:6)';
            wG = 1-rho*sum(pcoeff(:,1).*rho.^n);
            wL = 1-(1-rho)*sum(pcoeff(:,2).*rho.^n);
            wI = sum(pcoeff(:,3).*rho.^n);
            wP = sum(pcoeff(:,4).*rho.^n);
            etaL = rho*(1+(1-rho)*sum(pcoeff(:,5).*rho.^n));
            etaI = rho*(1-rho)*sum(pcoeff(:,6).*rho.^n);
            etaP = rho*(1-rho)*sum(pcoeff(:,7).*rho.^n);
            WG = (fG+fL)*wG;
            WL = (fG+fL)*wL;
            WI = (fG+fL)*wI;
            WP = (fG+fL)*wP;
            gammaG = WG/(2*sqrt(log(2)));
            gammaL = WL/2;
            gammaI = WI/(2*sqrt(2^(2/3)-1));
            gammaP = WP/(2*log(sqrt(2)+1));
            x = x-x0;
            G = 1/(sqrt(pi)*gammaG)*exp(-x.^2/gammaG^2);
            L = 1/(pi*gammaL)./(1+(x/gammaL).^2);
            I = 1/(2*gammaI)./(1+(x/gammaI).^2).^(3/2);
            P = 1/(2*gammaP)./cosh(x/gammaP).^2;
            pv = (1-etaL-etaI-etaP)*G + etaL*L + etaI*I + etaP*P;
            varargout{1} = amp*pv;
            if nargout == 2
                props.amp      = p(1);
                props.x0        = p(2);
                props.sigma     = p(3);
                props.gamma     = p(4);
                props.fwhm      = (fG.^5 + 2.69269*fG.^4.*fL + 2.42843*fG.^3.*fL.^2 ...
                + 4.47163*fG.^2.*fL.^3 + 0.07842*fG.*fL.^4 + fL.^5).^0.2;
                varargout{2} = props;
            end            
        end 
        function varargout = pseudovoigtLLHGD(x,p)
        % PSEUDOVOIGTLLHGD (Built-in Model) Pseudo-voigt type LLHGD  
        %   Y=(OBJ/LINEFIT).PSEUDOVOIGTLLHGD(X,[AMP,X0,SIGMA,GAMMA])
        %       AMP: Amplitude, X0: Location, SIGMA: Scale of Gaussian,
        %       GAMMA: Scale of Lorentzian     
        %
        %   Reference:
        %   (1) Yuyan Liu, et al., J. Opt. Soc. Am. B, 18, 666 (2001)
            amp    = p(1); 
            x0      = p(2);
            sigma   = p(3);
            gamma   = p(4);     
            fG  = 2*sigma*sqrt(2*log(2));
            fL  = 2*gamma;
            f = 0.5346*(2*gamma)+sqrt(0.2165975*(2*gamma).^2+(2*sigma*sqrt(2*log(2))).^2);
            d = (fL-fG)/(fL+fG);
            cG = 0.32460 - 0.61825*d + 0.17681*d^2 + 0.12109*d^3;
            cL = 0.68188 + 0.61293*d - 0.18384*d^2 - 0.11568*d^3;
            % new sigma and gamma for new Gaussian and Lorentzian
            sigma = f/2/sqrt(2*log(2));
            gamma = f/2;
            % new G and L with identical fwhm
            G = 1/sqrt(2*pi)./sigma.*exp(-(x-x0).^2/2./sigma.^2);
            L = 1/pi*gamma./((x-x0).^2+gamma.^2);
            % pv
            pv = cG*G + cL*L;
            varargout{1} = amp*pv;
            if nargout == 2
                props.amp      = p(1);
                props.x0        = p(2);
                props.sigma     = p(3);
                props.gamma     = p(4);
                props.fwhm      = f;
                varargout{2} = props;
            end            
        end 
        function varargout = skewpvSB(x,p)  
        % SKEWPVSB (Built-in Model) Asymmetric pseudo-Voigt type SB     
        %   Y=(OBJ/LINEFIT).SKEWPVSB(X,[AMP,X0,FWHM0,ETA,A]
        %       AMP: Amplitude, X0: Location, FWHM0: FWHM of symmetic
        %       Voigt, ETA: Fraction of Lorentzian contribution, A: Degreee
        %       of asymmetry 
        %
        %   Reference: 
        %   (1) A.L. Stancik and E.B. Brauns, Vib. Spectr. 47, 66 (2008)
            varargout{1} = linefit.skewpvSSG(x,[p(:);0]);
            if nargout == 2
                props.amp      = p(1);
                props.x0        = p(2);
                props.fwhm0     = p(3);
                props.eta       = p(4);
                props.a         = p(5);
                varargout{2} = props;
            end             
        end
        function varargout = skewpvSSG(x,p)
        % SKEWPVSSG  (Built-in Model) Asymmetric pseudo-Voigt type SSG     
        %   Y=(OBJ/LINEFIT).SKEWPVSSG(X,[AMP,X0,FWHM0,ETA,A,B]
        %       AMP: Amplitude, X0: Location, FWHM0: FWHM of symmetic
        %       Voigt, ETA: Fraction of Lorentzian contribution, A: Degreee
        %       of asymmetry, B: Location shift of sigmoidal
        %
        %   Reference: 
        %   (1) Martin Schmid et al., Surf. Interface Anal. 46, 505 (2014)
            amp    = p(1);
            x0      = p(2);
            fwhm0   = p(3);
            eta     = p(4);     % contribution of L 
            a       = p(5);
            b       = p(6);
            x = x-x0;
            fwhm = 2*fwhm0./(1+exp(-a*(x-b)));
            % new sigma and gamma for new Gaussian and Lorentzian
            sigma = fwhm/2/sqrt(2*log(2));
            gamma = fwhm/2;
            % new G and L with identical fwhm
            G = 1/sqrt(2*pi)./sigma.*exp(-x.^2/2./sigma.^2);
            L = 1/pi*gamma./(x.^2+gamma.^2);
            % pv 
            pv = (1-eta)*G + eta*L;
            varargout{1} = amp*pv;
            if nargout == 2
                props.amp      = p(1);
                props.x0        = p(2);
                props.fwhm0     = p(3);
                props.eta       = p(4);
                props.a         = p(5);
                props.b         = p(6);   
                varargout{2} = props;
            end  
        end      
        function varargout = skewnorm(x,p) 
        % SKEWNORM (Built-in Model) Skew-normal distribution
        %   Y=(OBJ/LINEFIT).SKEWNORM(X,[AMP,A,B,C])
        %       AMP: Amplitude, A: Location, B: Scale, C: Shape
        %
        %   Reference:
        %   (1) M.P. McLaughlin, Compendium of Common Probability
        %   Distributions, 2nd Edition,
        %   http://www.causascientia.org/math_stat/Dists/Compendium.pdf 
        %   (2) https://en.wikipedia.org/wiki/Skew_normal_distribution
            amp    = p(1);
            a       = p(2);
            b       = p(3);
            c       = p(4);
            varargout{1}    = amp/b*sqrt(2/pi)*exp(-(x-a).^2/(2*b^2))*(1/2).*(1+erf( (c*(x-a)/b)  /sqrt(2)));
            if nargout == 2
                delta = c/sqrt(1+c^2);
                props.amp      = amp;
                props.a         = a;
                props.b         = b;
                props.c         = c;
                props.mean      = a+b*delta*sqrt(2/pi);
                props.variance  = b^2*(1-2*delta^2/pi);
                props.skewness  = (4-pi)/2*(delta*sqrt(2/pi))^3/(1-2*delta^2/pi)^1.5;
                props.kurtosis  = 2*(pi-3)*(delta*sqrt(2/pi))^4/(1-2*delta^2/pi)^2;
                varargout{2} = props;   
            end
        end
        function varargout = logistic(x,p) 
        % LOGISTIC (Built-in Model) Logistic distribution
        %   Y=(OBJ/LINEFIT).LOGISTIC(X,[AMP,A,B])
        %       AMP: Amplitude, A: Location, B: Scale
        %
        %   Reference:
        %   (1) M.P. McLaughlin, Compendium of Common Probability
        %   Distributions, 2nd Edition,
        %   http://www.causascientia.org/math_stat/Dists/Compendium.pdf             
            amp    = p(1);
            a       = p(2);
            b       = p(3);
            varargout{1} = amp/b*exp(-(x-a)/b)./(1+exp(-(x-a)/b)).^2;
            if nargout == 2
                props.amp      = amp;
                props.a         = a;
                props.b         = b;
                props.mean      = a;
                props.variance  = pi^2*b^2/3;
                props.skewness  = 0;
                props.kurtosis  = 7/5;  
                props.mode      = a;
                props.median    = a;
                varargout{2} = props;                   
            end
        end
        function varargout = inversenorm(x,p) 
        % INVERSENORM (Built-in Model) Inverse normal distribution
        %   Y=(OBJ/LINEFIT).INVERSENORM(X,[AMP,X0,A,B])
        %       AMP: Amplitude, X0: Location, A: Location and scale, B: Scale
        %
        %   Reference:
        %   (1) M.P. McLaughlin, Compendium of Common Probability
        %   Distributions, 2nd Edition,
        %   http://www.causascientia.org/math_stat/Dists/Compendium.pdf             
        %   (2) https://en.wikipedia.org/wiki/Inverse_Gaussian_distribution
            amp    = p(1);
            x0      = p(2);
            a       = p(3);
            b       = p(4);
            y = zeros(size(x));
            ind = (x>x0);
            y(ind) = amp*sqrt(b./(2*pi*(x(ind)-x0).^3)).*exp(-b./(2*(x(ind)-x0)).*((x(ind)-x0-a)/a).^2);
            varargout{1} = y;
            if nargout == 2
                props.amp      = amp;
                props.x0         = x0;                
                props.a         = a;
                props.b         = b;
                props.mean      = x0+a;
                props.variance  = a^3/b;
                props.skewness  = 3*sqrt(a/b);
                props.kurtosis  = 15*a/b;
                props.mode      = x0+a/(2*b)*(sqrt(9*a^2+4*b^2)-3*a);
                varargout{2} = props;
            end
        end
        function varargout = expdist(x,p)  
        % EXPDIST (Built-in Model) Exponential distribution
        %   Y=(OBJ/LINEFIT).EXPDIST(X,[AMP,A,B])
        %       AMP: Amplitude, A: Location, B: Scale
        %
        %   Reference:
        %   (1) M.P. McLaughlin, Compendium of Common Probability
        %   Distributions, 2nd Edition,
        %   http://www.causascientia.org/math_stat/Dists/Compendium.pdf                
            amp    = p(1);
            a       = p(2);
            b       = p(3);
            y = zeros(size(x));
            ind = (x>=a);
            y(ind) = amp/b*exp(-(x(ind)-a)/b);
            varargout{1} = y;
            if nargout == 2
                props.amp      = amp;
                props.a         = a;
                props.b         = b;
                props.mean      = a+b;
                props.variance  = b^2;
                props.skewness  = 2;
                props.kurtosis  = 6;
                props.mode      = a;
                props.median    = a+b*log(2);
                varargout{2} = props;
            end
        end
        function varargout = gammadist(x,p) 
        % GAMMADIST (Built-in Model) Gamma distribution
        %   Y=(OBJ/LINEFIT).GAMMADIST(X,[AMP,A,B,C])
        %       AMP: Amplitude, A: Location, B: Scale, C: Shape
        %
        %   Reference:
        %   (1) M.P. McLaughlin, Compendium of Common Probability
        %   Distributions, 2nd Edition,
        %   http://www.causascientia.org/math_stat/Dists/Compendium.pdf              
            amp    = p(1);
            a       = p(2);
            b       = p(3);
            c       = p(4);
            y = zeros(size(x));
            ind = (x>a);
            y(ind) = amp/(gamma(c)*b)*((x(ind)-a)/b).^(c-1).*exp(-(x(ind)-a)/b);
            varargout{1} = y;
            if nargout == 2
                props.amp      = amp;
                props.a         = a;
                props.b         = b;
                props.c         = c;
                props.mean      = a+b*c;
                props.variance  = b^2*c;
                props.skewness  = 2/sqrt(c);
                props.kurtosis  = 6/c;
                if c>=1
                    props.mode      = a+b*(c-1);
                else
                    props.mode      = a;
                end
                varargout{2} = props;
            end
        end
        function varargout = inversegamma(x,p) 
        % INVERSEGAMMA (Built-in Model) Inverse gamma distribution
        %   Y=(OBJ/LINEFIT).INVERSEGAMMA(X,[AMP,X0,A,B])
        %       AMP: Amplitude, X0: Location, A: Scale, B: Shape
        %
        %   Reference:
        %   (1) M.P. McLaughlin, Compendium of Common Probability
        %   Distributions, 2nd Edition,
        %   http://www.causascientia.org/math_stat/Dists/Compendium.pdf  
        %   (2) https://en.wikipedia.org/wiki/Inverse-gamma_distribution
            amp     = p(1);
            x0      = p(2);
            a       = p(3);
            b       = p(4);
            y = zeros(size(x));
            ind = (x>x0);            
            y(ind) = amp*a^b/gamma(b)*(x(ind)-x0).^(-b-1).*exp(-a./(x(ind)-x0));
            varargout{1} = y;
            if nargout == 2
                props.amp       = amp;
                props.x0        = x0;
                props.a         = a;
                props.b         = b;
                if b>1
                    props.mean      = x0+a/(b-1);
                else
                    props.mean  = NaN;
                end
                if b>2
                    props.variance  = a^2/((b-2)*(b-1)^2);
                else
                    props.variance  = NaN;
                end
                if b>3
                    props.skewness  = 4*sqrt(b-2)/(b-3);
                else
                    props.skewness  = NaN;
                end
                if b>4
                    props.kurtosis  = 3*(b+5)*(b-2)/((b-4)*(b-3))-3;
                else
                    props.kurtosis  = NaN;
                end
                props.mode      = x0+a/(b+1);
                varargout{2} = props;
            end
        end
        function varargout = lognorm(x,p)   
        % LOGNORM (Built-in Model) Log-normal distribution
        %   Y=(OBJ/LINEFIT).LOGNORM(X,[AMP,X0,A,B])
        %       AMP: Amplitude, X0: Location, A: Log-location, B: Scale, both in log space
        %
        %   Reference:
        %   (1) M.P. McLaughlin, Compendium of Common Probability
        %   Distributions, 2nd Edition,
        %   http://www.causascientia.org/math_stat/Dists/Compendium.pdf              
        %   (2) https://en.wikipedia.org/wiki/Log-normal_distribution
            amp    = p(1);
            x0      = p(2);
            a       = p(3);
            b       = p(4);
            y = zeros(size(x));
            ind = (x>x0);
            y(ind) = amp/sqrt(2*pi)/b./(x(ind)-x0).*exp(-(log(x(ind)-x0)-a).^2/2/b^2);
            varargout{1} = y;
            if nargout == 2
                props.amp      = amp;
                props.a         = a;
                props.b         = b;
                props.mean      = x0+exp(a+b^2/2);
                props.variance  = (exp(b^2)-1)*exp(2*a+b^2);
                props.skewness  = (exp(b^2)+2)*sqrt(exp(b^2)-1);
                props.kurtosis  = exp(4*b^2)+2*exp(3*b^2)+3*exp(2*b^2)-6;
                props.mode      = x0+exp(a-b^2);
                props.median    = x0+exp(a);
                varargout{2} = props;
            end
        end
        function varargout = logcauchy(x,p) 
        % LOGCAUCHY (Built-in Model) Log-Cauchy distribution
        %   Y=(OBJ/LINEFIT).LOGCAUCHY(X,[AMP,X0,A,B])
        %       AMP: Amplitude, X0: Location, A: Log-location, B: Scale
        %
        %   Reference:
        %   (1) https://en.wikipedia.org/wiki/Log-Cauchy_distribution           
            amp    = p(1);
            x0      = p(2);
            a       = p(3);
            b       = p(4);
            y = zeros(size(x));
            ind = (x>x0);
            y(ind) = amp/pi./(x(ind)-x0)*b./((log(x(ind)-x0)-a).^2+b^2);
            varargout{1} = y; 
            if nargout == 2
                props.amp      = amp;
                props.x0        = x0;
                props.a         = a;
                props.b         = b;
                props.median    = x0+exp(a);
                props.variance  = Inf;
                varargout{2} = props;
            end
        end
        function varargout = levy(x,p) 
        % LEVY (Built-in Model) Levy distribution
        %   Y=(OBJ/LINEFIT).LEVY(X,[AMP,A,B])
        %       AMP: Amplitude, A: Location, B: Scale
        %
        %   Reference:
        %   (1) M.P. McLaughlin, Compendium of Common Probability
        %   Distributions, 2nd Edition,
        %   http://www.causascientia.org/math_stat/Dists/Compendium.pdf 
        %   (2) https://en.wikipedia.org/wiki/L%C3%A9vy_distribution
            amp    = p(1);
            a       = p(2);
            b       = p(3);
            y = zeros(size(x));
            ind = (x>a);
            y(ind) = amp*sqrt(b/(2*pi))*(x(ind)-a).^(-1.5).*exp(-b./(2*(x(ind)-a)));
            varargout{1} = y;
            if nargout == 2
                props.amp      = amp;
                props.a         = a;
                props.b         = b;
                props.mean      = Inf;
                props.variance  = Inf;
                props.mode      = a+b/3;
                props.median    = a+b/2/erfinv(0.5)^2;
                varargout{2} = props;
            end
        end
        function varargout = weibull(x,p)
        % WEIBULL (Built-in Model) Levy distribution
        %   Y=(OBJ/LINEFIT).WEIBULL(X,[AMP,A,B,C])
        %       AMP: Amplitude, A: Location, B: Scale, C: Shape
        %
        %   Reference:
        %   (1) M.P. McLaughlin, Compendium of Common Probability
        %   Distributions, 2nd Edition,
        %   http://www.causascientia.org/math_stat/Dists/Compendium.pdf    
        %   (2) https://en.wikipedia.org/wiki/Weibull_distribution
            amp    = p(1);
            a       = p(2);
            b       = p(3);
            c       = p(4);
            y = zeros(size(x));
            ind = (x>=a);
            y(ind) = amp*(c/b)*((x(ind)-a)/b).^(c-1).*exp(-((x(ind)-a)/b).^c);
            varargout{1} = y;
            if nargout == 2
                G = gamma(1+(1:4)/c);
                props.amp      = amp;
                props.a         = a;
                props.b         = b;
                props.c         = c;
                props.mean      = a+b*G(1);
                props.variance  = (-G(1)^2+G(2))*b^2;
                props.skewness  = (2*G(1)^3-3*G(1)*G(2)+G(3))*b^3/props.variance^1.5;
                props.kurtosis  = (-3*G(1)^4+6*G(1)^2*G(2)-4*G(1)*G(3)+G(4))*b^4/props.variance^2-3;
                if c<=1
                    props.mode  = a;
                else
                    props.mode  = a+b*(1-1/c)^(1/c);
                end
                props.median    = a+b*(log(2))^(1/c);
                varargout{2} = props;
            end
        end
        function varargout = paretoI(x,p) 
        % PARETOI (Built-in Model) Pareto distribution type I
        %   Y=(OBJ/LINEFIT).PARETOI(X,[AMP,A,B])
        %       AMP: Amplitude, A: Location and scale, B: Shape
        %
        %   Reference:
        %   (1) M.P. McLaughlin, Compendium of Common Probability
        %   Distributions, 2nd Edition,
        %   http://www.causascientia.org/math_stat/Dists/Compendium.pdf               
        %   (2) https://en.wikipedia.org/wiki/Pareto_distribution
            amp    = p(1);
            a       = p(2);
            b       = p(3);
            y = zeros(size(x));
            ind = (x>=a);   
            y(ind) = amp*b*a^b./x(ind).^(b+1);
            varargout{1} = y;
            if nargout == 2
                props.amp      = amp;
                props.a         = a;
                props.b         = b;
                if b>1
                    props.mean      = a*b/(b-1);
                else
                    props.mean      = Inf;
                end
                if b>2
                    props.variance  = a^2*b/(b-1)^2/(b-2);
                elseif b<=0
                    props.variance  = NaN;
                else
                    props.variance  = Inf;
                end
                if b>3
                    props.skewness  = 2*(1+b)/(b-3)*sqrt((b-2)/b);
                else
                    props.skewness  = NaN;
                end
                if b>4
                    props.kurtosis  = 6*(b^3+b^2-6*b-2)/(b*(b-3)*(b-4));
                else
                    props.kurtosis  = NaN;
                end
                props.mode      = a;
                props.median    = a*2^(1/b);
                varargout{2} = props;
            end
        end
        function varargout = burrXII(x,p)     
        % BURRXII (Built-in Model) Burr distribution type XII
        %   Y=(OBJ/LINEFIT).PARETOI(X,[AMP,X0,A,B,C])
        %       AMP: Amplitude, X0: Location, A: Scale, B,C: Shape
        %
        %   Reference:
        %   (1) M.P. McLaughlin, Compendium of Common Probability
        %   Distributions, 2nd Edition,
        %   http://www.causascientia.org/math_stat/Dists/Compendium.pdf               
            amp    = p(1);
            x0      = p(2);
            a       = p(3);
            b       = p(4);
            c       = p(5);
            y = zeros(size(x));
            ind = (x>x0);
            y(ind) = amp*b*c/a*((x(ind)-x0)/a).^(b-1)./(1+((x(ind)-x0)/a).^b).^(c+1);
            varargout{1} = y;
            if nargout == 2
                mu = nan(1,4);
                beta_A = c-(1:4)/b;
                beta_B = 1+(1:4)/b;
                ind = find(beta_A>=0);
                mu(ind) = a.^ind*c.*beta(beta_A(ind),beta_B(ind));
                props.amp      = amp;
                props.x0       = x0;
                props.a         = a;
                props.b         = b;
                props.c         = c;
                props.mean          = x0+mu(1);
                props.variance      = -mu(1)^2+mu(2);
                props.skewness      = (2*mu(1)^3-3*mu(1)*mu(2)+mu(3))/props.variance^1.5;
                props.kurtosis      = (-3*mu(1)^4 + 6*mu(1)^2*mu(2) - 4*mu(1)*mu(3) + mu(4))/props.variance^2-3;
                props.mode          = x0+a*((b-1)/(b*c+1))^(1/b);
                props.median        = x0+a*(2^(1/c)-1)^(1/b);                
                varargout{2} = props;
            end
        end
        function varargout = laplace(x,p) 
        % LAPLACE (Built-in Model) Laplace distribution
        %   Y=(OBJ/LINEFIT).LAPLACE(X,[AMP,A,B])
        %       AMP: Amplitude, A: Location, B: Scale
        %
        %   Reference:
        %   (1) M.P. McLaughlin, Compendium of Common Probability
        %   Distributions, 2nd Edition,
        %   http://www.causascientia.org/math_stat/Dists/Compendium.pdf             
            amp    = p(1);
            a       = p(2);
            b       = p(3);
            varargout{1} = amp/(2*b)*exp(-abs(x-a)/b);
            if nargout == 2
                props.amp      = amp;
                props.a         = a;
                props.b         = b;
                props.mean          = a;
                props.variance      = 2*b^2;
                props.skewness      = 0;
                props.kurtosis      = 3;
                props.mode          = a;                
                props.median        = a;                
                varargout{2} = props;                
            end
        end
        function varargout = uniform(x,p)
        % UNIFORM (Built-in Model) Uniform or window distribution
        %   Y=(OBJ/LINEFIT).UNIFORM(X,[AMP,A,B])
        %       AMP: Amplitude, A: Location, B: Scale (upper bound)
        %
        %   Reference:
        %   (1) https://en.wikipedia.org/wiki/Uniform_distribution_(continuous)           
            amp    = p(1);
            a       = p(2);
            b       = p(3);
            y = zeros(size(x));
            y(x>=a & x<=b)  = amp/(b-a);
            varargout{1} = y;
            if nargout == 2
                props.amp      = amp;
                props.a         = a;
                props.b         = b;
                props.mean          = (a+b)/2;
                props.variance      = (b-a)^2/12;
                props.skewness      = 0;
                props.kurtosis      = -6/5;
                props.median        = (a+b)/2;          
                varargout{2} = props;                
            end            
        end
        function varargout = skewlogistic(x,p)
        % SKEWLOGISTIC (Built-in Model) Skew-logistic or generalized logistic type I distribution
        %   Y=(OBJ/LINEFIT).SKEWLOGISTIC(X,[AMP,A,B,C])
        %       AMP: Amplitude, A: Location, B: Scale, C: Shape
        %
        %   Reference:
        %   (1) M.P. McLaughlin, Compendium of Common Probability
        %   Distributions, 2nd Edition,
        %   http://www.causascientia.org/math_stat/Dists/Compendium.pdf  
        %   (2) https://en.wikipedia.org/wiki/Generalized_logistic_distribution
            amp    = p(1);
            a       = p(2);
            b       = p(3);
            c       = p(4);
            varargout{1} = amp*c/b*exp(-(x-a)/b)./(1+exp(-(x-a)/b)).^(c+1);
            if nargout == 2
                props.amp      = amp;
                props.a         = a;
                props.b         = b;
                props.c         = c;
                props.mean          = a+(double(eulergamma)+psi(c))*b;
                props.variance      = (pi^2/6+psi(1,c))*b^2;
                props.skewness      = (psi(2,c)-psi(2,1))*b^3/props.variance^1.5;
                props.kurtosis      = ((psi(3,c)+psi(1,c)*(pi^2+3*psi(1,c))+3*pi^4/20)*b^4)/props.variance^2-3;
                props.mode          = a+b*log(c);
                props.median        = a-b*log(2^(1/c)-1);
                varargout{2} = props;
            end                        
        end
        function varargout = skewlaplace(x,p)
        % SKEWLAPLACE (Built-in Model) Skew-Laplace distribution
        %   Y=(OBJ/LINEFIT).SKEWLAPLACE(X,[AMP,A,B,C])
        %       AMP: Amplitude, A: Location, B,C: Scale
        %
        %   Reference:
        %   (1) M.P. McLaughlin, Compendium of Common Probability
        %   Distributions, 2nd Edition,
        %   http://www.causascientia.org/math_stat/Dists/Compendium.pdf  
        %   (2) https://en.wikipedia.org/wiki/Asymmetric_Laplace_distribution
            amp     = p(1);
            a       = p(2);
            b       = p(3);
            c       = p(4);
            y       = zeros(size(x));
            ind1    = x<=a;
            ind2    = x>a;
            y(ind1) = exp((x(ind1)-a)/b);
            y(ind2) = exp((a-x(ind2))/c);
            y = amp/(b+c)*y;
            varargout{1} = y;
            if nargout == 2
                props.amp      = amp;
                props.a         = a;
                props.b         = b;
                props.c         = c;
                props.mean      = a-b+c;
                props.variance  = b^2+c^2;
                props.skewnewss = 2*(c^3-b^3)/props.variance^1.5;
                props.kurtosis  = (9*b^4+6*b^2*c^2+9*c^4)/props.variance^2-3;
                props.mode      = a;
                props.median    = a+b*log((b+c)/(2*b));
                varargout{2} = props;
            end                       
        end
        function varargout = erf(x,p)
        % ERF (Built-in Model) Overloaded error function.
        %   Y=(OBJ/LINEFIT).ERF(X,[AMP,A,B])
        %       AMP: Amplitude, A: Location, B: Scale
            amp    = p(1);
            a       = p(2);
            b       = p(3);
            y = amp*erf((x-a)/b);
            varargout{1} = y;
            if nargout == 2
                props.amp      = amp;
                props.a         = a;
                props.b         = b;
                varargout{2} = props;
            end                       
        end
        function varargout = atan(x,p)
        % ATAN (Built-in Model) Overloaded inverse tangent function.
        %   Y=(OBJ/LINEFIT).ATAN(X,[AMP,A,B])
        %       AMP: Amplitude, A: Location, B: Scale
            amp    = p(1);
            a       = p(2);
            b       = p(3);
            y = amp*atan((x-a)/b);
            varargout{1} = y;
            if nargout == 2
                props.amp      = amp;
                props.a         = a;
                props.b         = b;
                varargout{2} = props;
            end                       
        end
        function varargout = powerlaw(x,p)
        % POWERLAW (Built-in Model) Power law function.
        %   Y=(OBJ/LINEFIT).POWERLAW(X,[AMP,A,B])
        %       AMP: Amplitude, A: Location, B: Shape
            amp     = p(1);
            a       = p(2);
            b       = p(3);
            y       = zeros(size(x));
            ind     = x>=a;
            y(ind) = amp*(x(ind)-a).^(-b);
            varargout{1} = y;
            if nargout == 2
                props.amp       = amp;
                props.a         = a;
                props.b         = b;
                varargout{2} = props;
            end                       
        end        
        function newobj = loadobj(obj)
        % LOADOBJ Called by matlab 'load' function to load and construct
        %   LineFit object
            newobj = linefit;
            if ~isempty(fieldnames(obj.CustomCurveModel))
                newobj.CustomCurveModel = obj.CustomCurveModel;
            end
            if ~isempty(fieldnames(obj.CustomBkgdModel))
                newobj.CustomBkgdModel = obj.CustomBkgdModel;
            end            
            flist = {...
                'Data','CurveModelIndex','BkgdModelIndex','Model',...
                'FlagNonNegBkgd','FlagLogScale','FitOptions','FitOutput','PeaksFound'};
            for ii=1:length(flist)
                newobj.(flist{ii}) = obj.(flist{ii});
            end
        end
    end % End of static method
end % End of classdef

% --- core fitting function
function [ydata_cal,iydata_cal,ydata_bkgd] = fcn_fit_core(x1,xdata,varargin)
FitParamsSet        = varargin{1};
Model               = varargin{2};
FlagLogScale        = varargin{3};
FlagNonNegBkgd      = varargin{4};
% --- Generate all the fitting parameters
fitflag = FitParamsSet.fitflag;
x=zeros(1,length(fitflag));
x(fitflag==1) = x1;
x(fitflag==0) = FitParamsSet.x2;
x_curve = reshape(x(1:prod(Model.CurveModel.NOfParams)),Model.CurveModel.NOfParams(1),[]);
x_bkgd = x(end+1-Model.BkgdModel.NOfParams(1):end);
% --- calculate curve
iydata_cal = nan(size(xdata,1),size(x_curve,2));
for ii=1:size(x_curve,2) % for multiple model curves
     iydata_cal(:,ii) = Model.CurveModel.ModelFcnHandle(xdata,x_curve(:,ii));
end
ydata_curve = sum(iydata_cal,2);
% --- calculate bkgd
ydata_bkgd = Model.BkgdModel.ModelFcnHandle(xdata,x_bkgd);
if FlagNonNegBkgd == 1
    ydata_bkgd(ydata_bkgd<0) = 0;
end
% --- assemble findal result
ydata_cal = ydata_curve + ydata_bkgd;
% --- check log scale
if FlagLogScale == 1
    ydata_cal = log10(ydata_cal);
end
end