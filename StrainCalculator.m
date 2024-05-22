classdef StrainCalculator < handle
    % S = StrainCalculator()
    %
    % StrainCalculator calculates strains of a flat 2D regular grid 
    % undergoing a series of deformations, resulting in displacements of
    % the original grid. A typical application is the post-processing of
    % `Digital Image Correlation (DIC) <https://en.wikipedia.org/wiki/Digital_image_correlation_and_tracking>`_ 
    % methods. In brief, a virtual or physical grid is somehow `applied to 
    % the sample <https://digitalimagecorrelation.org/#patterning>`_, and 
    % the sample is then loaded in a given number of
    % loading steps. The position of each grid element is retrieved at
    % each step, and their displacement is calculated. For example:
    %
    % .. figure:: images/pointcloud.gif
    %    :scale: 60 %
    %    :alt: Series of grid deformations
    %
    %    An example of a planar grid undergoing a series of deformations.
    %    Source: Source: `µDIC <https://mudic.readthedocs.io/en/latest/index.html>`_ 
	% 
    % StrainCalculator() can use the displacements of the grid at each step
    % to compute a strain map.
    %
    % Context:
    %
    % Several toolboxes are available online for the processing of DIC data 
    % and calculation of strains. However, the strain calculation is often 
    % well hidden or not well documentent. The aim of this toolbox is to 
    % propose different methods to calculate strains, with clear 
    % documentation and examples.
    % This program is aimed at engineers, engineering students, 
    % but also researchers using DIC and wishing to 
    % better understand what is calculated. Its aim is mostly educational.
    %
	% Example:
	%
	%    .. code:: matlab 
	%
	%       obj = StrainCalculator();
	%       obj.generateFakeData('HorizontalCompression');
	%       obj.calculateStrains;
	%       
	%       for k = 1 : obj.nSteps
	%           cla
	%           obj.plotStrain(k);
	%           pause(0.3)
	%       end
	%
    % Strains should be provided as an X and Y array of cells, each
    % containing the X and Y positions of the grid elements at each time
    % step.
    % For more detailed instructions look at the :ref:`examples`.
    % 
    %
    % .. todo::
    %
    %    Add documentation for strain functions returning [X,Y,XY] (and 
    %    adapt variable names)
    %
    % .. todo::
    %
    %    Add parfor to documentation (as well as the procedure to choose
    %    and dispatch strain calculation method)
    %

    properties (SetAccess = protected , GetAccess = protected)

    end
    properties (SetAccess = private , GetAccess = public)
        nWorkers = [];     % Parallel workers
        strain = struct(); % strain structure
    end

    properties (SetAccess = private , GetAccess = public, Dependent = true)
        nSteps % number of steps
		
		% Displacements are calculated dynamically from the provided displacements.
		% The returned variable is a structure with an X and Y field. X and Y are 
		% arrays of cells, containing displacements for each step. Displacements are
		% calculated relative to the initial positions.
		% For each loading step *n*, we have:
		%
		% .. code:: matlab
		%
		%     displacement.X{n} = = Xpositions{n} - Xpositions{1};
		%     displacement.Y{n} = = Ypositions{n} - Ypositions{1};
		%
		% From a mathematical point of view, let :math:`\textbf{X}` denote the original locations
		% of a material particle. At time :math:`\textbf{t}` (or step *n*), the vector pointing
		% to the new position of a point in the body, from its original location, is
		% :math:`\textbf{u}(\textbf{X},t)`.
		% Since the original coordinates are the independent variables, this is 
		% a Lagrangian formulation. Thus, the displacement provides the 
		% transformation from the material to the spatial frame,
		% :math:`\textbf{x}=\textbf{X}+\textbf{u}`.
        displacements
		
		% The size of the initial grid in the horizontal direction (a scalar). 
		% This is dynamically calculated from the first step of the positions.
		% It also assumes that the initial spacing is regular (all elements have the same size)
        XgridSpacing;
		
        % The size of the initial grid in the vertical direction (a scalar). 
		YgridSpacing;
    end
    properties (SetAccess = public , GetAccess = public , Dependent = true)
    
    end    
    properties (Constant, GetAccess = public)
        
    end %properties
    properties (SetAccess = public)
        Xpositions = cell(0);
        Ypositions = cell(0);
    end
    properties (SetAccess = private , GetAccess = private)
        displacements_ = [];
    end    

    
    methods
        
        function nSteps = get.nSteps(this)
            nSteps = length(this.Xpositions);
        end
        
        function XgridSpacing = get.XgridSpacing(this)
            % Calculates spacing from the provided position, assuming
            % regular grid
            XgridSpacing = this.Xpositions{1}(1,2) - this.Xpositions{1}(1,1);
        end
        function YgridSpacing = get.YgridSpacing(this)
            % Calculates spacing from the provided position, assuming
            % regular grid
            YgridSpacing = this.Ypositions{1}(2,1) - this.Ypositions{1}(1,1);
        end        
        
        function displacements = get.displacements(this)
            % Calculates displacements on demand the first time, then
            % return the stored value
            if isempty(this.displacements_)
                this.calcDisplacements();
            end
            
            displacements = this.displacements_;
        end

        function generateFakeData(this, varargin)
            % generateFakeData()
            %
            % ...generates Fake Data.
            %
            % Optional parameters:
            %
            % 'type': 
            %    Type of strain configuration
            %       * HorizontalTension (Default)
            %       * HorizontalCompression
            %       * VerticalCompression
            %       * VerticalTension
            %       * HorizontalShear
            %
            % 'maxStrain': maximum value of strain
            %       Default: 0.1 (10%)
            % 
            % Example:
            %
            %    obj = StrainCalculator();
            %    obj.generateFakeData('HorizontalCompression');
            %    obj.calculateStrains;
            %
            
            
            p = inputParser;
            validTypes = {'HorizontalTension', 'HorizontalCompression',...
                          'VerticalCompression', 'VerticalTension', 'HorizontalShear', 'PureShear'};
            checkType = @(x) any(validatestring(x,validTypes));
            addOptional(p,'type','',checkType)
            
            addOptional(p,'maxStrain',0.1 , @isnumeric)

            parse(p,varargin{:})
            
            maxStrain = p.Results.maxStrain;
            type = p.Results.type;
            % Print valid strain types if none were provided
            if isempty(type)
                disp('These are the available strain types:')
                disp(validTypes')
                disp('Example: StrainCalculator.generateFakeData(''HorizontalTension''')
                return
            end


            % Tension or compression

            nElements = 20;
            nSteps = 10;
            this.Xpositions = cell(nSteps,1);
            this.Ypositions = cell(nSteps,1);
            XinitialPosition = repmat(linspace(0,1,nElements), nElements, 1);
            YinitialPosition = flipud(XinitialPosition');
            % 10 steps
            for step = 0 : nSteps
                this.Xpositions{step+1} = zeros(nElements, nElements);
                this.Ypositions{step+1} = zeros(nElements, nElements);

                % Define deformation gradient
                switch type
                    case 'HorizontalTension'
                        lambda = 1 + maxStrain/(nSteps) * (step);
                        F = [  lambda 0 ; ...
                                 0 1   ];
                    case 'HorizontalCompression'
                        % Just need to change the sign
                        lambda = 1 - maxStrain/(nSteps) * (step);
                        F = [  lambda 0 ; ...
                                 0 1   ];

                    case 'VerticalTension'
                        lambda = 1 + maxStrain/(nSteps) * (step);
                        F = [  1 0 ; ...
                               0 lambda   ];

                    case 'VerticalCompression'
                        % Switch X and Y, and change Y sign
                        lambda = 1 - maxStrain/(nSteps) * (step);
                        F = [  1 0 ; ...
                               0 lambda   ];

                    case 'HorizontalShear'
                        gamma = maxStrain/(nSteps) * (step);
                        F = [  1  gamma; ...
                               0  1   ];

                    case 'PureShear'
                        gamma = maxStrain/(nSteps) * (step);
                        F = [  1  gamma; ...
                             gamma  1   ];

                end

                % Apply deformation gradient to all points
                for k = 1 : numel(XinitialPosition)
                    pos = F * [XinitialPosition(k) ; YinitialPosition(k)];
                    this.Xpositions{step+1}(k) = pos(1);
                    this.Ypositions{step+1}(k) = pos(2);
                end
            end
        end
        
        
        function h = plotStrain(this, varargin)
            
            p = inputParser;
            checkStep = @(x) isnumeric(x) && isscalar(x) && (x > 0) && (x <= this.nSteps);
            addOptional(p,'step',1, checkStep)
            
            validStrains = fieldnames(this.strain);
            checkStrains = @(x) any(validatestring(x,validStrains));
            addOptional(p,'direction','X',checkStrains)
            
            addOptional(p,'axis',gca,@ishandle)

            parse(p,varargin{:})
            
            step = p.Results.step;
            direction = p.Results.direction;
            ax = p.Results.axis;
            
            % Plot as a surface
            flatZCoord = this.Xpositions{1} .* 0;
            h = surf(ax, this.Xpositions{step}, this.Ypositions{step}, ...
                     flatZCoord, this.strain.(direction){step}, 'EdgeColor', 'None');
            view(0, 90)
            axis equal
            
            % Normalize with max and min strain
            if ~all(diff([this.strain.limits.(direction)])) == 0
                % Only do it if strain is not constant everywhere,
                % otherwise caxis will fail
                caxis([this.strain.limits.(direction)]);
            end
            
            
        end
        
        function plotGrid(this, varargin)
            if isempty(varargin)
                step = 1;
            else
                step = varargin{1};
            end
            
            plot(this.Xpositions{step}(:) , this.Ypositions{step}(:), 'r+');
            hold on;
        end
        
        
        function calcDisplacements(this)
            % Calculate displacements
            N  = size(this.Xpositions{1});
            for k = 1 : this.nSteps
                
                % Displacements
                this.displacements_.X{k} = ...
                    reshape(this.Xpositions{k}(:) - this.Xpositions{1}(:) , N(1) , N(2));
                this.displacements_.Y{k} = ...
                    reshape(this.Ypositions{k}(:) - this.Ypositions{1}(:) , N(1) , N(2));

                % Change sign of the Y displacement so that + corresponds to upwards
                this.displacements_.Y{k} = -this.displacements_.Y{k};                            
            end          
        end
        


        
        function strain = calcStrain(this, varargin)
            % strain = calcStrain(Method)
            % Calculates strains from provided positions
            % Strain calculation method can be one of the following:
            %
            % - EngStrain  (or Engineering)
            % - LogStrain  (or Log)
            % - Cauchy
            % - Green
            % - Almansi
            % - Infinitesimal
            
            % Check optional method
            if isempty(varargin)
                Method = 'LogStrain';
            else
                validMethods = {'Engineering', 'LogStrain', 'Log', 'Infinitesimal', ...
                                'Green', 'Almansi', 'True'};
                checkMethod = any(validatestring(varargin{1},validMethods));
                if checkMethod
                    Method = varargin{1};
                else
                    error(['Strain calculation method should be one of the following: ' ...
                        strjoin(validMethods,', ') '. Input was: ' varargin{1}])
                end
            end
            
            % Compute displacement and prepare strain structure for output
            this.calcDisplacements();
            this.strain = struct();
            
            % According to the method, set strain.method, the strain 
            % calculation function and also pre-calculate the 
            % relevant tensor (displacement gradient, deformation
            % gradient, or stretch tensor)
            
            switch Method
                case 'EngStrain'
                    strainFunction = @(varargin) this.StrainEng(varargin{:});
                    this.strain.Method = 'Engineering strain';   
                    T = this.DisplacementGradientTensor;
                case 'Engineering'
                    strainFunction = @(varargin) this.StrainEng(varargin{:});
                    this.strain.Method = 'Engineering strain';                    
                    T = this.DisplacementGradientTensor;
                case 'Log'
                    strainFunction = @(varargin) this.StrainLog(varargin{:});
                    this.strain.Method = 'Logarihtmic strain';                    
                    T = this.StretchTensor;  % This is needed for log strain, not F
                case 'LogStrain'
                    strainFunction = @(varargin) this.StrainLog(varargin{:});
                    this.strain.Method = 'Logarihtmic strain';                    
                    T = this.StretchTensor; % This is needed for log strain, not F
                case 'Infinitesimal'
                    strainFunction = @(varargin) this.StrainInfinitesimal(varargin{:});
                    this.strain.Method = 'Infinitesimal strain'; 
                    T = this.DeformationGradientTensor();
                case 'Green'
                    strainFunction = @(varargin) this.StrainGreen(varargin{:});
                    this.strain.Method = 'Green strain';
                    T = this.DeformationGradientTensor();
                case 'Almansi'
                    strainFunction = @(varargin) this.StrainAlmansi(varargin{:});
                    this.strain.Method = 'Almansi strain';
                    T = this.DeformationGradientTensor();
                case 'True'
                    strainFunction = @(varargin) this.StrainTrue(varargin{:});
                    this.strain.Method = 'True strain';
                    T = this.DeformationGradientTensor();
                otherwise
                    error(['Method to calculate deformation not recognized: ' Method]);
            end    

            % Calculate strains for each step and grid location
            nrows = size(this.Xpositions{1},1);
            ncols = size(this.Xpositions{1},2);
            [i, j, k] = meshgrid(2:this.nSteps, 1:nrows, 1:ncols);
            
            % indices for step number, horizontal and vertical location
            % of grid elements
            indices = [i(:), j(:) , k(:)];
            
            Xmat = zeros(length(indices), 1);
            Ymat = Xmat;
            XYmat = Xmat;
            for idx = 1 : length(indices)

                Tmat = T{indices(idx,1), indices(idx,2), indices(idx,3)};
                [Xmat(idx), Ymat(idx), XYmat(idx)] = ...
                        strainFunction(Tmat);  % Apply the chosen strain           

            end
            
            this.strain.X = cell(this.nSteps,1);
            this.strain.Y = this.strain.X;
            this.strain.XY = this.strain.X;
            zeroMatrix = zeros(nrows, ncols);
            for step = 2 : this.nSteps
                where = indices(:,1) == step;
                idx = sub2ind([nrows, ncols], indices(where,2), indices(where,3));
                
                this.strain.X{step}  = zeroMatrix;
                this.strain.Y{step}  = zeroMatrix;
                this.strain.XY{step} = zeroMatrix;
                
                this.strain.X{step}(idx)  = Xmat(where);
                this.strain.Y{step}(idx)  = Ymat(where);
                this.strain.XY{step}(idx) = XYmat(where);
            end


            
            %%% Calculate max and min strains
            this.StrainLimits();

            %%% Calculate principal strains
            this.StrainPrincipal();
            
            if nargout > 0
                strain = this.strain;
            end
        end
        
        
        function strain = StrainPrincipal(this)
            % strain = StrainPrincipal()
            %
            % Computes principal strain from eigenvalues and eigenvectors
            % of the strain tensor. Values will depend on the method
            % chosen to compute strains.
            %
            % For each element *i* of the grid, 
            % 
            % .. math::
            %
            %   \textbf{E}^{(i)} = \begin{bmatrix}
            %    E_{xx}^{(i)} & E_{xy}^{(i)} \\
            %    E_{xy}^{(i)} & E_{yy}^{(i)}\\
            %    \end{bmatrix}{}
            %
            % We calculate the eigenvalues (which correspond to the first
            % and second principal strains) and the eigenvectors (which
            % correspond to the directions of the principal strains).
            %
            % The maximum shear is calculated as:
            %
            % .. math::
            %
            %   \gamma_{max} = 0.5*(\max{([e_I, e_{II}])} - \min{([e_I, e_{II}])})
            %
            % A simple definition of Maximum shear is in:
            %
            % .. admonition:: References
            %   :class: admonition note
            %
            %   Mapping 3D Strains with Ultrasound Speckle Tracking: Method Validation and 
            %   Initial Results in Porcine Scleral Inflation.
            %   Cruz Perez et al., Annals of Biomedical Engineering 2015
            
            sz = size(this.Xpositions{1});
            nElements = numel(this.Xpositions{1});
            
            % Temporary variables for I and II principal strain
            I  = zeros(1,nElements);   
            II = zeros(1,nElements);    
            % Temporary variable for I principal strain vector
            Ivect = zeros(2,nElements); 
            % Temporary variables for max shear
            MaxShear = zeros(1,nElements); 
            
            this.strain.PrincipalStrainI  = cell(this.nSteps,1);
            this.strain.PrincipalStrainII = cell(this.nSteps,1);
            this.strain.PrincipalStrainVectors =  cell(this.nSteps,1);
            this.strain.PrincipalStrainShear   = cell(this.nSteps,1);
            
            % For each step and each element...
            for nStep = 2 : this.nSteps
                for nElement = 1 : nElements
                    % Calculate eigenvalues and eigenvectors of the strain tensor
                    [V , D] = eig([this.strain.X{nStep}(nElement)  , this.strain.XY{nStep}(nElement) ;...
                                   this.strain.XY{nStep}(nElement) , this.strain.Y{nStep}(nElement)]);
                    % Store them
                    I(nElement)  = D(1,1);
                    II(nElement) = D(2,2);
                    Ivect(:,nElement) = V(:,1);
                    % Compute max shear
                    MaxShear(nElement) = 0.5 * (max([D(1,1) , D(2,2)])...
                                                -min([D(1,1) , D(2,2)]));
                end

                % Store values as matrices
                this.strain.PrincipalStrainI{nStep}  = reshape(I,sz);
                this.strain.PrincipalStrainII{nStep} = reshape(II,sz);
                this.strain.PrincipalStrainVectors{nStep} = reshape(Ivect,[2,sz]);
                this.strain.PrincipalStrainShear{nStep}   = reshape(MaxShear,sz);
        
            end
        
            if nargout > 0
                strain = this.strain;
            end
            
        end
        
        
        
        function [X, Y, XY] = StrainInfinitesimal(~, Fmat)
            % strain = StrainInfinitesimal()
            %
            % Also called linear strain tensor, 
            % or small strain tensor. It is equivalent to the engineering
            % strain.
            %
            % The definition of the strain tensor is 
            % :math:`\epsilon=\frac{1}{2}(\textbf{F}^{T}+\textbf{F})-\textbf{I}`, 
            % where H is grad of displacement (TODO). 
            % 
            % \epsilon_cauchy = \lambda - 1
            % Si \lambda - 1 < 0.05
            % TOTO Since it is limited to small def, it is quite
            % independent from rotation, so it could be used for a material
            % model law.
            %
            %
            % TODO: should put to NaN if lambda-1 > 0.05
            
            Edef = 1/2 * (Fmat' + Fmat) - eye(2);

            X = Edef(1,1);
            Y = Edef(2,2);
            XY = Edef(1,2);

        end
        
        function [X, Y, XY] = StrainTrue(~, Fmat)            
            % 
            %
            % TODO
            % I - (F'F)^-0.5
            %

            % True strain tensor
            Edef = eye(2,2) - (Fmat' * Fmat)^-0.5;

            X = Edef(1,1);
            Y = Edef(2,2);
            XY = Edef(1,2);
        
        end
        
        
        function [X, Y, XY] = StrainAlmansi(~, Fmat)
            % strain = StrainEulerian()
            %
            % TODO Called Almansi
            % Calculates Eulerian strains, also called Almansi strains or spatial strains. 
            %
            % Definition:
            %
            %    The left Cauchy-Green tensor or Finger tensor is \textbf{B}=\textbf{F}\textbf{F}^T)
            %    The Almansi strain tensor (or spatial strain tensor) is:
            %
            % .. math::
            %
            %     \textbf{e}^{*}=\frac{1}{2}(\textbf{I}-(\textbf{F}\textbf{F}^T)^{-1})
            %

            % Left Cauchy-Green 
            B = Fmat * Fmat';

            % Strain tensor
            e = 0.5*(eye(2,2) - inv(B));

            X = e(1,1);
            Y = e(2,2);
            XY = e(1,2);

        end

        function [X, Y, XY] = StrainGreen(~, Fmat)
            % strain = StrainGreen()
            %
            % Calculates right Green strains.
            %
            % Definition:
            %
            %    The right Cauchy-Green deformation tensor can be
            %    calculated in terms of the deformation gradient **F** as
            %    :math:`\textbf{C}=\textbf{F}^{T}\textbf{F}`. The Green 
            %    strain tensor **E** (Lagrangian tensor, or material tensor) is:
            %
            % .. math::
            %
            %     \textbf{E}=1/2(\textbf{C}-\textbf{I})=1/2(\textbf{F}^{T}\textbf{F}-\textbf{I})
            %
            % Where :math:`\textbf{F}` is the 
            % :meth:`deformation gradient
            % <StrainCalculator.StrainCalculator.DeformationGradient>`.
            %

            % Cauchy-Green tensor
            C = Fmat' * Fmat;

            % Strain tensor Green-Lagrange
            % TODO: add green-Lagrange to doc
            Edef = 0.5*(C-eye(2));

            X = Edef(1,1);
            Y = Edef(2,2);
            XY = Edef(1,2);

        end


        function [X, Y, XY] = StrainEng(~, Hmat)
            % strain = calStraincEng()
            %
            % Calculates engineering strains for a given displacement graident 
            % tensor. Strains are directly updated as a property of the object (this.strain).
            %
            % Definition:
            %
            %   .. list-table:: Definitions
            %      :header-rows: 0
            %
            %      * - :math:`e_{x}=\frac{\partial u}{\partial x}`
            %        - Horizontal strain
            %      * - :math:`e_{y}=\frac{\partial v}{\partial y}`
            %        - Vertical strain
            %      * - :math:`e_{xy}=\frac{1}{2}(\frac{\partial u}{\partial y} + \frac{\partial v}{\partial x})`
            %        - Shear strain         
            %
            % Engineering strain is also known as Eulerian infinitesimal
            % strains.
            % It can be decomposed in a strain matrix and rotation matrix:
            %
            % .. math::
            %
            %    \textrm{strain matrix}=\begin{bmatrix}
            %    \frac{\partial{u}}{\partial{x}} & \frac{1}{2}(\frac{\partial{u}}{\partial{y}} + \frac{\partial{v}}{\partial{x}} )  \\
            %    \frac{1}{2}(\frac{\partial{u}}{\partial{y}} +\frac{\partial{v}}{\partial{x}}) & \frac{\partial{v}}{\partial{y}} 
            %    \end{bmatrix}
            %
            %    \textrm{rotation matrix}=\begin{bmatrix}
            %    0 & \frac{1}{2}(\frac{\partial{u}}{\partial{y}} - \frac{\partial{v}}{\partial{x}} )  \\
            %    -\frac{1}{2}(\frac{\partial{u}}{\partial{y}} +\frac{\partial{v}}{\partial{x}}) & 0 
            %    \end{bmatrix}
            %    
            %

            du_DX = Hmat(1,1);
            dv_DY = Hmat(2,2);
            dv_DX = Hmat(2,1);
            du_DY = Hmat(1,2);
            X = du_DX;
            Y = dv_DY;
            XY = 0.5 * (du_DY + dv_DX);
        end


        function [U, F] = StretchTensor(this)
            % [U, F] = StretchTensor()
            %
            % Returns:
            % 
            %    - *U*, a matrix of cells where U{i, j, k} is the 2x2 stretch  
            %       tensor calculated at the ith step and at the jth row and 
            %       kth column element.
            %
            %    - *F*, a matrix of cells where F{i, j, k} is the 2x2 deformation  
            %       tensor calculated at the ith step and at the jth row and 
            %       kth column element.
            %
            % The stretch tensor is defined by decomposing the deformation
            % gradient in a rotation tensor *R* and a stretch tensor *U*,
            % as :math:`F=RU`.
            % From the right Cauchy-Green tensor:
            % 
            % .. math::
            % 
            %    C = F^TF=U^TR^TRU
            %
            % Choosing U to be symmetric, it results :math:`C=U^2`.
            % Polar decomposition can be used to calculate the square root
            % of U (i.e., of a symmetric tensor). First we express the
            % tensor in its principal directions (so it will be diagonal),
            % then we take the square root of the diagonal values, and
            % finally rotate the tensor back in its initial orientation.
            %
            % In this implementation, we calculate the eigenvalues and
            % eigenvectors of U^2, square the eigenvalues and replace them
            % in the original orientation via:
            %
            % .. math::
            %
            %   C = \sum_{\alpha}\lambda_{\alpha}^2 N_{\alpha}\otimes  N_{\alpha}
            %
            %   U = \sum_{\alpha}\lambda_{\alpha} N_{\alpha}\otimes  N_{\alpha}
            % 
            % where :math:`\lambda_{\alpha}` is the square root of the
            % :math:`\alpha{th}` eigenvalue and :math:`N_{\alpha}` is the 
            % corresponding eigenvector.
            %
            % The function also returns *F* for convenience
            %


            F = this.DeformationGradientTensor();

            % U is a 3D matrix of cells, with a 2x2 tensor in each cell.
            U = cell(this.nSteps, size(F, 2), size(F,3));

            for step = 2 : this.nSteps
                for row = 1 : size(F,2)
                    for column = 1 : size(F,3)
                        % Deformation tensor for one element
                        Fmat = F{step, row, column};

                        % Cauchy-Green tensor
                        C = Fmat' * Fmat;
                        % Eigenvalues and eigenvectores
                        [eigVectors, eigValues] = eig(C);

                        % Square root of the eigenvectors
                        lambda1 = sqrt(eigValues(1,1));
                        lambda2 = sqrt(eigValues(2,2));
                        V1 = eigVectors(:,1);
                        V2 = eigVectors(:,2);

                        % [U,S,V] = svd(Fmat);
                        % TODO
                        % S = Stretch Tensor in the principal direcitons;
                        % V = transpose of rotation matrix
                        % V could be transposed....


                        % Calculate U for one element
                        U(step,row,column) = {lambda1*(V1*V1') + lambda2*(V2*V2')};
                    end
                end
            end

        end


        function [R, F] = RotationTensor(this)
            % [R, F] = RotationTensor(this)
            %
            % Returns:
            % 
            %    - *R*, a matrix of cells where R{i, j, k} is the 2x2 rotation  
            %       tensor calculated at the ith step and at the jth row and 
            %       kth column element.
            %
            %    - *F*, a matrix of cells where F{i, j, k} is the 2x2 deformation  
            %       tensor calculated at the ith step and at the jth row and 
            %       kth column element.
            %
            % The :meth:`~StrainCalculator.StrainCalculator.DeformationGradient`
            % can be decomposed into stretch and rotation components, and
            % therefore F can be expressed as :math:`F = RU`, where R is
            % the rotation tensor and U is the :meth:`~StrainCalculator.StrainCalculator.StretchTensor`.
            % Once the stretch tensor U is known, R can be evaluated as 
            % :math:`R = FU^{-1}`
            %
            % This method is not used in the toolbox, but it is kept 
            % for completeness


            [U, F] = this.StretchTensor();

            % R is a 3D matrix of cells, with a 2x2 tensor in each cell.
            R = cell(this.nSteps, size(F,2), size(F,3));

            for step = 2 : this.nSteps
                for row = 1 : size(F,2)
                    for column = 1 : size(F,3)
                        % R = F*inv(U)
                        R{step, row, column} = ...
                            F{step, row, column} / U{step, row, column};

                    end
                end
            end

        end
        
    function [X, Y, XY] = StrainLog(~, Umat)
        % strain = StrainLog()
        %
        % Calculates logarithmic strains. Log strain (or Hencky strain) 
        % is defined as :math:`H = ln(U)`, where *U* is the 
        % :meth:`~StrainCalculator.StrainCalculator.StretchTensor`.
        % 
        % In practice, a simple approach is to calculate *U* in the 
        % principal directions, and calculate the log in this reference:
        % 
        % .. math::
        %
        %    \epsilon_P=\begin{bmatrix}
        %    \ln{(\lambda_1)} & 0  \\
        %    0 & \ln{(\lambda_2)}
        %    \end{bmatrix}
        %
        % And then transform it back in the original coordinates. 
        %
        % This could be further simplified because :math:`C = U^T U`, 
        % this definition can be used:
        %
        % .. math::
        %
        %   H = \ln(U) = \ln(\sqrt{C}) = 1/2 \ln(C)
        %
        % In uniaxial or biaxial tests, where principal directions can
        % be assumed to be parallel to the grid; the calculation can be 
        % done directly:
        %
        %   .. list-table:: Definitions
        %      :header-rows: 0
        %
        %      * - :math:`\epsilon_{x}=\ln(1+\frac{\partial u}{\partial x})`
        %        - Horizontal strain
        %      * - :math:`\epsilon_{y}=\ln(1+\frac{\partial v}{\partial y})`
        %        - Vertical strain  
        %
        % The unidimensional case can be calculated by integrating the 
        % incremental strain :math:`\delta\epsilon=\delta l/ l`:
        %
        % .. math::
        %
        %   \epsilon=\int_{}^{}\delta l/ l=\ln(l/l_{0})=\ln(1+\Delta l/l_{0})=\ln(1+e)
        % 
        % Besides, the log strain tends to lose meaning when rotations
        % are present. Just like in engineering strain, imagine a 45°
        % rotation of an object, followed by a stretch along what was
        % the x direction (and is now at 45°): this will appear as
        % shear strain! This can be mitigated by calculating the log
        % strains in the principal directions, and then rotating back
        % to the original coordinate system.
        %
        % Keep in mind:
        %
        %    While logarithmic measures of strain are a favorite in one-dimensional
        %    or semi-qualitative treatment, they have never been successfully applied in
        %    general. Such simplicity for certain problems as may result from a particular
        %    strain measure is bought at the cost of complexity for other problems
        %
        %    -- Truesdell, Toupin: The Classical Field Theories
        %
        % Furthermore, log strains might not suitable to evaluate large
        % shear strains; its application is still controversial. See:
        %
        % .. admonition:: References
        %   :class: admonition note
        %
        %   Jonas et al., 2011. Problems with Using the Hencky Equivalent 
        %   Strain in Simple Shear. Materials transactions 52(9).
        %   `DOI <http://dx.doi.org/10.2320/matertrans.M2011086>`__
        %
        %   Onaka, 2015. Comment on “A comparison of the von
        %   Mises and Hencky equivalent strains for
        %   use in simple shear experiments”. Philosophical Magazine
        %   92(18): 2264-2271.
        %   `DOI <http://dx.doi.org/10.1080/14786435.2012.671551>`__.
        %
        % .. warning::
        % 
        %    Can I still use the elements outside the diagonal for
        %    strain?

        % Transform in principal directions
        % Eigenvalues and eigenvectores
        [eigVectors, eigValues] = eig(Umat);

        % Natural logarithm of the eigenvectors
        lambda1 = log(eigValues(1,1));
        lambda2 = log(eigValues(2,2));
        V1 = eigVectors(:,1);
        V2 = eigVectors(:,2);

        % H in original coordinates
        H = lambda1*(V1*V1') + lambda2*(V2*V2');

        X = H(1,1);
        Y = H(2,2);
        XY = H(1,2);
        
    end        
        
    end  % Methods
    
    
    
    
    
    methods (Access = private)
        
        function checkInput(this)
            % Checks the validity of the properties

            if isempty(this.Xpositions) || isempty(this.Ypositions)
                error('Please input Xpositions and Ypositions')
            end
            
            if length(this.Xpositions) ~= length(this.Ypositions)
                error('Xpositions and Ypositions should have the same length')
            end
            
            % Check that positions are numbers
            for k = 1 : length(this.Xpositions)
                if ~isnumeric(this.Xpositions{k})
                    error(['Xpositions should be an array of cells '..., 
                           'containing numbers only. Xpositions contains ',...
                           class(this.Xpositions{k})]);
                end
                if ~isnumeric(this.Ypositions{k})
                    error(['Ypositions should be an array of cells '..., 
                           'containing numbers only. Ypositions contains ',...
                           class(this.Ypositions{k})]);
                end                    
            end
            
            
        end
        
        
        function H = DisplacementGradientTensor(this)
            % H = DisplacementGradientTensor()
            %
            % .. Displacement Gradient Tensor:
            %
            % In matrix form:
            % 
            % .. math::
            %
            %   \textbf{H} = \begin{bmatrix}
            %    \frac{\partial u}{\partial X} & \frac{\partial u}{\partial Y} \\
            %    \frac{\partial v}{\partial X} & \frac{\partial v}{\partial Y}\\
            %    \end{bmatrix}{}
            %
            
            H = cell(this.nSteps, 1);
            for nStep = 2 : this.nSteps
                [du_DX, du_DY] = gradient(this.displacements.X{nStep}, ...
                                          this.XgridSpacing , this.YgridSpacing);
                [dv_DX, dv_DY] = gradient(this.displacements.Y{nStep} , ...
                                          this.XgridSpacing , this.YgridSpacing);

                % Store an H matrix for each element of the grid
                for row = 1 : size(du_DX,1)
                    for column = 1 : size(du_DX,2)
                        % Displacement gradient tensor H for one
                        % element.
                        % Y has a  minus sign to place origin bottom left, not
                        % top left as in Matlab's images.
                        H{nStep, row, column} = [du_DX(row,column), -du_DY(row,column);...
                                                 dv_DX(row,column), -dv_DY(row,column)];
                        
                        % Fill in the first matrix with NaN values
                        if nStep == 2
                            H{1, row, column} = nan(2,2);
                        end
                    end
                end
            end
            
        end
        
        
        function F = DeformationGradientTensor(this)
            % F = DeformationGradient()
            %
            % .. _Deformation Gradient:
            %
            % Calculates the deformation gradient F
            %
            % Definition;
            %    
            % .. math::
            %
            %   \textbf{F}=\frac{\textbf{x}}{\partial{\textbf{X}}}=
            %   \textbf{I}+\frac{\partial{\textbf{u}}}{\partial{\textbf{X}}}=
            %   \textbf{I}+\textbf{H}
            %
            % Where **H** is the displacement gradient tensor.
            % In matrix form **F** is:
            % 
            % .. math::
            %
            %   \textbf{F} = \begin{bmatrix}
            %    \frac{\partial x}{\partial X} & \frac{\partial x}{\partial Y} \\
            %    \frac{\partial y}{\partial X} & \frac{\partial y}{\partial Y}\\
            %    \end{bmatrix}{}
            %
            % .. note::
            %    F in general is not symmetric!
            %
            % In this code, the calculation of the gradient is implemented
            % using the `gradient()` function. An alternative would be
            % using an approach similar to the finite element method for
            % each element of the grid. This would have the advantage of
            % allowing to use higher-order polynomials to interpolate the
            % displacements. This approach would also be useful is the grid
            % was not horizontal/vertical, but at an angle.
            %

            H = this.DisplacementGradientTensor();
            F = cell(this.nSteps, 1);
            I = eye(2,2);
            % Store an F matrix for each loading step
            for nStep = 2 : this.nSteps
                % Store an F matrix for each element of the grid
                for row = 1 : size(H,2)
                    for column = 1 : size(H,3)
                        % Deformation gradient for one element
                        F{nStep, row, column} = H{nStep, row, column} + I;
                    end
                end
            end                
        end
        
        
        function StrainLimits(this)
            % Calculate max and min strain in each available direction
            
            fields = fieldnames(this.strain);
            fields(contains(fields, 'Method')) = [];
            for type = 1:length(fields)
                limits = [NaN, NaN];
                for step = 1 : length(this.strain.(fields{type}))
                    data = this.strain.(fields{type}){step}(:);
                    limits(1) = nanmin([limits(1), min(data)]);
                    limits(2) = nanmax([limits(2), max(data)]);
                end
                
                this.strain.limits.(fields{type}) = limits;
            end
            
            
        end
        

    end % Private Methods
    
end