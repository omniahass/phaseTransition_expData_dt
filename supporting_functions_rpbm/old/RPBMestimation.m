classdef RPBMestimation
    % RPBMestimation: RPBM class
    %
    % =====================================================================
    %
    %  Authors: Santiago Coelho (santiago.coelho@nyulangone.org)
    %  Copyright (c) 2026 New York University
    %              
    %   Permission is hereby granted, free of charge, to any non-commercial entity ('Recipient') obtaining a 
    %   copy of this software and associated documentation files (the 'Software'), to the Software solely for
    %   non-commercial research, including the rights to use, copy and modify the Software, subject to the 
    %   following conditions: 
    % 
    %     1. The above copyright notice and this permission notice shall be included by Recipient in all copies
    %     or substantial portions of the Software. 
    % 
    %     2. THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
    %     NOT LIMITED TO THE WARRANTIESOF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
    %     IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BELIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
    %     WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF ORIN CONNECTION WITH THE
    %     SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
    % 
    %     3. In no event shall NYU be liable for direct, indirect, special, incidental or consequential damages
    %     in connection with the Software. Recipient will defend, indemnify and hold NYU harmless from any 
    %     claims or liability resulting from the use of the Software by recipient. 
    % 
    %     4. Neither anything contained herein nor the delivery of the Software to recipient shall be deemed to
    %     grant the Recipient any right or licenses under any patents or patent application owned by NYU. 
    % 
    %     5. The Software may only be used for non-commercial research and may not be used for clinical care. 
    % 
    %     6. Any publication by Recipient of research involving the Software shall cite the references listed
    %     below.
    %
    %  REFERENCES:
    %  - Coelho, S., Fieremans, E., Novikov, D.S., 2026 (ArXiv)
    %

    methods     ( Static = true )
        % =================================================================
        % =================================================================
        % =================================================================
        % =================================================================
        % =================================================================
        % =================================================================
        % =================================================================
        % =================================================================
        % =================================================================
        % =================================================================
        % =================================================================
        % =================================================================
        % =================================================================
        function flag_outliers = detectVectorOutliers(x,threshold)
        % flag_outliers = detectVectorOutliers(x,threshold)
        %
        % if Gaussian dist, MAD*1.4826 = standard deviation
        % A good heuristic is threshold = 3 * 1.4826
        
        if ~exist('threshold', 'var') || isempty(threshold)
            threshold = 3;
        end
        
        x = x(:);  % ensure column vector
        
        % Compute median and MAD
        med = median(x);
        mad_val = mad(x, 1);  % '1' uses normalization by median, not mean
        
        % Identify non-outliers
        inliers = abs(x - med) < threshold * mad_val;
        
        % Filter out outliers
        flag_outliers = ~inliers;
        
        end
        % =================================================================
        function Dt = get_Dt_RPBM_vectorized(t, zeta)
        % Vectorized version of get_Dt_RPBM for paired (t, zeta) inputs.
        % t    - vector of size [N x 1] or [1 x N]
        % zeta - vector of size [N x 1] or [1 x N], same length as t
        % Dt   - vector of size [N x 1]
        
            % Flatten to column vectors [N x 1]
            sz   = size(t);
            t    = t(:);
            zeta = zeta(:);
        
            % Precompute per-pair scalars — all [N x 1]
            lambda = 1 + zeta;
            sq     = sqrt(lambda);
            A      = 2*(sq - 1) ./ lambda.^2;
            B      = 2*(sq - 1) .* (sq - 2) ./ lambda.^3;
            C      = (8*(sq - 1).^2 - zeta.*sq) ./ lambda.^4;
        
            % Integrand: y is scalar, all others are [N x 1] -> output is [N x 1]
            function val = integrand(y)
                sq_y  = sqrt(y);
                isy   = 1i * sq_y;
                inner = isy - 1;
                core  = imag(1 ./ (1 + zeta + 2*isy .* inner ...
                            .* (1 - sqrt(1 + zeta ./ inner.^2))));  % [N x 1]
                smooth = A * sq_y + C * y^(3/2);                    % [N x 1]
                val    = (smooth + core) .* exp(-t * y) / y^2;      % [N x 1]
            end
        
            result = integral(@integrand, 0, Inf, ...
                              'ArrayValued', true, ...
                              'RelTol', 1e-6, 'AbsTol', 1e-10);     % [N x 1]
        
            Dt = result ./ (pi * t) ...
               + 1 ./ lambda      ...
               + 2*A ./ sqrt(pi * t) ...
               + B ./ t            ...
               - C .* t.^(-3/2) / sqrt(pi);
        
            % Restore original shape
            Dt = reshape(Dt, sz);
        end
        % =================================================================
        function Dt=get_Dt_RPBM(t,zeta)
        % Finds D(t/tau) from D(omega) from Nature Physics 7:508 (2011), Eq (6)
        % by deforming the contour of integration from the Re omega axis 
        % to the one around the negative imaginary axis in the lower half plane.
        % A few terms are subtracted so that the integrand is regular at zero.
        %
        % (c) Dmitry S Novikov 2015 dima@alum.mit.edu
        
        lambda=1+zeta; % tortuosity
        A=2*(sqrt(lambda)-1)/lambda^2;
        B=2*(sqrt(lambda)-1)*(sqrt(lambda)-2)/lambda^3;
        C=(8*(sqrt(lambda)-1)^2-zeta*sqrt(lambda))/lambda^4;
        
        fint2=inline('(A*sqrt(y)+C*y.^(3/2)+imag(1./(1+zeta+2*i*sqrt(y).*(i*sqrt(y)-1).*(1 - sqrt(1 + zeta./(i*sqrt(y)-1).^2))))).*exp(-t*y)./y.^2', 'y', 't', 'zeta', 'A', 'C');
        
        %ymin=1e-10; ymax=1e10; tol=1e-10;
        Dt=integral(@(y)fint2(y,t,zeta,A,C), 0, Inf, 'ArrayValued', true)./(pi*t) + 1/lambda + 2*A./sqrt(pi*t) + B./t - C*t.^(-3/2)/sqrt(pi);
        end
        % =================================================================
        function [D0,a,kappa,tau,zeta] = sample_rpbm_param(bounds,M,diam_factor)
        % [D0,a,kappa,tau,zeta] = sample_rpbm_param(bounds,MN,diam_factor)
        % bounds is [3 x 2] lower + upper bounds for: D0, a, kappa
        %
        %
        % By: Santiago Coelho (06/03/2026)
        if ~exist('diam_factor', 'var') || isempty(diam_factor)
            diam_factor = 4;
        end
        D0 = rand(M,1) * ( bounds(1,2) - bounds(1,1) ) + bounds(1,1);       % [  1  -   3] D0 in um^2/ms
        a = rand(M,1) * ( bounds(2,2) - bounds(2,1) ) + bounds(2,1);     % [0.5  -   3] zeta dimensionless
        kappa = rand(M,1) * ( bounds(3,2) - bounds(3,1) ) + bounds(3,1); % [0.01 - 0.1] kappa in um/ms
        tau = D0./(2*kappa).^2;
        d = 2;
        zeta = (diam_factor * D0) ./ (a .* 2 .* kappa .* d);
        end
        % =================================================================
        function [fits] = noise_propagation(Dt_test,Dt_training,param_training,PRdegree)
        % [fits] = noise_propagation(Dt_test,Dt_training,param_training)
        %
        % Dt_training is [DIM x TrainingSize]
        % Dt_test is [1 x TrainingSize] or [DIM_output x TrainingSize]
        % param_training = [zeta ; kappa ; a ]  (3 x N);
        % (c) Santiago Coelho 2026 - santiago.Coelho@nyulangone.org

        if ~exist('PRdegree', 'var') || isempty(PRdegree)
            PRdegree = 3;
        end
        options.training_data = Dt_training;
        options.training_objective = param_training;
        options.test_data = Dt_test;
        % PR
        options.Degree = PRdegree;
        options.strategy = 'PR';
        % options.svdFLAG=1;
        % options.keep_sv=10;
        
        out = RPBMestimation.DataDriven_regression(options);
        fits.zeta = out.fits(:,1);
        fits.kappa = out.fits(:,2);
        fits.a = out.fits(:,3);
        
        end
        % =================================================================
        function out = DataDriven_regression(options)
        % [NET_TRAINED,TR] = DataDriven_regression(training_data,training_objective,Nneurons,NhiddenLayers,Nrounds)
        %
        % training_data is [DIM x TrainingSize]
        % training_objective is [1 x TrainingSize] or [DIM_output x TrainingSize]
        %
        % options.strategy = 'NN' or 'PR'
        %
        %
        %
        % Nneurons is the # of neurons per layer (NhiddenLayers is # of layers) 
        % Nrounds is the amount of times the NN is trained (random initializations)
        %
        % By: Santiago Coelho (26/02/2020) - santiago.Coelho@nyulangone.org

        training_data=options.training_data;
        training_objective=options.training_objective;
        test_data=options.test_data;
        DIM=size(training_data,1);
        if strcmp(options.strategy,'NN')||strcmp(options.strategy,'nn')
            if isfield(options, 'Nneurons')
                Nneurons=options.Nneurons;
            else
                Nneurons=ceil(1.5*DIM);
            end
            if isfield(options, 'NhiddenLayers')
                NhiddenLayers=options.NhiddenLayers;
            else
                NhiddenLayers=3;
            end
            if isfield(options, 'Nrounds')
                Nrounds=options.Nrounds;
            else
                Nrounds=5;
            end
            if isfield(options, 'useGPU')
                useGPU=options.useGPU;
            else
                useGPU='no';
            end
            trainFcn = 'trainscg';  % Levenberg-Marquardt backpropagation.
            % Create a Fitting Network
            hiddenLayerSize = Nneurons*ones(1,NhiddenLayers);
            net = fitnet(hiddenLayerSize,trainFcn);
            % Setup Division of Data for Training, Validation, Testing
            net.divideParam.trainRatio = 70/100;
            net.divideParam.valRatio = 15/100;
            net.divideParam.testRatio = 15/100;
            net.trainParam.showWindow = false;
            for round_id=1:Nrounds
                % =========================================================================
                % Train the Network
                [NET_trained{round_id},tr{round_id}] = train(net,training_data,training_objective,'useGPU',useGPU,'showResources','no');
                RMSE_best_perf(round_id)=tr{round_id}.best_tperf;
                % =========================================================================
            end
            [~,id_best_net]=min(RMSE_best_perf);
            out.NET_TRAINED=NET_trained{id_best_net};
            out.TR=tr{id_best_net};
            out.fits = out.NET_TRAINED(test_data,'useGPU',useGPU,'showResources','no');
            out.RMSE_training=RMSE_best_perf(id_best_net);
        elseif strcmp(options.strategy,'PR')||strcmp(options.strategy,'pr')
            if isfield(options, 'Degree')
                Degree=options.Degree;
            else
                Degree=2;
            end
            if isfield(options, 'svdFLAG')
                svdFLAG=options.svdFLAG;
                keep_sv=options.keep_sv;
            else
                svdFLAG=0;
            end
            if svdFLAG&&(DIM>1)
                % =========================================================================
                % Applying PCA for dimensionality reduction
                [coeff, score] = pca(training_data');
                reducedDimension = coeff(:,1:keep_sv);
                data_training = (training_data' * reducedDimension)';
                % =========================================================================
                % Applying training transformation to test dataset (optimal)
                data_testing = (test_data' * reducedDimension)';
            else
                data_training=training_data;
                data_testing=test_data;
            end
            % =========================================================================
            % Polynomial Regression FITTING
            % =========================================================================
            X_moments = RPBMestimation.Compute_extended_moments(data_training',Degree);
            X_moments_test = RPBMestimation.Compute_extended_moments(data_testing',Degree);        
            pinvX=pinv(X_moments);
            coeffs=pinvX*(training_objective');   
            training_obj=X_moments*coeffs;
            clearvars pinvX
            % Applying PR
            out.fits=X_moments_test*coeffs;
            out.RMSE_training=sqrt(mean((training_obj-training_objective').^2));
        end
        out.strategy=options.strategy;
        
        end
        % =================================================================
        function Moments_Extended = Compute_extended_moments(Moments,degree)
        % Moments_Extended = Compute_extended_moments(Moments,degree)
        %
        % This function computes the extended matrix of moments for polynomial
        % interpolation in multiple dimensions.
        %
        % Moments is [NxD] where N is the number of samples and D the dimension of
        % the problem (the number of different moments)
        %
        % e.g. in 1D this extended matrix is [1 X X^2 ... X^degree]
        %
        %
        % (c) Santiago Coelho 2026 - santiago.Coelho@nyulangone.org
        N=size(Moments,1);
        N_moments=size(Moments,2);
        % Preparing data points for polynomial fitting
        x_0=ones(N,1);
        if degree>=1
            x_1=Moments;
            X=[x_0 x_1];
        else
                X=x_0;
        end
        if degree>=2
            for current_degree=2:degree
                combinations=nchoosek(1:(N_moments+current_degree-1),current_degree);
                NtotalCombs=size(combinations,1);
                flags=zeros(NtotalCombs,N_moments+current_degree-1);
                sums=zeros(NtotalCombs,N_moments+current_degree-1);
                which=zeros(NtotalCombs,current_degree);
                for ii=1:NtotalCombs
                    flags(ii,combinations(ii,:))=1;
                    sums(ii,:)=cumsum(~flags(ii,:))+1;
                    which(ii,:)=sums(ii,combinations(ii,:));
                end
                Currentx=ones(N,NtotalCombs);
                for ii=1:NtotalCombs
                    for jj=1:current_degree
                        current_id=which(ii,jj);
                        Currentx(:,ii)=Currentx(:,ii).*Moments(:,current_id);
                    end
                end
                X=[X Currentx];
            end
        end
        Moments_Extended=X;
        end
        % =================================================================
        function [param] = fit4Ddata(Dperp,t,sigma,D0,mask,PRoptions)
        % [param] = fit4Ddata(Dperp,t,sigma,D0,mask,PRoptions)
        %
        % Dperp is a 4D array of size [ x y z Ntimes ] - Dperpendicular in [(\mu m)^2/ms]
        % t is [ 1 Ntimes] where Ntimes is the # of diffusion times in [ms]
        % sigma is the noise level of those Dperp estimates 
        % D0 is the value in which we are fixing our estimates in [(\mu m)^2/ms]
        % mask is a 3D logical array - if empty the fit is done in the
        % whole volume
        % PRoptions is a structure containing hyperparameters for the polynomial regression
        % estimation. Default values are:
        % bounds = [ 2 2 ; 10 50 ; 0.001 0.05 ]; % D0 in um^2/ms, a in um, kappa in um/ms
        % PRoptions.Ntrain = 1e4;
        % PRoptions.diam_factor = 4;
        % PRoptions.bounds = bounds;
        % PRoptions.PRdegree = 4;
        %
        % (c) Santiago Coelho 2026 - santiago.Coelho@nyulangone.org
        [x, y, z, ntimes] = size(Dperp);
        if ~exist('mask','var') || isempty(mask)
            mask = true(x, y, z);
        end
        Dperp = RPBMestimation.vectorize(Dperp, mask);
        
        % Compute training data
        if ~isfield(PRoptions,'Ntrain') 
            Ntrain = 1e4;
        else
            Ntrain = PRoptions.Ntrain;
        end
        if ~isfield(PRoptions,'diam_factor') 
            diam_factor = 4; % Nat Phys/NMR Biomed
        else
            diam_factor = PRoptions.diam_factor;
        end
        if ~isfield(PRoptions,'bounds') 
            bounds = [ 2 2 ; 10 50 ; 0.001 0.05 ]; % D0 in um^2/ms, a in um, kappa in um/ms
        else
            bounds = PRoptions.bounds;
        end
        if ~isfield(PRoptions,'PRdegree') 
            PRdegree = 4;
        else
            PRdegree = PRoptions.PRdegree;
        end
        t = t(:)';
        bounds(1,:) = [ D0 D0 ];
        [~,a_tr,kappa_tr,tau_tr,zeta_tr] = RPBMestimation.sample_rpbm_param(bounds,Ntrain,diam_factor);
        % Dt_train = zeros(Ntrain,length(t));
        % parfor ii = 1:Ntrain
        %     Dt_train(ii,:) = D0*RPBMestimation.get_Dt_RPBM(t./tau_tr(ii),zeta_tr(ii));
        % end
        tnorm_all = t./tau_tr;
        zeta_all = repmat(zeta_tr,1,length(t));
        Dt_train = RPBMestimation.get_Dt_RPBM_vectorized(tnorm_all, zeta_all);

        % PR
        if isscalar(sigma)
            sigma_sample = sigma * ones(1,Ntrain);
        else
            sigma_sample = RPBMestimation.vectorize(sigma, mask);
            flag_outliers = detectVectorOutliers(sigma_sample,5);
            sigma_sample = randsample(sigma_sample(~flag_outliers),Ntrain,true);
            % figure, histogram(sigma_sample,linspace(0,0.4,100))
        end
        options.Degree = PRdegree;
        options.strategy = 'PR';
        % options.svdFLAG=1;
        % options.keep_sv=10;
        options.training_data = (Dt_train + randn(Ntrain,length(t)).*sigma_sample' )';
        options.training_objective = [zeta_tr(:) kappa_tr(:) a_tr(:)]';
        options.test_data = Dperp;
        out = RPBMestimation.DataDriven_regression(options);
        param.zeta  = RPBMestimation.vectorize(out.fits(:,1)', mask);
        param.kappa = RPBMestimation.vectorize(out.fits(:,2)', mask);
        param.a     = RPBMestimation.vectorize(out.fits(:,3)', mask);

        end
        % =================================================================
        function [s, mask] = vectorize(S, mask)
            % This functions unrolls 4D data into a matrix according to a
            % 3D mask
            if nargin == 1
               mask = ~isnan(S(:,:,:,1));
            end
            if ismatrix(S)
                n = size(S, 1);
                [x, y, z] = size(mask);
                s = NaN([x, y, z, n], 'like', S);
                for i = 1:n
                    tmp = NaN(x, y, z, 'like', S);
                    tmp(mask(:)) = S(i, :);
                    s(:,:,:,i) = tmp;
                end
            else
                for i = 1:size(S, 4)
                   Si = S(:,:,:,i);
                   s(i, :) = Si(mask(:));
                end
            end
        end
        % =================================================================
        function [ sigma_Dpar, sigma_Dperp ] = LLS_noiseprop_DTI(Dref,bvals,dirs,sigma_b0)
        % [ sigma_Dpar, sigma_Dperp ] = LLS_noiseprop_DTI(Dref,bvals,dirs,sigma_b0)
        
        syms Dxx Dxy Dxz Dyy Dyz Dzz
        D = [Dxx Dxy Dxz ; Dxy Dyy Dyz ; Dxz Dyz Dzz];
        lambda = eig(D);
        % JD = jacobian(lambda,[Dxx, Dxy, Dxz, Dyy, Dyz, Dzz]);
        lambda_par = lambda(1);
        lambda_perp = 1/2 * (lambda(2) + lambda(3));
        lambda_fiber = [lambda_par ; lambda_perp];
        JD = jacobian(lambda_fiber,[Dxx, Dxy, Dxz, Dyy, Dyz, Dzz]);

        dirs_rots = RPBMestimation.get60dirs();
        Nrots = size(dirs_rots,2);
        for n_id = 1:Nrots
            [~,R_from_z_to_v] = RPBMestimation.rotate_vector(dirs_rots(:,n_id));
            D_rot_current = R_from_z_to_v * Dref * R_from_z_to_v';
            D_rot(:,:,n_id) = D_rot_current;
            Dxx = D_rot_current(1,1);
            Dxy = D_rot_current(1,2);
            Dxz = D_rot_current(1,3);
            Dyy = D_rot_current(2,2);
            Dyz = D_rot_current(2,3);
            Dzz = D_rot_current(3,3);
            Jac_current = double(subs(JD));
            Jac(:,:,n_id) = Jac_current;
        end
        Jac_real = real(Jac); % norm(Jac_real(:)-Jac(:))
        
        ind = [ 1     1;  1     2;   1     3;   2     2;   2     3;  3     3];
        cnt = [ 1     2     2     1     2     1 ];
        
        if size(dirs,1) ~= 3
            dirs = dirs';
        end
        clearvars grad
        grad = [dirs ; bvals(:)']';
        sz = size(grad);
        bmatrix = [ones([sz(1) 1]) -(grad(:,ones(1,6)*4).*grad(:,ind(1:6,1)).*grad(:,ind(1:6,2)))*diag(cnt)];
        BtB = bmatrix' * bmatrix;
        invBtB = inv(BtB);
        COND = cond(bmatrix);
        % LLS
        cov_D_LLS_matrix = sigma_b0^2 * invBtB ;
        cov_D_LLS_current = sqrt(diag(cov_D_LLS_matrix(2:7,2:7))');
        cov_lambda_matrix_current = 0;
        for n_id = 1:Nrots
            cov_lambda_matrix_current = cov_lambda_matrix_current + ...
            1/Nrots * Jac_real(:,:,n_id) * cov_D_LLS_matrix(2:7,2:7) * Jac_real(:,:,n_id)';
        end
        sigma_Dpar = sqrt(diag(cov_lambda_matrix_current(1,1)))';
        sigma_Dperp = sqrt(diag(cov_lambda_matrix_current(2,2)))';

        end
        % =================================================================
        function [ sigma_Dpar, sigma_Dperp ] = LLS_noiseprop_DTI_numerical(Dref,bvals,dirs,sigma_b0_norm,mask)
        % [ sigma_Dpar, sigma_Dperp ] = LLS_noiseprop_DTI_numerical(Dref,bvals,dirs,sigma_b0_norm,mask)
        if size(dirs,1) ~= 3
            dirs = dirs';
        end
        [x,y,z] = size(sigma_b0_norm);
        if ~exist('mask','var') || isempty(mask)
            mask = true(x,y,z);
        end
        if max(bvals(:))>100
            bvals = bvals/1e3;
        end
        Nb = length(bvals);
        sigma_b0_norm = RPBMestimation.vectorize(sigma_b0_norm,mask);
        flag_outliers = RPBMestimation.detectVectorOutliers(sigma_b0_norm,4);
        sigma_b0_ref = mean(sigma_b0_norm(~flag_outliers));
        M = length(sigma_b0_norm);
        Dpar_gt = Dref(1);
        Dperp_gt = Dref(2); 
        D = ones(6,M).*[Dperp_gt 0 0 Dperp_gt 0 Dpar_gt]';
        grad = [dirs ; bvals(:)']';
        ind = [ 1     1;  1     2;   1     3;   2     2;   2     3;  3     3];
        cnt = [ 1     2     2     1     2     1 ];
        s = size(grad);
        % b = [ones([s(1) 1],'like', grad) -(grad(:,ones(1,6)*4).*grad(:,ind(1:6,1)).*grad(:,ind(1:6,2)))*diag(cnt)];
        b = [ -(grad(:,ones(1,6)*4).*grad(:,ind(1:6,1)).*grad(:,ind(1:6,2)))*diag(cnt)];      
        S_dti = exp(b*D);
        dwi = RPBMestimation.vectorize(S_dti + randn(Nb,M)*sigma_b0_ref,mask);
        [b0, dt] = RPBMestimation.dti_fit(dwi, grad, mask);
        [~, ~, rd, ad, ~ , ~, ~] = RPBMestimation.dti_parameters(dt, mask);
        ad = RPBMestimation.vectorize(ad,mask);
        rd = RPBMestimation.vectorize(rd,mask);
        sigma_Dpar = RPBMestimation.vectorize(std(ad-Dpar_gt)/sigma_b0_ref * sigma_b0_norm,mask);
        sigma_Dperp = RPBMestimation.vectorize(std(rd-Dperp_gt)/sigma_b0_ref * sigma_b0_norm,mask);

        % figure('Position',[2571 673 1800 500])
        % subplot(131), h=histogram(sigma_b0_norm(:),linspace(0,0.1,25));
        % subplot(132), hold on, h=histogram(ad(:),linspace(1.2,2,25)); plot([1 1]*Dpar_gt,[0 1.1]*max(h.Values),'r-.','LineWidth',2)
        % subplot(133), hold on, h=histogram(rd(:),linspace(0.5,1.2,25)); plot([1 1]*Dperp_gt,[0 1.1]*max(h.Values),'r-.','LineWidth',2)

        end
        % =================================================================
        function dirs = get60dirs()
        dirs =  [-0.130068563916747	0.0932848110280071	0.301874114147997	-0.160588543132534	-0.963953191993350	0.978837928613104	-0.549522980273078	0.592328122242095	-0.729164145170031	0.238633693018993	0.434651955799002	0.652635346397925	-0.734000052801543	-0.403599416152604	0.853698627179037	0.179812898840564	-0.511230028936703	-0.347227091988192	-0.151596735206056	0.719922327180808	-0.822392486806688	0.704707871446281	0.829703699851419	-0.862550086857427	-0.335321484736507	-0.314225258918558	0.194088504613566	0.327542356254105	0.505947381273477	-0.834296964584067	0.309476571411950	-0.117340713967168	-0.363227137830263	0.521389418163999	0.945240558338114	-0.524417749263304	-0.763218546997936	-0.145525976291921	-0.0122722916686933	0.685740007952837	-0.956922291701663	0.484894981408397	-0.223076366317911	-0.0912265656525496	0.939489850665938	-0.548677930291736	0.0899830129379917	-0.529095236280458	0.857923466019652	0.0313910318301205	-0.503828318019078	0.488811191935837	0.294732743377875	-0.956802130789944	0.300000411642112	0.00897408722442505	-0.709091302055868	0.650561378894920	-0.711236989891623	0.789508073435427 ; ...
        0.282325584363997	-0.332500776697028	0.898485458170594	-0.968673668717596	0.261165366605327	-0.200111124452651	0.264272708506881	-0.285537557395265	-0.355775017966924	0.854283276620088	-0.848130577636177	0.307985185967729	-0.537419604056184	0.910143183919012	0.515134135617358	-0.483539314532076	0.478397264082337	-0.766740371788774	0.531267969479908	-0.684110754074832	-0.564282076688453	0.157821133257434	-0.133437282194303	0.257480038001757	-0.543060152023286	0.874112698528370	0.230046094506407	0.475172375772998	-0.530090908461732	0.190037107797886	-0.938443721190627	0.870791399069994	-0.419677486106270	0.848661761051783	0.130543357216993	-0.753119539980642	0.630507380213516	-0.179058930954968	-0.966700744188533	0.632133714249398	-0.162507387356043	-0.696844732289984	0.0147999107837224	0.687075735781306	0.218072671767995	-0.00591585331993446	-0.729496857445582	0.657056695425554	-0.330180891156839	0.998046624894397	-0.857385809697511	-0.101631079322322	0.0246017271241047	-0.166433335585371	0.635666172830116	-0.816956746605835	0.630753698155324	0.578870483786831	-0.154577580948680	-0.533415663960950 ; ...
        0.950460116519448	-0.938478117766921	-0.318741118563859	0.189427145224130	0.0508615270291595	-0.0428000861997289	-0.792580866341241	0.753402746821313	0.584588561288273	0.461794373984095	-0.302906257113693	-0.692251565441700	-0.415239800191922	-0.0935248418712461	0.0763902891366770	0.856654570239605	0.713988736067953	-0.539946801878329	-0.833530308074974	0.117065447446015	0.0725005901498665	0.691721985930035	-0.542020536671319	-0.435558696036320	0.769834510241357	0.370390978466020	-0.953629092875851	0.816656119894876	-0.680456373435304	0.517527267441655	0.153452383184467	-0.477444757172385	-0.831827418398286	0.0890297138886364	0.299129936247056	0.397236683548934	-0.141272407124139	0.973015975958025	-0.255615066149348	-0.360787484264772	-0.240605645609469	0.528473533949100	-0.974688615626135	0.720836074998898	-0.264202820430936	0.836012877586388	-0.678039373753192	-0.536968090245386	0.393634228145500	0.0540142356934065	-0.105102803442199	-0.866449503637717	0.955262982117209	0.238389654393286	-0.711286348620636	0.576629118180567	0.315181372565219	0.491608436960495	-0.685746101467561	-0.303553506686020 ] ;
        end
        % =================================================================
        function [R_from_v_to_z,R_from_z_to_v] = rotate_vector(target)
        % [R_from_v_to_z,R_from_z_to_v] = rotate_vector(target)
        %
        % This function provides the rotation matrices that rotate the target
        % vector to the z axis (R_from_v_to_z) and viceversa (R_from_z_to_v)
        %
        % By: Santiago Coelho (06/10/2021)
        [phi,theta,~]=cart2sph(target(1),target(2),target(3)); theta=pi/2-theta;
        R_from_z_to_v = eul2rotm([phi', theta', 0*phi'],'ZYZ'); % this rotates the z axis to the desired direction
        R_from_v_to_z = eul2rotm([0*phi', -theta', -phi'],'ZYZ'); %  this rotates a given direction to the z axis
        end
        % =================================================================
        function [b0, dt] = dti_fit(dwi, grad, mask)
        %
        % Modification of dki_fit for DTI only [b0, dt] = dti_fit(dwi, grad, mask)
        % Santiago Coelho but 99.99% is dki_fit
        %
            % Diffusion Kurtosis Imaging tensor estimation using 
            % (constrained) weighted linear least squares estimation 
            % -----------------------------------------------------------------------------------
            % please cite:  Veraart, J.; Sijbers, J.; Sunaert, S.; Leemans, A. & Jeurissen, B.,
            %               Weighted linear least squares estimation of diffusion MRI parameters: 
            %               strengths, limitations, and pitfalls. NeuroImage, 2013, 81, 335-346 
            %------------------------------------------------------------------------------------
            % 
            % Usage:
            % ------
            % [b0, dt] = dti_fit(dwi, grad, mask)
            %
            % Required input: 
            % ---------------
            %     1. dwi: diffusion-weighted images.
            %           [x, y, z, ndwis]
            %       
            %       Important: We recommend that you apply denoising, gibbs correction, motion-
            %       and eddy current correction to the diffusion-weighted image
            %       prior to tensor fitting. Thes steps are not includede in this
            %       tools, but we are happy to assist (Jelle.Veraart@nyumc.org).
            %
            %     2. grad: diffusion encoding information (gradient direction 'g = [gx, gy, gz]' and b-values 'b')
            %           [ndwis, 4] 
            %           format: [gx, gy, gx, b]
            %
            % Optional input:
            % ---------------
            %    3. mask (boolean; [x, y, x]), providing a mask limits the
            %       calculation to a user-defined region-of-interest.
            %       default: mask = full FOV
            %
            % 
            % Copyright: Vision Lab, University of Antwerp (2016)
            % 
            % This Source Code Form is subject to the terms of the Mozilla Public
            % License, v. 2.0. If a copy of the MPL was not distributed with this file,
            % You can obtain one at http://mozilla.org/MPL/2.0/
            % 
            % This code is distributed  WITHOUT ANY WARRANTY; without even the 
            % implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
            % 
            % For more details, contact: Jelle.Veraart@nyumc.org 
        
            % parameter checks
            dwi = double(dwi);
            dwi(dwi<=0)=eps;
            [x, y, z, ndwis] = size(dwi);
            if ~exist('grad','var') || size(grad,1) ~= ndwis || size(grad,2) ~= 4
                error('');
            end
            grad = double(grad);
            grad(:,1:3) = bsxfun(@rdivide,grad(:,1:3),sqrt(sum(grad(:,1:3).^2,2))); grad(isnan(grad)) = 0;
        
            if ~exist('mask','var') || isempty(mask)
                mask = true(x, y, z);
            end
            dwi = RPBMestimation.vectorize(dwi, mask);
        
            % tensor fit
        
            ind = [ 1     1;  1     2;   1     3;   2     2;   2     3;  3     3];
            cnt = [ 1     2     2     1     2     1 ];
            s = size(grad);
            b = [ones([s(1) 1],'like', grad) -(grad(:,ones(1,6)*4).*grad(:,ind(1:6,1)).*grad(:,ind(1:6,2)))*diag(cnt)];
        
            % unconstrained LLS fit
            dt = b\log(dwi);
            w = exp(b*dt);
        
            nvoxels = size(dwi,2);
        
            % WLLS fit initialized with LLS
            parfor i = 1:nvoxels
                wi = diag(w(:,i));
                logdwii = log(dwi(:,i));
                dt(:,i) = (wi*b)\(wi*logdwii);
            end
        
            b0 = exp(dt(1,:));
            dt = dt(2:end, :);
            b0 = RPBMestimation.vectorize(b0, mask);
            dt = RPBMestimation.vectorize(dt, mask);
        end
        % =================================================================
        function [fa, md, rd, ad, fe , lambdas, Vall] = dti_parameters(dt, mask)
        %
        % [fa, md, rd, ad, fe , lambdas, Vall] = dti_parameters(dt, mask) 
        % Modification of dki_parameters
        % Santiago Coelho but 99.99% is dki_fit
        %
            % diffusion tensor parameter calculation
            %
            % -----------------------------------------------------------------------------------
            % please cite: Veraart et al.  
            %              More Accurate Estimation of Diffusion Tensor Parameters Using Diffusion Kurtosis Imaging,
            %              MRM 65 (2011): 138-145.
            %------------------------------------------------------------------------------------
            % 
            % Usage:
            % ------
            % [fa, md, ad, rd, fe, mk, ak, rk] = dki_parameters(dt [, mask [, branch]])
            %
            % Required input: 
            % ---------------
            %     1. dt: diffusion kurtosis tensor (cf. order of tensor elements cf. dki_fit.m)
            %           [x, y, z, 21]
            %
            % Optional input:
            % ---------------
            %    2. mask (boolean; [x, y, x]), providing a mask limits the
            %       calculation to a user-defined region-of-interest.
            %       default: mask = full FOV
            %
            % 
            % output:
            % -------
            %  1. fa:                fractional anisitropy
            %  2. md:                mean diffusivity
            %  3. rd:                radial diffusivity
            %  4. ad:                axial diffusivity
            %  5. fe:                principal direction of diffusivity
            %
            % Important: The presence of outliers "black voxels" in the kurtosis maps 
            %            are we well-known, but inherent problem to DKI. Smoothing the 
            %            data in addition to the typical data preprocessing steps
            %            might minimize the impact of those voxels on the visual
            %            and statistical interpretation. However, smoothing comes
            %            with the cost of partial voluming. 
            %       
            % Copyright (c) 2016 New York University and University of Antwerp
            % 
            % This Source Code Form is subject to the terms of the Mozilla Public
            % License, v. 2.0. If a copy of the MPL was not distributed with this file,
            % You can obtain one at http://mozilla.org/MPL/2.0/
            % 
            % This code is distributed  WITHOUT ANY WARRANTY; without even the 
            % implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
            % 
            % For more details, contact: Jelle.Veraart@nyumc.org 
            
        [nx,ny,nz,n] = size(dt);
        if ndims(dt)~=4
            error('size of dt needs to be [x, y, z, 6]')
        end
        if n~=6
            error('dt needs to contain 6')
        end
        if ~exist('mask','var') || isempty(mask)     
                mask = ~isnan(dt(:,:,:,1));
        end
            
        dt = RPBMestimation.vectorize(dt, mask);
        nvoxels = size(dt, 2);
        
        
        % DTI parameters
        for i = 1:nvoxels
            DT = dt([1:3 2 4 5 3 5 6], i);
            DT = reshape(DT, [3 3]);
            [eigvec, eigval] = eigs(DT);
            eigval = diag(eigval);
            
            l1(i) = eigval(1,:);
            l2(i) = eigval(2,:);
            l3(i) = eigval(3,:);
            
            e1(:, i) = eigvec(:, 1); 
            Vall(:,i) = eigvec(:);
        end
        md = (l1+l2+l3)/3;
        rd = (l2+l3)/2;
        ad = l1;
        fa = sqrt(1/2).*sqrt((l1-l2).^2+(l2-l3).^2+(l3-l1).^2)./sqrt(l1.^2+l2.^2+l3.^2);
        
        % return maps
        fa = RPBMestimation.vectorize(fa, mask);
        md = RPBMestimation.vectorize(md, mask);
        ad = RPBMestimation.vectorize(ad, mask);
        rd = RPBMestimation.vectorize(rd, mask);
        fe = RPBMestimation.vectorize(e1, mask);
        
        nOutputs = nargout;
        if nOutputs >= 6
            lambdas = RPBMestimation.vectorize([l1;l2;l3], mask);
        end
        if nOutputs >= 7
            Vall = RPBMestimation.vectorize(Vall, mask);
            Vall = reshape(Vall,[nx,ny,nz 3 3]);
        end

        end
        % =================================================================
    end
end

