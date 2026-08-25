classdef DTI < handle
    %
    % Diffusion Tensor Imaging (DTI) class
    %
    % Usage:
    % dti = DTI(grad [,estimator='wlls'] [,options]);
    %
    % Required input: 
    % ---------------
    %
    %     grad: diffusion encoding information; [ndwis x 4] or [ndwis x 6]
    %
    %     Acceptable input formats: 
    %     
    %      1. [ndwis, 4], gradient direction and b-values, 
    %           format: [gx, gy, gx, bval]
    %
    %      2. [ndwis, 6], elements of the effective b-matrix (cf. Matiello et al., MRM 37:292-300, 1997), 
    %	        format: [bxx, bxy, bxz, byy, byz, bzz]
    %
    %      if the effective b-matrix is provided by the manufacturer, (2.)
    %      might be the preferred input format.
    %
    % Optional input:
    % ---------------
    %
    %    Estimator: selected parameter estimator
    %
    %    Available estimators:
    %	
    %       1. 'LLS': linear least squares estimator
    %       2. 'WLLS': Weighted linear least squares estimator        (cf. Veraart et al., NeuroImage 81:335-346, 2013)
    %       3. 'NLS': nonlinear least squares estimator         	  (cf. Koay et al., JMRI 182:115-125, 2006)
    %       4. 'CNLS': Constrained nonlinear least squares estimator  (cf. Koay et al., JMRI 182:115-125, 2006)
    %       5. 'CLS': Conditional least squares estimator             (cf. Veraart et al., MRM 70(4):972-984, 2013)
    %       6. 'CCLS': constrained conditinional least squares estimator  	  
    %    
    %    The default option is 'WLLS'.       
    %
    %    remark: CLS and CCLS require a noise level (or map) as an input and can only be used if your data is Rician distributed.        
    %
    % Available functions within the DTI class: 
    % -----------------------------------------
    %  1.  eddy_correct:      Correction of motion and eddy current distortion correction 
    %  2.  bmatrixrotation:   B-matrix rotation, only to be used in combination with (1)
    %  3.  skullstrip:        Skullstripping of diffusion-weighted data
    %  4.  smoothing:         Isotropic smoothing of diffusion-weighted data
    %  5.  dt:                Tensor estimation
    %  6.  outliers:          Automatic detection of outliers
    %  7.  getnoisemap        Estimation of 3D noisemap for Rician distributed data  
    %  8.  dwi:               Calculation of diffusion-weighted data based on tensors
    %  9.  eig:               Eigenvalue decomposition of the diffusion tensor
    %  10. adc:               Calculation of apparent diffusion coefficient along a particular direction
    %  11. fa:                Calculation of fractional anisitropy
    %  12. md:                Calculation of mean diffusivity
    %  13. rd:                Calculation of radial diffusivity
    %  14. ad:                Calculation of axial diffusivity
    %  15. fefa:              Calculation of color encoded fractional anisotropy map
    %  16. cl:                Calculation of linear diffusion
    %  17. cp:                Calculation of planar diffusion  
    %  18. cs:                Calculation of spherical diffusion  
    %  20. savenii:           Save data (diffusion-weigthed, tensor, or parameter) as NiFTi given a reference template
    %
    %  Remark 1: See functions for help on their use.
    %  Remark 2: (19) and (20) are currently specifically added for CBI
    %
    % Example:
    % --------
    %   dti = DTI(grad);
    %   dwi = dti.eddycorrect(dwi);
    %   [dwi, mask] = dti.skullstrip(dwi [, mask]);
    %   [dt, b0] = dti.dt(dwi, outliers);
    %   [eigval, eigvec] = DTI.eig(dt);
    %   md = DTI.md(eigval);
    %   fa = DTI.fa(eigval);
    %
    % Copyright 
    % Ben Jeurissen (ben.jeurissen@ua.ac.be)
    % Jelle Veraart (jelle.veraart@nyumc.org) 
    %
    properties (Constant = true, Access = protected)
        ind = [ 1     1
            1     2
            1     3
            2     2
            2     3
            3     3];
        cnt = [ 1     2     2     1     2     1 ];
    end
    
    properties (GetAccess = public, SetAccess = protected)
        b;
        estimator;
    end
    
    properties (GetAccess = private, SetAccess = protected)
        fun;
    end
    
    methods (Access = public)
        function this = DTI(grad, estimator, varargin)
            if nargin<1
               return 
            end
            if ~exist('grad','var') || isempty(grad)
                error('supply gradient table');
            end
            this.setGrad(grad);

            if ~exist('estimator','var') || isempty(estimator)
                estimator = 'wlls';
            end
            this.setEstimator(estimator,varargin{:});
        end
        
        function this = setGrad(this, grad)
            s = size(grad);
            if s(2) == 4 % gradient table
                
                
                normgrad = sqrt(sum(grad(:, 1:3).^2, 2)); normgrad(normgrad == 0) = 1;
                grad(:, 1:3) = grad(:, 1:3)./repmat(normgrad, [1 3]);
                this.b = DTI.grad2b(grad);
            elseif s(2) == 6 % b-matrix
                this.b = [ ones([s(1) 1],class(grad)) -grad ];
                this.b(:, [3 4 6]) = 2*this.b(:, [3 4 6]);
                disp(this.b)
            else
                error('gradient scheme not supported');
            end
        end
        
        function this = setEstimator(this, estimator, varargin)
            switch estimator
                case 'lls'
                    this.fun = @this.dt_lls;
                case 'wlls'
                    this.fun = @this.dt_wlls;
                    this.estimator = checkparams(varargin,...
                    {'wls_iter','integer',[1 10],1});
                case 'nls'
                    this.fun = @this.dt_nls;
                    this.estimator = checkparams(varargin,...
                    {'lambda','float',[],0.01,...
                    'tolX','float',[],1e-12,...
                    'tolFun','float',[],1e-12,...
                    'maxIter','integer',[],400,...
                    'wls_iter','integer',[1 10],1});
                case 'cls'
                    this.fun = @this.dt_cls;
                    this.estimator = checkparams(varargin,...
                    {'lambda','float',[],0.01,...
                    'tolX','float',[],1e-12,...
                    'tolFun','float',[],1e-12,...
                    'maxIter','integer',[],1000,...
                    'wls_iter','integer',[1 10],1});
                case 'ccls'
                    this.fun = @this.dt_ccls;
                    this.estimator = checkparams(varargin,...
                    {'lambda','float',[],0.01,...
                    'tolX','float',[],1e-12,...
                    'tolFun','float',[],1e-12,...
                    'maxIter','integer',[],1000,...
                    'wls_iter','integer',[1 10],1});
                case 'cnls'
                    this.fun = @this.dt_cnls;
                    this.estimator = checkparams(varargin,...
                    {'lambda','float',[],0.01,...
                    'tolX','float',[],1e-12,...
                    'tolFun','float',[],1e-12,...
                    'maxIter','integer',[],400,...
                    'wls_iter','integer',[1 10],1});
                otherwise
                    error('estimator not supported');
            end
            this.estimator.name = estimator;
        end
        
        function varargout = dt(this, varargin)           
            % [dt, b0] = *.dt(dwi[,outliers [, sigma]])
            %
            % Tensor estimation based on the predefined estimator. 
            % If data is corrupted with imaging artifacts, it might be
            % interesting to provide an outlier map (same size as the
            % diffusion weighted data). The 'outliers' will be excluded
            % from the fit. Check this class' function "outliers" for 
            % automatic outlier detection function. 
            % 
            % If the conditional least squares estimator is
            % chosen, the noise level or map must be provided.
            % 
            % The non-optional input variable dwi might be a vector [ndwis,
            % nvxls] or a 4D data stack [x, y, z, ndwis]. 
            % The estimated tensor dt is stored according to the
            % 'Antwerp'-convention (see changeTensorOrder). The dimensions
            % of dt is [nparams, nvxls] or [x, y, z, nparams] depending on
            % the size of dwi.
            
            if ndims(varargin{1}) ~= 2;  %#ok<*ISMAT>
                [varargin{1}, mask] = vec(varargin{1}); 
                if numel(varargin) > 1 && ~isempty(varargin{2})
                        varargin{2} = vec(varargin{2}, mask);
                end
                if numel(varargin) > 2 && ~isempty(varargin{3}) && ~isscalar(varargin{3})
                        varargin{3} = vec(varargin{3}, mask);
                end
            end
            
            varargin{1} = feval(class(this.b),varargin{1});
            fact = 1000/max(varargin{1}(:));
            varargin{1} = varargin{1}*fact;
            if numel(varargin) == 3
                varargin{3} = varargin{3}*fact;  
            end
%             varargin{1}(varargin{1}<=1)=1;
            varargin{1}(varargin{1}<=0)=eps;
            varargout{1} = this.fun(varargin{:});
            varargout{2} = exp(min(15,varargout{1}(1,:)))./fact;
            varargout{1} = varargout{1}(2:end,:);
            if exist('mask','var'); varargout{1} = unvec(varargout{1},mask); varargout{2} = unvec(varargout{2},mask); end;
        end
        
        function dwi_hat = fit(this, dwi)
            if ndims(dwi) ~= 2; [dwi, mask] = vec(dwi); end;
            [dt, b0] = this.dt(dwi);
            dwi_hat = this.dwi(dt, b0);
            if exist('mask','var'); dwi_hat = unvec(dwi_hat,mask); end;
        end
        
        function dwi = dwi(this, dt, b0)
            % dwi = *.dwi(dt, b0)
            %
            % Calculation of diffusion-weighted signals "dwi" based on diffusion
            % tensor(s) "dt" and nondiffusion-weighted signal "b0". dt can
            % be of size [nparams, nvxls] or [x, y, z, nparams]. b0 must be
            % size accoringly. Output is also sized accordingly.
            
            if ndims(dt) ~= 2; [dt, mask] = vec(dt); b0 = vec(b0,mask); end;
            dwi = exp(this.b*[log(b0); dt]);
            if exist('mask','var'); dwi = unvec(dwi,mask); end;
        end
        
        function dwi = sim_dwi(this, fa, diffusivity, el, az, b0)
            dwi = this.dwi(DTI.sim_dt(fa, diffusivity, el, az), repmat(b0, [1 size(el,1)]));
        end
        
        function [dwi, mask] = skullstrip(this, dwi, mask)
            % [dwi, mask] = *.skullstrip(dwi[, mask])
            % 
            % All diffusion-weighted data samples outside an anatomical
            % mask are set to NaN. If no mask is provided, a rudimentary 
            % brain mask will be calculated from the b0 image(s).
            dwi(dwi<=0)=eps;            
            dwi = feval('double', dwi);
            if ~exist('mask', 'var') || isempty(mask)
                grad = DTI.b2grad(this.b); bval = grad(:, 4); clear grad;
                mask = BrainMask(mean(dwi(:,:,:,bval==min(bval)), 4));
            end
            mask = repmat(mask, [1 1 1 size(dwi, 4)]);
            dwi(~mask) = NaN;
            mask = mask(:,:,:,1);
        end
        
        
        function [outliers] = outliers(this, dwi, sigma, method, t)
         %  outliers = *.outliers(dwi [, sigma])
         %  
         %  Automatic detection of outliers in diffusion weighted data
         %  "dwi". The algorithm is based on the RESTORE algorihtm, but 
         %  based on the  instead of the NLS to speed up to process  
         %  The non-optional input variable "dwi" might be a vector [ndwis,
         %  nvxls] or a 4D data stack [x, y, z, ndwis]. In the latter, the 
         %  number of false is reduced by removing isolated outliers.
         %  
         % The algorithm requires and estimate of the noise level (sigma),
         % if not provided, an estimate based on residuals is generated.
         
       
         if ndims(dwi) ~= 2; [dwi, mask] = vec(dwi); end;
         dwi(dwi <= 0) = eps;
         [ndwis, nvxls] = size(dwi);
         if ~exist('sigma', 'var') || isempty(sigma)   
            dt = pinv(this.b)*log(dwi);
            w = exp(this.b*dt);

            X = this.b; e = log(dwi) - X*dt; [n, N] = size(dwi); m = size(dt, 1); 
            sigma = zeros(1, N);
            whos e w
            class(e(:,1).*w(:,1))

            % tmp = e(:,i) .* w(:,i);
            % 
            % if ~isnumeric(tmp) || ~isreal(tmp)
            %     fprintf('Bad input at i=%d: class=%s, size=%s\n', ...
            %         i, class(tmp), mat2str(size(tmp)));
            %     keyboard % drop into debugger
            % end
            parfor i=1:size(dwi, 2)
                sigma(i) = sqrt(n/(n-m)) * 1.4826 * mad(double(e(:,i).*w(:,i)), 1);
            end
            sigma = nanmedian(sigma);
         end
         
         if ~exist('method', 'var') || isempty(method)         
            method = 'irwlls';
         end
         
         if ~exist('t', 'var') || isempty(t)         
           t = 0;
         end
         
         if ndims(sigma) == 3; sigma = vec(sigma, mask); end;
         if isscalar(sigma); sigma = sigma*ones([1, size(dwi, 2)]); end
         
         
         bval= -sum(this.b(:,[2 5 7]),2);
         b0_pos = bval==min(bval); 
         b0s = max(dwi(b0_pos, :), [], 1);
         outliers = dwi > repmat(b0s + 3*sigma, [ndwis, 1]);
         
         bmat = this.b;
         
         switch method
             case 'irwlls'
         
                 iters = ones([1, nvxls]);
    
                 warning off
                 parfor i = 1:nvxls
                     try
                         dwi_i = dwi(:,i);
                         out_i = outliers(:, i);


                         dwi_i = dwi_i(~out_i);
                         bmat_i = bmat(~out_i, :);
                         if size(bmat_i, 1) > size(bmat_i, 2)
                             sigma_i = sigma(i);
                             n_i = size(dwi_i, 1);

                             % goodness-of-fit
                             leverage = diag(bmat_i*inv(bmat_i'*bmat_i)*bmat_i');
                             dt = pinv(bmat_i)*log(dwi_i);
                             dwi_hat_i = exp(bmat_i*dt);                

                             e_l  = log(dwi_i) - log(dwi_hat_i);
                             e_nl = dwi_i - dwi_hat_i;

                             chi2res = sum(e_l.*e_l./(sigma_i./dwi_hat_i).^2)/(n_i-7) - 1; 
                             %gof = abs(chi2res) < 3*sqrt(2/(n_i-7)); 
                             %'b'
                             gof = abs(chi2res) < 25*sqrt(2/(n_i-7));


                             % iterative GMM
                             iter = 1;
                             while ~gof && iter < 25

                                 iter = iter+1;
                                 C = 1.4826*mad(e_nl(:),1)./dwi_hat_i;     
                                 GMM = C.^2./(e_l.^2 + C.^2).^2;  

                                 w = sqrt(GMM).*dwi_hat_i; 
                                 leverage = diag(bmat_i*inv(bmat_i'*diag(w.^2)*bmat_i)*bmat_i'*diag(w.^2));

                                 dt = mex_wls(double(bmat_i),double(log(dwi_i)), double(w));
                                 dwi_hat_i = exp(bmat_i*dt); dwi_hat_i(dwi_hat_i<1)=1;

                                 e_l = log(dwi_i) - log(dwi_hat_i);
                                 chi2res(iter) = sum(e_l.*e_l./(sigma_i./dwi_hat_i).^2)/(n_i-7) - 1;    
                                 gof = abs(chi2res(iter) - chi2res(iter-1)) < 1e-3;
                             end

                             % outlier selection

                             e_l = log(dwi_i) - log(dwi_hat_i);
                             lowerbound_linear = -3*sigma_i*sqrt(1-leverage)./dwi_hat_i;
                             e_nl = dwi_i - dwi_hat_i;
                             upperbound_nonlinear = 3*sigma_i*sqrt(1-leverage);

                             out_i(e_l < lowerbound_linear) = true;
                             out_i(e_nl > upperbound_nonlinear)  = true;

                             outliers(:, i) = out_i;
                         else
                             outliers(:, i) = 1;
                         end
                     catch
                         outliers(:, i) = 2;
                     end
                 end
                 
             case 'restore'
                 
                parfor i = 1:nvxls
                    i
                    
                    try
                    % initial estimate, goodness-of-fit criteria met?
                        dwi_i = dwi(:,i);
                        out_i = outliers(:, i);

                        bmat_i = bmat(~out_i, :);

                        if size(bmat_i, 1) > size(bmat_i, 2)
                            dwi_i = dwi_i(~out_i);

                            sigma_i = sigma(i);
                            n_i = size(dwi_i, 1);


                            dti = DTI(-bmat_i(:, 2:7));
                            dti.setEstimator('nls');

                            [dt_i, b0_i] = dti.dt(dwi_i);
                            dwi_hat = dti.dwi(dt_i, b0_i);
                            residu = dwi_i - dwi_hat;
                            chi2res = residu'*residu/sigma_i.^2/(n_i-7) - 1;
                            gof = abs(chi2res) < 3*sqrt(2/(n_i-7));

                            % if needed, iterative reweightning
                            iter = 1;
                            while ~gof && iter < 25
                                 iter = iter+1;

                                 C = 1.4826 * mad(residu, 1);
                                 GMM = 1./(residu.^2 + C^2).^2;  GMM = GMM ./ mean(GMM);    % aangepast met kwadraat!!!


                                 [dt_i, b0_i] = dti.dt(dwi_i, [], sqrt(GMM));
                                 dwi_hat = dti.dwi(dt_i, b0_i);
                                 residu = dwi_i - dwi_hat;
                                 chi2res = residu'*residu/sigma_i.^2/(n_i-7) - 1;
                                 gof = abs(chi2res) < 3*sqrt(2/(n_i-7));
                            end

                            % detect outliers + final fit
                            tmp = zeros(size(residu));
                            tmp(residu < -3*sigma_i) = 1;
                            tmp(residu > 3*sigma_i)  = 1;
                            tmp(b0_pos) = 0;

                            out_i(~out_i) = tmp;

                            if cond(bmat(~out_i,:)) < 1e15  % simple singularity criteria
                               outliers(:, i) = out_i;
                            end

                        else
                            outliers(:, i) = 1;
                        end
                    catch
                        outliers(:, i) = 2;
                    end
                end
                
                
             case 'irestore'
                 parfor i = 1:nvxls
                    try 
                        list = 1:ndwis;

                        dti = DTI(-bmat(:, 2:7));
                        dti.setEstimator('nls');

                        dwi_i = dwi(:,i);
                        [dt_i, b0_i] = dti.dt(dwi_i);
                        dwi_hat = dti.dwi(dt_i, b0_i);
                        residu = dwi_i - dwi_hat;

                        chi2res = (residu'*residu)/sigma(i).^2/(ndwis-7) - 1;
                        gof = abs(chi2res) < 2*sqrt(2/(ndwis-7)); %Change Criteria to 2 sigma

                        % if needed, iterative reweightning
                        iter = 1;
                        while ~gof && iter < ndwis*0.5

                            residu(b0_pos) = Inf;
                            [~, pos] = min(residu);
                            list(pos) = NaN;

                            bmati = bmat(~isnan(list), :);
                            if cond(bmati) < 1e15  % simple singularity criteria
                                dti.b = bmati;
                                [dt_i, b0_i] = dti.dt(dwi_i(~isnan(list)));
                                dwi_hat = dti.dwi(dt_i, b0_i);
                                residu = dwi_i(~isnan(list)) - dwi_hat;
                                chi2res = residu'*residu/sigma(i).^2/(ndwis-iter-7) - 1; 
                                gof = abs(chi2res) < 3*sqrt(2/((ndwis-iter)-7));
                                iter = iter+1;

                                tmp = Inf([1, ndwis]);
                                tmp(~isnan(list)) = residu;
                                residu = tmp;
                            else
                                list(isnan(list)) = 1;
                                gof = true;
                            end
                        end

                        tmp = zeros([1,ndwis]);
                        tmp(isnan(list)) = 1;
                        outliers(:, i) = tmp;
                    catch
                        outliers(:, i) = 2;
                    end
                end
                 
                 
             case 'acdc'
                grad = DTI.b2grad(this.b);
                bval = grad(:  , 4)./1000;
        
                bval_b0 = min(bval);
                list.b0 = find(bval < (bval_b0*1.01 + eps));
                b0 = mean(dwi(list.b0, :), 1);
       
        
                dADC = log(repmat(b0, [size(dwi, 1), 1])./dwi)./(bval(:, ones(1, size(dwi, 2)))+eps);
                dADC(list.b0, :) = 0; 
                outliers=dADC>t;
         end

         if exist('mask','var'); 
             disp(size(single(outliers)))
             disp(size(mask))
             % outliers = unvec(single(outliers),mask); 
             
             % outliers = unvec(single(outliers), mask, size(dwi));
             outliers(isnan(outliers)) = 0;
             %for j=1:size(outliers, 4)
             %   outliers(:,:,:,j) = imerode(outliers(:,:,:,j), ones([3 3 1]));
            %    outliers(:,:,:,j) = imdilate(outliers(:,:,:,j), ones([3 3 1]));
            %% end
         end;
        end
        
        
        function [dwi, params] = eddycorrect(this, dwi, mask, lowbval)
             
            % Eddy current distortion and motion correction algorithm. Unlike FSL's
            % eddy, there is no constraint on the gradient directions. The low b-valued 
            % DW images are used to predict the high b-valued DW images using the DTI 
            % model in order to create a stack of aligned images with similar constrast 
            % as the acquired images. Registration is expected to be much more robust 
            % in this case. Signal modulation is applied, transformation parameters are
            % provided for b-matrix rotation. 
            %
            % remark: if the gradient directions are distributed over the
            % entire unit sphere, you might consider FSL's eddy tool
            % instead (http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/EDDY)
            %
            % The algorithm makes use of elastix (http://elastix.isi.uu.nl/index.php) for the registration. 
            % Please download and install before using this function.
            % 
            % usage [dwi, params] = *.eddy_correct(dwi, [[mask], lowbval])
            % 
            % Input: 
            %  dwi: 4D diffusion-weighted data 
            %  mask (optional): 3D mask. If empty or not provided, mask is created temporatly using
            %                            BrainMask.m to define exclude the image background 
            %                            during registration
            %  lowbval (optional): The maximal bvalue used for DTI based signal prediction. Default 1000 s/mm2 
            %
            % Output:
            %   dwi: Corrected 4D diffusion-weighted data
            %   params: affine correction parameters (rotation, shearing, scaling, translation)
            
            [s, ~] = system('elastix --help');
            if s
                error('elastix must be properly installed and on path (see installation documentation)')
            end
            root = cd;
            mkdir(fullfile(root, 'input'));
            mkdir(fullfile(root, 'output'));
            DTI.createparameterfile(fullfile(root, 'input', 'parameters.txt'));
            
            grad = DTI.b2grad(this.b); bval = grad(:, 4); clear grad;
            if ~exist('lowbval', 'var') || isempty(lowbval)
                lowbval = 1000;
            end
            list = find(floor(bval)<=lowbval);
            b0_pos = bval==min(bval);
            
            
            if ~exist('mask', 'file') || isempty(mask)
                mask = BrainMask(dwi(:,:,:,b0_pos));
            end
            
            % step 1: registration of low b-valued diffusion-weighted data
            disp('step 1: registration of low b-valued diffusion-weighted data') 
    
            dwi_ = dwi;
            params = zeros(12, size(dwi, 4));     
            for i=1:numel(list)
                flt = dwi(:,:,:,list(i));
                ref = dwi(:,:,:,b0_pos(1));

                nii = make_nii(flt, [2.5 2.5 2.5]); save_nii(nii, fullfile(root, 'input', 'flt.nii'))
                nii = make_nii(ref, [2.5 2.5 2.5]); save_nii(nii, fullfile(root, 'input', 'ref.nii'))
                nii = make_nii(single(mask), [2.5 2.5 2.5]); save_nii(nii, fullfile(root, 'input', 'mask.nii'))

                f = sprintf('elastix -f %s -fMask %s -m %s -out %s -p %s', ...
                            fullfile(root, 'input', 'ref.nii'), ... 
                            fullfile(root, 'input', 'mask.nii'), ...
                            fullfile(root, 'input', 'flt.nii'), ...
                            fullfile(root, 'output'), ...
                            fullfile(root, 'input', 'parameters.txt'));

                [~,~] = system(f);
                nii = load_untouch_nii(fullfile(root, 'output', 'result.0.nii'));  dwi_(:,:,:,list(i)) = nii.img;
                
               f = sprintf('%s',  fullfile(root, 'output', 'TransformParameters.0.txt'));                                    
               fid = fopen(f);
               C = textscan(fid, '%s');
               fclose(fid);

               params(1:11, list(i)) = str2double(C{1}(6:16));
               tmp = cell2mat(C{1}(17));
               params(12, list(i)) = str2double(tmp(1:end-1));
            end

            % step 2a: If high b-values are included in the DW stack
            if nnz(bval > lowbval) > 0
                  % step 2a: generate DW contrast based on DTI fit
                  disp('step 2: generate DW contrast based on DTI fit') 
                  dwi_ = vec(dwi_, mask); dwi_(dwi_<1)=1;
                  
                  b_ = this.b(list,:);
                  dt = b_(:, 1:7)\log(dwi_(list, :));
                  dwi_ = exp(this.b(:, 1:7)*dt); 
                  dwi_(~isfinite(dwi_))=0;

%                   
                  
                   b0s = max(dwi_(b0_pos, :), [], 1);
                   remove = dwi_ > repmat(b0s, [size(dwi_, 1), 1]);
                   dwi_(remove) = 0;
                   dwi_ = unvec(dwi_, mask); dwi_(isnan(dwi_)) = 0;
                   
                  % step 2b:
                  disp('step 3: registration of all DW images to the predicted images');
                  params = zeros(12, size(dwi, 4));
                  for i=1:size(dwi, 4)
                    flt = dwi(:,:,:,i);
                    ref = dwi_(:,:,:,i);

                    nii = make_nii(flt, [2.5 2.5 2.5]); save_nii(nii, fullfile(root, 'input', 'flt.nii'))
                    nii = make_nii(ref, [2.5 2.5 2.5]); save_nii(nii, fullfile(root, 'input', 'ref.nii'))
                    nii = make_nii(single(mask), [2.5 2.5 2.5]); save_nii(nii, fullfile(root, 'input', 'mask.nii'))

                    f = sprintf('elastix -f %s -fMask %s -m %s -out %s -p %s', ...
                                fullfile(root, 'input', 'ref.nii'), ... 
                                fullfile(root, 'input', 'mask.nii'), ...
                                fullfile(root, 'input', 'flt.nii'), ...
                                fullfile(root, 'output'), ...
                                fullfile(root, 'input', 'parameters.txt'));

                    [~,~] = system(f);
                    nii = load_untouch_nii(fullfile(root, 'output', 'result.0.nii'));  dwi(:,:,:,i) = nii.img;

                    f = sprintf('%s',  fullfile(root, 'output', 'TransformParameters.0.txt'));                                    
                    fid = fopen(f);
                    C = textscan(fid, '%s');
                    fclose(fid);

                    params(1:11, i) = str2double(C{1}(6:16));
                    tmp = cell2mat(C{1}(17));
                    params(12, i) = str2double(tmp(1:end-1));     
                  end
            else
                dwi = dwi_;
            end
            
            % signal modulation
            for i=1:size(dwi, 4)
                dwi(:,:,:,i) = dwi(:,:,:,i)*prod(params([7 8 9], i));  
            end
            
           
            % clean up tmp files
            files = dir(fullfile(root, 'input'));
            for i = 1:numel(files)
                if ~isdir(fullfile(root, 'input', files(i).name))
                    delete(fullfile(root, 'input', files(i).name))
                end
            end
            rmdir(fullfile(root, 'input'))

            files = dir(fullfile(root, 'output'));
            for i = 1:numel(files)
                if ~isdir(fullfile(root, 'output', files(i).name))
                    delete(fullfile(root, 'output', files(i).name))
                end
            end
            rmdir(fullfile(root, 'output')); 
         end
         function dwi_ = gibbsremoval(this, dwi, nc, maxbval, method)
             
             
            if ~exist('method', 'var') || isempty(method)
                method = 'TGV';
            end
            
            grad = DTI.b2grad(this.b); bval = grad(:, 4); clear grad;
            if ~exist('maxbval', 'var') || isempty(maxbval)
                maxbval = max(bval); % usually there is no need to change this
            end
            
            [nx, ny, nz, ndwis] = size(dwi);
   
            if isscalar(nc)
                nc = [nc nc];
            end
            
            kmask =  zpad(ones(nx,ny), nc(1), nc(2)); 
            FT = p2DFT(kmask, nc, 1, 2);
            
       
            switch method
                case 'TGV'
                    
                    imgSens = ones(nc(1), nc(2));
                    w = ones(nc(1), nc(2)); 
                    reduction = 2^(-8);     % usually there is no need to change this
                    maxit = 1000;           % use 500 Iterations for optimal image quality
                    alpha = 1e-4; % 2.5e-3; %1e-4; 
                    
                case 'TV'
                    param = init;
                    param.XFM = 1;
                    param.xfmWeight = 0;
                    param.Itnlim = 1000;
                    param.FT = FT;
                    param.TV = TVOP;
                    param.TVWeight = 1e-4;
            end
            
            
            dwi_ = zeros(nc(1),nc(2), nz, ndwis);

            list.do = find(bval<=maxbval);
            list.dont = find(bval>maxbval);
            
            for i=1:numel(list.do)
                disp(['processing ', num2str(i), ' out of', num2str(numel(list.do))]);
                for j  = 1:nz
                    k = zpad(fft2c(dwi(:,:,j,list.do(i))), nc(1), nc(2));
                    maxk = max(k(:));                  
                    switch method
                        case 'TGV'              
                            dwi_(:,:,j, list.do(i)) = tgv2_l2_2D_pd(imgSens, k/maxk, FT, w, 2*alpha, alpha, maxit, reduction)*maxk;
                        case 'TV'
                            param.data =  k / maxk;
                            dwi_(:,:,j, list.do(i)) = FT'*param.data;
                            for n=1:1
                              dwi_(:,:,j, list.do(i)) = fnlCg(dwi_(:,:,j, list.do(i)),param);
                            end
                            dwi_(:,:,j, list.do(i)) = dwi_(:,:,j, list.do(i))*maxk;
                    end
                end
            end
            for i=1:numel(list.dont)
               for j  = 1:nz
                    k = zpad(fft2c(dwi(:,:,j,list.dont(i))), nc(1), nc(2));
                    dwi_(:,:,j, list.dont(i)) = abs(ifft2c(k));  
                   
                end
            end
            
         end
         function dwi = smoothing(this, dwi, maxbval, kernelsize, width, csfmask, dim)
            % dwi = *.smoothing(dwi [, maxbval[, kernelsize[, width[, csfmask]]]])
            %
            % Isotropic smoothing of dwi data, mainly in order to reduce
            % the Gibbs ringing. It might be recommended to only smooth the
            % high SNR (or low b-valued) data in order not to alter the
            % Rice distribution of the low SNR data. This important to maintain 
            % the high accuracy of the WLLS The max bval can be set (default 
            % 1000 s/mm^2). The kernelsize is by default [5 5]. The width, 
            % i.e. the FWHM in voxels, is default 1.25.
            % if a csf mask is given as an additional argument, CSF
            % infiltration in microstructural signal is avoided during
            % smoothing.
            
            grad = DTI.b2grad(this.b); bval = grad(:, 4); clear grad;
            
            dwi = feval('double', dwi);
            dwi(dwi<=0)=eps;   
            if ~exist('maxbval', 'var') || isempty(maxbval)
                maxbval = 1000;
            end
           
            if ~exist('width', 'var') || isempty(width)
                width = 1.25;
            end
            if ~exist('dim', 'var') || isempty(dim)
                dim = '2d';
            end
            
            if ~exist('kernelsize', 'var') || isempty(kernelsize)
                if strcmpi(dim, '2d')
                    kernelsize = [5 5];
                else
                    kernelsize = [5 5 5];
                end
            elseif numel(kernelsize)==1
                if strcmpi(dim, '2d')
                    kernelsize = [kernselsize kernelsize];
                else
                    kernelsize = [kernselsize kernelsize kernelsize];
                end
            end 
            
            if strcmpi(dim, '2d') && numel(kernelsize)~=2
                    kernelsize = [5 5];
                    warning('incorrect kernelsize definition, default values used')
            end

            if strcmpi(dim, '3d') && numel(kernelsize)~=3
                    kernelsize = [5 5 5];
                    warning('incorrect kernelsize definition, default values used')
            end
            % smoothing performed in 2D!!!
            
            
            if exist('csfmask', 'var') && ~isempty(csfmask)
                
                bgmask = isnan(dwi(:,:,:,1));
                list = find(bval<=maxbval);
                for i= 1:numel(list)
                    
                    wmgm = dwi(:,:,:,list(i)); wmgm(csfmask) = NaN;
                    
                    if strcmpi(dim, '2d')
                        wmgm = nansmoothing(wmgm, kernelsize, width);
                    else
                        wmgm = nansmoothing3d(wmgm, kernelsize, width);
                    end
                    csf = dwi(:,:,:,list(i)); csf(~csfmask) = NaN;
                    if strcmpi(dim, '2d')
                        csf = nansmoothing(csf, kernelsize, width);
                    else
                        csf = nansmoothing3d(csf, kernelsize, width);
                    end
                    
                    total = nansum(cat(4, wmgm, csf), 4); 
                    total(bgmask) = NaN;
                    
                    dwi(:,:,:,list(i)) = total;    
                end  

                
            else
                list = find(bval<=maxbval);
                for i= 1:numel(list)
                    for j=1:size(dwi, 3)
                        if strcmpi(dim, '2d')
                            h = fspecial('gaussian', kernelsize, width/(2*sqrt(2*log(2)))); 
                            dwi(:,:,j,list(i)) = filter2(h, dwi(:,:,j,list(i)), 'same');  
                        else
                            sigma = width/(2*sqrt(2*log(2)));  
                            siz   = (kernelsize-1)/2;
                            [x,y, z] = ndgrid(-siz(1):siz(1),-siz(2):siz(2),-siz(3):siz(3));
                            h = exp(-(x.*x/2/sigma^2 + y.*y/2/sigma^2 + z.*z/2/sigma^2));
                            h = h/sum(h(:));
     
                            dwi(:,:,j,list(i)) = imfilter(dwi(:,:,j,list(i)),h);
                        end
                    end
                end  
            end 
         end
         
         function bmatrixrotation(this, params)
            % *.bmatrixrotation(params)
            % Rotation of the b-matrix (cf. Leemans and Jones, MRM 61(6):1336-49, 2009)
            % The function should only be used in combination with this
            % class' eddycorrect function. eddycorrect's "second output variable, that is, 
            % "params" is the (only) input variable of bmatrixrotation. The class property "bmatrix" 
            % is changed, no explicit output generated.
             
            sz = size(params);
            if sz(1)~= 12 || sz(2)~= size(this.b, 1)
                warning('incorrect variable params ... no bmatrix rotation done')
            end
            b_ = zeros(sz(2), 6);
            for i=1:sz(2)
                alpha = params(1, i); beta = params(2, i); gamma = params(3, i);
                Rx = [1 0 0; 0 cos(alpha) -sin(alpha); 0 sin(alpha) cos(alpha)];
                Ry = [cos(beta) 0 sin(beta); 0 1 0; -sin(beta) 0 cos(beta)];
                Rz = [cos(gamma) -sin(gamma) 0; sin(gamma) cos(gamma) 0; 0 0 1];
                R = Rx*Ry*Rz;
                
                B(1, 2) =  -this.b(i, 3)/2;
                B(1, 3) =  -this.b(i, 4)/2;
                B(2, 3) =  -this.b(i, 6)/2;
                
                B(2, 1) =  -this.b(i, 3)/2;
                B(3, 1) =  -this.b(i, 4)/2;
                B(3, 2) =  -this.b(i, 6)/2;
                
                B(1, 1) =  -this.b(i, 2);
                B(2, 2) =  -this.b(i, 5);
                B(3, 3) =  -this.b(i, 7);

                
                B = R*B*R';
                
                b_(i, 1) = B(1, 1);
                b_(i, 2) = B(1, 2);
                b_(i, 3) = B(1, 3);
                b_(i, 4) = B(2, 2);
                b_(i, 5) = B(2, 3);
                b_(i, 6) = B(3, 3);            
            end
            this.setGrad(b_);
         end
         
         function tdi(this, root, target, shell, ntracts, resolution, lmax) 
             
             
             if ~exist('resolution', 'var') || isempty(resolution)         
                 resolution = .5;
             end
             if ~exist('ntracts', 'var') || isempty(ntracts)         
                 ntracts = 1000000;
             end
             
             if ~exist('target', 'var') || isempty(target)         
                 error('target directory must be provided');
             end
             
             nii = load_untouch_nii(root); dwi = nii.img;
             
             currentdir = pwd;  mkdir(target); cd(target)
             
             
             grad = DTI.b2grad(this.b); 
             bval = grad(:, 4); bval = 100*round((bval+1)/100);          % it's multishell. Effective b-values will not work
             bvec = grad(:, 1:3); bvec(:, [2 3]) = -bvec(:, [2 3]);  % different software packages ... different conventions
             
             if ~exist('shell', 'var') || isempty(shell)
                 shell = max(bval);
             end
             
             mask = BrainMask(dwi(:,:,:,bval==min(bval)), 75);
             DTI.savenii(mask, fullfile(target, ['mask.nii']), root);
             
             list = [find(bval==0); find(bval==shell)];
             dwi = dwi(:,:,:,list);    DTI.savenii(dwi, fullfile(target, ['dwi.nii']), root)
             b = [bvec(list, :), bval(list)];  save(fullfile(target, 'grad.txt'), 'b', '-ascii');
             
             n = nnz(bval==shell);             
             if ~exist('lmax', 'var') || isempty(lmax)         
                 lmax = floor(2.*(floor(sqrt(1+8.*n)-3)./4));
                 if mod(lmax, 2)==1
                     lmax = lmax - 1;
                 end
             end
             
             
             
             system('dwi2tensor dwi.nii -grad grad.txt dt.nii')
             system('tensor2FA dt.nii - | mrmult - mask.nii fa.nii');
             system('tensor2vector dt.nii - | mrmult - fa.nii ev.nii');
             system('erode mask.nii -npass 3 - | mrmult fa.nii - - | threshold - -abs 0.7 sf.nii');
             system('estimate_response dwi.nii -grad grad.txt sf.nii response.txt')
             system(['csdeconv dwi.nii -grad grad.txt response.txt -lmax ', num2str(lmax),' -mask mask.nii CSD.nii'])
             system(['streamtrack SD_PROB CSD.nii -seed mask.nii -mask mask.nii whole_brain.tck -num ', num2str(ntracts)])
             disp('done')
             system(['tracks2prob whole_brain.tck -vox ', num2str(resolution), ' tdi.nii'])
             system(['tracks2prob whole_brain.tck -vox ', num2str(resolution), ' -colour dectdi.nii'])
             cd(currentdir)
         end
         
         
         function noisemap = getnoisemap(this, dwi, bshell, brainmask)
            % noisemap = *.getnoisemap(dwi[, bshell[, brainmask]])
            % 
            % Calculation of 3D noise map of a Rician distributed
            % single-shell diffusion weighted data set. If a multi-shell
            % data set is provided ("dwi"), the outer shell is by default
            % selected, although the shell might also be user defined
            % ("bshell"). The method requires a brainmask, if not provided,
            % a rudimentary mask will be computed.
            %
            % The algorithm is based on Veraart et al., MRM 70(4):972-984, 2013 
            dwi = feval('double', dwi); 
           
            grad = DTI.b2grad(this.b); bval = grad(:, 4); grad = grad(:, 1:3);
            if ~exist('brainmask', 'var') || isempty(brainmask)
                [~, brainmask] = this.skullstrip(dwi); 
            end
            if ~exist('bshell', 'var') || isempty(bshell)
                bshell = max(bval);
            end
            
            % calculate mean based on spherical harmonics
            [x, y, z, ~] = size(dwi);
            mask = true(x, y, z);

            mn_ = [];
            List = [];
            for bvals = bshell
                list = abs(bval-bvals) < 10;
                dwi = vec(dwi, mask);
                
                lmax = min(6, SH.maxlmax(nnz(list)));
                lmax  =lmax - mod(lmax, 2);
                
                sh = SH(grad(list, :), lmax, 'lls');
                coef = sh.coef(dwi(list, :));
                
                mn = sh.amp(coef);
                mn = unvec(mn, mask);
                
                mn_ = cat(4, mn_, mn); clear mn
                List = cat(1, List, list); 
            end;
            
            dwi = dwi(List==1, :); dwi = unvec(dwi, mask);
            
            [x1, y1, z1] = meshgrid(0.5:y+0.5, 0.5:x+0.5, 0.5:z+0.5);
            [xi, yi, zi] = meshgrid(1:y,1:x,1:z);
            
            n = size(mn_, 4);
            mn = zeros([size(x1), n]);

            for i=1:n
                mn(:,:,:,i) = interpn(yi, xi, zi, mn_(:,:,:,i), y1, x1, z1, 'linear');
            end;
            clear mn_

            
            % calculate deviations based on wavelet decomposition
            hhh=zeros(x+1,y+1,z+1, n);
            for k=1:n
                level=2;
                wt1 = dwt3J(dwi(:,:,:,k),'db1',level);
                hhh(:,:,:,k)=wt1.dec{2,2,2};
            end;
            sd = zeros(x+1, y+1, z+1, n);
            for k=1:n
                tmp = medfilt3(abs(hhh(:,:,:,k)), [1 1 1])/0.6745;
                sd(:,:,:,k) = tmp;
            end;
            
            % apply Koay's correction term (and some filtering to smooth the map)
            snr = mn(2:x, 2:y, 2:z, :)./sd(2:x, 2:y, 2:z, :);
            mask = true([x-1, y-1, z-1]);
            snr = vec(snr, mask);
            snr(snr<1.9131) = 1.9131;
            
            [corrTerm, ~] = rootFinder(double(snr), 1);
            sd = vec(sd(2:x, 2:y, 2:z, :), mask);
            
            noisemap =sd./sqrt(corrTerm);
            noisemap = nanmedian(unvec(noisemap, mask), 4);
           
            [x1, y1, z1] = meshgrid(1.5:y, 1.5:x, 1.5:z);
            [xi, yi, zi] = meshgrid(1:y,1:x,1:z);

            noisemap = interpn(y1, x1, z1, noisemap, yi, xi, zi, 'linear');

      
            
 %           outliers = prctile(noisemap(brainmask(:)), 95);
 %           noisemap(noisemap>outliers)=NaN;
            
%            noisemap(~brainmask) = NaN;
%            noisemap = nanmedian(noisemap, 3);
%            noisemap(isnan(noisemap)) = 0;
%             h = fspecial('average', [7 7]);
%             noisemap = imfilter(noisemap, h, 'replicate');
%             noisemap = repmat(noisemap, [1 1 z]);  
         end
    end
    
    methods (Access = private)
        function dt = dt_lls(this, dwi, outliers)
            
             %dt = this.b\log(dwi);
             if ~exist('outliers', 'var') || isempty(outliers)
               outliers = isnan(dwi);
             end
             outliers_ = sum(outliers, 1);
             list.out_ = outliers_ >0;
             list.in_ =  outliers_ ==0;
            
             npars = size(this.b, 2);
             nvxls = size(dwi, 2);
             dt = zeros(npars, nvxls);
             
             dt(:, list.in_) = pinv(this.b)*log(dwi(:, list.in_));
             
             if nnz(list.out_) > 0
                for j = find(list.out_)
                        in_ = outliers(:, j) == 0;
                        dt(:, j) = pinv(this.b(in_, :))*log(dwi(in_, j));
               end
            end
        end
        
        function dt = dt_wlls(this, dwi, outliers)
            if ~exist('outliers', 'var') || isempty(outliers)
               outliers = isnan(dwi);
            end
            
            outliers_ = sum(outliers, 1);
            list.out_ = outliers_ >0;
            list.in_ =  outliers_ ==0;
                                    
            dt = this.dt_lls(dwi, outliers);
            
            for i = 1:this.estimator.wls_iter
               dt(:, list.in_) = mex_wls(this.b,log(dwi(:, list.in_)),exp(this.b*dt(:, list.in_)));
            end
            if nnz(list.out_) > 0
                for i = 1:this.estimator.wls_iter
                    for j=find(list.out_)
                        in_ = outliers(:, j) == 0;
                         if ~isempty(this.b(in_, :)) && cond(this.b(in_, :))<1e10 
                            dt(:, j) = mex_wls(this.b(in_, :),log(dwi(in_, j)),exp(this.b(in_, :)*dt(:, j)));
                        else
                            dt(:,j) = NaN;
                        end
                    end
                end
            end
        end
        
        function dt = dt_nls(this, dwi, outliers, weights)
            if ~exist('outliers', 'var') || isempty(outliers)
               outliers = isnan(dwi);
            end
            if ~exist('weights', 'var') || isempty(weights)
               weights = ones([size(this.b, 1),1]);
            end
            dt = this.dt_wlls(dwi, outliers);
            bmat = double(this.b);
            options = optimset('lsqnonlin'); options = optimset(options,'Jacobian','on','TolFun',this.estimator.tolFun,'TolX',this.estimator.tolX,'MaxIter',this.estimator.maxIter,'Display','off');
            parfor i = 1:size(dwi,2)
                in_ = outliers(:, i) == 0;
                dt(:,i) = lsqnonlin(@(x)DTI.residual(x,bmat(in_, :),double(dwi(in_,i)),double(weights(in_))),double(dt(:,i)),[],[],options);
            end
        end
        
        function dt = dt_cls(this, dwi, outliers, sigma)
            
            if ~exist('outliers', 'var') || isempty(outliers)
               outliers = isnan(dwi);
            end
            if ~exist('sigma', 'var') || isempty(sigma)
                error('noise lelvels must be provided in order to use conditional least squares estimator')
            end
                          
            dt = this.dt_wlls(dwi, outliers);
            bmat = double(this.b);
            %options = optimset('lsqnonlin'); options = optimset(options, 'Algorithm', {'levenberg-marquardt',.01}, 'Jacobian','on','TolFun',this.estimator.tolFun,'TolX',this.estimator.tolX,'MaxIter',this.estimator.maxIter, 'TypicalX', [1; ones(6,1)*1e-4], 'Display','off');
            
            TypX = [1; ones(6,1)*1e-4];
            if size(dt, 1) == 22
                TypX = [1; ones(6,1)*1e-6; ones(15, 1)*1e-6];
            end
            options = optimset('Algorithm', {'levenberg-marquardt',.01}, 'Jacobian', 'on', 'Display','off', 'DerivativeCheck', 'off', 'Diagnostics', 'off', 'TolFun', this.estimator.tolFun, 'MaxFunEvals' , 7000, 'MaxIter' , this.estimator.maxIter, 'TolX', this.estimator.tolX, 'TypicalX', TypX);
            nvxls = size(dwi, 2);
            if isscalar(sigma)
               sigma = sigma*ones([nvxls, 1]);
            end
            
            dt_ = dt;
            parfor i = 1:nvxls
                in_ = outliers(:, i) == 0;
                [dt_(:,i)] = lsqnonlin(@(x)DTI.condresidual(x,bmat(in_, :),double(dwi(in_,i)), double(sigma(i))),double(dt(:,i)),[],[],options);
            end
            dt = dt_;
        end
        
        function dt = dt_ccls(this, dwi, outliers, sigma)
            
            if ~exist('outliers', 'var') || isempty(outliers)
               outliers = isnan(dwi);
            end
            if ~exist('sigma', 'var') || isempty(sigma)
                error('noise lelvels must be provided in order to use conditional least squares estimator')
            end
                          
            dt = this.dt_wlls(dwi, outliers);
            [eigval, eigvec] = DTI.eig(dt(2:7,:)); eigval(eigval<eps)=eps;
            for i = 1:size(dwi,2)
                c = reshape(eigvec(:,i),[3 3])*diag(eigval(:,i))*reshape(eigvec(:,i),[3 3])';
                c = chol(c);
                dt(2:7,i) = c([1; 4; 7; 5; 8; 9]);
            end
            bmat = double(this.b);
            
            TypX = [1; ones(6,1)*1e-4];
            if size(dt, 1) == 22
                TypX = [1; ones(6,1)*1e-6; ones(15, 1)*1e-6];
            end
            
            
            %options = optimset('lsqnonlin'); options = optimset(options, 'Algorithm', {'levenberg-marquardt',.01}, 'Jacobian','on','TolFun',this.estimator.tolFun,'TolX',this.estimator.tolX,'MaxIter',this.estimator.maxIter, 'TypicalX', [1; ones(6,1)*1e-4], 'Display','off');
            options = optimset('Algorithm', {'levenberg-marquardt',.01}, 'Jacobian', 'on', 'Display','off', 'DerivativeCheck', 'off', 'Diagnostics', 'off', 'TolFun', this.estimator.tolFun, 'MaxFunEvals' , 7000, 'MaxIter' , this.estimator.maxIter, 'TolX', this.estimator.tolX, 'TypicalX', TypX);
            nvxls = size(dwi, 2);
            if isscalar(sigma)
               sigma = sigma*ones([nvxls, 1]);
            end
            
            dt_ = dt;
            parfor i = 1:nvxls
                in_ = outliers(:, i) == 0;
                [dt_(:,i)] = lsqnonlin(@(x)DTI.condcholresidual(x,bmat(in_, :), double(dwi(in_,i)), double(sigma(i))),double(dt(:,i)),[],[],options);
            end
            dt_(1:7,:) = [dt_(1,:); dt_(2,:).^2; dt_(2,:).*dt_(3,:); dt_(2,:).*dt_(4,:); dt_(3,:).^2+dt_(5,:).^2; dt_(3,:).*dt_(4,:)+dt_(5,:).*dt_(6,:); dt_(4,:).^2+dt_(6,:).^2+dt_(7,:).^2];
            dt = dt_; 
        end
        
        function dt = dt_cnls(this, dwi, outliers)
            if ~exist('outliers', 'var') || isempty(outliers)
               outliers = isnan(dwi);
            end
            dt = this.dt_wlls(dwi, outliers);
            [eigval, eigvec] = DTI.eig(dt(2:7,:)); eigval(eigval<eps)=eps;
            for i = 1:size(dwi,2)
                c = reshape(eigvec(:,i),[3 3])*diag(eigval(:,i))*reshape(eigvec(:,i),[3 3])';
                c = chol(c);
                dt(2:7,i) = c([1; 4; 7; 5; 8; 9]);
            end
            bmat = double(this.b);
            options = optimset('lsqnonlin'); options = optimset(options,'Jacobian','on','TolFun',this.estimator.tolFun,'TolX',this.estimator.tolX,'MaxIter',this.estimator.maxIter,'Display','off');
            parfor i = 1:size(dwi,2)
                in_ = outliers(:, i) == 0;
                dt(:,i) = lsqnonlin(@(x)DTI.chol_residual(x,bmat(in_, :),double(dwi(in_,i))),double(dt(:,i)),[],[],options);
            end
            dt(1:7,:) = [dt(1,:); dt(2,:).^2; dt(2,:).*dt(3,:); dt(2,:).*dt(4,:); dt(3,:).^2+dt(5,:).^2; dt(3,:).*dt(4,:)+dt(5,:).*dt(6,:); dt(4,:).^2+dt(6,:).^2+dt(7,:).^2];
        end
    end
    
    methods (Access = public, Static = true)

        function [eigval, eigvec] = eig(dt)
            % [eigval, eigvec] = DTI.eig(dt)
            %
            % Eigenvalue decomposition of the diffuiosn tensor "dt". dt can
            % be of size [nparams, nvxls] or [x, y, z, nparams]. Output arguments 
            % are sized accordingly.
            
            if ndims(dt) ~= 2; [dt, mask] = vec(dt); end;
            dt = dt(1:6,:);
            [eigval, eigvec] = mex_dti_eig(dt);
            eigvec = eigvec([7 8 9 4 5 6 1 2 3],:);
            eigval = eigval([3 2 1],:);
            %eigval(eigval<0) = 0;

            if exist('mask','var'); eigval = unvec(eigval, mask); eigvec = unvec(eigvec, mask); end;
        end
        
        function fa = fa(eigval)
            % fa = DTI.fa(eigval)
            %
            % Calculation of the fractional anisotropy "fa" using the eigenvalues "eigval". eigval can
            % be of size [3, nvxls] or [x, y, z, 3]. Output arguments 
            % are sized accordingly.
            
            if ndims(eigval) ~= 2; [eigval, mask] = vec(eigval); end;
            l1 = eigval(1,:); l2 = eigval(2,:); l3 = eigval(3,:);
            fa = sqrt(1/2).*sqrt((l1-l2).^2+(l2-l3).^2+(l3-l1).^2)./sqrt(l1.^2+l2.^2+l3.^2);
            if exist('mask','var'); fa = unvec(fa, mask); end;
        end
        
        function md = md(eigval)
            % md = DTI.md(eigval)
            %
            % Calculation of the mean diffusivity "md" using the eigenvalues "eigval". eigval can
            % be of size [3, nvxls] or [x, y, z, 3]. Output arguments 
            % are sized accordingly.
            if ndims(eigval) ~= 2; [eigval, mask] = vec(eigval); end;
            md = mean(eigval,1);
            if exist('mask','var'); md = unvec(md, mask); end;
        end
        
        function ad = ad(eigval)
            % ad = DTI.ad(eigval)
            %
            % Calculation of the axial diffusivity "ad" using the eigenvalues "eigval". eigval can
            % be of size [3, nvxls] or [x, y, z, 3]. Output arguments 
            % are sized accordingly.
            
            if ndims(eigval) ~= 2; [eigval, mask] = vec(eigval); end;
            ad = eigval(1,:);
            if exist('mask','var'); ad = unvec(ad, mask); end;
        end

        function l2 = l2(eigval)
            % l2 = DTI.l2(eigval)
            %
            % "L2" using the eigenvalues "eigval". eigval can
            % be of size [3, nvxls] or [x, y, z, 3]. Output arguments 
            % are sized accordingly.
            
            if ndims(eigval) ~= 2; [eigval, mask] = vec(eigval); end;
            l2 = eigval(2,:);
            if exist('mask','var'); l2 = unvec(l2, mask); end;
        end

        function l3 = l3(eigval)
            % l3 = DTI.l3(eigval)
            %
            % "L3" using the eigenvalues "eigval". eigval can
            % be of size [3, nvxls] or [x, y, z, 3]. Output arguments 
            % are sized accordingly.
            
            if ndims(eigval) ~= 2; [eigval, mask] = vec(eigval); end;
            l3 = eigval(3,:);
            if exist('mask','var'); l3 = unvec(l3, mask); end;
        end
        
        function rd = rd(eigval)
            % rd = DTI.rd(eigval)
            %
            % Calculation of the radial diffusivity "rd" using the eigenvalues "eigval". eigval can
            % be of size [3, nvxls] or [x, y, z, 3]. Output arguments 
            % are sized accordingly.
            if ndims(eigval) ~= 2; [eigval, mask] = vec(eigval); end;
            rd = mean(eigval(2:3,:),1);
            if exist('mask','var'); rd = unvec(rd, mask); end;
        end
        
        function fefa = fefa(fe, fa)
            % fefa = DTI.fefa(fe, fa)
            %
            % Calculation of the color encoded fractional anisotropy "fefa" using the first eigenvalue "fe" and the fractional anisotropy "fa". fe can
            % be of size [3 (or 9), nvxls] or [x, y, z, 3 (or 9)]. fa must be sized accordingly. The output argument is also 
            % sized accordingly.
            if ndims(fe) ~= 2; [fe, mask] = vec(fe); fa = vec(fa,mask); end;
            fefa = abs(fe(1:3, :)).*fa([1 1 1],:);
            if exist('mask','var'); fefa = unvec(fefa, mask); end;
        end
        
        function adc = adc(dt, dir)
            % adc = DTI.adc(dt, dir)
            %
            % Calculation of the apparent diffusion coeffient "adc" along user-defined directions "dir"  based on the diffusion tensor "dt". dt can
            % be of size [npars, nvxls] or [x, y, z, npars]. The output argument is 
            % sized accordingly.
            if ndims(dt) ~= 2; [dt, mask] = vec(dt); end; 
            adc = (dir(:,DTI.ind(1:6,1)).*dir(:,DTI.ind(1:6,2))) * diag(DTI.cnt) * dt;
            if exist('mask','var'); adc = unvec(adc, mask); end;
        end
       
        function cl = cl(eigval)
            % cl = DTI.cl(eigval)
            %
            % Calculation of the linear diffusion "cl" using the eigenvalues "eigval". eigval can
            % be of size [3, nvxls] or [x, y, z, 3]. Output argument is sized accordingly.
            % Linear diffusion is one of the Westin metrics (Westin et al. ISMRM 1997)
            
            if ndims(eigval) ~= 2; [eigval, mask] = vec(eigval); end;
            cl = (eigval(1,:)-eigval(2,:))./eigval(1,:);
            if exist('mask','var'); cl = unvec(cl, mask); end;
        end
        
        function cp = cp(eigval)
            % cp = DTI.cp(eigval)
            %
            % Calculation of the planar diffusion "cl" using the eigenvalues "eigval". eigval can
            % be of size [3, nvxls] or [x, y, z, 3]. Output argument is sized accordingly.
            % Planar diffusion is one of the Westin metrics (Westin et al. ISMRM 1997)
            
            if ndims(eigval) ~= 2; [eigval, mask] = vec(eigval); end;
            cp = (eigval(2,:)-eigval(3,:))./eigval(1,:);
            if exist('mask','var'); cp = unvec(cp, mask); end;
        end
        
        function cs = cs(eigval)
            % cs = DTI.cs(eigval)
            %
            % Calculation of the spherical diffusion "cs" using the eigenvalues "eigval". eigval can
            % be of size [3, nvxls] or [x, y, z, 3]. Output argument is sized accordingly.
            % Spherical diffusion is one of the Westin metrics (Westin et al. ISMRM 1997)
            
            if ndims(eigval) ~= 2; [eigval, mask] = vec(eigval); end;
            cs = eigval(3,:)./eigval(1,:);
            if exist('mask','var'); cs = unvec(cs, mask); end;
        end
        

        
        function dt = sim_dt(fa, diffusivity, el, az)
            
            a = fa/sqrt(3-2*fa^2);
            dt_ = diffusivity*[ 1-a 0 0; 0 1-a 0; 0 0 1+2*a ];
            s = size(el);
            sinaz = sin(az);
            cosaz = cos(az);
            sinel = sin(el);
            cosel = cos(el);
            zeross = zeros(s);
            oness = ones(s);
            
            R_az = reshape([cosaz -sinaz zeross; sinaz cosaz zeross; zeross zeross oness],[s(1) 3 3]);
            R_el = reshape([cosel zeross sinel; zeross oness zeross; -sinel zeross cosel],[s(1) 3 3]);
            
            dt = zeros(6, s(1));
            for i = 1:s(1)
                R_az_i = squeeze(R_az(i,:,:));
                R_el_i = squeeze(R_el(i,:,:));
                dt_i = R_az_i*R_el_i*dt_*R_el_i'*R_az_i';
                dt(:,i) = dt_i([1 2 3 5 6 9]);
            end
        end
        
        function b = grad2b(grad)
            s = size(grad);
            b = [ones([s(1) 1],class(grad)) -(grad(:,ones(1,6)*4).*grad(:,DTI.ind(1:6,1)).*grad(:,DTI.ind(1:6,2)))*diag(DTI.cnt)];
        end
        
        function grad = b2grad(b)
            if all(b(:,1)==1)      %% interal b-matrix format [1 -b]
                grad(:,4) = sum(-b(:,[2 5 7]),2);
                grad(:,1:3) = sqrt(-b(:,[2 5 7])./repmat(grad(:,4),[1 3]));
                
                sx = ones(size(b,1), 1);
                sy = sign(-b(:, 3)); sy(sy==0) = 1;
                sz = sign(-b(:, 4)); sz(sz==0) = sign(-b(sz==0, 6)).*sy(sz==0); sz(sz==0) = 1;
                
                grad(:, 1:3) = grad(:, 1:3)./[sx, sy, sz];
                
                grad(isnan(grad)) = 0;
            else                  %% original b-matrix format [b]
                grad(:,4) = sum(b(:,[1 4 6]),2);
                grad(:,1:3) = sqrt(b(:,[1 4 6])./repmat(grad(:,4),[1 3]));
                
                sx = ones(size(b,1), 1);
                sy = sign(b(:, 2)); sy(sy==0) = 1;
                sz = sign(b(:, 3)); sz(sz==0) = sign(b(sz==0, 5)).*sy(sz==0); sz(sz==0) = 1;
                
                grad(:, 1:3) = grad(:, 1:3)./[sx, sy, sz];
                
                
                grad(isnan(grad)) = 0;
            end
        end
        
%         function [dt, flag] = changeTensorOrder(dt, type)
%             % dt = DTI.changeTensorOrder(dt, type)
%             % This static function can be used to switch between two
%             % different ways to store the diffusion tensor.
%             % Currently, tensors are stored differently in Antwerp and NYU,
%             % both in order of tensor elements and units. The DTI and DKI
%             % class internally use the Antwerp format.
%             %
%             % input
%             % dt: diffusion tensor
%             % type: 'Antwerp' or 'NYU'
%             %   
%             % Antwerp: [Dxx Dxy Dxz Dyy Dyx Dzz], D in mm^2/s
%             %
%             % NYU:     [Dxx Dyy Dzz Dxy Dxz Dyz], D in um^2/ms 
%             if ndims(dt) ~= 2; [dt, mask] = vec(dt); end;
% 
%             vals = sum(abs(dt(1:6, :)))/6;
%             OoM = floor(log(abs(vals))./log(10));
% 
%             if median(OoM) <= -2
%                 currenttype = 'Antwerp';
%             else
%                 currenttype = 'NYU';
%             end
% 
%             [~, nvxls] = size(dt);
%             flag = false;
%             if ~strcmpi(currenttype, type)
%                 switch type
%                     case 'Antwerp'
%                       neworder = [1 4 5 2 6 3]';
%                       scaling  = [1e-3 1e-3 1e-3 1e-3 1e-3 1e-3]';  
%                     case 'NYU'
%                       neworder = [1 4 6 2 3 5]';
%                       scaling  = [1e3 1e3 1e3 1e3 1e3 1e3]';
%                 end
%                 flag = true;
%                 dt = dt(neworder(1:6), :).*repmat(scaling(1:6), [1 nvxls]); 
%             end
%             
%             if exist('mask','var'); dt = unvec(dt,mask); end;
%         end
        
        function savenii(var, fname, example, bgvalue)
            % DTI.savenii(var, fname, example)
            %
            % Saves variable "var", sized [x, y, z, m], as a NiFTi file.
            % NiFTi filename "fname" must be provided by the user. A NiFTi
            % template must be given as well "example". Typically, this is
            % the original diffusion-weigthed data to ensure all maps are
            % aligned with the original data. variables are converted to
            % singles.
            
            if ~exist('bgvalue', 'var') || isempty(bgvalue)
               bgvalue = 0;
            end
            
            var(isnan(var)) = bgvalue;
            
            nii = load_untouch_nii(example);
            nii.hdr.dime.dim(2) = size(var, 1);
            nii.hdr.dime.dim(3) = size(var, 2);
            nii.hdr.dime.dim(4) = size(var, 3);
            nii.hdr.dime.dim(5) = size(var, 4);
            nii.hdr.dime.datatype = 16;
            nii.hdr.dime.bitpix =  32;
            nii.img = feval('single',var);
            
            save_untouch_nii(nii, fname);
            
        end
    end
        
    methods (Access = public, Static = true)
        function [F,J] = residual(dt,b,dwi, weights)
            if ~exist('weights', 'var') || isempty(weights)
                weights = ones([size(b, 1),1]);
            end;
            
            dwi_hat = exp(b*dt);
            F = diag(weights)*(dwi_hat-dwi);
            if nargout > 1
                J = diag(weights)*dwi_hat(:,ones(1, size(b,2))).*b;
            end
        end
        
       
        
        function [F,J] = chol_residual(x,b,dwi)
            dt = x;
            dt(1:7) = [x(1);
                       x(2)^2;
                       x(2)*x(3);
                       x(2)*x(4);
                       x(3)^2+x(5)^2;
                       x(3)*x(4)+x(5)*x(6);
                       x(4)^2+x(6)^2+x(7)^2];
            dwi_hat = exp(b*dt);
            F = dwi_hat-dwi;
            if nargout > 1
                d_dt = eye(numel(dt));
                d_dt(1:7, 1:7) =...
                    [1      0      0      0      0      0      0;
                     0      2*x(2) x(3)   x(4)   0      0      0;
                     0      0      x(2)   0      2*x(3) x(4)   0;
                     0      0      0      x(2)   0      x(3)   2*x(4);
                     0      0      0      0      2*x(5) x(6)   0;
                     0      0      0      0      0      x(5)   2*x(6);
                     0      0      0      0      0      0      2*x(7)];
                J = dwi_hat(:,ones(1,size(b,2))).*b*d_dt;
            end

        end

        
        function [F,J] = condcholresidual(x,b,dwi, sigma)
            
            dt = x;
            dt(1:7) = [x(1);
                       x(2)^2;
                       x(2)*x(3);
                       x(2)*x(4);
                       x(3)^2+x(5)^2;
                       x(3)*x(4)+x(5)*x(6);
                       x(4)^2+x(6)^2+x(7)^2];
                   
            dwi_hat = exp(b*dt);
            [mu, dmu] = chimean(dwi_hat, sigma, 1, [true false]); 
                nanlocs = find(isnan(mu) | ~isfinite(mu));
                mu(nanlocs) = dwi_hat(nanlocs);
                dmu(nanlocs,1) = 1;
            F = double(mu-dwi);
            
            if nargout > 1
                
                
                d_dt = eye(numel(dt));
                d_dt(1:7, 1:7) =...
                    [1      0      0      0      0      0      0;
                     0      2*x(2) x(3)   x(4)   0      0      0;
                     0      0      x(2)   0      2*x(3) x(4)   0;
                     0      0      0      x(2)   0      x(3)   2*x(4);
                     0      0      0      0      2*x(5) x(6)   0;
                     0      0      0      0      0      x(5)   2*x(6);
                     0      0      0      0      0      0      2*x(7)];
                 
                 
                npars = numel(dt);
                ddwi_dt = dwi_hat(:, ones(1, npars)).*b*d_dt;                    
                dmu_ddwi = dmu(:,1);
                dF_ddwi =  dmu_ddwi;
                J =    dF_ddwi(:, ones(1, npars)) .* ddwi_dt;
            end
        end
    
        
                
        function [F,J] = condresidual(dt,b,dwi, sigma)
            
            dwi_hat = exp(b*dt);
            [mu, dmu] = chimean(dwi_hat, sigma, 1, [true false]); 
                nanlocs = find(isnan(mu) | ~isfinite(mu));
                mu(nanlocs) = dwi_hat(nanlocs);
                dmu(nanlocs,1) = 1;
            F = double(mu-dwi);
            
            if nargout > 1
                npars = numel(dt);
                ddwi_dt = dwi_hat(:, ones(1, npars)).*b;                    
                dmu_ddwi = dmu(:,1);
                dF_ddwi =  dmu_ddwi;
                J =    dF_ddwi(:, ones(1, npars)) .* ddwi_dt;
            end
        end
    
        function createparameterfile(fname)
            fid = fopen(fname, 'w+');
            fprintf(fid, '%s\n', '// Parameter file for Eddy Current distortion correction');
            fprintf(fid, '%s\n', '// Read Elastix documantation for definition of parameters');
            fprintf(fid, '%s\n', '(FixedInternalImagePixelType "float")');
            fprintf(fid, '%s\n', '(MovingInternalImagePixelType "float")');
            fprintf(fid, '%s\n', '(FixedImageDimension 3)');
            fprintf(fid, '%s\n', '(MovingImageDimension 3)');
            fprintf(fid, '%s\n', '(Registration "MultiResolutionRegistration")');
            fprintf(fid, '%s\n', '(FixedImagePyramid "FixedSmoothingImagePyramid")');
            fprintf(fid, '%s\n', '(MovingImagePyramid "MovingSmoothingImagePyramid")');
            fprintf(fid, '%s\n', '(Interpolator "BSplineInterpolator")');
            fprintf(fid, '%s\n', '(ResampleInterpolator "FinalBSplineInterpolator")');
            fprintf(fid, '%s\n', '(Resampler "DefaultResampler")');
            fprintf(fid, '%s\n', '(Optimizer "AdaptiveStochasticGradientDescent")');
            fprintf(fid, '%s\n', '(Transform "AffineDTITransform")'); %AffineDTITransform
            fprintf(fid, '%s\n', '(Scales -1 -1 -1 1-24 -1 1e24 1e24 -1 1e24 -1 -1 -1)');
            fprintf(fid, '%s\n', '(Metric "AdvancedMattesMutualInformation" )');
            fprintf(fid, '%s\n', '(AutomaticScalesEstimation "true")');
            fprintf(fid, '%s\n', '(AutomaticTransformInitialization "false")');
            fprintf(fid, '%s\n', '(NumberOfResolutions 2)');
            fprintf(fid, '%s\n', '(ImagePyramidSchedule 1 1 1 1 0 0)');
            fprintf(fid, '%s\n', '(WriteResultImage "true")');
            fprintf(fid, '%s\n', '(ResultImagePixelType "float")');
            fprintf(fid, '%s\n', '(ResultImageFormat "nii")');
            fprintf(fid, '%s\n', '(ErodeMask "false")');
            fprintf(fid, '%s\n', '(HowToCombineTransforms "Compose")');
            fprintf(fid, '%s\n', '(NumberOfSpatialSamples 8192)');
            fprintf(fid, '%s\n', '(NewSamplesEveryIteration "true")');
            fprintf(fid, '%s\n', '(ImageSampler "Random")');
            fprintf(fid, '%s\n', '(NumberOfHistogramBins 256)');
            fprintf(fid, '%s\n', '(FixedKernelBSplineOrder 3)');
            fprintf(fid, '%s\n', '(MovingKernelBSplineOrder 3)');
            fprintf(fid, '%s\n', '(BSplineInterpolationOrder 1)');
            fprintf(fid, '%s\n', '(FinalBSplineInterpolationOrder 3)');
            fprintf(fid, '%s\n', '(DefaultPixelValue 0)');
            fprintf(fid, '%s\n', '(MaximumNumberOfIterations 100 500 500 500)');
            fprintf(fid, '%s\n', '(SigmoidInitialTime 4.0)');
            fclose(fid);
        end    
   end
end