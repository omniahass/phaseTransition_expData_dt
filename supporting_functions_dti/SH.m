classdef SH < handle
    properties (GetAccess = public)
        dir;
        lmax;
        estimator;
        options;
        sh
    end
    properties (Access = private)
        fun;
    end
            
    methods (Access = public)
        function this = SH(dir, lmax, estimator, options)
            sdir = size(dir);
            if ~isfloat(dir) && sdir(2) ~= 3
                error('dir should be an nx3 matrix with floating point numbers');
            end
            this.dir = dir;
            this.lmax = SH.maxlmax(sdir(1));
            if exist('lmax','var') && ~isempty(lmax)
                if ~isscalar(lmax)  || lmax < 0 || lmax > this.lmax
                    error(['lmax should be even and in [' int2str(0) '-' int2str(this.lmax) ']']);
                end
                this.lmax = lmax;
            end
            this.sh = SH.eval(this.dir, this.lmax);
            if exist('estimator','var') && ~isempty(estimator)
                if exist('options','var') && ~isempty(options)
                    this.setEstimator(estimator, options);
                else
                    this.setEstimator(estimator);
                end
            else
                this.setEstimator('lls');
            end
        end
        
        function this = setEstimator(this, estimator, options)
            switch estimator
                case 'lls'
                    this.fun = @this.coef_lls;
                case 'mle'
                    this.fun = @this.coef_mle;
                    if nargin < 3
                        options = optimset('fminunc'); options = optimset(options,'GradObj','on','Hessian','off','Display','off', 'MaxFunEvals', 20000, 'MaxIter', 2000);
                    end
                    this.options = options;
                otherwise
                    error('estimator not supported');
            end
            this.estimator = estimator;
        end
        
        function varargout = coef(this, varargin)
            if ndims(varargin{1}) ~= 2; [varargin{1}, mask] = vec(varargin{1}); end;
            varargin{1} = feval(class(this.sh),varargin{1});
            if strcmp(this.estimator,'mle') && size(varargin,2)>1 && ~isempty(varargin{2}) &&~isscalar(varargin{2}) && exist('mask','var'); varargin{2} = vec(varargin{2},mask); end
            varargout{1} = this.fun(varargin{:});
            if exist('mask','var'); varargout{1} = unvec(varargout{1},mask); end;
        end
        
        function coef = coef_lls(this, amp, lambda)
            if nargin<3 
                coef = this.sh\amp;
            else
                B = this.sh;
                ord = [0:2:this.lmax];
                nOrd = ord*2 + 1;
                d = [];
                for i = 1:numel(nOrd)
                    d = [d; (ord(i)^2*(ord(i)+1)^2)*ones(nOrd(i), 1)];
                end
                L = diag(d);
                coef = (B'*B + lambda*L)\(B')*amp;
            end;
        end
        function coef = coef_wlls(this, amp, var)
            
            % Cx = d
            
%             for i = 1:size(amp, 2)
%                    y = amp(:,i);
%                    X = this.sh;
%                    W = diag(1./var);cond(W)
%                    C = X'*W*X; cond(C)
%                    d = X'*W*y; cond(d) 
%                    coef(1:size(this.sh,2), i) = lsqlin(C, d, -C, zeros([size(C, 1), 1]));
%             end;
            
           y = amp;
           X = this.sh;
           W = diag(1./var);
           C = X'*W*X;
           d = X'*W*y;  
           coef = C\d;
        end
        
        
        function coef = coef_wnls(this, amp, alpha, sigma)         
                    
           y = amp;
           X = this.sh;
           W = diag(1./alpha);
           C = X'*W*X;
           d = X'*W*y;  
           coef = C\d;
           
           sh_ = double(this.sh);
           
           for i=1:size(coef, 2)
            coef(:,i) = fminunc(@(x)SH.weightedresidual(x,sh_,double(amp(:,i)),double(alpha), double(sigma)),double(coef(:,i)),this.options);
           end;      
        end
        
                
        function coef = coef_mle(this, amp, sigma)
            coef = this.coef_lls(amp);
            sh_ = double(this.sh);
            if exist('sigma','var') && ~isempty(sigma)
                if isscalar(sigma)
                    sigma = repmat(sigma,[1 size(amp,2)]);
                end
                parfor i = 1:size(amp,2)
                    coef(:,i) = fminunc(@(x)SH.LogLikelihood(x,sh_,double(amp(:,i)),double(sigma(:,i))),double(coef(:,i)),this.options);
                end
            else
                amp_hat = this.sh*coef;
                resnrm = sum((amp-amp_hat).^2);
                sigma = sqrt(median(resnrm)/(size(amp,1)-size(this.sh,2)));
                coef(size(coef,1)+1,:) = sigma;
                parfor i = 1:size(amp,2)
                    coef(:,i) = fminunc(@(x)SH.LogLikelihood(x,sh_,double(amp(:,i))),double(coef(:,i)),this.options);
                end
            end
        end
        
        function amp = amp(this, coef)
            if ndims(coef) ~= 2; [coef, mask] = vec(coef); end;
            coef = feval(class(this.sh),coef);
            amp = this.sh*coef;
            if exist('mask','var'); amp = unvec(amp,mask); end;
        end
    end
    
    methods (Access = public, Static = true)
        function n = lmax2n(lmax)
            n = (lmax+1).*(lmax+2)./2;
        end
        
        function lmax = n2lmax(n)
            lmax = 2.*(floor(sqrt(1+8.*n)-3)./4);
        end
        
        function lmax = maxlmax(n)
            lmax = floor(SH.n2lmax(n));
%             if ~iseven(lmax)
%                 lmax = lmax - 1;
%             end
        end
        
        function sqresid = weightedresidual(x,sh_,amp,alpha,sigma)
                
                amp_hat = sh_*x;
                y = -((amp_hat).^2)/(2*sigma^2);
                I0 = besselmx(double('I'),0,-y/2,0);  
                I1 = besselmx(double('I'),1,-y/2,0);
                mu = (sigma*sqrt(pi/2)*exp(y/2).*((1-y).*I0 - y.*I1))';
                var = 2*sigma^2 + amp_hat'.^2 - mu.^2;
                var = alpha.*var';
                
                y = amp;
                X = sh_;
                W = diag(1./var);
                C = X'*W*X;
                d = X'*W*y;  
                coef = C\d;
                sqresid = sum((sh_*coef - amp).^2, 1);              
        end;
        
        
        function [f, df_del, df_daz, d2f_del2, d2f_daz2, d2f_deldaz] = eval(dir, lmax)
            classdir = class(dir);
            if size(dir,2) == 3
                dir = c2s(dir);
            end
            el = dir(:,1)';
            az = dir(:,2);
            
            num = size(dir,1);
            
            f = zeros(num,SH.lmax2n(lmax),classdir);
            if nargout > 1
                df_del = zeros(size(f),classdir);
                df_daz = zeros(size(f),classdir);
                if nargout > 3
                    d2f_del2 = zeros(size(f),classdir);
                    d2f_daz2 = zeros(size(f),classdir);
                    d2f_deldaz = zeros(size(f),classdir);
                end
            end
            sign = (-1).^(0:lmax);
            for l=0:2:lmax
                q = sign(ones(1,num),1:l+1);
                q = q.*legendre(l,cos(el),'sch')';
                q = q.*sqrt((2*l+1)/(4*pi));
                loff = l*(l+1)/2 + 1;
                for m=0:l
                    if m
                        cosmaz = cos(m*az); sinmaz = sin(m*az);
                        f(:,loff+m) = q(:,m+1).*cosmaz;
                        f(:,loff-m) = q(:,m+1).*sinmaz;
                        if nargout > 1
                            if m == 1
                                tmp = sqrt((l+m)*(l-m+1)*2).*q(:,m);
                            else
                                tmp = sqrt((l+m)*(l-m+1)  ).*q(:,m);
                            end
                            if m < l
                                tmp = tmp - sqrt((l-m)*(l+m+1)).*q(:,m+2);
                            end
                            tmp = -tmp/2;
                            if nargout > 3
                                tmp2 = -((l+m)*(l-m+1) + (l-m)*(l+m+1)).*q(:,m+1);
                                if (m == 1)
                                    tmp2 = tmp2 - ((l+1)*l).*q(:,2);
                                else
                                    if (m == 2)
                                        tmp2 = tmp2 + sqrt((l+m)*(l-m+1)*(l+m-1)*(l-m+2)*2).*q(:,m-1);
                                    else
                                        tmp2 = tmp2 + sqrt((l+m)*(l-m+1)*(l+m-1)*(l-m+2)).*q(:,m-1);
                                    end
                                end
                                if (l > m+1)
                                    tmp2 = tmp2 + sqrt((l-m)*(l+m+1)*(l-m-1)*(l+m+2)).*q(:,m+3);
                                end
                                tmp2 = tmp2/4.0;
                                
                                d2f_del2(:,loff+m) = tmp2.*cosmaz;
                                d2f_del2(:,loff-m) = tmp2.*sinmaz;
                                d2f_daz2(:,loff+m) = -q(:,m+1).*cosmaz.*m.*m;
                                d2f_daz2(:,loff-m) = -q(:,m+1).*sinmaz.*m.*m;
                                d2f_deldaz(:,loff+m) = -tmp.*sinmaz.*m;
                                d2f_deldaz(:,loff-m) =  tmp.*cosmaz.*m;
                            end
                            df_del(:,loff+m) = tmp.*cosmaz;
                            df_del(:,loff-m) = tmp.*sinmaz;
                            
                            df_daz(:,loff+m) = -q(:,m+1).*sinmaz.*m;
                            df_daz(:,loff-m) =  q(:,m+1).*cosmaz.*m;
                        end
                    else
                        f(:,loff) = q(:,1);
                        if l > 1
                            if nargout > 1
                                df_del(:,loff) = q(:,2).*sqrt(l*(l+1)/2);
                                if nargout > 3
                                    d2f_del2(:,loff) = (sqrt(l*(l+1)*(l-1)*(l+2)/2) * q(:,3) - l*(l+1) * q(:,1))/2;
                                end
                            end
                        end
                    end
                end
            end
        end
        
        function [f, grad, hess] = LogLikelihood(coef, X, amp, sigma)
            p = size(X,2);
            amp_hat = X*coef(1:p);
            sigma_not_supplied = nargin < 4;
            if (sigma_not_supplied); sigma = coef(end); end;

            if nargout > 1
                if nargout > 2
                    [f, grad, hess__] = logricepdf(amp, amp_hat, sigma, [false true sigma_not_supplied]);
                    if sigma_not_supplied
                        [j,i] = find(triu(ones(p))');
                        [j1,i1] = find(triu(ones(p+1))');
                        hess_ = zeros(size(X,1),size(j1,1));
                        m = (j1 < p+1) & (i1 < p+1);
                        hess_(:,m) = hess__(:,ones(sum(1:p),1)).*X(:,i).*X(:,j);
                        hess_(:,end) = hess__(:,3);
                        m = ~m; m(end) = 0;
                        hess_(:,m) = hess__(:,ones(p,1)*2).*X;
                        % sum over n
                        hess_ = -sum(hess_,1)';
                        % reorder in p x p
                        hess = zeros([p+1 p+1]);
                        liniu = sub2ind([p+1 p+1], i1, j1);
                        linil = sub2ind([p+1 p+1], j1, i1);
                        hess(liniu) = hess_;
                        hess(linil) = hess_; % p x p
                    else
                        [j,i] = find(triu(ones(p))');
                        hess_ = hess__(:,ones(sum(1:p),1)).*X(:,i).*X(:,j);
                        % sum over n
                        hess_ = -sum(hess_,1)';
                        % reorder in p x p
                        hess = zeros([p p]);
                        liniu = sub2ind([p p], i, j);
                        linil = sub2ind([p p], j, i);
                        hess(liniu) = hess_;
                        hess(linil) = hess_; % p x p
                    end
                else
                    [f, grad] = logricepdf(amp, amp_hat, sigma, [false true sigma_not_supplied]);
                end
                if sigma_not_supplied
                    grad = cat(2,grad(:,ones(p,1)).*X,grad(:,2));
                else
                    grad = grad(:,ones(p,1)).*X;
                end
                grad = -sum(grad,1)'; % p x 1
            else
                f = logricepdf(amp, amp_hat, sigma);
            end
            f = -sum(f,1)'; % 1 x 1
        end
        
        
    end
end
