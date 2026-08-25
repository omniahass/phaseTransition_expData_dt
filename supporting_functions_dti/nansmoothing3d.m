function B = nansmoothing3d(A, kernelsize, width, version)
    
     if ~exist('kernelsize', 'var') || isempty(kernelsize)
                kernelsize = [5 5 5];
     elseif numel(kernelsize)==1
                kernelsize = [kernelsize, kernelsize kernelsize];
     elseif numel(kernelsize)~=3
          warning('incorrect kernel size. Must contain 3 elements. [5 5 5] selected instead')
          kernelsize = [5 5 5];
     end 
     if ~exist('width', 'var') || isempty(width)
          width = 1.25;
     end
     if ~exist('version', 'var') || isempty(version)
          version = 2;
     end
     A_orig = A;
     bordersize = floor(kernelsize/2);
     [mx,my, mz] = size(A);
     sx = mx + 2*bordersize(1);
     sy = my + 2*bordersize(2);
     sz = mz + 2*bordersize(3);
     
     res = NaN(sx,sy, sy);
    
     idxx = sx/2+1-mx/2 : sx/2+mx/2;
     idxy = sy/2+1-my/2 : sy/2+my/2;
     idxz = sz/2+1-mz/2 : sz/2+mz/2;
    
     res(idxx,idxy, idxz) = A;
     A = res;

            
     sigma = width/(2*sqrt(2*log(2)));  
     siz   = (kernelsize-1)/2;
     [x,y, z] = ndgrid(-siz(1):siz(1),-siz(2):siz(2),-siz(3):siz(3));
     h = exp(-(x.*x/2/sigma^2 + y.*y/2/sigma^2 + z.*z/2/sigma^2));
     h = h/sum(h(:));
     
     
     switch version
         case 1
             x_ = 0; 
             for x = idxx
                 x_ = x_+1;
                 y_ = 0;
                 for y = idxy
                     y_ = y_+1;
                     z_ = 0;
                     for z = idxz
                         z_ = z_+1;
                         block = A(x-bordersize(1):x+bordersize(1), y-bordersize(2):y+bordersize(2), z-bordersize(3):z+bordersize(3));
                         results = nansum(block(:).*h(:));

                         h_ = h;
                         h_(isnan(block)) = 0;
                         n = sum(h_(:), 1);
                         B(x_, y_, z_) = results./(n+eps);

                     end
                 end
             end
         case 2
             k= 0;
             for x_ = 1:kernelsize(1)
                for y_ = 1:kernelsize(2)
                    for z_ = 1:kernelsize(3)
                        k=k+1;
                        tmp(:,:,:,k) = A(x_:x_+mx-1, y_:y_+my-1, z_:z_+mz-1);
                        
                    end
                end
             end
            
             segment = vec(tmp, true(mx, my, mz)); 

            H = repmat(h(:), [1 size(segment, 2)]);
            results = nansum(segment.*H);

            H_ = H;
            H_(isnan(segment)) = 0;
            n = sum(H_, 1);

            results = results./(n+eps);
            B = unvec(results, true(mx, my, mz));
     end
     B(isnan(A_orig)) = NaN;
end

