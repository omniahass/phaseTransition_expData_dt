function B = nansmoothing(A, kernelsize, width)
    
     if ~exist('kernelsize', 'var') || isempty(kernelsize)
                kernelsize = [5 5];
     elseif numel(kernelsize)==1
                kernelsize = [kernelsize, kernelsize];
     elseif numel(kernelsize)~=2
          warning('incorrect kernel size. Must contain 2 elements. [5 5] selected instead')
          kernelsize = [5 5];
     end 
     if ~exist('width', 'var') || isempty(width)
          width = 1.25;
     end
     

     for i=1:size(A, 3)
        B(:,:,i) = colfilt(A(:,:,i),[5 5],'sliding',@(x)nansmoothing_f(x, kernelsize, width));
     end
     
     B(isnan(A)) = NaN;
end

function results = nansmoothing_f(segment, kernelsize, width)
    
    % normal filtering ... NaNs are ignored.
    h = fspecial('gaussian', kernelsize, width/(2*sqrt(2*log(2)))); 
    H = repmat(h(:), [1 size(segment, 2)]);
    results = nansum(segment.*H);
    
    % by just ignoring the NaNs the filter no longer adds up to 1! i
    % correct for this by applying a normalization term. i counted the
    % percentage of the filter that has not been applied because of NaNs 
    % in the data.

    H_ = H;
    H_(isnan(segment)) = 0;
    n = sum(H_, 1);
    
    results = results./(n+eps);
end

