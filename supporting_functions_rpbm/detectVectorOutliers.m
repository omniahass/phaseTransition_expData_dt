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