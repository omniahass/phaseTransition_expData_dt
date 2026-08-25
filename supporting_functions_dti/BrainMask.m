function [mask, t] = BrainMask(b0,t)
%
% Calculates a rudimentary brain mask from the b0 image(s).
% The threshold t can be calculated automatically (value is returned)
% or manually (set it as input t).
%
% Copyright Ben Jeurissen (ben.jeurissen@ua.ac.be)
%

verbose = false; % set this to true if you want to see the histogram
b0 = mean(double(b0),4);
if nargin < 2
    min_ = min(b0(:)); max_ = max(b0(:));
    bins = min_:((max_-min_)/100):max_; bins = bins(2:end);
    histogram = histc(b0(:), bins);
    [peaks, loc] = findpeaks(-histogram,'npeaks',1);
    t = bins(loc(1));
    if verbose
        figure; bar(bins, histogram); hold all;
        line([t t], [0 max(histogram)],'color','red');
    end
end    
% threshold image
mask = logical(b0 > t);
% erode image
mask = imerode(mask,ones(3,3,3));
% keep the biggest connected component
CC = bwconncomp(mask,26);
[~,idx] = max(cellfun(@numel,CC.PixelIdxList));
mask = false(size(mask)); mask(CC.PixelIdxList{idx}) = true;
% dilate image
mask = imdilate(mask,ones(3,3,3));
% dilate image
mask = imdilate(mask,ones(3,3,3));
% fill the holes
mask = imfill(mask,'holes');
% erode image
mask = imerode(mask,ones(3,3,3));
% mask = imerode(mask,ones(3,3,3));
% mask = imerode(mask,ones(3,3,3));
end