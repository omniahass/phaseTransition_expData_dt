function S = c2s(C)  
%
% Converts cartesian coordinates to spherical coordinates.
%
% Copyright Ben Jeurissen (ben.jeurissen@ua.ac.be)
%
norm = sqrt(sum(C.^2, 2));
S(:,1) = acos(C(:,3)./norm);
S(:,2) = atan2(C(:,2), C(:,1));
end