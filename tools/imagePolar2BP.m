function BP = imagePolar2BP(image,varargin) 
%
%       creates beam pattern from image defined in polar coordinates 
%

[mC indC] = max(max(image)); 
if nargin > 1 
    BP = image(:,indC) ./ mC ; 
else
    BP = image(:,indC) ; 
end