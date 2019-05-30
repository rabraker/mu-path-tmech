% [ I_fit ] = detrend_plane(I)
% 
% Fits a plane to a (raster) data set and removes it. If the data was CS
% sub-sampled, use detrend_sampled_plane.m instead. 
% inputs
%  -----
%  I : pix x pix image data.
function [ I_fit ] = detrend_plane(I)


[ypix, xpix] = size(I);

% Fit a plane to the data and remove. 
xs = [1:1:xpix];
X = reshape(repmat(xs, ypix,1),[],1);
ys = [1:1:ypix]';
Y = reshape(repmat(ys, xpix,1),[],1);

PHI = [X, Y, X*0 + 1;];


Z = reshape(I, [],1);

coeffs = PHI\Z;


mx = coeffs(1);
my = coeffs(2);
b = coeffs(3);

z_fit = X*mx + Y*my + b;

I_fit = reshape(Z-z_fit, ypix, xpix);

end

