function Mpix = pixvec2mat(vpix, nrows, ncol)
  % Mpix = pixvec2mat(vpix, nrows, ncol)
  %
  % Converts the (forced to be) column vector vpix into a matrix such that
  % The rows of Mpix are taken as contiguous chunks of vpix. I.e.,
  %
  % Mpix = [ vpix(1:ncol)';
  %          vpix(ncol+1:2*ncol)'
  %              :              ]
  % where ncol = length(vpix)/nrows;
  %
  % See Also: PixelVectorToMatrix, which is a less efficient implementation of
  % the same operation.

  % assure vpix is a column vector to start
  vpix = vpix(:);
  
  if nargin == 2
    Mpix = reshape(vpix, [], nrows)';
  else
    Mpix = reshape(vpix', nrows, ncol)';
  end
  
end