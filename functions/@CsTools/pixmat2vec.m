function vpix = pixmat2vec(Mat)
% Convert the matrix X to a single vector of pixels, going along the rows of X.
% 
% In other words,
% v = [X(1,:), X(2,:), ...]'

  [n,m] = size(Mat);
    
  vpix = reshape(Mat', n*m,1);
  
end
