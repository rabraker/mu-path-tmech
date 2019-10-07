

function Im_fit = norm_align(Im_sub, Im)
  
  D = normxcorr2(Im_sub, Im);

  [~, idx] = max(abs(D(:)));
  
  [i, j] = ind2sub(size(D), idx);
  
  Im_fit = Im(i-size(Im_sub,1)+1:i, j-size(Im_sub,1)+1:j);
  
end