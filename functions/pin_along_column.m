% pixmat = pin_along_column(pixmat, x1, x2)
function pixmat = pin_along_column(pixmat, x1, x2)

  xs = 1:size(pixmat,2);
  for k=1:size(pixmat,1)
      row = detrend(pixmat(k, :));
      y1 = row(x1);
      y2 = row(x2);
%     y1 = pixmat(k, x1);
%     y2 = pixmat(k, x2);
    
    m = (y2-y1)/(x2-x1);
    b = y1 - m*x1;
    line = m*xs + b;
%     pixmat(k,:) = pixmat(k,:) - line;
    pixmat(k, :) = row - line;
  end
  

end
