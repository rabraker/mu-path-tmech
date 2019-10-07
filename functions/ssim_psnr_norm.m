function [psn, ssm] = ssim_psnr_norm(Im_master, Im_compare, varargin)
%   Im_master = detrend_lines_local(Im_master);
%   Im_compare = detrend_lines_local(Im_compare);
  
  if length(varargin) > 1
    DRng = varargin{1};
  else
    Im_master = Im_master - mean(Im_master(:));
    Im_compare = Im_compare - mean(Im_compare(:));
    
    mx1 = max((Im_master(:)));
    mx2 = max((Im_compare(:)));
    
    mn1 = min((Im_master(:)));
    mn2 = min((Im_compare(:)));
    
    mx = max(mx1, mx2);
    mn = min(mn1, mn2);
    DRng = mx - mn;
  end
  psn = psnr(Im_master, Im_compare, DRng);
  ssm = ssim(Im_master, Im_compare, 'DynamicRange', DRng);
  
end

function X=detrend_lines_local(X)
  
  for k=1:size(X,1)
    X(k,:) = detrend(X(k,:));
  end
  
end