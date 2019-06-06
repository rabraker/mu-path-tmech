function [idx_left, idx_right] = find_pin_idxs(IM)

IM = detrend_plane(IM);

im_var = var(IM);

xpix = 512;
edge = 15;
pin_width = 100;


[~, idx_left] = min(im_var(edge:edge+pin_width)) ;
idx_left = idx_left + edge + 1;

[~, idx_right] = min(im_var(end-edge-pin_width:end-edge)) ;
idx_right = (xpix - edge-pin_width) + idx_right+1;


% pixmat_ = pin_along_column(IM, idx_left, idx_right);
% 
% figure(1)
% imagesc(pixmat_-mean(pixmat_(:)), [-thresh, thresh])
% colormap('gray')

end