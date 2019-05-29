clear
addpath('functions/')
clc


im_og = double(rgb2gray(imread('dnn_dna_gt.png')));
im_og = im_og / max(im_og(:));
n = 256;
im_og = im_og(1:n, 1:n);

mu_len = 25;
samp_frac = 0.25;

[mask, pix_idx] = CsTools.mu_path_mask(mu_len, n, n, samp_frac);

U_fun = @(x) idct(x);
Ut_fun = @(x) dct(x);

E_fun = @(x) CsTools.E_fun1(x, pix_idx);
Et_fun = @(x) CsTools.Et_fun1(x, pix_idx, n, n);


opts = NESTA_opts('U', U_fun, 'Ut', Ut_fun, 'alpha_v', 0.25, 'alpha_h', 1,...
    'verbose', 5, 'TolVar', 1e-5);

b = CsTools.pixmat2vec(im_og);
b = b(pix_idx);

mu = 1e-3;
delta = .001

x_bptv = NESTA_mine(E_fun, Et_fun, b, mu, delta, opts);
X_bptv = CsTools.pixvec2mat(x_bptv, n);


opts = NESTA_opts('U', U_fun, 'Ut', Ut_fun, 'alpha_v', 0, 'alpha_h', 0,...
    'verbose', 5, 'TolVar', 1e-5);


x_bp = NESTA_mine(E_fun, Et_fun, b, mu, delta, opts);
X_bp = CsTools.pixvec2mat(x_bp, n);


[Fig_dna_og, ax_dna_og] = imshow_local(im_og, 1);
[Fig_dna_bp, ax_dna_bp] = imshow_local(X_bp, 2);
[Fig_dna_bptv, ax_dna_bptv] = imshow_local(X_bptv, 3);

% bottom, left corner of box
x0 = 55;
y0 = 66;

y1 = 24;
x1 = 96;

slice_x = x0:x1;
slice_y = y1:y0;

size(slice_x)
size(slice_y)

plot_box(ax_dna_og, x0, y0, x1, y1);

[Fig_dna_og_zoom, ax_dna_og_zoom] = imshow_local(im_og(slice_y, slice_x), 4);
[Fig_dna_bp_zoom, ax_dna_bp_zoom] = imshow_local(X_bp(slice_y, slice_x), 5);
[Fig_dna_bptv_zoom, ax_dna_bptv_zoom] = imshow_local(X_bptv(slice_y, slice_x), 6);

save_fig(Fig_dna_og, fullfile(PATHS.tmech_fig(), 'dna_og'), false)
save_fig(Fig_dna_bp, fullfile(PATHS.tmech_fig(), 'dna_bp'), false)
save_fig(Fig_dna_bptv, fullfile(PATHS.tmech_fig(), 'dna_bptv'), false)

save_fig(Fig_dna_og_zoom, fullfile(PATHS.tmech_fig(), 'dna_og_zoom'), false)
save_fig(Fig_dna_bp_zoom, fullfile(PATHS.tmech_fig(), 'dna_bp_zoom'), false)
save_fig(Fig_dna_bptv_zoom, fullfile(PATHS.tmech_fig(), 'dna_bptv_zoom'), false)


function plot_box(ax, x0, y0, x1, y1)
    
   hold(ax, 'on')
   lw = 1.2
   %left edge
   plot(ax, [x0, x0], [y0, y1], 'r', 'LineWidth', lw);
   %top
   plot(ax, [x0, x1], [y1, y1], 'r', 'LineWidth', lw);
   %right edge
   plot(ax, [x1, x1], [y1, y0], 'r', 'LineWidth', lw);
   %bottom
   plot(ax, [x0, x1], [y0, y0], 'r', 'LineWidth', lw);
    
end

function [fig, ax] = imshow_local(im, fignum)
    fig = mkfig(fignum, 2, 2);
    ax = tight_subplot(1, 1, 0.01, 0.01, 0.01, true);
    
    imagesc(ax, im, [0, 1])
    colormap(ax, 'gray')
end


