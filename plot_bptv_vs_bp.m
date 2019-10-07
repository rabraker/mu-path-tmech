clear
addpath('functions/')
clc

% ------------ Row 1, CS20NG -------------------------
im_og = double(imread('cs20ng.png'));
% im_og = im_og(1:2:512, 1:2:512);
mu_len = 50;
samp_frac = 0.15;
n = size(im_og, 1);

[mask, pix_idx] = CsTools.mu_path_mask(mu_len, n, n, samp_frac);

figbase = 1;
t_recon_cs20ng = cp_bptv_bp(im_og, pix_idx, 'cs20ng', figbase);

U_fun = @(x) idct(x);
Ut_fun = @(x) dct(x);
    
E_fun = @(x) CsTools.E_fun1(x, pix_idx);
Et_fun = @(x) CsTools.Et_fun1(x, pix_idx, n, n);
%        
% clc
% im_og = im_og / max(im_og(:));
% b = CsTools.pixmat2vec(im_og);
% b = b(pix_idx);
% alpv = .01;
% alph = .91;
% opts = NESTA_opts('U', U_fun, 'Ut', Ut_fun, 'verbose', 0, 'TolVar', 1e-5);

% [ssm, X] = fun(im_og, b, E_fun, Et_fun, opts, alpv, alph);
% alps_0 = [alpv; alph];
% %%
% 
% 
% f = @(alps) fun(im_og, b, E_fun, Et_fun, opts, alps(1), alps(2));
% 
% fmopts = optimoptions('fmincon', 'MaxIterations', 10, 'Disp', 'iter');
% alps = fmincon(f, alps_0, [], [], [], [], [0;0], [], [], fmopts)

%

% ------------ Row 2, DNA -------------------------
im_og = double(rgb2gray(imread('dnn_dna_gt.png')));

n = 256;
im_og = im_og(1:n, 1:n);

mu_len = 25;
samp_frac = 0.25;
figbase = 100;

[mask, pix_idx] = CsTools.mu_path_mask(mu_len, n, n, samp_frac, figbase);

t_recon_dna = cp_bptv_bp(im_og, pix_idx, 'dna', figbase);

fprintf('CS20NG reconstruction time: %.2f\n', t_recon_cs20ng);
fprintf('DNA reconstruction time: %.2f\n', t_recon_dna);


function t_recon = cp_bptv_bp(im_og, pix_idx, prefix, figbase)
    
    im_og = im_og / max(im_og(:));
    n = size(im_og, 1);
    
    U_fun = @(x) idct(x);
    Ut_fun = @(x) dct(x);
    
    E_fun = @(x) CsTools.E_fun1(x, pix_idx);
    Et_fun = @(x) CsTools.Et_fun1(x, pix_idx, n, n);
    
    
    opts = NESTA_opts('U', U_fun, 'Ut', Ut_fun, 'alpha_v', 0.1, 'alpha_h', .75,...
        'verbose', 5, 'TolVar', 1e-5, 'mu', 1e-3, 'sigma', 0.001);
    
    b = CsTools.pixmat2vec(im_og);
    b = b(pix_idx);
    
    tic 
    x_bptv = NESTA_mine(E_fun, Et_fun, b, opts);
    t_recon = toc;
    X_bptv = CsTools.pixvec2mat(x_bptv, n);
    
    
    opts = NESTA_opts('U', U_fun, 'Ut', Ut_fun, 'alpha_v', 0, 'alpha_h', 0,...
        'verbose', 5, 'TolVar', 1e-5);
    
    
    x_bp = NESTA_mine(E_fun, Et_fun, b, opts);
    X_bp = CsTools.pixvec2mat(x_bp, n);
    
    
    [Fig_dna_og, ax_og] = imshow_local(im_og, 1 + figbase);
    [Fig_dna_bp, ax_bp] = imshow_local(X_bp, 2 + figbase);
    [Fig_dna_bptv, ax_bptv] = imshow_local(X_bptv, 3 + figbase);
    
    % bottom, left corner of box
    
    x0 = 58;
    y0 = 66;
    
    y1 = 24;
    x1 = 96;
    
    slice_x = x0:x1;
    slice_y = y1:y0;
    
    size(slice_x)
    size(slice_y)
    
    plot_box(ax_og, x0, y0, x1, y1);
    
    [Fig_og_zoom] = imshow_local(im_og(slice_y, slice_x), 4 + figbase);
    [Fig_bp_zoom] = imshow_local(X_bp(slice_y, slice_x), 5 + figbase);
    [Fig_bptv_zoom] = imshow_local(X_bptv(slice_y, slice_x), 6 + figbase);
    
    save_fig(Fig_dna_og, fullfile(PATHS.tmech_fig(), [prefix, '_og']), false)
    save_fig(Fig_dna_bp, fullfile(PATHS.tmech_fig(), [prefix, '_bp']), false)
    save_fig(Fig_dna_bptv, fullfile(PATHS.tmech_fig(), [prefix, '_bptv']), false)
    
    save_fig(Fig_og_zoom, fullfile(PATHS.tmech_fig(), [prefix, '_og_zoom']), false)
    save_fig(Fig_bp_zoom, fullfile(PATHS.tmech_fig(), [prefix, '_bp_zoom']), false)
    save_fig(Fig_bptv_zoom, fullfile(PATHS.tmech_fig(), [prefix, '_bptv_zoom']), false)

end

function plot_box(ax, x0, y0, x1, y1)
    
   hold(ax, 'on')
   lw = 1.2;
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


function [ssm, X_bptv] = fun(im_og, b, E_fun, Et_fun, opts, alpv, alph)
    n = size(im_og, 1);
    
%     U_fun = @(x) idct(x);
%     Ut_fun = @(x) dct(x);
%     
%     E_fun = @(x) CsTools.E_fun1(x, pix_idx);
%     Et_fun = @(x) CsTools.Et_fun1(x, pix_idx, n, n);
%         
%     opts = NESTA_opts('U', U_fun, 'Ut', Ut_fun, 'alpha_v', alpv, 'alpha_h', alph,...
%         'verbose', 5, 'TolVar', 1e-5);

    opts.alpha_v = alpv;
    opts.alpha_h = alph;
    opts.mu = 1e-5;
    opts.sigma = 0.01;
%     mu = 1e-5;
%     delta = 0.01;

    x_bptv = NESTA_mine(E_fun, Et_fun, b, opts);

    X_bptv = CsTools.pixvec2mat(x_bptv, n);
    figure(7)
    imagesc(X_bptv, [0, 1])
    colormap('gray')
    
    ssm = -psnr(im_og, X_bptv, 1);
    s= sprintf('alph = %f, alpv = %f, ssm = %f\n', alph, alpv, ssm);
    title(s)
    fprintf('%s', s);
end
