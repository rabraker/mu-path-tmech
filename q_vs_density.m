% ------------ Row 1, CS20NG -------------------------


im_og = double(imread('cs20ng.png'));
% im_og = double(imread('cameraman.tif'));
im_og = im_og - min(im_og(:));
im_og = (im_og / max(im_og(:)))*256;
% im_og = im_og(1:2:512, 1:2:512);

samp_frac = 0.15;
npix = size(im_og, 1);

mx = max(im_og(:)) - min(im_og(:))
skip = 4;
h = skip:skip:npix;
im_128_sub = im_og(h, :);

im_128 = interp_y(im_128_sub, npix);

ssm_128 = ssim(im_128, im_og, 'DynamicRange', mx);
psn_128 = psnr(im_128, im_og, mx);

U_fun = @(x) idct(x);
Ut_fun = @(x) dct(x);


%%
fprintf('128: ssm: %.3f, psn: %.3f\n', ssm_128, psn_128);
fprintf('===============================================\n');
figure(11); clf
subplot(1,2,1),
imagesc(im_128), colormap('gray')
drawnow()

mu = 1e-7;
delta = 0.01;
% alpv = 0.01;
alph = 0.75;
opts = NESTA_opts('U', U_fun, 'Ut', Ut_fun, 'verbose', 0, 'TolVar', 1e-5,...
    'alpha_v', alpv, 'alpha_h', alph);


qs = (5:3:50);
fracs = (5:2:50)/100;
% fracs = 50/100;
dat = zeros(length(qs), 6);
map = zeros(length(qs), length(fracs));
gam = npix/50;
for j = 1:length(qs)
    store = true;
    opts.alpha_v = 0.5;
    opts.alpha_h = 2*(gam*qs(j))/npix;
    for k=1:length(fracs)
        mu_len = qs(j);
        samp_frac = fracs(k);
        [mask, pix_idx] = CsTools.mu_path_mask(mu_len, npix, npix, samp_frac);
        E_fun = @(x) CsTools.E_fun1(x, pix_idx);
        Et_fun = @(x) CsTools.Et_fun1(x, pix_idx, npix, npix);
        
        b = CsTools.pixmat2vec(im_og);
        b = b(pix_idx);
        x_bptv = NESTA_mine(E_fun, Et_fun, b, mu, delta, opts);
        
        X_bptv = CsTools.pixvec2mat(x_bptv, npix);
        subplot(1,2,2), imagesc(X_bptv), colormap('gray')
        drawnow()
        ssm_bp_k = ssim(X_bptv, im_og, 'DynamicRange', mx);
        psn_bp_k = psnr(X_bptv, im_og, mx);
        fprintf('q=%d, frac=%.3f, ssm_bp: %.3f, psn_bp: %.3f\n', mu_len,...
            samp_frac, ssm_bp_k, psn_bp_k)
        map(j, k) = ssm_bp_k;
        if (ssm_bp_k > ssm_128 && store)
            dat(j, :) = [ssm_bp_k, psn_bp_k, qs(j), fracs(k), k, j];
%             break
           store = false;
        end
    end
end
%%
figure(10);clf
hold on
surf(map)
xlabel('frac')
ylabel('q')
%%
figure(11); clf
hold on

subplot(3,1,1)
plot(dat(:, 3), dat(:,4))

function img_interp = interp_y(im_og, ypix)

    ypix_old = size(im_og, 1);
    xpix = size(im_og, 2);
    img_interp = zeros(ypix, xpix);
    
    % Now interpolate the missing rows
    vq = 1:ypix;
    skip = floor(ypix/ypix_old);
    h = skip:skip:ypix;
    
    for k=1:xpix
        img_interp(:, k) = interp1(h, im_og(:, k), vq, 'linear', 'extrap');
    end


end
