% ------------ Row 1, CS20NG -------------------------

addpath functions
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

skip = 8;
h = skip:skip:npix;
im_64_sub = im_og(h, :);
im_64 = interp_y(im_64_sub, npix);
ssm_64 = ssim(im_64, im_og, 'DynamicRange', mx);
psn_64 = psnr(im_64, im_og, mx);


U_fun = @(x) idct(x);
Ut_fun = @(x) dct(x);

mu = 1e-7;
delta = 0.01;
alpv = 0.01;
alph = 0.75;
opts = NESTA_opts('U', U_fun, 'Ut', Ut_fun, 'verbose', 0, 'TolVar', 1e-5,...
    'alpha_v', alpv, 'alpha_h', alph);

%%
qs = [(1:2:10), (15:5:50), (50:15:512)];
fracs = (5:2:50)/100;
% fracs = 50/100;
dat128 = zeros(length(qs), 7);
dat64 = zeros(length(qs), 7);
map_ssm = zeros(length(qs), length(fracs));
map_psn = zeros(length(qs), length(fracs));

Y = map_ssm;
X = map_ssm;
gam = npix/50;
fprintf('128: ssm: %.3f, psn: %.3f\n', ssm_128, psn_128);
fprintf('64: ssm: %.3f, psn: %.3f\n', ssm_64, psn_64);
fprintf('===============================================\n');
figure(11); clf
subplot(1,3,1),
imagesc(im_128, [1, 255]), colormap('gray')
title('128 lines')
subplot(1,3,2)
imagesc(im_128, [1, 255]), colormap('gray')
title('64 lines')
drawnow()


% qs = 150
% fracs = .4
for yidx = 1:length(qs)
    store128 = true;
    store64 = true;
    opts.alpha_v = 1;
    opts.alpha_h = 1; %2*(gam*qs(yidx))/npix;
    for xidx=1:length(fracs)
        mu_len = qs(yidx);
        samp_frac = fracs(xidx);
        if 1
        [mask, pix_idx, npaths] = CsTools.mu_path_mask(mu_len, npix, npix, samp_frac);
        E_fun = @(x) CsTools.E_fun1(x, pix_idx);
        Et_fun = @(x) CsTools.Et_fun1(x, pix_idx, npix, npix);
        
        b = CsTools.pixmat2vec(im_og);
        b = b(pix_idx);
        x_bptv = NESTA_mine(E_fun, Et_fun, b, mu, delta, opts);
        
        X_bptv = CsTools.pixvec2mat(x_bptv, npix);
%         figure(11);
        subplot(1,3,3)
        imagesc(X_bptv, [1, 255])
        colormap('gray')
        drawnow()
        
        ssm_bp_k = ssim(X_bptv, im_og, 'DynamicRange', mx);
        psn_bp_k = psnr(X_bptv, im_og, mx);
        fprintf('q=%d, frac=%.3f %%, ssm_bp: %.3f, psn_bp: %.3f\n', mu_len,...
            samp_frac*100, ssm_bp_k, psn_bp_k)
        map_ssm(yidx, xidx) = ssm_bp_k;
        map_psn(yidx, xidx) = ssm_bp_k;
        end
        Y(yidx, xidx) = qs(yidx);
        X(yidx, xidx) = fracs(xidx);
        if (ssm_bp_k > ssm_128 && store128)
            dat128(yidx, :) = [ssm_bp_k, psn_bp_k, qs(yidx), fracs(xidx), npaths, xidx, yidx];
           store128 = false;
        end
        if (ssm_bp_k > ssm_64 && store64)
            dat64(yidx, :) = [ssm_bp_k, psn_bp_k, qs(yidx), fracs(xidx), npaths, xidx, yidx];
           store64 = false;
        end        
    end
end
%%
idx = find(qs > 365, 1, 'first')-1;
%%
dat128 = dat128(1:idx, :);
dat64 = dat64(1:idx, :);
map_ssm = map_ssm(1:idx, :);
map_psn = map_psn(1:idx, :);
%%
clc
% for j=1:length(qs)
%    for k=1:length(fracs)
%        dat64(
%    end
% end
y = qs(1:idx); %X(1, :);
x = fracs; %Y(1, :);
figure(10);clf
ax = gca();
hold on

surfl(x, y, map_ssm, 'light'), colormap('parula')

ylabel('q')
xlabel('frac')
zlabel('SSIM')

q128_opt = dat128(:, 3);
idx128_keep = 1:length(q128_opt); %find(q128_opt ~= 0);
q128_opt = q128_opt(idx128_keep);
f128_opt = dat128(idx128_keep, 4);

q64_opt = dat64(:, 3);
idx64_keep = 1:length(q64_opt); %find(q64_opt ~= 0);
q64_opt = q64_opt(idx64_keep);
idx64_keep = find(q64_opt ~= 0);
f64_opt = dat64(idx64_keep, 4);


h1 = plot3(f128_opt, q128_opt,  dat128(idx128_keep,1), 'r', 'LineWidth', 2);
h2 = plot3(dat64(idx64_keep,4), dat64(idx64_keep,3),  dat64(idx64_keep,1), 'b', 'LineWidth', 2);

ax.CameraPosition = [-2.4088  359.7179    1.4406];
ax.CameraViewAngle = 9.2102;
ax.CameraTarget = [0.2500   25.0000    0.7500];

h1.DisplayName = '128 line contour';
h2.DisplayName = '64 line contour';

legend([h1, h2])
title('cameraman')

% ylim([1, 100])
%%
figure(11); clf
hold on

% subplot(3,1,1)
plot(dat128(:, 3), dat128(:,4))
xlabel('q')
ylabel('frac')

hold on
plot(dat64(:, 3), dat64(:,4), 'r')




%%
tmu_s = load('tmu_s.mat')

tmu_ = (mean(tmu_s.tmu_12) + mean(tmu_s.tmu_15) + mean(tmu_s.tmu_25))/3

% tmu_ = tmu_ + 150*AFM.Ts;
if 0
q_opt = dat128(:, 3);
f_opt = dat128(:, 4);
npath_s = dat128(:, 5);
else
q_opt = dat64(:, 3);
f_opt = dat64(:, 4);
npath_s = dat64(:, 5);
end

idx_drop = find(q_opt == 0);
q_opt(idx_drop) = [];
f_opt(idx_drop) = [];
npath_s(idx_drop) = [];


npix = 512;
rate_1hz = 1; %Hz
rate_5hz = 5; %Hz
rate_10hz = 10; %Hz
time_128_1hz = q_opt*0;
time_128_5hz = q_opt*0;
time_128_10hz = q_opt*0;
for k=1:size(f_opt, 1)
    samp_frac = f_opt(k);
    q = q_opt(k);
    %[pmask, ~, npaths] = CsTools.mu_path_mask(q, npix, npix, samp_frac);
    npaths = npath_s(k);
    time_128_1hz(k) = (npaths * tmu_) + npaths * (q/npix)*(1/(2*rate_1hz));
    time_128_5hz(k) = (npaths * tmu_) + npaths * (q/npix)*(1/(2*rate_5hz));
    time_128_10hz(k) = (npaths * tmu_) + npaths * (q/npix)*(1/(2*rate_10hz));
end

figure(12);clf
hold on
t_512_1 = 512;
t_512_5 = 512/5;
t_512_10 = 512/10;
h1hz = plot(q_opt, time_128_1hz/t_512_1);
h5hz = plot(q_opt, time_128_5hz/t_512_5);
h10hz = plot(q_opt, time_128_10hz/t_512_10);

if 0
    h128 = plot([q_opt(1), q_opt(end)], ([128, 128])/512);
else
    h128 = plot([q_opt(1), q_opt(end)], ([64, 64])/512);
end
xlabel('q')
ylabel('fractional improvement')
ylim([0, 2])
grid on

h1hz.DisplayName = 'CS at 1Hz';
h5hz.DisplayName = 'CS at 5Hz';
h10hz.DisplayName = 'CS at 10Hz';

legend([h1hz, h5hz, h10hz])
% xlim([1, 150])

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
