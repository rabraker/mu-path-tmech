% This file compares 8 slow raster scans all taken back to back at 0.5 Hz overa
% 5 micron area. The goal is to establish a baseline of what we should expect
% for the performance metrics.

clc
clear
close all
addpath functions
addpath ~/matlab/dependencies/l1c/
addpath ~/matlab/dependencies/l1c/lib

% initialize paths.

figbase = 20;

size_dir = '5microns';
fname = fullfile(PATHS.sysid, 'x-axis_sines_info_intsamps_zaxisFourierCoef_10-29-2018-01.mat');
models = load(fname);
G = -models.modelFit.G_zdir;
p = pole(G);
z = zero(G);
gg = zpk(z(end-1:end), p(1:2), 1, G.Ts);
hole_depth = (20);

chan_map = ChannelMap([1:5]);
exp_date = '3-21-2019'
% ----------------------- Load and Process CS-data -----------------------------
dat_root = PATHS.raster_image_data(size_dir, exp_date);

% ----------------------- Load and Process raster-data -------------------------
raster_files = {...
'raster_scan_512pix_5mic_01Hz_out_3-21-2019-01.csv',...
'raster_scan_512pix_5mic_01Hz_out_3-21-2019-02.csv',...
'raster_scan_512pix_5mic_01Hz_out_3-21-2019-06.csv',...
};
% 'raster_scan_512pix_5mic_01Hz_out_3-21-2019-03.csv',...
% 'raster_scan_512pix_5mic_01Hz_out_3-21-2019-04.csv',...
% 'raster_scan_512pix_5mic_01Hz_out_3-21-2019-05.csv',...


pixmaps = cell(length(raster_files), 1);
for k=1:length(raster_files)
  raster_paths = get_raster_paths(dat_root, raster_files{k});
  rast_exps{k} = RasterExp(raster_paths, 'reload_raw', false);
end

%%
% npix = rast_exps{1}.npix_x;
% Ts = rast_exps{1}.Ts;
% width = rast_exps{1}.width;
% stit = sprintf('Scan %d', k)

%%
x1s = [44,   49,  52,  57, 55];
x2s = [459, 459, 467. 474, 474];
figbase = 10;

for k=1:length(rast_exps)
  rast_exps{k}.bin_raster_really_slow(@detrend);
  
  pixmats_raw{k} = rast_exps{k}.pix_mat(1:end, 1:end);
%   rast_exps{k}.pix_mat_pinned = pixmats_raw{k};
  pixmat_ = pin_along_column(rast_exps{k}.pix_mat, x1s(k), x2s(k));
  
  rast_exps{k}.pix_mat_pinned = pixmat_ - mean(pixmat_(:));
  rast_exps{k}.pin_idx_s = [x1s(k), x2s(k)];
  stit = sprintf('(raster) %.2f Hz', rast_exps{k}.meta_in.raster_freq);
  plot_raster_data(rast_exps{k}.pix_mat_pinned, figbase*k, stit)
  pixmats{k} = rast_exps{k}.pix_mat_pinned;
end

for k=1:length(rast_exps)
  rast_exps{k}.save()
end

%%
slice = 25:511-25;

imm_idx = 1;

im_master = pixmats{imm_idx};
im_master = im_master - mean(im_master(:));
% figure, imagesc(im_master), colormap('gray')

modes = [false, true];
F = mkfig(3000, 3.5, 3.5); clf
[ha, pos] = tight_subplot(2, 2, [.02, .01 ], [.027, 0.055], .005);
ha = reshape(ha', 2, [])';

F6 = mkfig(3002, 3.5, 3.5/2); clf
[ha6, pos] = tight_subplot(1, 2, [.02, .01 ], [.027, 0.055], .005);


clc
suffix = {'NoAlgn', 'Algn'}
names = {'one2two', 'one2six'};
ST = struct();
for j=1:2
  mode = modes(j);

  thresh = (20/7)*(1/1000)*20;
  
  titles = {'image 1 \& 2', 'image 1 \& 6'};
  for k=1:2 
    imk = pixmats{k+1};
    imk = imk - mean(imk(:));

    if mode == true
      imk_slice = imk(slice, slice);
      im1_ontok_fit = norm_align(imk_slice, im_master);
    else
      imk_slice = imk;
      im1_ontok_fit = im_master;
    end
    [psn_1k, ssm_1k] = ssim_psnr_norm(im1_ontok_fit, imk_slice, 2*thresh);
    name = [names{k}, suffix{j}]
    psn_name = [name, 'Psnr'];
    ssm_name = [name, 'Ssim'];
    ST.(psn_name) = psn_1k;
    ST.(ssm_name) = ssm_1k;
    
    dmg = rast_exps{k+1}.damage_metric();
    mx = max(imk_slice(:));
    mn = min(imk_slice(:));
    stit = titles{k}; 
    
%     fprintf("%s: 'PSNR=%.2f, SSIM=%.2f'\n", stit, psn_1k, ssm_1k);
    ime=imk_slice - im1_ontok_fit; %, 'parent', h(k), 'Scaling', 'joint')
    ime = ime-mean(ime(:));
    imagesc(ha(j, k), ime, [-thresh, thresh]); %, [-0.5*thresh, 0.5*thresh]);
    colormap('gray')
    
    set(ha(j, k), 'YTick', [], 'XTick', [])
    if j == 1
        title(ha(j, k), stit)
    end
    
    if k==2
        imagesc(ha6(j), ime, [-thresh, thresh]);
        colormap('gray')
        set(ha6(j), 'YTick', [], 'XTick', [])
    end
  end

%   set(ha(end), 'Visible', 'off');

end

data_writer({ST}, 'latex/baseline_errors.txt', '%.2f')
save_fig(F, 'latex/figures/baseline_errors_aligned_1Hz', false)
save_fig(F6, 'latex/figures/baseline_errors_1and6', false)


%%
function plot_raster_data(pixmat2, figbase, stit)

  thresh = (20/7)*(1/1000)*20;
  pixmat2 = pixmat2 - mean(pixmat2(:));
  F10 = figure(figbase+2); clf
  ax1 = subplot(3,1,[1,2]);
  ax2 =subplot(3,1,3);
    

  lo = min(min(pixmat2));
  hi = max(max(pixmat2));
  
  
  imshow_dataview(flipud(pixmat2 - mean(pixmat2(:))), [-thresh, thresh], ax1, ax2)
  try
    grid(ax1, 'on')
  catch
    keyboard
  end
  
  colormap(ax1, 'gray')
  grid(ax2, 'on')
  ax1.GridAlpha = 1;
  ax2.GridAlpha = 1;
  title(ax1, stit)
  title(ax2, stit)
  
%   figure(figbase+4)
%   mesh(pixmat2)
%   xlabel('x')
%   ylabel('y')
%   title(stit)
end
