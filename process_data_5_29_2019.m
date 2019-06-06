clc
clear
%%
close all
addpath functions

% initialize paths.
% init_paths();
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
exp_date = '5-30-2019'
% ----------------------- Load and Process CS-data -----------------------------
dat_root = PATHS.raster_image_data(size_dir, exp_date);

% ----------------------- Load and Process raster-data -------------------------
raster_files = {...
'raster_scan_512pix_5mic_1.00e+00Hz_fourier_out_5-30-2019-01.csv',...
'raster_scan_512pix_5mic_4.00e+00Hz_fourier_out_5-30-2019-01.csv',...
'raster_scan_512pix_5mic_05Hz_fourier_out_5-30-2019-01.csv',...
'raster_scan_512pix_5mic_8.00e+00Hz_fourier_out_5-30-2019-01.csv',...
'raster_scan_512pix_5mic_10Hz_fourier_out_5-30-2019-01.csv',...
'raster_scan_512pix_5mic_1.50e+01Hz_fourier_out_5-30-2019-01.csv',...
'raster_scan_128pix_5mic_05Hz_fourier_out_5-30-2019-01.csv',...
};

cs_files = {...
'cs-traj-512pix-10perc-500nm-5mic-01Hz-150prescan-notconnect_out_5-30-2019-01.csv',...
'cs-traj-512pix-10perc-500nm-5mic-02Hz-150prescan-notconnect_out_5-30-2019-01.csv',...
'cs-traj-512pix-10perc-500nm-5mic-04Hz-150prescan-notconnect_out_5-30-2019-01.csv',...
'cs-traj-512pix-10perc-500nm-5mic-05Hz-150prescan-notconnect_out_5-30-2019-01.csv',...
'cs-traj-512pix-10perc-500nm-5mic-08Hz-150prescan-notconnect_out_5-30-2019-01.csv',...
'cs-traj-512pix-15perc-500nm-5mic-01Hz-150prescan-notconnect_out_5-30-2019-01.csv',...
'cs-traj-512pix-15perc-500nm-5mic-02Hz-150prescan-notconnect_out_5-30-2019-01.csv',...
'cs-traj-512pix-15perc-500nm-5mic-04Hz-150prescan-notconnect_out_5-30-2019-01.csv',...
'cs-traj-512pix-15perc-500nm-5mic-05Hz-150prescan-notconnect_out_5-30-2019-01.csv',...
'cs-traj-512pix-15perc-500nm-5mic-08Hz-150prescan-notconnect_out_5-30-2019-01.csv',...
};


%%
rast_exps = {};
for k=1:length(raster_files)
  raster_paths = get_raster_paths(dat_root, raster_files{k});
  rast_exps{k} = RasterExp(raster_paths, 'load_full', true, 'reload_raw', false);
  if rast_exps{k}.time_total == 0
      rast_exps{k}.time_total = rast_exps{k}.samps_per_period*rast_exps{k}.npix*AFM.Ts;
  end
end


% dat = loadjson('/media/labserver/afm-cs/imaging/raster/5microns/parents/raster_scan_512pix_5mic_10Hz.json')

%%
use_ze = false;
% last image is 128 pixels
x1s = [25,   27,  27,  27,  27, 27, 6*4];
x2s = [493, 493, 493, 495, 494, 494, 122*4]; 
figbase = 10;
for k=1:length(rast_exps)
    if k==7
        rast_exps{k}.bin_raster_really_slow(@detrend, use_ze, 512);
        rast_exps{k}.pix_mat = rast_exps{k}.interp_y(512);
        stit = sprintf('(raster %.1f sec) interpolated from %d pix %.2f Hz',...
            length(rast_exps{k}.ze)*AFM.Ts, rast_exps{k}.npix_y,...
            rast_exps{k}.meta_in.raster_freq);
    else
        rast_exps{k}.bin_raster_really_slow(@detrend, use_ze);
        stit = sprintf('(raster %.1f sec) %.2f Hz', length(rast_exps{k}.ze)*AFM.Ts,...
            rast_exps{k}.meta_in.raster_freq);
    end
  
  pixmats_raw{k} = rast_exps{k}.pix_mat(1:end, 1:end);
%   rast_exps{k}.pix_mat_pinned = pixmats_raw{k};
  pixmat_ = pin_along_column(rast_exps{k}.pix_mat, x1s(k), x2s(k));
%   pixmat_ = pixmats_raw{k};
  rast_exps{k}.pix_mat_pinned = pixmat_ - mean(pixmat_(:));
  rast_exps{k}.pin_idx_s = [x1s(k), x2s(k)];
  
  plot_raster_data(rast_exps{k}.pix_mat_pinned, figbase*k, stit)

end

% save the first image so we can simulate reconstruction in plot_bptv_vs_bp.m
%%
img = rast_exps{1}.pix_mat_pinned;
img = img - min(img(:));
img = (img/max(img(:)) ) * 255/1.7;
img = uint8(img);

imwrite(img, 'cs20ng.png')

%%
for k=1:length(rast_exps)
  rast_exps{k}.save()
end
%%
use_ze = false;
data_root = PATHS.cs_image_data(size_dir, exp_date);
cs_exps = {};
%%
for k=1:length(cs_files)
  cs_paths = get_cs_paths(data_root, cs_files{k});

  gg = @(u, idx_state_s) log_creep_detrend(u, idx_state_s);
  
  cs_exps{k} = CsExp(cs_paths, 'feature_height', hole_depth, 'gg', gg,...
      'load_full', true, 'reload_raw', true);

  cs_exps{k}.print_state_times();
  cs_exps{k}.sub_sample_frac()
end

%%
bp = true;
recalc = true;
use_dct2 = false;
thresh = (20/7)*(1/1000)*20;


register_uzk = true;
for k=5:5 %length(cs_exps)
  cs_exps{k}.process_cs_data(false, [], use_ze, register_uzk);
  fprintf('finished processing raw CS data...\n');
  fprintf('nperc=%.3f\n', sum(cs_exps{k}.pix_mask(:))/cs_exps{k}.npix^2);
  ht = cs_exps{k}.feature_height;
  if bp
      % opts = l1qc_dct_opts('l1_tol', 0.00001, 'epsilon', 0.1);
      % cs_exps{k}.solve_bp(recalc, use_dct2, opts);
      U_fun = @(x) idct(x);
      Ut_fun = @(x) dct(x);

      opts = NESTA_opts('U', U_fun, 'Ut', Ut_fun, 'alpha_v', 0.1, 'alpha_h', .5,...
          'verbose', 5, 'TolVar', 1e-5);
      
      cs_exps{k}.solve_nesta(recalc, use_dct2, opts);


    fprintf('Finished solving bp problem #%d\n', k);
    stit = sprintf('(CS) %.2f Hz equiv, \\%% %.2f sampling\nTotal time: %.2f',...
      cs_exps{k}.meta_in.tip_velocity/(cs_exps{k}.meta_in.width * 2),...
      cs_exps{k}.meta_in.actual_sub_samble_perc,...
      length(cs_exps{k}.x)*cs_exps{k}.Ts);
    
    f10=mkfig(1001 + 2*k, 6, 7.5); clf
    ax4 = subplot(3,1,[1,2]);
    ax4_2 =subplot(3,1,3);
    
    ImshowDataView.setup(f10);
    cb_exp =  @(event_obj)cs_exps{k}.dataview_callback(event_obj, ax4, ax4_2);
    bar_bp = mean(cs_exps{k}.Img_bp(:));
    
    cs_exps{k}.Img_bp = cs_exps{k}.Img_bp - bar_bp;
    cs_exps{k}.Img_raw = cs_exps{k}.Img_raw - bar_bp;
    
    im_tmp = cs_exps{k}.Img_bp;
    im_tmp = detrend_plane(im_tmp);
    im_tmp = im_tmp - mean(im_tmp(:));
%     im_tmp = SplitBregmanROF(im_tmp, 100, 0.001);
    ImshowDataView.imshow(im_tmp,...
      [-thresh, thresh], ax4, ax4_2, cb_exp)
    title(ax4, stit)
  end
% keyboard
%   cs_exps{k}.save();
end
%%
[~, axs1] = make_cs_traj_figs(figbase, 4);
cs_exps{1}.plot_all_cycles(axs1{1:4});

[~, axs2] = make_cs_traj_figs(figbase+100, 4);
cs_exps{2}.plot_all_cycles(axs2{1:4});

%%
for k=1:length(cs_exps)
  cs_exps{k}.save()
end
%%
Fig1 = mkfig(3000, 7, 5.5); clf;
ha1 = tight_subplot(3, 4, [0.032, 0.015], [.02, .04], [.02, .02], true);

Fig_err = mkfig(3001, 9, 9); clf
ha_err = tight_subplot(4, 4, [0.025, 0.015], [.01, .02], [.05, .02], true);
%
Fig_rows = mkfig(3002, 7, 4.5); clf
ha_row = tight_subplot(2, 1, [0.1, 0.015], [.1, .05], [.085, .02], false);
xlabel(ha_row(2), 'x-direction pixel', 'FontSize', 14)
ylabel(ha_row(1), 'height [nm]', 'FontSize', 14)
title(ha_row(1), 'raster', 'FontSize', 14)
title(ha_row(2), 'CS', 'FontSize', 14)
ylabel(ha_row(2), 'height [nm]', 'FontSize', 14)

grid(ha_row(1), 'on')
grid(ha_row(2), 'on')

mu = Inf;
Img_filts = {};
mxs = [];
thresh = (20/7)*(1/1000)*20;
DRng = 2*thresh;

for k=1:length(cs_exps)
  Img_filts{k} = cs_exps{k}.Img_bp - min(cs_exps{k}.Img_bp(:));
  if ~isinf(mu)
    Img_filts{k} = SplitBregmanROFAn(Img_filts{k}, mu, 0.001);
  end
%   figure(101+k)
  mx = max(Img_filts{k}(:)) - min(Img_filts{k}(:));
  mxs = [mxs; mx];
%   mesh(Img_filts{k}), colormap('gray')
end

slice = 30:512-30;
master_idx = 1;
im_master = rast_exps{master_idx}.pix_mat_pinned - mean(rast_exps{master_idx}.pix_mat_pinned(:));
if ~isinf(mu)
  %im_master = SplitBregmanROF(im_master, mu, 0.001);
end
fprintf('---------------------------------------------------\n\n');
scan_metrics = {};
row_idx = 166;


excludes = [];
do_rows = [1, 3, 5];
err_idx = 1;
plt_idx = 1;
row_plt_idx = 1;
for k=1:length(rast_exps)
  imk = rast_exps{k}.pix_mat_pinned - mean(rast_exps{k}.pix_mat_pinned(:));
  if ~isinf(mu)
    %imk = SplitBregmanROF(imk, mu, 0.001);
  end

  im1_ontok_fit = im_master(slice, slice);
  imk_slice = norm_align(im1_ontok_fit, imk);

%   imk_slice = imk(slice, slice);
%   im1_ontok_fit = norm_align(imk_slice, im_master);
  
  
%   [im1_ontok_fit, imk_slice] = align_by_metric(im_master, imk, [], 'psnr');
  [psn_1k, ssm_1k] = ssim_psnr_norm(im1_ontok_fit, imk_slice, DRng);
  damage = rast_exps{k}.damage_metric();
  quality = rast_exps{k}.quality_metric();
  
  csm = ScanMetrics('ssim', ssm_1k, 'psnr', psn_1k,...
    'quality', quality, 'damage', damage,...
    'rate', rast_exps{k}.meta_in.raster_freq,...
    'time', rast_exps{k}.time_total,...
    'coverage', 100,...
    'type', 'raster');
    
  scan_metrics{end+1} = csm;
  
  if k ~= 1
      im_err = im1_ontok_fit - imk_slice;
      imagesc(ha_err(err_idx), im_err, [-thresh, thresh]);
      colormap(ha_err(err_idx), 'gray');
      stit_err = sprintf('raster (%.1f Hz)', csm.rate);
      title(ha_err(err_idx), stit_err);
      err_idx = err_idx + 1;
  end
  
  if ismember(k, [1, 3, 5, 6])
      imagesc(ha1(plt_idx), imk, [-thresh, thresh]);
      colormap(ha1(plt_idx), 'gray');
      stit = sprintf('raster (%.1f Hz)', csm.rate);
      title(ha1(plt_idx), stit);
      
      
      if ~isempty(intersect(k, do_rows))
          hold(ha1(plt_idx), 'on')
          hold(ha_row(1), 'on')
          
          plot(ha1(plt_idx), [1, 512], [row_idx, row_idx], 'r');
          hl_row = plot(ha_row(1), imk_slice(row_idx, :)*AFM.volts2nm_z());
          hl_row.DisplayName = sprintf('%.1f Hz', csm.rate);
          hold(ha1(plt_idx), 'on')
          hold(ha_row(1), 'on')
          
      end
      plt_idx = plt_idx + 1;
  end

  fprintf('(raster %.1f Hz) psnr: %.4f, ssm: %.4f, damage: %.4f, time: %.4f\n',...
    rast_exps{k}.meta_in.raster_freq, psn_1k, ssm_1k, damage, csm.time);

end



remove_ticks(ha1)
remove_ticks(ha_err)
%

set(ha_row(1), 'XLim', [1, 512])
axes(ha_row(1))
leg1 = legend();
set(leg1, 'NumColumns', 3, 'FontSize', 11, 'Position', [0.5574 0.8968 0.4223 0.0500])


j = length(rast_exps);
fprintf('---------------------------------------------------\n');
cs_stats = {};
do_rows = [1, 4, 6, 10];
for k=1:length(cs_exps)

  imk = cs_exps{k}.Img_bp - mean(cs_exps{k}.Img_bp(:));
  if ~isinf(mu)
    imk = SplitBregmanROF(imk, mu, 0.001);
  end
%   imk_slice = imk(slice, slice);
%   im1_ontok_fit = norm_align(imk_slice, im_master);
  
  im1_ontok_fit = im_master(slice, slice);
  imk_slice = norm_align(im1_ontok_fit, imk);  
  
  [psn_1k, ssm_1k] = ssim_psnr_norm(im1_ontok_fit, imk_slice, DRng);

  damage = cs_exps{k}.damage_metric();
  quality = cs_exps{k}.quality_metric();

  csm = ScanMetrics('ssim', ssm_1k, 'psnr', psn_1k,...
    'quality', quality, 'damage', damage,...
    'rate', cs_exps{k}.meta_in.tip_velocity/(cs_exps{2}.meta_in.width * 2),...
    'time', length(cs_exps{k}.x)*cs_exps{k}.Ts,...
    'coverage', cs_exps{k}.meta_in.actual_sub_samble_perc,...
    'type', 'CS');
  
  fprintf('(CS %.1f Hz, %% %.1f, %.1f) psnr: %.4f, ssm: %.4f, damage: %.4f, time: %.4f\n',...
    csm.rate, csm.coverage, cs_exps{k}.sub_sample_frac*100, psn_1k, ssm_1k, damage, csm.time);

    scan_metrics{end+1} = csm;
  [t_cycle_avg, t_connect] = cs_exps{k}.estimate_mpt_connect_savings();
  cs_stats{k, 1} = csm;
  cs_stats{k, 2} = t_cycle_avg;
  cs_stats{k, 3} = t_connect;
  

  if csm.rate ~= 4
      im_err = im1_ontok_fit - imk_slice;
      imagesc(ha_err(err_idx), im_err, [-thresh, thresh]);
      colormap(ha_err(err_idx), 'gray');
      title(ha_err(err_idx), sprintf('CS (%.1f \\%%, %.1f Hz)', csm.coverage, csm.rate));
      
      err_idx = err_idx + 1;
  end
  
  if csm.rate ~= 4
      imagesc(ha1(plt_idx), imk, [-thresh, thresh]);
      colormap(ha1(plt_idx), 'gray');
      title(ha1(plt_idx), sprintf('CS (%.1f \\%%, %.1f Hz)', csm.coverage, csm.rate));

      
      if ~isempty(intersect(k, do_rows))
          hold(ha1(plt_idx), 'on')
          hold(ha_row(2), 'on')
          
          plot(ha1(plt_idx), [1, 512], [row_idx, row_idx], 'r');
          hl_row = plot(ha_row(2), imk_slice(row_idx, :)*AFM.volts2nm_z());
          hl_row.DisplayName = sprintf('%.1f Hz, %.1f \\%%', csm.rate, csm.coverage);
      end
      plt_idx = plt_idx + 1;
 end
  
  j=j+1;
end

remove_ticks(ha1)
remove_ticks(ha_err)

set(ha_row(2), 'XLim', [1, 512])
axes(ha_row(2))
leg2 = legend();
set(leg2, 'NumColumns', 4, 'FontSize', 11, 'Position', [0.1032 0.4253 0.8589 0.0500])

%%
csm_r = {};
csm_cs10 = {};
csm_cs15 = {};

tr = [];
Hzr = [];
dmg_r = [];
psnr_r = [];
ssim_r = [];

t_cs10 = [];
Hz_cs10 = [];
dmg_cs10 = [];
psnr_cs10 = [];
ssim_cs10 = [];

t_cs15 = [];
Hz_cs15 = [];
dmg_cs15 = [];
psnr_cs15 = [];
ssim_cs15 = [];
for k = 1:length(scan_metrics)
    if strcmp(scan_metrics{k}.type, 'raster')
        csm_r{end+1} = scan_metrics{k};
        tr(end+1) = scan_metrics{k}.time;
        Hzr(end+1) = scan_metrics{k}.rate;
        dmg_r(end+1) = scan_metrics{k}.damage;
        ssim_r(end+1) = scan_metrics{k}.ssim;
        psnr_r(end+1) = scan_metrics{k}.psnr;
    elseif scan_metrics{k}.coverage < 14
        t_cs10(end+1) = scan_metrics{k}.time;
        Hz_cs10(end+1) = scan_metrics{k}.rate;
        dmg_cs10(end+1) = scan_metrics{k}.damage;
        ssim_cs10(end+1) = scan_metrics{k}.ssim;
        psnr_cs10(end+1) = scan_metrics{k}.psnr;
        
        csm_cs10{end+k} = scan_metrics{k};
    else
        t_cs15(end+1) = scan_metrics{k}.time;
        Hz_cs15(end+1) = scan_metrics{k}.rate;
        dmg_cs15(end+1) = scan_metrics{k}.damage;
        ssim_cs15(end+1) = scan_metrics{k}.ssim;
        psnr_cs15(end+1) = scan_metrics{k}.psnr;
        
        csm_cs15{end+k} = scan_metrics{k};
    end
end

% save('damage_data.mat', 'Hzr', 'tr', 'dmg_r', 't_cs10', 'Hz_cs10', 'dmg_cs10',...
%     't_cs15', 'Hz_cs15', 'dmg_cs15')

idx_512 = 1:6;
idx_128 = 7;

F_rdi = mkfig(3, 3., 3); clf
ha = tight_subplot(1, 1, .1, [0.155, 0.05], [0.162, 0.04], false);
hold(ha(1), 'on')
grid(ha(1), 'on')
set(ha(1), 'XScale', 'log');

h3 = semilogx(ha(1), tr(idx_512), dmg_r(idx_512), 'x', 'MarkerSize', 8);
h3a = semilogx(ha(1), tr(idx_128), dmg_r(idx_128), '+', 'MarkerSize', 8);

xlabel(ha(1), 'total time [s]')
ylabel(ha(1), 'RDI')



h4 = semilogx(ha(1), t_cs10, dmg_cs10, 'or');
h5 = semilogx(ha(1), t_cs15, dmg_cs15, '*b');


h3.DisplayName = 'raster 512';
h3a.DisplayName = 'raster 128';

h4.DisplayName = 'CS: 10\%';
h5.DisplayName = 'CS: 15\%';
legend([h3, h3a, h4, h5], 'location', 'northeast')

offsets = [50, 10, 10, 5, 5];
for k=1:5
    st = sprintf('%.0f Hz', Hzr(k));
    t1 = text(tr(k)+offsets(k), dmg_r(k), st, 'HorizontalAlignment', 'left');
end

offsets = [10, 10, 10, 5, 5];
for k=1:5
    st = sprintf('%d Hz', Hz_cs15(k));
    t1 = text(t_cs15(k)+offsets(k), dmg_cs15(k), st, 'HorizontalAlignment', 'left');
end


%%

save_fig(F_rdi, fullfile(PATHS.tmech_fig(), 'cs_rast_damage'), false)
%%

% ---------------------------------------------------------------------- %
 
idx_512 = 2:6;
idx_128 = 7;

F_ssm = mkfig(4, 3, 3.); clf
ha_ssm = tight_subplot(1, 1, .1, [0.155, 0.05], [0.162, 0.02], false);
F_psn = mkfig(5, 3, 3.); clf
ha_psn = tight_subplot(1, 1, .1, [0.155, 0.05], [0.162, 0.02], false);

xlabel(ha_ssm, 'total time [s]')
ylabel(ha_ssm, 'SSIM')
hold(ha_ssm, 'on')
grid(ha_ssm, 'on')

set(ha_ssm, 'XScale', 'log')

ylabel(ha_psn, 'PSNR')
xlabel(ha_psn, 'total time [s]')

hold(ha_psn, 'on')
grid(ha_psn, 'on')
set(ha_psn, 'XScale', 'log')

h4 = semilogx(ha_ssm, tr(idx_512), ssim_r(idx_512), 'x', 'MarkerSize', 8);
h4a = semilogx(ha_ssm, tr(idx_128), ssim_r(idx_128), '+', 'MarkerSize', 8);

h5 = semilogx(ha_ssm, t_cs10, ssim_cs10, 'or');
h6 = semilogx(ha_ssm, t_cs15, ssim_cs15, '*b');


h7 = semilogx(ha_psn, tr(idx_512), psnr_r(idx_512), 'x', 'MarkerSize', 8);
h7a = semilogx(ha_psn, tr(idx_128), psnr_r(idx_128), '+', 'MarkerSize', 8);

hold(ha_ssm, 'on')
grid(ha_ssm, 'on')

h8 = semilogx(ha_psn, t_cs10, psnr_cs10, 'or');
h9 = semilogx(ha_psn, t_cs15, psnr_cs15, '*b');


h7.DisplayName = 'raster 512';
h7a.DisplayName = 'raster 512';

h8.DisplayName = 'CS: 10\%';
h9.DisplayName = 'CS: 15\%';
leg = legend([h7, h7a, h8, h9], 'location', 'Northwest');

%%
save_fig(F_ssm, fullfile(PATHS.tmech_fig(), 'cs_rast_time_vs_ssim'), false)
save_fig(F_psn, fullfile(PATHS.tmech_fig(), 'cs_rast_time_vs_psnr'), false)

%%


% set(ha1(end), 'Visible', 'off')
% set(ha_err(end), 'Visible', 'off')


save_fig(Fig1, fullfile(PATHS.tmech_fig(), 'cs_raster_images_4-26-2019'), false)
save_fig(Fig_err, fullfile(PATHS.tmech_fig(), 'cs_raster_images_err_4-26-2019'),false)
save_fig(Fig_rows, fullfile(PATHS.tmech_fig(), 'cs_raster_pixel_rows_4-26-2019'), false)

%%

% skip the 0.5Hz scan
sm = scan_metrics(1:end)
S1 = scan_metrics_table(sm)

S2 = state_times_table(cs_exps)



fid = fopen(fullfile(PATHS.tmech_table(), 'cs_raster_table_4-26-2019_muInf_dct2.tex'), 'w+');
fprintf(fid, '%s', S1);
fclose(fid);

fid = fopen(fullfile(PATHS.tmech_table(), 'cs_state_times_table_4-26-2019_muInf_dct2.tex'), 'w+');
fprintf(fid, '%s', S2);
fclose(fid);
%%
% Compare the 5-hz, 128 line raster to 4Hz, 15\% CS

Fig_subl = mkfig(1010, 3.5, 2); clf
ha = tight_subplot(1, 2, [0.02, 0.02], [0.01, 0.16], [0.01, 0.01], true);

imagesc(ha(1), rast_exps{end}.pix_mat_pinned, [-thresh, thresh]);
colormap(ha(1), 'gray')
imagesc(ha(2), cs_exps{end-2}.Img_bp, [-thresh, thresh]);
colormap(ha(2), 'gray')
remove_ticks(ha)

stit1 = sprintf( '(raster 25.6 sec.) \n5~Hz, 128 lines interpolated')
title(ha(1),stit1)
stit2 = sprintf('(CS %.1f sec.) \n%.1f \\%%, %.0f~Hz',...
    cs_exps{8}.time_total(),...
    cs_exps{8}.sub_sample_frac()*100,...
    cs_exps{8}.meta_in.tip_velocity*0.5/cs_exps{8}.meta_in.width)
title(ha(2), stit2)

save_fig(Fig_subl, fullfile(PATHS.tmech_fig, 'subline_vs_cs'), false)


function S = scan_metrics_table(csm_s)
  S = 'type &  rate (Hz) & PSNR & SSIM & time [s] & damage\\';
  S = sprintf('%s\n\\toprule\n', S);
  
  for k=1:length(csm_s)
    csm = csm_s{k};
    if strcmp(csm.type, 'CS')
      description = sprintf('CS (%.2f~\\%%)', csm.coverage);
    else
      description = 'raster';
    end
    S = sprintf('%s%s & %.2f & %.2f & %.2f & %.1f & %.2f\\\\\n', S,...
      description, csm.rate, csm.psnr, csm.ssim, csm.time, csm.damage);
  end
end

function S = state_times_table(cs_exps)
  S = 'description & move & engage & pre-scan & scan & tip-up & total\\';
  S = sprintf('%s\n\\toprule\n', S);
  for k=1:length(cs_exps)
    cst = cs_exps{k}.get_state_times();
    description = sprintf('%.1f Hz, %.2f~\\%%',...
      cs_exps{k}.meta_in.tip_velocity/(cs_exps{2}.meta_in.width * 2),...
      cs_exps{k}.meta_in.actual_sub_samble_perc);
  
    S = sprintf('%s%s &%.2f & %.2f & %.2f & %.2f & %.2f & %.2f\\\\\n', S,...
      description, cst.move, cst.tdown, cst.tsettle, cst.scan,...
      cst.tup, cst.total);
  end
end




%%
function plot_raster_data(pixmat2, figbase, stit, plot_mesh)
  if nargin < 4
    plot_mesh = false;
  end
  thresh = (20/7)*(1/1000)*20;
  pixmat2 = pixmat2 - mean(pixmat2(:));
%   F10 = figure(figbase+2); clf
  F10 = mkfig(figbase+2, 6, 7.5, false);
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
  if plot_mesh
    figure(figbase+4)
    mesh(pixmat2)
    xlabel('x')
    ylabel('y')
    title(stit)
  end
end
