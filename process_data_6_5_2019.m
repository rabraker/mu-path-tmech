clc
clear
%%
close all
addpath functions

% initialize paths.

size_dir = '5microns';
fname = fullfile(PATHS.sysid, 'x-axis_sines_info_intsamps_zaxisFourierCoef_10-29-2018-01.mat');
models = load(fname);
G = -models.modelFit.G_zdir;
p = pole(G);
z = zero(G);
gg = zpk(z(end-1:end), p(1:2), 1, G.Ts);
hole_depth = (20);

chan_map = ChannelMap(1:5);

% ----------------------- Load and Process raster-data -------------------------
fls = ls('/media/labserver/afm-cs/imaging/raster/5microns/6-5-2019/*.csv', '-t');
raster_files = strsplit(fls);
% The last element is always an empty string.
raster_files = flipud(raster_files(1:end-1)')


fls = ls('/media/labserver/afm-cs/imaging/cs-imaging/5microns/6-5-2019/*.csv', '-t');
cs_files = strsplit(fls);
% The last element is always an empty string.
cs_files = flipud(cs_files(1:end-1)')



rast_exps = cell(length(raster_files), 1);
for k=1:length(raster_files)
    [dat_root, dat_name, ext] = fileparts(raster_files{k});
    
    raster_paths = get_raster_paths(dat_root, [dat_name, ext], 'parents-loop');
    rast_exps{k} = RasterExp(raster_paths, 'load_full', true, 'reload_raw', false);
    if rast_exps{k}.time_total == 0
        rast_exps{k}.time_total = rast_exps{k}.samps_per_period*rast_exps{k}.npix*AFM.Ts;
    end
end
fprintf('finished loading raster data\n')


target_pix = 512;
use_ze = false;
% last image is 128 pixels
x1s = [36,   27,  27,  27,  27, 27, 6*4];
x2s = [493, 493, 493, 495, 494, 494, 122*4];
figbase = 10;
for k=1:1 %length(rast_exps)
    rast_exps{k}.uz = log_creep_detrend(rast_exps{k}.uz, []);
    
    if rast_exps{k}.npix_y_og ~= target_pix
        rast_exps{k}.bin_raster_really_slow([], use_ze, target_pix);
        rast_exps{k}.bin_raster_really_slow([], true, target_pix);
        rast_exps{k}.interp_y(target_pix, true);
        stit = sprintf('(raster %.1f sec) interpolated from %d pix %.2f Hz',...
            length(rast_exps{k}.ze)*AFM.Ts, rast_exps{k}.npix_y_og,...
            rast_exps{k}.meta_in.raster_freq);
        
    else
        rast_exps{k}.bin_raster_really_slow([], use_ze);
        rast_exps{k}.bin_raster_really_slow([], true);
        %rast_exps{k}.pix_mat = rast_exps{k}.interp_missing(rast_exps{k}.pix_mat);
        
        stit = sprintf('(raster %.1f sec) %.2f Hz', length(rast_exps{k}.ze)*AFM.Ts,...
            rast_exps{k}.meta_in.raster_freq);
    end
    
    %[idx_left, idx_right] = find_pin_idxs(rast_exps{k}.pix_mat);
    idx_left = 61;
    idx_right = 480;
    %   pixmat_ = rast_exps{k}.pix_mat;
    pixmat_ = pin_along_column(rast_exps{k}.pix_mat, idx_left, idx_right);
    %   pixmat_ = pixmats_raw{k};
    rast_exps{k}.pix_mat_pinned = pixmat_ - mean(pixmat_(:)); %mean(vec(pixmat_(:, 12:500)));
    rast_exps{k}.pin_idx_s = [idx_left, idx_right];
    if false
        [~, ax1] = plot_raster_data(rast_exps{k}.pix_mat_pinned, figbase*k, stit);
        hold(ax1, 'on')
        plot(ax1, [idx_left, idx_left], [1, 512], 'r')
        plot(ax1, [idx_right, idx_right], [1, 512], 'r')
    end
end

fprintf('Finished processing raster data\n');
if false
    save the first image so we can simulate reconstruction in plot_bptv_vs_bp.m
    
    img = rast_exps{1}.pix_mat_pinned;
    img = img - min(img(:));
    img = (img/max(img(:)) ) * 255/1.7;
    img = uint8(img);
    
    imwrite(img, 'cs20ng.png')
end

for k=1:length(rast_exps)
%     rast_exps{k}.save(true)
end


%%
cs_exps = cell(length(cs_files), 1);
for k=1:length(cs_files)
    [dat_root, dat_name, ext] = fileparts(cs_files{k});
    cs_paths = get_cs_paths(dat_root, [dat_name, ext], 'parents-loop');
    
    gg = @(u, idx_state_s) log_creep_detrend(u, idx_state_s);
    
    cs_exps{k} = CsExp(cs_paths, 'feature_height', hole_depth, 'gg', gg,...
        'load_full', true, 'reload_raw', false);
%     cs_exps{k}.uz = fft_notch(cs_exps{k}.uz, AFM.Ts, 212, 218);
    cs_exps{k}.print_state_times();
end
fprintf('finished loading cs data\n')

bp = true;
recalc = false;
use_dct2 = false;
thresh = (20/7)*(1/1000)*20;

use_ze=false;

register_uzk = true;
% fprintf('  Hz  |  actual perc  |  nom. perc  |  mu-path overhead |\n')
% fprintf('---------------------------------------------------------\n')
for k=1:length(cs_exps)
    cs_exps{k}.process_cs_data(use_ze, register_uzk);
    cs_exps{k}.process_cs_data(true, false);
    %     fprintf('finished processing raw CS data..\n');
    fprintf('  Hz  |  actual perc  |  nom. perc  |  mu-path overhead |\n')
    fprintf('---------------------------------------------------------\n')
    fprintf('  %.1f     %.3f          %.2f            %.4f\n',...
        cs_exps{k}.equiv_raster_rate,...
        cs_exps{k}.sub_sample_frac()*100,...
        cs_exps{k}.meta_in.actual_sub_samble_perc,...
        cs_exps{k}.get_mean_mu_overhead());
    
    ht = cs_exps{k}.feature_height;
    U_fun = @(x) idct(x);
    Ut_fun = @(x) dct(x);
    
    % U_fun = @(x) CsTools.Ufun_dct2(x, 512);
    % Ut_fun = @(x) CsTools.Utfun_dct2(x, 512);
    opts = NESTA_opts('U', U_fun, 'Ut', Ut_fun, 'alpha_v', 0.1, 'alpha_h', .75,...
        'verbose', 0, 'TolVar', 1e-5, 'mu', 1e-5, 'sigma', 1e-2);
    if 1
        cs_exps{k}.solve_nesta(recalc, use_dct2, use_ze, opts);
        cs_exps{k}.solve_nesta(recalc, use_dct2, true, opts);
        %[idx_left, idx_right] = find_pin_idxs(cs_exps{k}.pix_mat_uz);
        % idx_left = idx_left - 5;
        % cs_exps{k}.pix_mat_uz = pin_along_column(cs_exps{k}.pix_mat_uz, idx_left, idx_right);
    end
    fprintf('Finished solving bp problem #%d\n', k);
    stit = sprintf('(CS) %.2f Hz equiv, \\%% %.2f sampling\nTotal time: %.2f',...
        cs_exps{k}.meta_in.tip_velocity/(cs_exps{k}.meta_in.width * 2),...
        cs_exps{k}.meta_in.actual_sub_samble_perc,...
        length(cs_exps{k}.x)*cs_exps{k}.Ts);
    
    bar_bp = mean(cs_exps{k}.pix_mat_uz(:));
    cs_exps{k}.pix_mat_uz = cs_exps{k}.pix_mat_uz - bar_bp;
    cs_exps{k}.pix_mat_raw_uz = cs_exps{k}.pix_mat_raw_uz - bar_bp;
       
    if false
        f10=mkfig(1001 + 2*k, 6, 7.5); clf
        ax4 = subplot(3,1,[1,2]);
        ax4_2 =subplot(3,1,3);
        ImshowDataView.setup(f10);
        cb_exp =  @(event_obj)cs_exps{k}.dataview_callback(event_obj, ax4, ax4_2);
        im_tmp = cs_exps{k}.pix_mat_uz;
        % im_tmp = detrend_plane(im_tmp);
        im_tmp = im_tmp - mean(im_tmp(:));
        % im_tmp = SplitBregmanROF(im_tmp, 100, 0.001);
        ImshowDataView.imshow(im_tmp, [-thresh, .25*thresh], ax4, ax4_2, cb_exp)
        title(ax4, stit)
        hold(ax4, 'on')
        drawnow();
        
        %plot(ax4, [idx_left, idx_left], [1, 512], 'r')
        %plot(ax4, [idx_right, idx_right], [1, 512], 'r')
    end
end    
write_cs_meta_data(cs_exps, opts, 'latex/cs_data.txt');
fprintf('Finished processing CS data\n');
% 
for k=1:length(cs_exps)
%     cs_exps{k}.save()
end


% [~, axs1] = make_cs_traj_figs(figbase, 3);
% cs_exps{1}.plot_all_cycles(axs1{:}, [],[], 512);
% %%
% [~, axs2] = make_cs_traj_figs(figbase+100, 4);
% % cs_exps{end}.plot_all_cycles(axs2{1:4});
% cs_exps{end}.plot_all_cycles([], [], axs2{1}, [],[], 512);
%%

thresh = (20/7)*(1/1000)*20;
DRng = 2*thresh;

slice = 30:512-30;
master_rast_idx = 1;
assert(rast_exps{master_rast_idx}.npix_y == 512)
assert(abs(rast_exps{master_rast_idx}.meta_in.raster_freq - 1) < 1e-3)
im_rast_master = rast_exps{master_rast_idx}.pix_mat_pinned; % - mean(rast_exps{master_idx}.pix_mat_pinned(:));

master_cs_idx = 1;
assert(rast_exps{master_cs_idx}.npix_y == 512)
assert(abs(rast_exps{master_cs_idx}.meta_in.raster_freq - 1) < 1e-3)
im_cs_master = rast_exps{master_cs_idx}.pix_mat_pinned; % - mean(rast_exps{master_idx}.pix_mat_pinned(:));


fprintf('---------------------------------------------------\n\n');
raster_metrics = cell(length(rast_exps), 1);
cs_metrics = cell(length(cs_exps), 1);

excludes = [];

err_idx = 1;
plt_idx = 1;
row_plt_idx = 1;
rastm_512 = ScanMetrics();
rastm_128 = ScanMetrics();
rastm_64 = ScanMetrics();
rast_imk_fit = {};

dct_ = @(x) dct(x);
for k=1:length(rast_exps)
    if k== master_rast_idx
        continue
    end
    imk = rast_exps{k}.pix_mat_pinned; % - mean(rast_exps{k}.pix_mat_pinned(:));
   
    im_master_slice = im_rast_master(slice, slice);
    imk_slice = norm_align(im_master_slice, imk);
    rast_imk_fit{k} = imk_slice;
    %   [im1_ontok_fit, imk_slice] = align_by_metric(im_master, imk, [], 'psnr');
    [psn_1k, ssm_1k] = ssim_psnr_norm(im_master_slice, imk_slice, DRng);
    % damage = rast_exps{k}.damage_metric();
    quality = rast_exps{k}.quality_metric();
    damage = quality;
    time_k = length(rast_exps{k}.x)*AFM.Ts;
    rast_exps{k}.pix_mat_ze = rast_exps{k}.pix_mat_ze - mean(rast_exps{k}.pix_mat_ze(:));
    rate_k = rast_exps{k}.meta_in.raster_freq;
    
    args = {'ssim', ssm_1k, 'psnr', psn_1k,...
        'quality', quality, 'damage', damage,...
        'rate', rate_k,...
        'time', rast_exps{k}.time_total,...
        'coverage', 100,...
        'type', 'raster',...
        'ypix', rast_exps{k}.npix_y_og};

    if abs(rate_k - 2.5) < .5
        %continue
    end
    switch rast_exps{k}.npix_y_og
        case 512
            rastm_512.append_metrics(args{:});
        case 128
            rastm_128.append_metrics(args{:});
        case 64
            rastm_64.append_metrics(args{:});
        otherwise
            error('unknown ypix')
    end
    
    fprintf('(raster %.1f Hz) psnr: %.4f, ssm: %.4f, damage: %.4f, time: %.4f\n',...
        rate_k, psn_1k, ssm_1k, damage, time_k);    
end
%

csm_12 = ScanMetrics();
csm_15 = ScanMetrics();
csm_25 = ScanMetrics();

csm_12_full = ScanMetrics();
csm_15_full = ScanMetrics();
csm_25_full = ScanMetrics();

cs_imk_fit = {};
for k=1:length(cs_exps)
    
    imk = cs_exps{k}.pix_mat_uz; % - mean(cs_exps{k}.pix_mat_uz(:));
    
    im_master_slice = im_cs_master(slice, slice);
    imk_slice = norm_align(im_master_slice, imk);
    cs_imk_fit{k} = imk_slice;
    
    [psn_1k, ssm_1k] = ssim_psnr_norm(im_master_slice, imk_slice, DRng);
    
    % damage = cs_exps{k}.damage_metric();
    quality = cs_exps{k}.quality_metric();
    damage = quality;
    time_k = length(cs_exps{k}.x)*cs_exps{k}.Ts;
    rate_k = cs_exps{k}.meta_in.tip_velocity/(cs_exps{2}.meta_in.width * 2);
    coverage_k = cs_exps{k}.meta_in.actual_sub_samble_perc;
    
    cs_exps{k}.pix_mat_ze = cs_exps{k}.pix_mat_ze - mean(cs_exps{k}.pix_mat_ze(:));
    % sum(cs_exps{k}.pix_mat_ze(:).^2/512^2)
    args = {'ssim', ssm_1k, 'psnr', psn_1k,...
        'quality', quality, 'damage', damage,...
        'rate', rate_k,...
        'time', time_k,...
        'coverage', coverage_k,...
        'type', 'CS'};
    if rate_k == 4 || rate_k== 2
        %continue
    end
    if abs((cs_exps{k}.meta_in.actual_sub_samble_perc - 25)) < 2
        csm_25.append_metrics(args{:});
        if psn_1k > 25
            %keyboard
        end
    elseif abs((cs_exps{k}.meta_in.actual_sub_samble_perc - 15)) < 1
        csm_15.append_metrics(args{:});
    elseif abs((cs_exps{k}.meta_in.actual_sub_samble_perc - 12.5)) < 2
        csm_12.append_metrics(args{:});
    else
        error('no known percent')
    end
    
    fprintf('(CS %.1f Hz, %% %.1f, %.1f) psnr: %.4f, ssm: %.4f, damage: %.4f, time: %.4f\n',...
        rate_k, coverage_k, cs_exps{k}.sub_sample_frac*100, psn_1k, ssm_1k, damage, time_k);
    
    [t_cycle_avg, t_connect] = cs_exps{k}.estimate_mpt_connect_savings();
    
end

%
%%
% Fig1 = mkfig(3000, 5.5, (5.5/3)*5); clf;
% ha1 = tight_subplot(5, 3, [0.022, 0.01], [.02, .02], [.01, .01], true);
% ha1 = reshape(ha1', 3, [])';
Fig1 = mkfig(3000, 5.5, 9.1*(3/5)); clf;
ha1 = tight_subplot(3, 3, [0.03, 0.01], [.01, .03], [.01, .01], true);
ha1 = reshape(ha1', 3, [])';

Fig_subl = mkfig(3001, 5.6, 9.1*(2/5)); clf
ha_subl = tight_subplot(2, 3, [0.04, 0.01], [.01, .04], [.01, .01], true);
ha_subl = reshape(ha_subl', 3, [])';

Fig_rows = mkfig(3002, 7, 4.5); clf
ha_row = tight_subplot(2, 1, [0.1, 0.015], [.1, .05], [.085, .02], false);
xlabel(ha_row(2), '$X$-direction pixel', 'FontSize', 14)
ylabel(ha_row(1), 'height [nm]', 'FontSize', 14)
title(ha_row(1), 'raster', 'FontSize', 14)
title(ha_row(2), 'CS', 'FontSize', 14)
ylabel(ha_row(2), 'height [nm]', 'FontSize', 14)

grid(ha_row(1), 'on')
grid(ha_row(2), 'on')
row_idx = 169;

row_hands = gobjects(3, 1);

clr_range = [-thresh, thresh];
for k=1:length(rast_exps)
    switch rast_exps{k}.npix_y_og
        case 512
            fig_row = 1;
            ax_im = ha1(1, :);
        case 128
            fig_row = NaN;
            ax_im = ha_subl(1, :);
        case 64
            fig_row = NaN;
            ax_im = ha_subl(2, :);
        otherwise
            fprintf('Dony know what to do with npix = %d\n', rast_exps{k}.npix_y_og)
            continue
    end
    switch round(rast_exps{k}.meta_in.raster_freq, 1)
        case 1.0
            fig_col = 1;
        case 5.0 
            fig_col = 2;
        case 8
            fig_col = 3;
        case 15 
            fig_col = 4;
        otherwise
            fprintf('Dony know what to do with raster freq = %f\n',...
                rast_exps{k}.meta_in.raster_freq)
            continue
    end

    rate_k = rast_exps{k}.meta_in.raster_freq;
    

    
    ax_im = ax_im(fig_col);
    
    imk = rast_exps{k}.pix_mat_pinned;% - mean(rast_exps{k}.pix_mat_pinned(:));
    
    imagesc(ax_im, imk, clr_range);
    colormap(ax_im, 'gray');
    stit = sprintf('raster %d (%.1f Hz)', rast_exps{k}.npix_y_og, rate_k);
    title(ax_im, stit, 'FontSize', 7.5)
        



    if fig_row == 1
        hold(ha1(fig_row, fig_col), 'on')
        hold(ha_row(1), 'on')
        ha_row(1).ColorOrderIndex = fig_col;
        plot(ha1(fig_row, fig_col), [1, 512], [row_idx, row_idx], 'r');
        if abs(rate_k - 1.0)<0.1 && k~=1
            continue
        end
        row_hands(fig_col) = plot(ha_row(1), imk(row_idx, :)*AFM.volts2nm_z());
        row_hands(fig_col).DisplayName = sprintf('%.1f Hz', rate_k);
        hold(ha_row(1), 'on')
        %keyboard
    end
%     pause

end
remove_ticks(ha1)
remove_ticks(ha_subl)

set(ha_row(1), 'XLim', [1, 512])
set(ha_row(1), 'YLim', [-20, 10])
% axes(ha_row(1))
leg1 = legend(row_hands);
set(leg1, 'NumColumns', 4, 'FontSize', 11, 'Position', [0.4149 0.8945 0.5703 0.0500])
%


fprintf('---------------------------------------------------\n');

do_rows = [1, 4, 6, 10];
row_hands = gobjects(4, 1);

row_iter = 1;
for k=1:length(cs_exps)
    switch round(cs_exps{k}.meta_in.actual_sub_samble_perc, 0)
%         round(cs_exps{k}.sub_sample_frac*100, 0) 
        case 12
            fig_row = 2;
        case 13
            fig_row = 2;            
        case 25
            fig_row = 3;
        otherwise
            fprintf('Dony know what to do with sample frac = %.0f\n', cs_exps{k}.sub_sample_frac*100)
            continue
    end
    switch round(cs_exps{k}.equiv_raster_rate, 1)
        case 1.0
            fig_col = 1;
        case 5
            fig_col = 2;
        case 8 
            fig_col = 3;
        otherwise
            fprintf('Dony know what to do with rate = %f\n', cs_exps{k}.meta_in.raster_freq)
            continue
    end    

    imk = cs_exps{k}.pix_mat_uz;
    imagesc(ha1(fig_row, fig_col), imk, clr_range);
    colormap(ha1(fig_row, fig_col), 'gray');
    stit = sprintf('CS (%.1f Hz, %.1f \\%%)',...
        cs_exps{k}.equiv_raster_rate(), cs_exps{k}.sub_sample_frac()*100);
    title(ha1(fig_row, fig_col), stit, 'FontSize', 7.5);
    
    if abs(cs_exps{k}.equiv_raster_rate()-2)< .1 ||...
            (abs(cs_exps{k}.equiv_raster_rate()- 5)<.1 && abs(cs_exps{k}.sub_sample_frac()*100 - 12)< 1) ||...
        (abs(cs_exps{k}.equiv_raster_rate()-1)<.1 && abs(cs_exps{k}.sub_sample_frac()*100-25)<1.5)
        continue
    end
    
    hold(ha1(fig_row, fig_col), 'on')
    hold(ha_row(2), 'on')
    plot(ha1(fig_row, fig_col), [1, 512], [row_idx, row_idx], 'r');    
    row_hands(row_iter) = plot(ha_row(2), imk(row_idx, :)*AFM.volts2nm_z());
    hold(ha_row(2), 'on')
    
    row_hands(row_iter).DisplayName = sprintf('%.1f Hz, %.1f \\%%',...
        cs_exps{k}.equiv_raster_rate(), cs_exps{k}.sub_sample_frac()*100);
        %end
        plt_idx = plt_idx + 1;

    
    row_iter = row_iter + 1;
end
%
remove_ticks(ha1)
% remove_ticks(ha_err)

set(ha_row(2), 'YLim', [-20, 10])
set(ha_row(2), 'XLim', [1, 512])
axes(ha_row(2))
leg2 = legend();
set(leg2, 'NumColumns', 4, 'FontSize', 11, 'Position', [0.1032 0.4253 0.8589 0.0500])

%%
color_ord_map = [1, 2.0, 2.5, 4.0, 5.0, 8.0, 10.0;
                 1, 2,   3,   4,   5,   6,   7];
             
F_rdi = mkfig(3, 4., 3); clf
ha = tight_subplot(1, 1, .01, [0.1, 0.1], [0.1, 0.02], false);
ha.FontSize = 8;
hold(ha(1), 'on')
grid(ha(1), 'on')
set(ha(1), 'XScale', 'log');

h3 = semilogx_color_points2(ha(1), rastm_512, 'damage', 'x', 1, color_ord_map);
h3a = semilogx_color_points2(ha(1), rastm_128, 'damage',  '+', 1, color_ord_map);
h3b = semilogx_color_points2(ha(1), rastm_64, 'damage',  's', 1, color_ord_map);

h4 = semilogx_color_points2(ha(1), csm_12, 'damage', 'o', 1, color_ord_map);
h5 = semilogx_color_points2(ha(1), csm_25, 'damage', '*', 1, color_ord_map);

plot(ha(1), rastm_512.time(1:end), rastm_512.damage(1:end), ':k')
plot(ha(1), rastm_128.time, rastm_128.damage, ':k')
plot(ha(1), rastm_64.time, rastm_64.damage, ':k')
plot(ha(1), csm_12.time, csm_12.damage, ':k')
plot(ha(1), csm_25.time, csm_25.damage, ':k')

xlabel(ha(1), 'total time [s]')
ylabel(ha(1), 'RDI')

ha_none = axes('Position', [0.1473 0.8912 0.1788 0.0358]);
hold(ha_none, 'on');
h3_ = plot(1,1, 'xk', 'DisplayName', 'raster 512 lines');
h3a_ = plot(1,1, '+k', 'DisplayName', 'raster 128 lines');
h3b_ = plot(1,1, 'sk', 'DisplayName', 'raster 64 lines');
h4_ = plot(1,1, 'ok', 'DisplayName', 'CS: 12.5\%');
h5_ = plot(1,1, '*k', 'DisplayName','CS: 25\%');
leg0 = legend([h3_, h3a_, h3b_, h4_, h5_], 'Position', [0.0998 0.9108 0.8724 0.0639], 'NumColumns', 3);
ha_none.Visible = 'off';
xlim(ha_none, [0, 0.01])

xs = logspace(log10(6.25), log10(600), 8); %[7, 15, 40, 100, 200];
wo = 7;
yo = 35;
h = 3;
fac = (xs(1) + wo)/xs(1);
text_color = [0.97, 0.97, 0.97];
for k=1:7
% rate_k = rastm_512.rate(k);
rate_k = color_ord_map(1, k);
clr = ha.ColorOrder;

x1 = xs(k);

x2 = x1*fac;
w = (xs(k+1) - xs(k))*.92;

X = [x1, x1+w, x1+w, x1];
Y = [yo, yo, yo+h, yo+h];
C= clr(k,:);
pt = patch(ha, X, Y, C);

st = sprintf('%.1f Hz', rate_k);
tx = text(ha, mean(X(1:2))-0.5*x1/wo, mean(Y([1,3])), st, 'VerticalAlignment', 'middle',...
    'HorizontalAlignment', 'center', 'Color', text_color);

end
ylim(ha, [0, 38])
xlim(ha, [6, 600])
%%

F_ssm = mkfig(4, 4, 3.); clf
ha_ssm = tight_subplot(1, 1, .01, [0.1, 0.1], [0.1, 0.02], false);
ha_ssm.FontSize = 8;
xlabel(ha_ssm, 'total time [s]')
ylabel(ha_ssm, 'SSIM')
hold(ha_ssm, 'on')
grid(ha_ssm, 'on')

set(ha_ssm, 'XScale', 'log')

h3 = semilogx_color_points2(ha_ssm, rastm_512, 'ssim', 'x', 1, color_ord_map);
h3a = semilogx_color_points2(ha_ssm, rastm_128, 'ssim',  '+', 1, color_ord_map);
h3b = semilogx_color_points2(ha_ssm, rastm_64, 'ssim',  's', 1, color_ord_map);

h4 = semilogx_color_points2(ha_ssm, csm_12, 'ssim', 'o', 1, color_ord_map);
h5 = semilogx_color_points2(ha_ssm, csm_25, 'ssim', '*', 1, color_ord_map);

plot(ha_ssm, rastm_512.time(2:end), rastm_512.ssim(2:end), ':k')
plot(ha_ssm, rastm_128.time, rastm_128.ssim, ':k')
plot(ha_ssm, rastm_64.time, rastm_64.ssim, ':k')
plot(ha_ssm, csm_12.time, csm_12.ssim, ':k')
plot(ha_ssm, csm_25.time, csm_25.ssim, ':k')

% I want all the marker legends to be the same color. The idea here is to plot
% dummy markers on a second axes, produce the legend for that, then make the
% other axis not visible, and set its xlim to be smaller than the x-value of the
% markers, which will make them invisible.
ha_none = axes('Position', [0.1473 0.8912 0.1788 0.0358]);
hold(ha_none, 'on');
h4_ = plot(1,1, 'xk', 'DisplayName', 'raster 512 lines');
h4a_ = plot(1,1, '+k', 'DisplayName', 'raster 128 lines');
h4b_ = plot(1,1, 'sk', 'DisplayName', 'raster 64 lines');
h5_ = plot(1,1, 'ok', 'DisplayName', 'CS: 12.5\%');
h6_ = plot(1,1, '*k', 'DisplayName','CS: 25\%');
leg1 = legend([h4_, h4a_, h4b_, h5_, h6_], 'Position', [0.0998 0.8949 0.8724 0.0956], 'NumColumns', 3);
ha_none.Visible = 'off';
xlim(ha_none, [0, 0.01])


ylim(ha_ssm, [0.65, 0.83])
% xs = logspace(log10(7), log10(600), 8); %[7, 15, 40, 100, 200];
xs = logspace(log10(6.25), log10(230), 8); %[7, 15, 40, 100, 200];
wo = 7;
yo = 0.815;
h = 0.015;
fac = (xs(1) + wo)/xs(1);
for k=1:7
% rate_k = rastm_512.rate(k);
rate_k = color_ord_map(1, k);
clr = ha_ssm.ColorOrder;

x1 = xs(k);

x2 = x1*fac;

w = (xs(k+1) - xs(k))*.9;

X = [x1, x1+w, x1+w, x1];
Y = [yo, yo, yo+h, yo+h];
C= clr(k,:);
pt = patch(ha_ssm, X, Y, C);

st = sprintf('%.1f Hz', rate_k);
tx = text(ha_ssm, mean(X(1:2))-0.25*x1/wo, mean(Y([1,3])), st,...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'Color', text_color);

end
xlim(ha_ssm, [6, 230])

%%
% -------------- PSNR
F_psn = mkfig(5, 4, 3.); clf
ha_psn = tight_subplot(1, 1, .01, [0.1, 0.1], [0.1, 0.02], false);
ha_psn.FontSize = 8;
ylabel(ha_psn, 'PSNR')
xlabel(ha_psn, 'total time [s]')

hold(ha_psn, 'on')
grid(ha_psn, 'on')
set(ha_psn, 'XScale', 'log')


h3 = semilogx_color_points2(ha_psn, rastm_512, 'psnr', 'x', 1, color_ord_map);
h3a = semilogx_color_points2(ha_psn, rastm_128, 'psnr',  '+', 1, color_ord_map);
h3b = semilogx_color_points2(ha_psn, rastm_64, 'psnr',  's', 1, color_ord_map);
h4 = semilogx_color_points2(ha_psn, csm_12, 'psnr', 'o', 1, color_ord_map);
h5 = semilogx_color_points2(ha_psn, csm_25, 'psnr', '*', 1, color_ord_map);

plot(ha_psn, rastm_512.time(2:end), rastm_512.psnr(2:end), ':k')
plot(ha_psn, rastm_128.time, rastm_128.psnr, ':k')
plot(ha_psn, rastm_64.time, rastm_64.psnr, ':k')
plot(ha_psn, csm_12.time, csm_12.psnr, ':k')
plot(ha_psn, csm_25.time, csm_25.psnr, ':k')

ha_none = axes('Position', [0.1473 0.8912 0.1788 0.0358]);
hold(ha_none, 'on');
h7_ = plot(1,1, 'xk', 'DisplayName', 'raster 512 lines');
h7a_ = plot(1,1, '+k', 'DisplayName', 'raster 128 lines');
h7b_ = plot(1,1, 'sk', 'DisplayName', 'raster 64 lines');
h8_ = plot(1,1, 'ok', 'DisplayName', 'CS: 12.5\%');
h9_ = plot(1,1, '*k', 'DisplayName','CS: 25\%');
leg2 = legend([h7_, h7a_, h7b_, h8_, h9_], 'Position', [0.0998 0.8914 0.8724 0.0956], 'NumColumns', 3);
ha_none.Visible = 'off';
xlim(ha_none, [0, 0.01])

xs = logspace(log10(6.25), log10(230), 8); %[7, 15, 40, 100, 200];
wo = 5;
fac = (xs(1) + wo)/xs(1);
h = 0.6;
yo = 24.4;
for k=1:7
% rate_k = rastm_512.rate(k);
rate_k = color_ord_map(1, k);
clr = ha_psn.ColorOrder;

x1 = xs(k);
x2 = x1*fac;

w = (xs(k+1) - xs(k))*.9;

X = [x1, x1+w, x1+w, x1];
Y = [yo, yo, yo+h, yo+h];
C= clr(k,:);
pt = patch(ha_psn, X, Y, C);

st = sprintf('%.1f Hz', rate_k);
tx = text(ha_psn, mean(X(1:2))-0.25*x1/wo, mean(Y([1,3])), st,...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center', 'Color', text_color);

end
ylim(ha_psn, [18, 25])
xlim(ha_psn, [6, 230])
%%
% set(ha1(end), 'Visible', 'off')
% set(ha_err(end), 'Visible', 'off')
save_fig(Fig1, fullfile(PATHS.tmech_fig(), 'cs_raster_images_6-5-2019'), false)
%%
save_fig(Fig_subl, fullfile(PATHS.tmech_fig(), 'cs_raster_images_subl_6-5-2019'),false)
save_fig(Fig_rows, fullfile(PATHS.tmech_fig(), 'cs_raster_pixel_rows_6-5-2019'), false)


%%
save_fig(F_rdi, fullfile(PATHS.tmech_fig(), 'cs_rast_damage'), false)
save_fig(F_ssm, fullfile(PATHS.tmech_fig(), 'cs_rast_time_vs_ssim'), false)
save_fig(F_psn, fullfile(PATHS.tmech_fig(), 'cs_rast_time_vs_psnr'), false)
%%

F_trade = mkfig(10, 3.5, 2.25); clf
ha_t = tight_subplot(1, 1, 0.01, [0.15, 0.05], [0.15, 0.02], false);
% plot(rastm_512.rate, rastm_512.time)
hold on

nrm_s = rastm_512.time;
h1 = plot(rastm_64.rate, rastm_64.time./nrm_s, '-o');
hold on
nrm_s_cs12 = 512./csm_12.rate;
h2 = plot(csm_12.rate, csm_12.time./nrm_s_cs12, '--o');


h3 = plot(rastm_128.rate, rastm_128.time./nrm_s, '-o');

% h4a = plot(csm_15.rate, csm_15.time./nrm_s, '--o');
nrm_s_cs25 = 512./csm_25.rate;
h4 = plot(csm_25.rate, csm_25.time./nrm_s_cs25, '--o');


% title('25~\% sampling')
% plot(csm_15.rate, csm_15.time)
h1.DisplayName = 'raster 64 lines';
h2.DisplayName = 'CS 12.5 \%';

h3.DisplayName = 'raster 128 lines';
% h4a.DisplayName = 'CS 15 \%';
h4.DisplayName = 'CS 25 \%';
leg_t = legend([h1, h3, h2, h4], 'Position', [0.1706 0.6866 0.3648 0.2655],...
    'NumColumns', 1, 'Box', 'off');

xlabel('rate [Hz]')
ylabel('fraction of 512 line scan time')
grid on
%%
save_fig(F_trade, 'latex/figures/improvements', false)
%%
clc
tmu_12 = [];
tmu_15 = [];
tmu_25 = [];
rates_12 = [];
rates_15 = [];
rates_25 = [];
for k=1:length(cs_exps)
   tmu = cs_exps{k}.get_mean_mu_overhead();
   cs_exps{k}.sub_sample_frac()*100;
   if abs(cs_exps{k}.sub_sample_frac()*100 - 12) < 1
       tmu_12(end+1) = tmu;
       rates_12(end+1) = cs_exps{k}.equiv_raster_rate();
   elseif abs(cs_exps{k}.sub_sample_frac()*100 - 15) <  1
       tmu_15(end+1) = tmu;
       rates_15(end+1) = cs_exps{k}.equiv_raster_rate();
   elseif abs(cs_exps{k}.sub_sample_frac()*100 - 25) < 1.5
       tmu_25(end+1) = tmu;
       rates_25(end+1) = cs_exps{k}.equiv_raster_rate();
   end
end

figure(15)
x12 = 12*ones(length(tmu_12), 1);
x15 = 15*ones(length(tmu_15), 1);
x25 = 25*ones(length(tmu_25), 1);
plot(x12, tmu_12, 'bx')
hold on
plot(x15, tmu_15, 'rx')
plot(x25, tmu_25, 'mx')

save('tmu_s.mat', 'tmu_12', 'tmu_15', 'tmu_25', 'rates_12', 'rates_15', 'rates_25');

%%
% skip the 0.5Hz scan
% sm = scan_metrics(1:end)
% S1 = scan_metrics_table(sm)
% fid = fopen(fullfile(PATHS.tmech_table(), 'cs_raster_table_6-5-2019_muInf_dct2.tex'), 'w+');
% fprintf(fid, '%s', S1);
% fclose(fid);

hz_skip = [2.0, 4.0];
frac_skip = 15.02;
S2 = state_times_table(cs_exps, hz_skip, frac_skip)
fid = fopen(fullfile(PATHS.tmech_table(), 'cs_state_times_table_6-5-2019_muInf_dct2.tex'), 'w+');
fprintf(fid, '%s', S2);
fclose(fid);
%%



function write_cs_meta_data(cs_exps, nesta_opts, fname)
    
%     fracs = zeros(length(cs_exps), 1);
%     rates = zeros(length(cs_exps), 1);
%     
%     for k=1:length(cs_exps)
%        fracs(k) = cs_exps{k}.meta_in.actual_sub_samble_perc; 
%        rates(k) = cs_exps{k}.equiv_raster_rate(); 
%     end
%     fracs = unique(round(fracs, 1));
%     rates = unique(round(rates, 1));

% latex catchfilebetweentags doesnt seem to like snake_case, so use CamelCase.
    mu_len_mic = cs_exps{1}.meta_in.mu_length;
    mu_len_pix = cs_exps{1}.meta_in.npix * (mu_len_mic / cs_exps{1}.meta_in.width);
    scan_meta = struct('preScanSamples', cs_exps{1}.meta_in.pre_scan_samples,...
        'xyBoundaryMic', cs_exps{1}.meta_exp.state_machine_params.xy_error_threshold*AFM.volts2mic_xy,...
        'xySetSamps', cs_exps{1}.meta_exp.state_machine_params.xy_settled_samples_threshold,...
        'muLenNM', mu_len_mic*1000,...
        'muLenPix', round(mu_len_pix, 0),...
        'spScan', cs_exps{1}.meta_exp.z_axis_params.setpoint_scan,...
        'spMv', cs_exps{1}.meta_exp.z_axis_params.setpoint_up,...
        'alph', nesta_opts.alpha_h,...
        'alpv', nesta_opts.alpha_v,...
        'muf', nesta_opts.mu,...
        'sigma', nesta_opts.sigma);
    
    data_writer({scan_meta}, fname);
end

function ha = semilogx_color_points2(ax, SM, fld, marker, ord_start, ord_map)
    if nargin < 5
        ord_start = 1;
    end
    ax.ColorOrderIndex = ord_start;
    colors = ax.ColorOrder;
    for k =1:length(SM.time)
      idx = find(abs(ord_map(1, :) - SM.rate(k)) < .1, 1, 'first') + ord_start-1;
      clr = colors(idx, :);
       ha = semilogx(ax, SM.time(k), SM.(fld)(k),  marker, 'MarkerSize', 8,...
                     'color', clr);
    end
end

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

function S = state_times_table(cs_exps, hz_skip, frac_skip)
    S = 'description & move & engage & pre-scan & scan & tip-up & total\\';
    S = sprintf('%s\n\\toprule\n', S);
    for k=1:length(cs_exps)
        cst = cs_exps{k}.get_state_times();
        rate = cs_exps{k}.equiv_raster_rate();
%         frac = cs_exps{k}.meta_in.actual_sub_samble_perc;
        frac = cs_exps{k}.sub_sample_frac()*100;
        if any(abs(frac_skip - frac) < 0.5)
            continue
        end
        if any(abs(hz_skip - rate) < 0.2)
            continue
        end

        description = sprintf('%.1f Hz, %.1f~\\%%',rate, frac);
        
        S = sprintf('%s%s &%.2f & %.2f & %.2f & %.2f & %.2f & %.2f\\\\\n', S,...
            description, cst.move, cst.tdown, cst.tsettle, cst.scan,...
            cst.tup, cst.total);
    end
end


%%
function [F10, ax1, ax2] = plot_raster_data(pixmat2, figbase, stit, plot_mesh)
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
