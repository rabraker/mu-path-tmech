% 
% This script plots the z-bounce like experiment over a flat portion of the
% CS20ng grating. Both controllers are based on the loop-shaping transfer
% function control. The -03 experiment has the feedforward filter, the 02 does
% not.

% close all
clear PATHS
addpath functions
addpath ~/matlab/dependencies/l1c/
addpath ~/matlab/dependencies/l1c/lib

addpath ~/matlab/afm-cs/matlab-code/functions/
% initialize paths.

figbase = 20;


size_dir = '5microns';
hole_depth = (20);

chan_map = ChannelMap([1:5]);
% exp_date = '3-20-2019'
% ----------------------- Load and Process CS-data -----------------------------


data_root = '/media/labserver/afm-cs/z-bounce/6-12-2019/';
% ----------------------- Load and Process raster-data -------------------------
cs_files = {...
'cs-traj-512pix-3perc-500nm-5mic-01Hz-150prescan-notconnect_out_6-12-2019-03.csv',...
'cs-traj-512pix-3perc-500nm-5mic-01Hz-150prescan-notconnect_out_6-12-2019-02.csv',...
};

cs_exps = {};
for k=1:length(cs_files)
  cs_paths = get_cs_paths(data_root, cs_files{k});
  gg = @(u, idx_state_s) log_creep_detrend(u, idx_state_s);
  cs_exps{k} = CsExp(cs_paths, 'feature_height', hole_depth, 'gg', gg, 'load_full', true);
  cs_exps{k}.print_state_times();
end

name_s = {'CR w/o notch', 'CR-w/notch', 'CZ w/o notch', 'CZ w/ notch'};
for k=1:length(cs_exps)
   cs_exps{k}.uz  = cs_exps{k}.uz*AFM.volts2nm_z;
end

%%
if 0
    [~, axs1] = make_cs_traj_figs(figbase+100, 3);
    cs_exps{1}.plot_all_cycles(axs1{1:3}, [], [], 512);
    
    [~, axs2] = make_cs_traj_figs(figbase+200, 3);
    cs_exps{2}.plot_all_cycles(axs2{1:3}, [], [], 512);
    
    linkaxes([axs1{:}, axs2{:}], 'x')
end

width = 3.5;
height = 2.5;
F5 = mkfig(25, width, height); clf
margh = [0.1, 0.02];
margw = [0.105, 0.05];
axs2 = tight_subplot(2, 2, .06, margh, margw, false);
axs2 = reshape(axs2', 2, [])';

% title(axs1(1), '$xy$ fixed', 'FontSize', 14)
tstart1 = 2.103;
tstart2 = 2;
tend1 = 2.25;
tend2 = 2.15;
zbounce = false;
plot_uz_time_range(cs_exps{1}, axs2(1,2), tstart1, tend1, zbounce, 'uz');
plot_uz_time_range(cs_exps{1}, axs2(2,2), tstart1, tend1, zbounce, 'x');
zbounce = false;

plot_uz_time_range(cs_exps{2}, axs2(1, 1), tstart2, tend2, zbounce, 'uz');
h = plot_uz_time_range(cs_exps{2}, axs2(2, 1), tstart2, tend2, zbounce, 'x');

grid(axs2(1), 'on')
grid(axs2(2), 'on')
grid(axs2(3), 'on')
grid(axs2(4), 'on')

ylim(axs2(1,1), [-10, 50])
ylim(axs2(1,2), [-10, 50])

ylim(axs2(2,1), [-0.25, 6])
ylim(axs2(2,2), [-0.25, 6])

xlim(axs2(:,2), [tstart1, tend1])
xlim(axs2(:,1), [tstart2, tend2])

set(axs2(:, 2), 'YTickLabel', '')
set(axs2(1, :), 'XTickLabel', '')

ylabel(axs2(1), '$u_Z$~[nm]', 'FontSize', 10)
ylabel(axs2(2), '$x$~[$\mu$m]', 'FontSize', 10)
xlabel(axs2(2), 'time [s]', 'FontSize', 10)
xlabel(axs2(2,2), 'time [s]', 'FontSize', 10)

leg = legend(h, 'location', 'northeast','box', 'off', 'Position', [0.2360 0.3185 0.3085 0.1841]);

save_fig(F5, 'latex/figures/with_wo_ff', false)


function [UZ, UZ_pre_jump, UZ_jump, freqs] = psd_around_jumps(cs_exp, skips)
    N_mu = length(cs_exp.idx_state_s.scan);
    
    UZ = 0;
    UZ_pre_jump = 0;
    UZ_jump = 0;
    
    n_pre_jump=0;
    n_jump=0;

    for k = 1:N_mu
        idxk = cs_exp.idx_state_s.scan{k};
        uzk = cs_exp.uz(idxk);
        
        [UZK, freqs] = power_spectrum_local(uzk); %, AFM.Ts);
        
        UZ = UZK + UZ;
        
        if ~isempty(intersect(k, skips-1))
            UZ_pre_jump = UZK + UZ_pre_jump;
            n_pre_jump = n_pre_jump + 1;
        end
        
        if ~isempty(intersect(k, skips))
            UZ_jump = UZK + UZ_jump;
            n_jump = n_jump + 1;
        end
    end
    UZ = UZ/N_mu;
    UZ_pre_jump = UZ_pre_jump/n_pre_jump;
    UZ_jump = UZ_jump/n_jump;
end

function hands = plot_uz_time_range(cse, ax, tstart, tend, zbounce, signal)
indc = {'k',        'r',   [0, .75, .75],       'b',        [.93 .69 .13], 'm';
        'xy-move', 'tip down', 'tip settle',  '$\mu$-path scan', 'tip up', 'connect'};
    if signal == 'x'
        scl = 5;
    else
        scl = 1;
    end
    hold(ax, 'on')
    state_s = {'move', 'tdown', 'tsettle', 'scan', 'tup'};
    hands = gobjects(length(state_s), 1);
    for k=1:length(state_s)
        if k==4 && zbounce
            clr = indc{1,3};
        else
            clr = indc{1,k};
        end
        state = state_s{k};
        
        idxs = cse.get_idx_by_state_in_time_range(state, tstart, tend);
        
        for j=1:length(idxs)
            t = cse.t(idxs{j});
            uzk = cse.(signal)(idxs{j});
            h = plot(ax, t, uzk*scl, 'Color', clr);
        end
        h.DisplayName = indc{2, k};
        hands(k) = h;
    end
end


function [Pxx, freqs] = power_spectrum_local(X)

  N = length(X);

  [Pxx, freqs] = periodogram(detrend(X), [], N, 1/AFM.Ts);
end



