clc
clear
%%
% close all
clear
clc
addpath functions/scanning_v1/
addpath ~/matlab/dependencies/l1c/

% initialize paths.

figbase = 20;

size_dir = '5microns';
hole_depth = (20);

chan_map = ChannelMap([1:5]);
exp_date = '6-5-2019'
% ----------------------- Load and Process CS-data -----------------------------
dat_root = PATHS.raster_image_data(size_dir, exp_date);

% ----------------------- Load and Process raster-data -------------------------
cs_files = {...
'cs-traj-512pix-12perc-500nm-5mic-01Hz-150prescan-notconnect_out_6-5-2019-01.csv',...
};


data_root = PATHS.cs_image_data(size_dir, exp_date);
cs_exps = {};
for k=1:length(cs_files)
  cs_paths = get_cs_paths(data_root, cs_files{k});
  gg = @(u, idx_state_s) log_creep_detrend(u, idx_state_s);
  cs_exps{k} = CsExp(cs_paths, 'feature_height', hole_depth, 'gg', gg, 'load_full', true);

  cs_exps{k}.print_state_times();
  cs_exps{k}.sub_sample_frac()
  cs_exps{k}.t = (0:length(cs_exps{k}.x)-1)'*AFM.Ts;
end

% taken from cycle 3, where tip detaches.
cs_exp1 = cs_exps{1};
zdfl_free = -0.637-.05;
cs_exp1.ze = cs_exp1.ze - zdfl_free;


cs_exp1.x = cs_exp1.x*AFM.volts2mic_xy;
cs_exp1.y = cs_exp1.y*AFM.volts2mic_xy;
cs_exp1.uz = cs_exp1.uz*AFM.volts2nm_z;

%%
[~, axs] = make_cs_traj_figs(figbase, 4);
cs_exps{1}.plot_all_cycles(axs{1:4});



%%

indc = {   'k',        'r',   [0, .75, .75],       'b',        [.93 .69 .13], 'm';
          'xy-move', 'tip down', 'tip settle',  '$\mu$-path scan', 'tip up', 'connect'};
        
Fig = mkfig(100, 3.5, 4); clf
ha = tight_subplot(4, 1, [.05, .05], [0.075, .03], [.15, .03], false);


tstart = 7.25;
CS_idx1 = cs_exps{1}.find_cycle_idx(tstart);
CS_idx2 = CS_idx1 + 2;

%
% The Bad Case

state_names = {'move', 'tdown', 'tsettle', 'scan', 'tup'};

% The Good Case
traj_names = {'y', 'x', 'uz', 'ze'};
names = {'$Y$ [$\mu$m]', '$X$ [$\mu$m]', '$u_Z$ [nm]', 'deflection [v]'};



xlm1 = [7.286, 7.434];
hands = gobjects(length(state_names), 1);
for j=1:length(traj_names)
  h_j = ha(j);
  for k=1:length(state_names)
    state_name = state_names{k};
    h = cs_exp1.plot_traj_from_csidx_by_state(CS_idx1, CS_idx2,...
        state_name, traj_names{j}, h_j, 0, 'color', indc{1, k});
    h.DisplayName = indc{2, k};
    if strcmp(traj_names{j}, 'x')
        hands(k) = h;
    end
    
  end
  
  xlim(h_j, xlm1);
  grid(h_j, 'on')
  ylabel(h_j, [names{j} ], 'FontSize', 12);
  if j<3
    set(h_j, 'XTickLabel', [])
  else
    xlabel(h_j, 'time [s]')
  end
end

ylim(ha(2), [0, 5]);
ylim(ha(3), [-20, 90])
ylim(ha(4), [-0.15, 0.45])
hands(3).DisplayName = 'tip settle/pre-scan';

leg = legend(hands, 'NumColumns', 1, 'Position', [0.6064 0.6017 0.3808 0.1656]);
%%
save_fig(Fig, 'latex/figures/cs_cycle', false)
