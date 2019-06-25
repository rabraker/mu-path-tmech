% Implement augmented state space integral control with shift register to
% 'estimate' the delay states.


clear all
%%
% clc
% Options
figbase  = 50;
verbose = 0;
rmpath('/home/arnold/matlab/afm-cs/matlab-code/functions/')
addpath('functions');
addpath(fullfile('functions', 'state_space_x'));


% ---- Paths for shuffling data to labview and back. ------


TOL = 14/512; % max volts by pixels
% TOL = .01;
tol_mode = 'abs';
% which simulations to run

do_sim_hyst = false;
do_inv_hyst = false;
do_drift = false;
do_invdrift = false;

plotstate = false;
plotpoles = false;
plotdriftpoles = false;
plotdriftbode = false;
saveon = true;

cntrl_type = 'const-sig';


% ------- Load Plants -----
[plants, frf_dataX] = CanonPlants.plants_ns14(9, '5micron');
Gxyz_frd = get_H();

Ts  = plants.SYS.Ts;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                  Design reference "trajectory"                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get a ref trajectory to track.
N  = 800;
r1 =1.37;

step_ref = StepRef([r1], N);
yref = step_ref.yref;
yref.Data = yref.Data*1;
step_descr = 'single_step';

dist_traj = yref;
dist_traj.Data = dist_traj.Data*0;
thenoise = dist_traj;


F_yudu = figure(60); clf
subplot(3,1,1)
hold on, grid on;
step_ref.plot(F_yudu, '-k', 'LineWidth', 0.5);

F_y = figure(61); clf
hold on, grid on
if max(abs(yref.Data)) > 0
  step_ref.plot(F_y);
  step_ref.plot_settle_boundary(F_y, TOL, tol_mode);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Design control/estimator gains. This gains are what we actually         %
% use in both the simulation and experimentant system.                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 3). LQR generation gain.
% -------------------------------------------------------------------------
% Adjust the du_max to account for the gain of gdrift_inv.
du_max_orig = StageParams.du_max;
du_max = du_max_orig/norm(plants.gdrift_inv, Inf);



xdir_tf_cntrl = get_xdir_loop_shaped_control();
%

ydir_cntrl = get_ydir_standard_control();


% Debug: (why is bw different between frd and model??) try replacing the xdir frf.
Gvibx_ = frd(frf_dataX.G_uz2stage, frf_dataX.freqs_Hz, AFM.Ts, 'FrequencyUnit', 'Hz');
Gvibx_ = frd(freqresp(Gvibx_, Gxyz_frd.Frequency*2*pi), Gxyz_frd.Frequency, AFM.Ts, 'FrequencyUnit', 'Hz');
Gxyz_frd(1,1) = Gvibx_(1,1);
gdi = blkdiag(plants.gdrift_inv, 1, 1);
Gxyz_frd = Gxyz_frd*gdi;
[DD, Dzki, Dz_inv] = get_gz_dz(Gxyz_frd(3,3));
%
g_static = zpk([], [], 1, AFM.Ts);
Dx_ff = xdir_tf_cntrl.M;
Dy_ff = ydir_cntrl.M;

xdir_tf_cntrl.M = g_static;
ydir_cntrl.M = g_static;
[HH1] = close_three_axis(Gxyz_frd, xdir_tf_cntrl, ydir_cntrl, Dzki, Dz_inv);

xdir_tf_cntrl.M = Dx_ff;
ydir_cntrl.M = Dy_ff;
[HH2] = close_three_axis(Gxyz_frd, xdir_tf_cntrl, ydir_cntrl, Dzki, Dz_inv);

%%
F = bode_local(HH1, HH2, Gxyz_frd, Dx_ff, Dy_ff);

save_fig(F, 'latex/figures/MIMO_CL_uxuy', false)
%
return
%%
if 0
    % xdir_cntrl = get_xdir_standard_control(cntrl_type);
    sys_obsDist = xdir_cntrl.sys_obsDist;
    K_lqr = xdir_cntrl.K_lqr;
    Nx = xdir_cntrl.Nx;
    L_dist = xdir_cntrl.L_dist;

    % -------------------------------------------------------------------
    % -------------------- Setup Fixed Point stuff ----------------------
    % For now, keep populating the state space stuff, until we are make it work
    % reliably without.
    A_obs_cl = sys_obsDist.a - L_dist*sys_obsDist.c;
    
    nw = 32;
    nf = 26;
    
    du_max_fxp = fi(du_max, 1, 32, 26);
    Nx_fxp = fi(Nx, 1, 32, 30);
    L_fxp = fi(L_dist, 1, 32, 30);
    
    sys_obs_fxp.a = fi(sys_obsDist.a -L_dist*sys_obsDist.c, 1, nw, nw-7);
    sys_obs_fxp.b = fi(sys_obsDist.b, 1, nw, 29);
    sys_obs_fxp.c = fi(sys_obsDist.c, 1, nw, 28);
    
    K_fxp = fi(K_lqr, 1, nw,32-10);
    sims_fxpl = SimAFM(plants.PLANT, K_fxp, Nx_fxp, sys_obs_fxp, L_fxp, du_max_fxp,...
        true, 'nw', nw, 'nf', nf, 'thenoise', thenoise);
    
   
    if 1
        sims_fxpl.r = plants.hyst_sat.r;
        sims_fxpl.w = plants.hyst_sat.w;
        sims_fxpl.rp = fi(plants.hyst_sat.rp, 1, 16, 11);
        sims_fxpl.wp = fi(plants.hyst_sat.wp, 1, 16, 11);
        sims_fxpl.d = plants.hyst_sat.d;
        sims_fxpl.ws = plants.hyst_sat.ws;
        sims_fxpl.dp = fi(plants.hyst_sat.dp, 1, 16, 11);
        sims_fxpl.wsp = fi(plants.hyst_sat.wsp, 1, 16, 11);
    end
    if 1
        sims_fxpl.gdrift_inv = plants.gdrift_inv;
        sims_fxpl.gdrift = plants.gdrift;
        sims_fxpl.Dx_ff = Dx_ff;
        sims_fxpl.Dy_ff = Dy_ff;
        
    end
    % ---------- Transfer function compensators. ---------------
    sims_fxpl.Dy = ydir_cntrl.D;
    sims_fxpl.Dx = xdir_tf_cntrl.D;
    sims_fxpl.Ki_y = ydir_cntrl.Ki;
    sims_fxpl.Ki_x = xdir_tf_cntrl.Ki;    
    
    fprintf('===========================================================\n');
    fprintf('Writing control data...\n');
    fprintf('===========================================================\n');
    
    
    sims_fxpl.sys_obs_fp = sys_obsDist;
    sims_fxpl.sys_obs_fp.a = sys_obsDist.a - L_dist*sys_obsDist.c;
  
    % sims_fxpl.write_control_data(controlDataPath, yref, traj_path)
    control_path = fullfile(PATHS.step_exp, sprintf('LinControls-%s_5micron_xyff_DyDx.json', cntrl_type))
    sims_fxpl.write_control_data_json(control_path);
    
    %%
    
    control_path_mat = fullfile(PATHS.step_exp, sprintf('LinControls-%s_5micron_xyff_DyDx.mat', cntrl_type))

    save(control_path_mat, 'xdir_tf_cntrl', 'xdir_cntrl', 'ydir_cntrl')
    
    
end


function F11 = bode_local(HH1, HH2, Gxyz_frd, Dx_ff, Dy_ff)
    F11 = mkfig(11, 7, 5.45, true); clf
    % [ha, pos] = tight_subplot(3, 3, .01, [.061, 0.03], [.065, .02]);
    [ha] = tight_subplot(3, 2, .01, [.061, 0.03], [.065, .02]);
    ha = reshape(ha', [], 3)';
    freq_s = logspace(log10(1), log10(12500), 350)';

    h4 = frf_bode_mag(Dx_ff, freq_s, ha(1,1), 'Hz', '--m');
    h5 = frf_bode_mag(Dy_ff, freq_s, ha(2,2), 'Hz', '--m');
    
    h1 = mimo_bode_mag(HH1(:,1:2), freq_s, ha, '-r');
    h2 = mimo_bode_mag(HH2(:, 1:2), freq_s, ha, '-b');
    h3 = mimo_bode_mag(Gxyz_frd(:, 1:2), freq_s, ha, '-k');
    
%     ylabel(ha(2,2), '')
    
    
    fprintf('Bandwidth HH1: %f\n', bandwidth(HH1(1,1))/2/pi)
    fprintf('Bandwidth HH2: %f\n', bandwidth(HH2(1,1))/2/pi)
  
    h1(1,1).DisplayName = 'C.L w/o F.F';
    h2(1,1).DisplayName = 'C.L w/ F.F';
    h3(1,1).DisplayName = 'O.L. Plant';

    h4.DisplayName = '$F_X$';

    h5.DisplayName = '$F_Y$';
    
    leg = legend([h1(1,1), h2(1,1), h3(1,1), h4], 'FontSize', 10);
    legy = legend(h5, 'FontSize', 12);
    set(leg, 'location', 'southwest')
    set(legy, 'location', 'southwest')
    
    xlim(ha, [1, 4000])
    ylim(ha(1:2, :), [-75, 15])
    ylim(ha(3, :), [-75, 35])
end

function H_frd = get_H()
  %
  sysx_fn = fullfile(PATHS.sysid(), 'ALL-axis_sines_info_intsamps_quickFourierCoef_1-17-2019x-drive-01.json');
  sysy_fn = fullfile(PATHS.sysid(), 'ALL-axis_sines_info_intsamps_quickFourierCoef_1-17-2019y-drive-01.json');
  sysz_fn = fullfile(PATHS.sysid(), 'ALL-axis_sines_info_intsamps_quickFourierCoef_1-17-2019z-drive-03.json');
  sysx_ = SweptSinesOnline(sysx_fn);
  sysy_ = SweptSinesOnline(sysy_fn);
  sysz_ = SweptSinesOnline(sysz_fn);
  Gux = sysx_.FRF_from_FC(1, [2,3,4]);
  Guy = sysy_.FRF_from_FC(1, [2,3,4]);
  Guz = sysz_.FRF_from_FC(1, [2,3,4]);
  
  H_frd = [Gux, Guy, Guz];
  
end

function [DD, D_ki, D_inv] = get_gz_dz(Gzz)
    % Load initial guess.
    root = fullfile(PATHS.sysid, 'z_axis_evolution', 'batch-1');
    file = 'first_res_fit-3-17-2019-1.json';
    dat = loadjson(fullfile(root, file));
    
    idx1 = find(Gzz.Frequency < 180, 1, 'last');
    idx2 = find(Gzz.Frequency > 230, 1, 'first');
    % Initial guess
    gd0 = tf(dat.Dinv_den, dat.Dinv_Num, AFM.Ts)/dat.K;
    sos_fos = SosFos(gd0);
    lg = LogCostZPK(Gzz.ResponseData(idx1:idx2), Gzz.Frequency(idx1:idx2)*2*pi, sos_fos);
    lg.solve_lsq(1)
    D = lg.sos_fos.realize();
    
    D_inv = 1/D;
    
    KI = -0.07;
    
    D_ki = zpk([], [1], KI, AFM.Ts);
    
    DD = D_ki * D_inv;
    
end

function analyze_margins(plants, sys_obsDist, K_lqr, L_dist, verbose)
  [Sens, Hyd, Hyr, Hyeta, Loop] = ss_loops_delta_dist(plants.SYS, plants.sys_recyc,...
    sys_obsDist, K_lqr, L_dist);
  
  F_clbode = figure(25);clf; hold on, grid on
  Hbode_sens = bodeplot(Sens);
  setoptions(Hbode_sens, 'FreqUnits', 'Hz')
  legend('S: linear')
  grid on, hold on;
  
  [Gm_lin, Pm_lin] = margin(Loop);
  
  fprintf('-------- MARGINS ------------------\n')
  fprintf('Linear: GM = %.2f [dB], PM = %.2f [deg]\n', 20*log10(Gm_lin), Pm_lin)
  fprintf('Bandwidth: %f [Hz]\n', bandwidth(minreal(Hyr))/2/pi);
  
  if verbose >=2  
    figure(101)
    rlocus(Loop);
    title('Klin')
  end
  
  if verbose >=1
    figure(104)
    nyquist(Loop)
    title('Klin')
    xlim([-2, 0.1])
    ylim([-1, 1])
  end
  
  
end


