addpath functions
clear
% close all
root = fullfile(PATHS.sysid, 'z_axis_evolution', 'batch-1');
file = 'first_res_fit-3-17-2019-1.json';

dat = loadjson(fullfile(root, file));

gz_k = dat.SOS_frf_data.resp_real + dat.SOS_frf_data.resp_imag*1i;

ss_fname = fullfile(root, 'ALL-axis_sines_info_intsamps_quickFourierCoef_3-17-2019-01.json');
ss_data = SweptSinesOnline(ss_fname);

G_frf = frd(ss_data.FC_s(:,2)./ss_data.FC_s(:,1), ss_data.freq_s, AFM.Ts, 'FrequencyUnit', 'Hz');


KI = -0.06;
D_inv = dat.K*frd(tf(dat.Dinv_Num, dat.Dinv_den, AFM.Ts),  ss_data.freq_s,  'FrequencyUnit', 'Hz');

KI_norm = KI;

D_ki = zpk([], [1], KI, AFM.Ts);
D_ki_old = zpk([], [1], 0.005, AFM.Ts);

DD = D_ki * D_inv;

Loop = DD * G_frf;
Loop_noinv = D_ki * G_frf;

Huz_d = feedback(D_ki, D_inv*G_frf);
H2 = D_ki/(1+Loop);

Hyr = feedback(D_ki*D_inv*G_frf, 1);
Hyr_old = feedback(D_ki_old*G_frf, 1);
Hyr_noinv = feedback(D_ki*dat.K*G_frf, 1);


F2 = mkfig(2, 3.48, 2); clf
[ha, pos] = tight_subplot(1, 1, [.02, .01 ], [.15, 0.08], [.1, .02]);


h1 = frf_bode_mag(Loop, ss_data.freq_s, ha, 'Hz', ':b');
h2 = frf_bode_mag(G_frf, ss_data.freq_s, ha, 'Hz', '-k'); %, 'LineWidth', 1.1);
h4 = frf_bode_mag(Hyr, ss_data.freq_s, ha, 'Hz', '-.r'); % 'LineWidth', 1.1);

ylim(ha(1), [-40, 20]);

h1.DisplayName = 'Loop gain';
h2.DisplayName = '$G_{Z_d,u_Z}$';
h4.DisplayName = '$H_{Z_d, r_{Z}}$ ';
leg = legend([h1, h2, h4]);
set(ha, 'FontSize', 8)
set(leg, 'Location', 'southwest', 'FontSize', 8);
ylabel(ha, 'Output $Z$ (Mag [dB])')
title(ha, 'Input $Z$')
%%
% [pm, gm] = margin(Loop)
save_fig(F2, fullfile(PATHS.tmech_fig(), 'z_control_design'), false);




