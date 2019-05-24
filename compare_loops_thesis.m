
clear
clc
Ts = 40e-6;

LPF = zpk([], [0.7, 0.7], 1, Ts);
LPF = LPF/dcgain(LPF);
models = load(fullfile(PATHS.sysid, 'x-axis_sines_infoFourierCoef_10-21-2018-03.mat'));
G = -models.modelFit.G_zdir;
Ki_s = [0.005, 0.01, 0.05];


root = '/media/labserver/afm-cs/z-scope';

fnames = {'data-out_KI_p005.csv',...
  'data-out_KI_p05.csv'};

idx = Inf;


% F1 = mkfig(1, 5, 2.5, true); clf
% % ax1 = gca();
% ax1 =axes('Position', [0.13 0.15 0.775 0.815]);
% ax1.Color = fig_color;
% hold on, grid on,
% 
% % figure(2); clf
% F2 = mkfig(2, 5, 2.5, true); clf
% % ax2 = gca();
% ax2 =axes('Position', [0.13 0.15 0.775 0.815]);
% ax2.Color = fig_color;
% hold on, grid on,
% 
% h_ze = gobjects(length(fnames),1);
% h_uz = gobjects(length(fnames),1);
% clrs = {'b', 'r', 'k'};

% KI_s = [0.005, 0.05]; % for legend.
% for k=1:length(fnames)
%   Dz = zpk([0], [1], Ki_s(k), Ts);
%   
%   H_ud = minreal((1 + G*Dz)/Dz)*LPF;
%   
%   dat_k = csvread(fullfile(root, fnames{k}));
%   
%   ze_k = dat_k(:,1);
%   uz_k = dat_k(:,2);
%   
%   idx = min(idx, length(ze_k));
%     
%   t_k = (0:length(ze_k)-1)'*Ts;
% 
%   idx1 = 480;
%   idx2 = 56532;
%   h_ze_k = plot(ax1, t_k(idx1: idx2), ze_k(idx1:idx2), 'Color', clrs{k});
%   h_ze_k.DisplayName = sprintf('$K_I=%.3f$', KI_s(k));
%   uu = detrend(uz_k);
%   h_uz_k = plot(ax2, t_k(idx1:idx2), uu(idx1:idx2), 'Color', clrs{k});
%   h_uz_k.DisplayName = sprintf('$K_I=%.3f$', KI_s(k));
%   
%   h_ze(k) = h_ze_k;
%   h_uz(k) = h_uz_k;
% 
% end
% 
% linkaxes([ax1, ax2], 'x')
% xlim(ax1, [0, 1.5])


% title(ax1, 'z-error', 'FontSize', 14)
% xlabel(ax1, 'time [s]')
% xlabel(ax2, 'time [s]')
% leg1 = legend(ax1, h_ze);
% set(leg1, 'FontSize', 14)
% 
% % title(ax2, 'Uz', 'FontSize', 14)
% leg2 = legend(ax2, h_uz);
% set(leg2, 'FontSize', 14);
% linkaxes([ax1, ax2], 'x')
% ylabel(ax1, '$e_z$', 'FontSize', 16);
% ylabel(ax2, '$u_z$', 'FontSize', 16);
% 
% print(F2, '-dpdf', 'figures/uz_noise_example_Kis.pdf')
% print(F1, '-dpdf', 'figures/ze_noise_example_Kis.pdf')


freqs = logspace(log10(.1), log10(12.5e3), 1000)';

rand_fname = fullfile(PATHS.sysid, 'rand_noise_zaxis_10-30-2018_01.mat');
models = load(rand_fname);

Ki_s = [0.001, 0.01, 0.05];

G = minreal(ss(models.modelFit.G_zdir)/models.modelFit.gdrift)*2;

KK = abs(dcgain(G));

Dinv = models.modelFit.Dinv;
p = pole(G);

D2 = zpk(p(end-1:end), [0,0], 1, Ts);
D2 = D2/dcgain(D2);

F5 = mkfig(5, 5, 4, true); clf
set(gca(), 'Color', fig_color)
% Dinv = 1;
F6 = mkfig(6, 5, 4, true); clf
set(gca(), 'Color', fig_color)

Hu = gobjects(length(Ki_s),1);
clrs = {'b', 'r', 'k'};

for k=1:length(Ki_s)
  
  Dz = zpk([0], [1], -Ki_s(k)/KK, Ts);
  
  H_de = -minreal(1/(1+Dz*Dinv*G));
  %   H_du = -minreal(Dz/(1+Dz*Dinv*G));
  H_du = feedback(Dz, Dinv*G);
  frf_bode_mag(H_de, freqs, F5, 'Hz', clrs{k}, 'LineWidth', 1.5);
  hold on, grid on
  
  Hu(k) = frf_bode_mag(H_du, freqs, F6, 'Hz', clrs{k}, 'LineWidth', 1.5);
  Hu(k).DisplayName = sprintf('$K_I = %.3f$', Ki_s(k));
end

figure(F5);
title('Disturbance to error')
ylim([-50, 5])
lower = [-75, -75];
upper = [5,5];
h = ciplot(lower, upper, [7, 30],'g');
alpha(h, '.25')

figure(F6)
ylim([-50, 5])
% leg = legend(Hu);
lower = [-75, -75];
upper = [5,5];
h = ciplot(lower, upper, [7, 30],'g');
alpha(h, '.15')


lower = [-75, -75];
upper = [5,5];
h2 = ciplot(lower, upper, [.1, 1],'b');
alpha(h2, '.15')
%
st1 = text(.28, -48, 'Sample Frequency Content', 'rot', 90, 'FontSize', 15)
st2 = text(15, -48, 'Building Vibration', 'rot', 90, 'FontSize', 15)

st3 = text(102, -30, '$K_I=0.001$', 'rot', -59, 'FontSize', 15, 'interpreter', 'latex')
st4 = text(454, -21, '$K_I=0.01$', 'rot', -59, 'FontSize', 15, 'interpreter', 'latex')
st5 = text(909, -12, '$K_I=0.05$', 'rot', -62, 'FontSize', 15, 'interpreter', 'latex')

print(F6, '-dpdf', 'figures/Hdu_1.pdf')

Hu(2).Visible = 'off';
Hu(3).Visible = 'off';
st4.Visible = 'off';
st5.Visible = 'off';

save_fig(F6, fullfile(PATHS.defense_fig, 'Hdu_01'), true)




Hu(3).Visible = 'on';
lower = [-75, -75];
upper = [5,5];
h2 = ciplot(lower, upper, [.1, 100],'b');
alpha(h2, '.15')

h.Visible = 'off';
st2.Visible = 'off';
st5.Visible = 'on';
%%
save_fig(F6, fullfile(PATHS.defense_fig, 'Hdu_02'), true)



%%
% h.DisplayName = 'Building Vibration Band';
set(leg, 'FontSize', 14)
set(leg, 'Position', [0.3743 0.1111 0.5298 0.2477])
title('Disturbance to control')


%%
print(F5, '-dpdf', 'figures/Hde.pdf')
print(F6, '-dpdf', 'figures/Hdu.pdf')


