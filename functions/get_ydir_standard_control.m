function ydirControl = get_ydir_standard_control(verbose)
    
    if nargin < 1
        verbose = true;
    end
    
    % ------- Load Plants -----
    sysy_fn = fullfile(PATHS.sysid(), 'ALL-axis_sines_info_intsamps_quickFourierCoef_1-17-2019y-drive-01.json');
    sysy_ = SweptSinesOnline(sysy_fn);
    Gy_frd = sysy_.FRF_from_FC(1, [3]); % idx=2 is xdir coupling
    
    [Gy, gy_eigen] = fit_gy_prox(Gy_frd, true);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                                                                         %
    % Design control/estimator gains. This gains are what we actually         %
    % use in both the simulation and experimentant system.                    %
    %                                                                         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    Ki = 0.025;
    Dki = zpk([0], [1], Ki, AFM.Ts);
    Dy_ff = make_dy_ff();
    
    Dy_notch = zpk([0.978+1i*0.164, 0.978-1i*0.164], [0.7, 0.7], 1, AFM.Ts);
    Dy_notch = Dy_notch/dcgain(Dy_notch);
    % Dy = Dy_notch;
    
    [pw, zw, pl_zet, zr_zet] = pole_zero_freq(Gy);
    [~, idx_plsort] = sort(pw);
    [~, idx_zrsort] = sort(zw);
    pl = pole(Gy);
    zr = zero(Gy);
    
    z1 = pl(idx_plsort(1:2));
    p1 = zr(idx_zrsort(1:2));
    z2 = pl(idx_plsort(3:4));
    p2 = zr(idx_zrsort(3:4));
    
    
    g1 = balreal(zpk(z1, p1, 1, Gy.Ts));
    g2 = balreal(zpk(z2, p2, 1, Gy.Ts));
    
    Dy = balreal(Dy_notch) * (g1/dcgain(g1)) * (g2 / dcgain(g2));
    
    [Sens, Hy_rprime, Hyr, Loop] = tf_loops(Gy, Dki * Dy, 1, Dy_ff);
    

    
    ydirControl = AfmSisoController(Gy,...
        'Ki', Ki,...
        'Dki', Dki,...
        'D', Dy,...
        'M', Dy_ff);
    
    if verbose
        figure(1),
        step(ydirControl.Hyr, ydirControl.Hy_rprime)
        title('y-dir step response (loop shaped)')
    end    
    
    [gm, pm] = margin(ydirControl.Loop);
    fprintf('----------- y-direction-------------\n')
    fprintf('(no ff) bandwidth: %.2f Hz\n', bandwidth(ydirControl.Hy_rprime)/2/pi);
    fprintf('(with ff) bandwidth: %.2f Hz\n', bandwidth(ydirControl.Hyr)/2/pi);
    fprintf('Phase Margin: %.2f [deg]\n', pm);
    fprintf('Gain Margin: %.2f [deg]\n', 20*log10(gm));
end

function [Sens, Hyr, Hyr_prime, Loop] = tf_loops(G, D1, D2, M)
%            rprime    
%              ^
%              |
%    R-->[ M ]--->O---->[ D1 ]-->[ G ]--+---> y
%                 |                     |           
%                 +-----[ D2 ]<---------+
%
%                G*D1
% y =      --------------- * M * R
%            I + D2*G*G1

   Loop =  G*D1*D2;
   Hyr_prime = minreal(balreal(G)*D1 / (1 + Loop));
   Sens = 1 / (1 + Loop);
   Hyr = Hyr_prime * M;
    
end

function Dy_ff = make_dy_ff()
    wz1 = 214 * 2 * pi;
    wz2 = 505 * 2 * pi;
    wz3 = 382 * 2 * pi;
    zz1 = 0.03;
    zz2 = 0.04;
    
    g = zpk(tf([1, 2*zz1*wz1, wz1^2], conv([1, 300*2*pi],  [1, 300*2*pi])));
    g2 = zpk(tf([1, 2*zz2*wz2, wz2^2], conv([1, 450*2*pi],  [1, 450*2*pi])));
    g3 = zpk(tf([1, 2*zz2*wz3, wz3^2], conv([1, 450*2*pi],  [1, 450*2*pi])));
    
    g = c2d(g, AFM.Ts, 'matched') * c2d(g2, AFM.Ts, 'matched') * c2d(g3, AFM.Ts, 'matched');
    Dy_ff = g * (1/dcgain(g));
end

function [gy_lg, gy_eigen] = fit_gy_prox(Gy_frd, verbose)

    ssopts = frf2ss_opts('Ts', AFM.Ts, 'Nd', 9);
    gy_ = frf2ss(squeeze(Gy_frd.Response(1:end-1)), Gy_frd.Frequency(1:end-1)*2*pi, 9, ssopts);
    gy_eigen = gy_.realize(13);
    
    sos_fos_y = SosFos(gy_eigen, 'iodelay', 9);
    enforce_stable = true;
    enforce_nmp = true;
    enforce_pos_k = true;
    lg = LogCostZPK(Gy_frd.Response(1:end-1), Gy_frd.Frequency(1:end-1)*2*pi,...
        sos_fos_y, enforce_stable, enforce_nmp, enforce_pos_k);
    lg.solve_lsq(3);
    gy_lg = lg.sos_fos.realize();
    nd = round(lg.sos_fos.theta(end), 0);
    gy_lg.IODelay = nd;
    
    if verbose
        figure(14); clf
        bode(Gy_frd, gy_eigen, gy_lg, '--')
        legend('FRF', 'Eigenspace', 'Log-fit')
    end
end