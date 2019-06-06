
function xdirControl = get_xdir_loop_shaped_control(verbose)
    if nargin < 1
        verbose = true;
    end

    addpath(fullfile(getCsRoot(), 'matlab-code', 'functions'));
    
    % ------- Load Plants -----
    [plants, frf_data] = CanonPlants.plants_ns14(9, '5micron');
    
    % It is important to id the mode with a very small (0.1 volt) amplitude.
    newest_xfrf_path = fullfile(PATHS.sysid, 'x-axis_sines_info_first_resFourierCoef_6-1-2019-01.json');
    sysx_new = SweptSinesOnline(newest_xfrf_path);
    
    
    [z, p, k] = zpkdata(plants.Gvib, 'v');
    p1 = p(end-1:end);
    z1 = z(end-3:end-2);
    p(end-1:end) = [];
    z(end-3:end-2) = [];
    
    g_tmp = zpk(z, p, k, plants.Gvib.Ts)*plants.gdrift;
    g_tmp = g_tmp*dcgain(plants.Gvib*plants.gdrift)/dcgain(g_tmp);
    
    Gx_bend_ = sysx_new.FRF_from_FC(1, [2]);
    Gx_bend_frd = Gx_bend_/g_tmp;
    
    % Initial guess.
    gx_bend0 = zpk(z1, p1, 1, plants.Gvib.Ts);
    
    sos_fos_xbend = SosFos(gx_bend0, 'iodelay', 9);
    lg = LogCostZPK(squeeze(Gx_bend_frd.Response), Gx_bend_frd.Frequency*2*pi, sos_fos_xbend);
    lg.solve_lsq(1);
    gx_bend_lg = lg.sos_fos.realize();
    
    if verbose
        figure(100)
        opts = bodeoptions();
        opts.FreqUnits = 'Hz';
        bodeplot(Gx_bend_frd, opts)
        grid
        hold on
        bodeplot(gx_bend_lg, Gx_bend_frd.Frequency*2*pi);
    end

    Dinv = 1 / gx_bend_lg;
    Dinv = Dinv / dcgain(Dinv);
    
    D_notch= zpk([0.979 + 1i*0.163, 0.979  - 1i*0.163], [0.8, 0.8], 1, AFM.Ts);
    D_notch = D_notch/dcgain(D_notch);
    

   
    Ki = 0.035;
    Dki = zpk(0, 1, Ki, AFM.Ts);
    % THe gain crossover is nominally at 615 rad/s. The real pole at 400 hz and the
    % integrator collide and exit unit circle.
    wz = 250 * 2*pi;
    wp = 2.2*wz;

    kk = wp/wz;
    D_lead_ = kk*tf([1, wz], [1, wp]); 
    D_lead = c2d(D_lead_, AFM.Ts, 'matched');
    
    Dx = balreal(Dinv)*balreal(D_notch)*balreal(D_lead);
    
    Dx_ff = make_dx_ff();
    
    xdirControl = AfmSisoController(plants.Gvib,...
                                   'Ki', Ki,...
                                   'Dki', Dki,...
                                   'D', Dx,...
                                   'M', make_dx_ff());
    if verbose
        figure(2)
        step(xdirControl.Hy_rprime, xdirControl.Hyr)
        title('x-dir step response (loop shaped)')
        
        figure(3); clf
        opts = bodeoptions();
        opts.FreqUnits = 'Hz';
        h = bodeplot(xdirControl.Loop, xdirControl.Hy_rprime, D_lead, opts);
        axs = h.getaxes();
        set(axs(1,1,2), 'YLim', [-360, 45])
        set(axs(1,1,1), 'YLim', [-80, 20])
        grid on
        
        figure(4)
        rlocus(xdirControl.Loop)
        xlim([0.7, 1.1])
        ylim([-0.6, 0.6])
    end
    
    [gm, pm] = margin(xdirControl.Loop);
    fprintf('----------- x-direction-------------\n')
    fprintf('(no ff) bandwidth: %.2f Hz\n', bandwidth(xdirControl.Hy_rprime)/2/pi);
    fprintf('(with ff) bandwidth: %.2f Hz\n', bandwidth(xdirControl.Hyr)/2/pi);
    fprintf('Phase Margin: %.2f [deg]\n', pm);
    fprintf('Gain Margin: %.2f [deg]\n', 20*log10(gm));
end


function Dx_ff = make_dx_ff()
    wz1 = 214 * 2 * pi;
    wz2 = 505 * 2 * pi;
    wz3 = 372 * 2 * pi;
    zz = 0.01;
    
    
    g = zpk(tf([1, 2*zz*wz1, wz1^2], conv([1, 350*2*pi],  [1, 350*2*pi])));
    g2 = zpk(tf([1, 2*zz*wz2, wz2^2], conv([1, 450*2*pi],  [1, 450*2*pi])));
    g3 = zpk(tf([1, 2*zz*wz3, wz3^2], conv([1, 450*2*pi],  [1, 450*2*pi])));
    w_lpf = 2 * pi * 500;
    g4 = zpk([], [-w_lpf], w_lpf);
    
    g = c2d(g, AFM.Ts, 'matched') * c2d(g2, AFM.Ts, 'matched')...
        * c2d(g3, AFM.Ts, 'matched') * c2d(g4, AFM.Ts, 'matched');
    Dx_ff = g * (1/dcgain(g));
end
