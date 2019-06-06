function xdirControl = get_xdir_standard_control(type, gam)
addpath(fullfile(getCsRoot(), 'matlab-code', 'functions'));
addpath(fullfile(getCsRoot(), 'matlab-code', 'functions', 'state_space_x'));

md = 1;



% ------- Load Plants -----
[plants, frf_data] = CanonPlants.plants_ns14(9);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Design control/estimator gains. This gains are what we actually         %
% use in both the simulation and experimentant system.                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 3). LQR generation gain.
% -------------------------------------------------------------------------
if strcmp(type, 'const-sig')
  cmplx_rad = 0.9;
  [Q1, R0, S1] = build_control_constsigma(plants.sys_recyc, cmplx_rad);
    gam_rob = 46.4;
elseif strcmp(type, 'choose-zeta')
  can_cntrl = CanonCntrlParamsChoozeZeta();
  [Q1, R0, S1] = build_control_choosezeta(plants.sys_recyc, can_cntrl);
  gam_rob = 25;
else
    error('unrecognized control type.');
end

if exist('gam', 'var')
    gam_rob = gam;
end


% ------------------------- Observer Gain ---------------------------------
can_obs_params = CanonObsParams_01();
can_obs_params.beta = 50;
[sys_obsDist, L_dist] = build_obs(plants.SYS, can_obs_params);

% 2). Design FeedForward gains.
K_lqr = dlqr(plants.sys_recyc.a, plants.sys_recyc.b, Q1, R0+gam_rob, S1);

if 0
  verbose = 0;
  analyze_margins(plants, sys_obsDist, K_lqr, L_dist, verbose);
end

[Sens, Hyd, Hyr, Hyeta, Loop, D1, D2] = ss_loops_delta_dist(plants.SYS, plants.sys_recyc,...
    sys_obsDist, K_lqr, L_dist);



  Dx_ff = make_dx_ff();
  
  xdirControl = AfmSisoSSController(plants.SYS, plants.sys_recyc,...
                                    sys_obsDist, L_dist, K_lqr, 'M', Dx_ff);


end


function Dx_ff = make_dx_ff()
    wz1 = 214 * 2 * pi;
    wz2 = 505 * 2 * pi;
    zz = 0.01;
    
    
    g = zpk(tf([1, 2*zz*wz1, wz1^2], conv([1, 350*2*pi],  [1, 350*2*pi])));
    g2 = zpk(tf([1, 2*zz*wz2, wz2^2], conv([1, 450*2*pi],  [1, 450*2*pi])));
    w_lpf = 2 * pi * 500;
    g3 = zpk([], [-w_lpf], w_lpf);
    
    g = c2d(g, AFM.Ts, 'matched') * c2d(g2, AFM.Ts, 'matched') * c2d(g3, AFM.Ts, 'matched');
    Dx_ff = g * (1/dcgain(g));
end
