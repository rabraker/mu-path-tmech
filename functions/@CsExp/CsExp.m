classdef CsExp < handle
  properties
    x;
    x_positive;
    y;
    y_positive;
    uz;
    ze;
    z_friction;
    t;
    Ts;
    met_ind;
    idx_state_s;
    npix;
    width;
    channel_map;
   
    Gz;
    meta_exp;
    meta_in;
    state_times;
    time_total;
    feature_height
    gg;
    cs_paths;

    pix_mat_raw_uz;
    pix_mat_uz;
    pix_mask_uz;

    pix_mat_raw_ze;
    pix_mask_ze;
    pix_mat_ze;
    
    UserData = [];
  end

  methods
    function self = CsExp(cs_paths, varargin)
    % channel_map, Ts, feature_height_nm, gg)
    % Obtain cs_paths (which is a struct) from, e.g.,
    % cs_exp_paths(data_root, data_name)
      default_chan_map = ChannelMap([1:6]);
      optP = inputParser();
      optP.addParameter('reload_raw', false, @(s)islogical(s));
      optP.addParameter('channel_map', default_chan_map);
      optP.addParameter('feature_height', Inf);
      optP.addParameter('gg', []);
      optP.addParameter('Ts', AFM.Ts);
      optP.addParameter('load_full', false, @(s)islogical(s));
      optP.parse(varargin{:});
      opts = optP.Results;
      % Check to see if an already processed mat file exists.
      mat_exists = exist(cs_paths.data_path_mat, 'file');

      if ~opts.reload_raw && mat_exists == 2
        fprintf('loading from mat file...\n')
        self = self.load_mat(cs_paths.data_path_mat, opts.load_full);
        fprintf('\bdone\n')
      else
        fprintf('loading from raw data...\n')
        self = self.load_raw_data(cs_paths, opts);
        fprintf('\bdone\n')
      end

      self.time_total = length(self.x)*self.Ts;
    end
    function rate = equiv_raster_rate(self)
       rate = self.meta_in.tip_velocity * 0.5 / self.meta_in.width; 
    end
    function N = N_scans(self)
        N = length(self.idx_state_s.scan);
    end
    
    function [tmu_overhead_avg, tscan_avg] = get_mean_mu_overhead(self)
    % [tscan_avg, tmu_overhead_avg] = get_mean_state_times(self)
       N = self.N_scans();
       times = self.get_state_times();
       % compute the average move time:
       tscan_avg = times.scan/N;
       tmu_overhead_avg = (times.move + times.tdown + times.tsettle + times.tup)/N;
    end
    
    function frac = sub_sample_frac(self)
      N = self.npix^2;
      n_samps = sum(self.pix_mask_uz(:));
      frac = n_samps/N;
    end
    
    function [t_cycle_est, t_connect] = estimate_mpt_connect_savings(self)
      N_up = length(self.idx_state_s.tup);
      N_scan = length(self.idx_state_s.scan);
      N_move = length(self.idx_state_s.move);
      N_down = length(self.idx_state_s.tdown);
      N_settle = length(self.idx_state_s.tsettle);
      N_connect = length(self.idx_state_s.connect);
      
      times = self.get_state_times();
      % compute the average move time:
      % scan_avg = times.scan/N_scan;
      % move_avg = times.move/N_move;
      % tdown_avg = times.tdown/N_down;
      % tsettle_avg = times.tdown/N_settle;
      
      % Average time to do CS cycle, excluding scanning:
      t_cycle = (times.move + times.tdown + times.tsettle + times.tup);
      t_cycle_avg = t_cycle/(N_move);
      
      % Estimated time to do a cycle, rather than than connect:
      t_cycle_est = N_connect*t_cycle_avg;
      t_connect = times.connect;      
    end
    
    function times = get_state_times(self)
      flds = fields(self.idx_state_s)';
      total = 0;
      for fld_ = flds
        fld = fld_{1};
        N_state = 0;
        for k=1:length(self.idx_state_s.(fld))
          N_state = N_state + length(self.idx_state_s.(fld){k});
        end
        times.(fld) = N_state * self.Ts;
        total = total + times.(fld);
      end
      times.total = total;
    end
    
    function damage = damage_metric(self)
    % Compute a damage metric based on the deflection signals positivity.
    % This is computed as the power of the positive values of the negative
    % error signal. The motivation is that, for a given setpoint, we do not
    % care, from a damage perspective, if the error dips below the setpoint
    % (though that will affect image quality), because this corresponds to the
    % tip parachiting off a ledge. Rather from a damage perspective, what we
    % care about is events where (ze - ref) signal becomes positive.
    
      ref = self.meta_exp.z_axis_params.setpoint_scan;
      % rather than subtracting mean, subtract the reference value.
      err = self.ze - ref;  % shift to zero.
      
      % This should *should* drop most everything not part of scan/settling.
      % If it doesnt, well, our control isnt doing well and it should part of
      % the metric.
      err_pos = err(err>0);
      damage = sum(err_pos.^2)/length(err_pos)/self.Ts;
    end

    function quality = quality_metric(self)
    % Compute a quality metric based on the deflection signal's power.

      ref = self.meta_exp.z_axis_params.setpoint_scan;
      % rather than subtracting mean, subtract the reference value.
      ze_cs_scan = [];
      for k=1:length(self.idx_state_s.scan)
        ze_cs_scan = [ze_cs_scan; self.ze(self.idx_state_s.scan{k})]; %#ok<AGROW>
      end
      err = ze_cs_scan - ref;
      quality = sum(abs(err).^2)/length(err)/self.Ts;
    end
    
    function save(self, force_save)
    % Serialize to a .mat file to the location contained
    % in raster_paths.data_path_mat.
    
      if nargin <2
        force_save = false;
      end

      % Remove empty fields, so we don't overwrite data we potentially didn't
      % load with empty.
      if ~force_save && ~self.raw_data_loaded()
        warning(['Not saving data because the raw data is not loaded',...
          'and force_save flag is false. Saving as an append operation',...
          'is very time consuming so is disabled by default.'])
        return
      end
      % Go ahead and save it.
      warning('off', 'MATLAB:structOnObject');
      self_struct = struct(self);
      if ~self.raw_data_loaded() % weve been instructed to save anyway
        for fld=fieldnames(self_struct)'
          
          if isempty(self_struct.(fld{1}))
            self_struct = rmfield(self_struct, fld{1});
          end
        end
        save(self.cs_paths.data_path_mat, '-struct', 'self_struct', '-append');
      else
        save(self.cs_paths.data_path_mat, '-struct', 'self_struct');
      end
      warning('on', 'MATLAB:structOnObject');
    end    
    
    function self = load_mat(self, data_path_mat, load_full)
    % Load ourself from the location contained in data_path_mat. If
    % load_full=false (the default), then the original time-series
    % data, x,y,uz,ze, x_positive,y_positive will not be loaded.
    % This is to speed things up when we just want to work with the
    % already processed images.

    % I tried use the matfile. This saves us zero time. I cant
    % figure out how to do this without creating an extra structure
    % and passing that into or out of self. This is wastful...
      self_mat = matfile(data_path_mat);
      to_load_list = properties(self);
      if ~load_full
        no_loads = {'x', 'y', 'x_positive', 'y_positive', 'uz', 'ze'};
        to_load_list = setdiff(to_load_list, no_loads);
      end

      data = load(data_path_mat, to_load_list{:});

      for prop = to_load_list'
        prop = prop{1};
        self.(prop) = data.(prop);
      end

    end

    function h = plot_traj_in_time_interval_by_state(self, t0, t1, state_name,...
        traj_name, ax, t_offset, varargin)
      
      idx_s = self.get_idx_by_state_in_time_range(state_name, t0, t1);
      
      for k=1:length(idx_s)
        h = plot(ax, self.t(idx_s{k})-t_offset, self.(traj_name)(idx_s{k}), varargin{:});
        hold(ax, 'on')
      end
    end
    
    function h = plot_traj_from_csidx_by_state(self, idx_start, idx_end, state_name,...
        traj_name, ax, t_offset, varargin)
      
      for k=idx_start:idx_end
        idx = self.idx_state_s.(state_name){k};
        h = plot(ax, self.t(idx)-t_offset, self.(traj_name)(idx), varargin{:});
        hold(ax, 'on')
      end
    end    
    % ************************************************************************ %
    % ******************* methods in other files ***************************** %
    
    [idx_cell, N_min]= get_idx_by_state_in_time_range(self,...
        state_name, t_start, t_end)
      
    % Defined in psd_from_intervals.m
    [sig_psd, freqs, k] = psd_from_intervals(self, signal, state, ...
                                             starts, ends)
    % Defined in process_cs_data.m
    [pix_mask] = process_cs_data(self, use_ze, register_uzk)
      
    function [CS_idx, start_idx, end_idx] = find_cycle_idx(self, time)
    % find the CS-cycle index corresponding to time. A single cycle is defined
    % as:
    % xymove --> tip-down--> tip-settle-->scan-->tip-up.
    %
    % Very niave search through everything. Would be faster to bisect.

      for CS_idx=1:length(self.idx_state_s.tup)
        % first state is xy-move, last is tip-up. Thus, it is sufficient to check if
        % time is between the first time of move and last time of tup.
        idx_mov = self.idx_state_s.move{CS_idx};
        idx_tup = self.idx_state_s.tup{CS_idx};
        t_mov = self.t(idx_mov);
        t_up = self.t(idx_tup);
        if t_mov(1) <= time && time <= t_up(end)
          start_idx = idx_mov(1);
          end_idx = idx_tup(end);
          return
        end
      end
      % if we get here, we didnt find it.
      CS_idx = [];
      start_idx = [];
      end_idx = [];
      warning('time not found\n');
    end

    function xy_positive(self)
      % ---------- We really need something more sophisticated here. 
      %
      %
      % Move x-y data to the positive othant. Need to do this based on
      % MEASUREMENT DATA, because with pre-scan, we may be purposefully outside
      % e.g., to the left of the y-axis. Moving everything to positive orthant
      % will leave us with a strip of "unmeasured" data in the final image.
      xmin_meas = [];
      ymin_meas = []; % n.b. min([ [], 9]) = 9.
      for k=1:length(self.idx_state_s.scan)
        idx_k = self.idx_state_s.scan{k};
        xmin_k = min(self.x(idx_k));
        ymin_k = min(self.y(idx_k));
        xmin_meas = min([xmin_meas, xmin_k]); % brackets necessary
        ymin_meas = min([ymin_meas, ymin_k]); % brackets necessary
      end

      self.x_positive = self.x - xmin_meas;
      self.y_positive = self.y - ymin_meas;
    end

    function print_state_times(self)
      if self.state_times(1) == 0
        idx_shift = 1;
      else
        idx_shift = 0;
      end
      s_times = self.get_state_times();
      tmove = s_times.move;
      tlower = s_times.tdown;
      tsettle = s_times.tsettle;
      tscan = s_times.scan;
      tup = s_times.tup;
      t_con = s_times.connect;
      total = s_times.total;
      fprintf(['Total times\n--------\n',...
               'move   |  lower  |  settle  | scan   |connect |   up  |  total   |\n']);
      fprintf('%.3f  | %.3f   | %.3f    | %.3f | %.3f | %.3f |  %.3f   |\n', tmove,...
        tlower, tsettle, tscan, t_con, tup, total);
    end

    function [x_k, y_k, uz_k, ze_k, t_k] = get_state_cycle_k(self, k, state)
      idx_k = self.idx_state_s.(state){k};
      x_k = self.x(idx_k);
      y_k = self.y(idx_k);
      uz_k = self.uz(idx_k);
      ze_k = self.ze(idx_k);
      t_k = idx_k*self.Ts;

    end

    function [x_k, y_k, uz_k, ze_k, t_k] = get_tup_k(self, k)
      [x_k, y_k, uz_k, ze_k, t_k] =self.get_state_cycle_k(k, 'tup');
    end
    function [x_k, y_k, uz_k, ze_k, t_k] = get_settle_k(self, k)
      [x_k, y_k, uz_k, ze_k, t_k] =self.get_state_cycle_k(k, 'tsettle');
    end

    function [x_k, y_k, uz_k, ze_k, t_k] = get_down_k(self, k)
      [x_k, y_k, uz_k, ze_k, t_k] =self.get_state_cycle_k(k, 'tdown');
    end

    function [x_k, y_k, uz_k, ze_k, t_k] = get_scan_k(self, k, positive_xy)
      if exist('positive_xy', 'var') && positive_xy
        idx_k = self.idx_state_s.scan{k};
        x_k = self.x_positive(idx_k);
        y_k = self.y_positive(idx_k);
        [~, ~, uz_k, ze_k, t_k] =self.get_state_cycle_k(k, 'scan');
      else
        [x_k, y_k, uz_k, ze_k, t_k] =self.get_state_cycle_k(k, 'scan');
      end
    end

    function [U_scan, U_z, U_orig] = dynamic_detrend(self, idx)
        [~, ~, U_scan] = self.get_scan_k(idx);
        [~, ~, U_down] = self.get_down_k(idx);
        [~, ~, U_settle] = self.get_settle_k(idx);
        U_orig = [U_down(:); U_settle(:); U_scan(:)];
        scan_start_idx = length(U_down) + length(U_settle) + 1;
        t_k = (0:length(U_orig)-1)'*self.Ts;
        U_z = lsim(self.Gz, U_orig(:)-U_orig(1), t_k);

        U_scan = U_z(scan_start_idx:end);
    end
    function [U_scan, U_z, U_orig] = dynamic_detrend_ze(self, idx)
        [~, ~, U_scan] = self.get_scan_k(idx);
        [~, ~, U_down] = self.get_down_k(idx);
        [~, ~, U_settle] = self.get_settle_k(idx);
        U_orig = [U_down(:); U_settle(:); U_scan(:)];

        scan_start_idx = length(U_down) + length(U_settle) + 1;

        t_k = (0:length(U_orig)-1)'*self.Ts;
        U_z = lsim(self.Gz, U_orig(:)-U_orig(1), t_k);
        %U_z = U_orig;
        U_scan = U_z(scan_start_idx:end);
    end

    function [y_idx, x_idx, U_k] = mu_data2pix_xy(self, X_raw, Y_raw, U_ks)
%       y_idx = min(self.npix-1, floor(mean(Y_raw))) + 1;

      x_spread = max(X_raw) - min(X_raw);
      y_spread = max(Y_raw) - min(Y_raw);
      xpix_start = floor(min(X_raw));
      ypix_start = floor(min(Y_raw));
      
      npix_x_path_k = ceil(x_spread); % number of bins for this path.
      npix_y_path_k = ceil(y_spread); % number of bins for this path.

      % Now, define a set of bins for the x-direction. Each bin will have a
      % different number of data points.
      % If we include 0, then there are three points for two pixel bins.
      xbins = linspace(min(X_raw), max(X_raw), npix_x_path_k+1);
      ybins = linspace(min(Y_raw), max(Y_raw), npix_y_path_k+1);
      U_k = [];
      x_idx = [];
      y_idx = [];
      
%       for kk=1:npix_y_path_k
%         ind_y = find(Y_raw >= ybins(kk) & Y_raw < ybins(kk+1));
%         if ypix_start+kk > self.npix || isempty(ind_y)
%             continue
%         end
        for jj = 1:npix_x_path_k
          if xpix_start+jj > self.npix
            continue
          end
          % Get the indeces correspondiing to the current x-data bin.
          ind_x = find(X_raw >= xbins(jj) & X_raw < xbins(jj+1));
          %idx_kj = intersect(ind_y, ind_x);
          if ~isempty(ind_x) % Avoid errors if it is empty.
            % Slice out the corresponding height data, and average it
            % together since this is a single pixel.
            u_pix_jj = mean(U_ks(ind_x));
            
            % Collect the uz data into a string of pixels
            U_k = [U_k; u_pix_jj];
            x_idx = [x_idx, xpix_start+jj];
            y_idx = [y_idx, max(1, int32(mean(Y_raw(ind_x)))) ];
          end
        end
      %end

    end % bin_data_into_pix

    function [y_idx, x_idx, U_k] = mu_data2pix(self, X_raw, Y_raw, U_ks)
      % Make the assumption that the y-data for each path is constant enough.
      % Since we start at the (0,0) corner the xplane, we'll take the floor,
      % and add 1 for 1-based indexing.
      y_idx = min(self.npix-1, floor(mean(Y_raw))) + 1;

      x_spread = max(X_raw) - min(X_raw);
      xpix_start = floor(min(X_raw));
      npix_path_k = ceil(x_spread); % number of bins for this path.

      % Now, define a set of bins for the x-direction. Each bin will have a
      % different number of data points.
      % If we include 0, then there are three points for two pixel bins.
      xbins = linspace(min(X_raw), max(X_raw), npix_path_k+1);
      U_k = [];
      x_idx = [];
      for jj = 1:npix_path_k
        if xpix_start+jj > self.npix
          continue
        end
        % Get the indeces correspondiing to the current x-data bin.
        ind_x = find(X_raw >= xbins(jj) & X_raw < xbins(jj+1));
        if ~isempty(ind_x) % Avoid errors if it is empty.
                           % Slice out the corresponding height data, and average it
                           % together since this is a single pixel.
          u_pix_jj = mean(U_ks(ind_x));

          % Collect the uz data into a string of pixels
          U_k = [U_k; u_pix_jj];
          x_idx = [x_idx, xpix_start+jj];
        end
      end
      % Register the control data to zero. We can do this because we are
      % scanning long enough that we are always guaranteed to exit a hole.
      U_k = U_k - max(U_k);

    end % bin_data_into_pix


    function solve_bp(self, recalc, use_2d, use_ze, opts)
    % Solve the Basis Pursuit problem in either 1d or 2d. If in 1D, use the mex
    % function. 
    % Options
    % -------
    % recalc : (true|false), default false. Do not optimize if self.Img_bp is
    %          non-empty and non-zero.
    % use_2d : (true|false), default false. If true, compute using 2D-dct.
      if nargin <2
        recalc = false;
      end
      if nargin <3
        use_2d = false;
      end
      if nargin < 4
          use_ze = false;
      end
      assert(islogical(use_ze));
      
      if use_ze
          error('solve_bp cannot yet solve with self.pix_mat_ze')
      end
      if ~recalc && ~isempty(self.Img_bp) && sum(self.Img_bp(:)) ~= 0
        warning(['BP solution already calculated, so skipping optimization.',...
          'Pass recalc flag to recompute']);
        return;
      end
      
      [n m] = size(self.pix_mat_raw_uz);
      
      if ~exist('opts', 'var')
        opts = l1qc_dct_opts('verbose', 2, 'l1_tol', 0); %'epsilon', 0.01);
      end
      pix_idx = find(CsTools.pixmat2vec(self.pix_mask_uz) > 0.5);
      % b, set of measurements. have to remove all the spots we didn't sample.
      b = CsTools.pixmat2vec(self.pix_mat_raw_uz);
      b = b(pix_idx);
      min_b = min(b);
      max_b = max(b);
      max_diff_b = max_b - min_b;
      b = b/max_diff_b;
      tic
      if use_2d
        % A = @(x) CsTools.Afun_dct2(x, pix_idx, n, m);
        % At = @(x) CsTools.Atfun_dct2(x, pix_idx, n, m);
        % x0 = At(b);
        % eta_vec = CsTools.l1qc_logbarrier(x0, A, At, b, opts);
        % self.Img_bp = idct2(CsTools.pixvec2mat(eta_vec, n))*max_diff_b;
        [x_est, LBRes] = l1qc_dct(n, m, b, pix_idx, opts);
        self.pix_mat_uz = CsTools.pixvec2mat(x_est*max_diff_b, n);
      else
        [x_est, LBRes] = l1qc_dct(n*m, 1, b, pix_idx, opts);
        self.pix_mat_uz = CsTools.pixvec2mat(x_est*max_diff_b, n);
      end
      
      time_bp = toc;
      
      fprintf('BP Time: %f\n', time_bp);
      
    end
function solve_nesta(self, recalc, use_2d, use_ze, opts)
    % Solve the Basis Pursuit problem in either 1d or 2d, using NESTA (my version) 
    % To set opts, use NESTA_opts.m
    % Options
    % -------
    % recalc : (true|false), default false. Do not optimize if self.Img_bp is
    %          non-empty and non-zero.
    % use_2d : (true|false), default false. If true, compute using 2D-dct.
      if nargin <2
        recalc = false;
      end
      if nargin <3
        use_2d = false;
      end
      if nargin < 4
          use_ze = false;
      end
      assert(islogical(use_ze));
      
      if use_ze
          img_raw = self.pix_mat_raw_ze;
          img = self.pix_mat_ze;
          pix_mask = self.pix_mask_ze;
      else
          img_raw = self.pix_mat_raw_uz;
          img = self.pix_mat_uz;
          pix_mask = self.pix_mask_uz;
      end
      
      
      if ~recalc && ~isempty(self.pix_mat_uz) && sum(self.pix_mat_uz(:)) ~= 0
        warning(['BP solution already calculated, so skipping optimization.',...
          ' Pass recalc flag to recompute.']);
        return;
      end
      
      [n, m] = size(img_raw);
      
      pix_idx = find(CsTools.pixmat2vec(pix_mask) > 0.5);
      
      if use_2d
          M_fun = @(x) CsTools.pixmat2vec(dct2(CsTools.pixvec2mat(x, n)));
          Mt_fun = @(x) CsTools.pixmat2vec(idct2(CsTools.pixvec2mat(x, n)));
      else
          M_fun = @(x) dct(x);
          Mt_fun = @(x) idct(x);
      end
      E_fun = @(x) CsTools.E_fun1(x, pix_idx);
      Et_fun = @(x) CsTools.Et_fun1(x, pix_idx, n, m);
      
      % b, set of measurements. have to remove all the spots we didn't sample.
      b = CsTools.pixmat2vec(img_raw);
      b = b(pix_idx);
      min_b = min(b);
      max_b = max(b);
      max_diff_b = max_b - min_b;
      b = b/max_diff_b;

      if ~exist('opts', 'var')
        opts = NESTA_opts('Verbose', 10, 'errFcn', @(x)norm(x),...
            'U', M_fun, 'Ut', Mt_fun, 'mu', 1e-5, 'sigma', 1e-2);
      end
      
      tic
      [x_est] = NESTA_mine(E_fun, Et_fun, b, opts);

      img = CsTools.pixvec2mat(x_est*max_diff_b, n);

      if use_ze
          self.pix_mat_ze = img;
      else
          self.pix_mat_uz = img;
      end
      time_nesta = toc;
      
      fprintf('NESTA Time: %f\n', time_nesta);
      
    end
  end
  
  methods (Access = 'private')
    function flag = raw_data_loaded(self)
      flag = ~isempty(self.x) || ~isempty(self.y)...
        || ~isempty(self.uz) || ~isempty(self.ze);
    end
  end
  
  methods (Static)

    function [figs, axs] = make_traj_figs(figbase)
      Fig_uz = figure(20+figbase); clf

      ax1 = gca();
      Fig_ze = figure(30+figbase); clf
      ax2 = gca();

      Fig_x = figure(40+figbase); clf
      ax3 = gca();
      Fig_y = figure(50+figbase); clf
      ax4 = gca();

      figs = {Fig_uz, Fig_ze, Fig_x, Fig_y};
      axs = {ax1, ax2, ax3, ax4};
    end

    function idx_state_s = divide_by_state(met_ind)
      names = {'move', 'tdown', 'tsettle', 'scan', 'tup', 'connect'};
      names_idx = [1:6];
      idx_state_s.scan = {};
      idx_state_s.move = {};
      idx_state_s.tdown = {};
      idx_state_s.tsettle = {};
      idx_state_s.tup = {};
      idx_state_s.connect = {};
      met_ind(met_ind <=-6 ) = -6;
      
      met_ind_temp = abs(met_ind);
      idx_end = 0;
      
      break_out = false;
      k=1;
      
      master_idx = 1;
      while ~isempty(met_ind_temp)
        
        state_cur = met_ind_temp(1);
        idx_ = find(met_ind_temp ~= state_cur, 1, 'first')-1;
        slice = 1:idx_;
        if isempty(slice)
          break;
        end
        name = names{state_cur};
        
        n = length(idx_state_s.(name));
        
        idx_set = slice + master_idx - 1;
        idx_state_s.(name){n+1} = idx_set;
      
        master_idx = idx_set(end)+1;

        met_ind_temp(slice) = [];
      end
    end


  end

end




% %     function self = fit_gdrift_per_cycle(self, ax1, ax2, go, gvib)
% %       indc = {'k',        'r', [0, .75, .75], 'b', [.93 .69 .13], ;
% %        'xy-move', 'tip down', 'tip settle',  '$\mu$-path scan', 'tip up',};
% %
% %       state_seq = {'move', 'tdown', 'tsettle', 'scan', 'tup'};
% %       hold(ax1, 'on')
% %       hold(ax2, 'on')
% % %       title(ax1, 'uz')
% % %       title(ax2, 'z-err')
% %       grid(ax1, 'on')
% %       grid(ax2, 'on')
% %       [z, p, k] = zpkdata(go, 'v');
% %       theta0 = [z;p;k*2];
% %       np = length(p);
% %       ub = ones(2*np+1);
% %       ub(end) = Inf;
% %
% %       for idx_cs_seq = 1:length(self.idx_state_s.move);
% %
% %           idx_down = self.idx_state_s.tdown{idx_cs_seq};
% %           idx_settle = self.idx_state_s.tsettle{idx_cs_seq};
% %           idx_scan = self.idx_state_s.scan{idx_cs_seq};
% %           scan_start_idx = length(idx_down) + length(idx_settle) + 1;
% %
% %           uz_ = self.uz([idx_down, idx_settle]);
% %           ze_ = self.ze([idx_down, idx_settle]);
% %           t_  = [idx_down, idx_settle;]*self.Ts;
% %           
% %           uz_whole = self.uz([idx_down, idx_settle, idx_scan]);
% %           ze_whole = self.ze([idx_down, idx_settle, idx_scan]);
% %           t_whole  = [idx_down, idx_settle, idx_scan]*self.Ts;
% %           
% %           u0 = uz_(1);
% %           ze0 = ze_(1);
% %           uz_ = uz_ - u0;
% %           ze_ = ze_ - ze0;
% %           
% %           fun = @(theta) fit_gdrift(theta, gvib, ze_, uz_, t_, np);
% %           
% %           theta = lsqnonlin(fun, theta0, 0.9*ub, ub);
% %           gdrift = zpk(theta(np+1:end-1), theta(1:np), theta(end), self.Ts);
% %           z_fit1 = lsim(gdrift*gvib, uz_whole-u0, t_whole)
% %           z_fit2 = -lsim(gdrift*gvib, uz_whole-u0, t_whole) + u0;
% %           
% %           k = idx_down(1);
% %           plot(ax1, t_whole(idx_down-k+1), z_fit2(idx_down-k+1), 'color', indc{1,2}, 'LineStyle', '--')          
% %           plot(ax1, t_whole(idx_settle-k+1), z_fit2(idx_settle-k+1), 'color', indc{1,3}, 'LineStyle', '-')
% %           plot(ax1, t_whole(idx_scan-k+1), z_fit2(idx_scan-k+1), 'color', indc{1,4}, 'LineStyle', '--')
% %           
% %           plot(ax1, t_whole(idx_down-k+1), uz_whole(idx_down-k+1), 'color', indc{1,2}, 'LineStyle', '-')          
% %           plot(ax1, t_whole(idx_settle-k+1), uz_whole(idx_settle-k+1), 'color', indc{1,3}, 'LineStyle', '-')
% %           plot(ax1, t_whole(idx_scan-k+1), uz_whole(idx_scan-k+1), 'color', indc{1,4}, 'LineStyle', '-')
% % 
% %           plot(ax2, t_whole(idx_down-k+1), z_fit1(idx_down-k+1), 'color', indc{1,2}, 'LineStyle', '--')          
% %           plot(ax2, t_whole(idx_settle-k+1), z_fit1(idx_settle-k+1), 'color', indc{1,3}, 'LineStyle', '-')
% %           plot(ax2, t_whole(idx_scan-k+1), z_fit1(idx_scan-k+1), 'color', indc{1,4}, 'LineStyle', '--')
% % 
% %           plot(ax2, t_whole(idx_down-k+1), ze_whole(idx_down-k+1)-ze0, 'color', indc{1,2}, 'LineStyle', '-')          
% %           plot(ax2, t_whole(idx_settle-k+1), ze_whole(idx_settle-k+1)-ze0, 'color', indc{1,3}, 'LineStyle', '-')
% %           plot(ax2, t_whole(idx_scan-k+1), ze_whole(idx_scan-k+1)-ze0, 'color', indc{1,4}, 'LineStyle', '-')
% %           
% %         keyboard
% %       end
% % 
% % %       linkaxes([ax1, ax2], 'x')
% % %       plot(ax2, [self.t(1), self.t(end)], [.05, .05], '--k')
% % %       plot(ax2, [self.t(1), self.t(end)], -[.05, .05], '--k')
% % 
% %     end      
