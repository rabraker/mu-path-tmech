function process_cs_data(self, verbose, figs, use_ze, register_uzk)
  %  process_cs_data(self, verbose, figs, use_ze, register_uzk)
  %
  %  Process the raw CS data. This will populate self.Img_raw and 
  %  self.pix_mask.
  %
  %  Arguments
  %  -----------
  %     verbose: integer. Plot some stuff if greater than 0 (broken at the moment).
  %     figs: only used with verbose
  %     use_ze: [true|false] if true, process deflection data, rather than
  %     control data.
  %     register_uzk: [true|false] Only used if use_ze==false. If true, subtract
  %     the max value of each mu-path off.
  %
  if isempty(self.x) || isempty(self.y) || isempty(self.uz) || isempty(self.ze)
    warning(['You have empty raw data fields. Reload with ',...
      '"load_full=true".', 'Skipping']);
    return
  end
  if ~exist('use_ze', 'var') || isempty(use_ze)
      use_ze = false;
  end
  if ~exist('register_uzk', 'var') || isempty(register_uzk)
      register_uzk = true;
  end
  
  % bin all the data into pixels.
  tend_last = 0;
  microns_per_volt = AFM.volts2mic_xy;
  pix_per_volt = (self.npix/self.width)*microns_per_volt;
  if isempty(self.x_positive) || isempty(self.y_positive)
    self.xy_positive();
  end
  if nargin <2
    verbose = false;
  end
  if verbose && ~exist('figs', 'var')
    figs{1} = figure;
    figs{2} = figure;
    figs{3} = figure;
  end
  if verbose
    [~, ax1] = parse_fig_ax(figs{1});
    [~, ax2] = parse_fig_ax(figs{2});
    [~, ax3] = parse_fig_ax(figs{3});
  end
  
  % Reset the pix_mask
  Ki=self.meta_exp.z_axis_params.Ki_scan;
  DI = zpk(0, 1, Ki, AFM.Ts);
  ref = self.meta_exp.z_axis_params.setpoint_scan;
  self.pix_mask = self.pix_mask*0;
  for k = 1:length(self.idx_state_s.scan)
    
    % Get the data for the current mu-path.
    if use_ze
        [X_raw, Y_raw, ~, U_scan] = self.get_scan_k(k, true);
    else
        [X_raw, Y_raw, U_scan] = self.get_scan_k(k, true);
    end
    Y_raw = Y_raw*pix_per_volt;
    X_raw = X_raw*pix_per_volt;
    
    if max(U_scan) - min(U_scan) > 0.4 % throw out rediculous data.
      fprintf('skipping cycle %d\n', k)
      continue
    end

    [y_idx, x_idx, U_k] = self.mu_data2pix_xy(X_raw, Y_raw, U_scan);

    if ~use_ze && register_uzk
        % Register the control data to zero. We can do this because we are
        % scanning long enough that we are always guaranteed to exit a hole.
        U_k = U_k - max(U_k);
    end
    
    assert(length(y_idx) == length(x_idx))
    for k=1:length(x_idx)
        self.Img_raw(y_idx(k), x_idx(k)) = U_k(k);
        self.pix_mask(y_idx(k), x_idx(k)) = 1;
    end
    
%  THIS IS BASICALLY BROKEN WITH THE MU-PATH-CONNECT SCHEME    
%     if verbose
%       t_k = (0:length(U_z)-1)'*self.Ts;
%       idx = length(U_z) - length(U_scan);
%       % figure(Fig1);
%       plot(ax1, t_k(1:idx)+tend_last, U_orig(1:idx)-U_orig(end), '--g');
%       h_f1_uog = plot(ax1, t_k(idx:end)+tend_last, U_orig(idx:end)-U_orig(end), 'k');
%       
%       plot(ax1, t_k(1:idx)+tend_last, U_z(1:idx)-U_z(end), '--b');
%       h_f1_uz = plot(ax1, t_k(idx:end)+tend_last, U_z(idx:end)-U_z(end), 'r');
%       tend_last = t_k(end) + tend_last;
%       U_ = U_orig(idx:end);
%       t_ = (0:length(U_)-1)'*self.Ts;
%       % figure(Fig2);
%       h_f2_uog = plot(ax2, t_, U_ - max(U_), 'k');
%       h_f2_uz = plot(ax2, t_, U_z(idx:end) - max(U_z(idx:end)), 'r');
%       
%       % -------------------
%       % ---- visualize ------
%       if abs(max(U_k) - min(U_k))> self.feature_height
%         % Then we have an edge.
%         cs = 'r';
%       else
%         cs = 'b';
%       end
%       plot(ax3, U_k, 'color', cs)
%     end
    drawnow();
  end % main loop
  
  if verbose % draw legends
    h_f1_uog.DisplayName = 'original';
    h_f1_uz.DisplayName = 'Dynamic Detrend';
    h_f2_uog.DisplayName = 'original';
    h_f2_uz.DisplayName = 'Dynamic Detrend';
    legend(ax1, [h_f1_uog, h_f1_uz])
    legend(ax2, [h_f2_uog, h_f2_uz])
  end
  
  [n, m] = size(self.Img_raw);
  if n > self.npix || m > self.npix
     warning(['Raw image is oversize: expected %d x %d, but have %d x %d.\n',... 
         'Resizing...'], self.npix, self.npix, n, m) 
  end
  self.Img_raw = self.Img_raw(1:self.npix, 1:self.npix);
  self.pix_mask = self.pix_mask(1:self.npix, 1:self.npix);
end % process_cs_data()