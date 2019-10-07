function process_cs_data(self, use_ze, register_uzk)
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
  
  
  pix_mask = zeros(self.npix, self.npix);
  pix_mat_raw = zeros(self.npix, self.npix);
  
  % bin all the data into pixels.

  microns_per_volt = AFM.volts2mic_xy;
  pix_per_volt = (self.npix/self.width)*microns_per_volt;
  if isempty(self.x_positive) || isempty(self.y_positive)
    self.xy_positive();
  end

  
  for k = 1:length(self.idx_state_s.scan)
    
    % Get the data for the current mu-path.
    if use_ze
        [X_raw, Y_raw, ~, U_scan] = self.get_scan_k(k, true);
    else
        [X_raw, Y_raw, U_scan] = self.get_scan_k(k, true);
    end
    Y_raw = Y_raw*pix_per_volt;
    X_raw = X_raw*pix_per_volt;
%     figure(1)
%     subplot(2,1,1);
%     cla
%     plot(detrend(X_raw))
%     hold on
%     plot([0, length(X_raw)], [0.5, 0.5], 'r')
%     plot([0, length(X_raw)], -[0.5, 0.5], 'r')
%     subplot(2,1,2)
%     plot(Y_raw)
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
    for j=1:length(x_idx)
        pix_mat_raw(y_idx(j), x_idx(j)) = U_k(j);
        pix_mask(y_idx(j), x_idx(j)) = 1;
    end

  end % main loop
  
  
  [n, m] = size(pix_mat_raw);
  if n > self.npix || m > self.npix
     warning(['Raw image is oversize: expected %d x %d, but have %d x %d.\n',... 
         'Resizing...'], self.npix, self.npix, n, m) 
  end
  if use_ze
      self.pix_mat_raw_ze = pix_mat_raw(1:self.npix, 1:self.npix);
      self.pix_mask_ze = pix_mask(1:self.npix, 1:self.npix);
  else
      self.pix_mat_raw_uz = pix_mat_raw(1:self.npix, 1:self.npix);
      self.pix_mask_uz = pix_mask(1:self.npix, 1:self.npix);
  end
end % process_cs_data()