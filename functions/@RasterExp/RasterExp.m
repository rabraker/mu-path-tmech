classdef RasterExp <matlab.mixin.Copyable
% rast_exp = RasterExp(raster_paths, varargin)
% 
% 'rel
  properties
    raster_paths;
    channel_map;
    % parent_data;
    xref;
    yref;
%     meta_data;
    Ts;
    met_ind;
    npix_x;
    npix_y;
    npix_y_og;
    width;
    volts2pix_y;
    volts2pix_x;
    micron2pix_y;
    micron2pix_x;

    samps_per_period;
    samps_per_line;
    gg;
    meta_exp;
    meta_in;
    time_total;
    uz;
    ze;
    z_friction;
    x;
    y;
    pix_mat;
    pix_mat_pinned;
    pin_idx_s;
    pix_mask;
    
    pix_mask_ze;
    pix_mat_ze;
%     UserData;
  end
  
  methods
    function self= RasterExp(raster_paths, varargin)
    % self= RasterExp(raster_paths, varargin)
    
    % [datmat, samps_period, samps_line] 
      self.raster_paths = raster_paths;
      default_chan_map = ChannelMap([1:6]);
      optP = inputParser();
      optP.addParameter('reload_raw', false, @(s)islogical(s));
      optP.addParameter('channel_map', default_chan_map);
      optP.addParameter('gg', []);
      optP.addParameter('load_full', false, @(s)islogical(s));
      
      optP.parse(varargin{:});
      opts = optP.Results;
      
      % Check to see if an already processed mat file exists.
      mat_exists = exist(raster_paths.data_path_mat, 'file');
      if ~opts.reload_raw && mat_exists
        fprintf('loading from mat file...\n')
        self = self.load_mat(raster_paths.data_path_mat, opts.load_full);
        fprintf('\bdone\n')
      else
        fprintf('loading from raw data...\n')
        self = self.load_raw_data(raster_paths, opts);
        self.npix_y_og = self.npix_y;
        fprintf('\bdone\n')
      end      
      self.time_total = length(self.x)*self.Ts;
      
    end
    
    % Methods defined in other files
    self = load_raw_data(self, raster_paths, npix, width, opts)
    [ self] = bin_raster_really_slow(self, line_detrender, use_error, npix_x)
    trace_inds = get_trace_indeces(self)
    
    im = interp_missing(self, im, mask);
    
    function meta_data(self)
      error(['For consistency with CsExp, the meta_data property has been',...
        ' removed. That data is now located in self.meta_in.'])
    end
    
    function interp_y(self, ypix, fill_missing)
        
        if size(self.pix_mat, 1) == ypix
            return
        end
        
        if nargin < 3
            fill_missing = false;
        end
        
        if fill_missing
            self.pix_mat = self.interp_missing(self.pix_mat);
        end
        ypix_old = self.npix_y_og;
        xpix = size(self.pix_mat, 2);
        img_uz = zeros(ypix, xpix);
        img_ze = zeros(ypix, xpix);
        
        % Now interpolate the missing rows
        vq = 1:ypix;
        skip = floor(ypix/ypix_old);
        h = skip:skip:ypix;
        
        for k=1:xpix
           img_uz(:, k) = interp1(h, self.pix_mat(:, k), vq, 'linear', 'extrap');
           img_ze(:, k) = interp1(h, self.pix_mat_ze(:, k), vq, 'linear', 'extrap');
        end
        self.npix_y = ypix;
        self.pix_mat = img_uz;
        self.pix_mat_ze = img_ze;
        % Create a new pix_mask
        %  pix_mask_old = self.pix_mask;
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
      err_pos = err(err>0);
      damage = sum(err_pos.^2)/length(err_pos)/self.Ts;
    end

    function quality = quality_metric(self)
    % Compute a quality metric based on the deflection signal's power.

      ref = self.meta_exp.z_axis_params.setpoint_scan;
      % rather than subtracting mean, subtract the reference value.
      err = self.ze - ref;  % shift to zero.
      quality = sum(abs(err).^2)/length(err)/self.Ts;
    end
    
    function offset = find_lag(self)
    % Use cross correlation between the reference and the second period to find
    % the number of samples of offset between the two waveforms.
        xr = self.xref/AFM.volts2mic_xy;
        
        % take the second period, after we have steady state.
        idx_start = self.samps_per_period;
        idx_end = idx_start+self.samps_per_period;
        x_test = self.x(idx_start+1:idx_end);
        
        [xc, lags] = xcorr(xr, x_test);
        [~, idx] = max(abs(xc));
        
        offset = -lags(idx);
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
        save(self.raster_paths.data_path_mat, '-struct', 'self_struct', '-append');
      else
        save(self.raster_paths.data_path_mat, '-struct', 'self_struct');
      end
      warning('on', 'MATLAB:structOnObject');
    end

    function self = load_mat(self, data_path_mat, load_full)
    % Load ourself from the location contained in data_path_mat. If
    % load_full=false (the default), then the original time-series
    % data, x,y,uz,ze, x_positive,y_positive will not be loaded.
    % This is to speed things up when we just want to work with the
    % already processed images.

    % I tried use the matfile() function. This saves us zero time. I cant
    % figure out how to do this without creating an extra structure
    % and passing that into or out of self. This is wastful...
    %%%       self_mat = matfile(data_path_mat);
    
      to_load_list = properties(self);
      if ~load_full
        no_loads = {'x', 'y', 'uz', 'ze'};
        to_load_list = setdiff(to_load_list, no_loads);
      end

      data = load(data_path_mat, to_load_list{:});

      for prop = to_load_list'
        self.(prop{1}) = data.(prop{1});
      end

    end
    
    
    function plot_n_periods(self, signal, ax, N0, N1)
    % plot_n_periods(self, signal, ax, N0, N1)
    % 
    % Plot the signal ('x', 'y', 'ze', 'uz') to axis ax between
    % raster periods N0 and N1.
      sig = self.(signal)(N0 * self.samps_per_line + 1: N1* ...
                          self.samps_per_line);
      
      plot(ax, sig);
    end


    function [uz_row, x_row] = get_row(self, row)
    % [uz_row, x_row] = get_row(row)
    % Extract from the raw data a single row of control and
    % x-direction data, which should
    % correspond to an actual row in the post-processed image.      
      start_idx = self.samps_per_period*(row-1) + 1;
      end_idx = start_idx + self.samps_per_line;

      x_row = self.x(start_idx:end_idx);
      uz_row = detrend(self.uz(start_idx:end_idx) );
      x_row = x_row - x_row(1);
      
      x_row = x_row * (self.npix/x_row(end) );
    end
  end
  
  methods (Access = 'private')
    function flag = raw_data_loaded(self)
      flag = ~isempty(self.x) || ~isempty(self.y)...
        || ~isempty(self.uz) || ~isempty(self.ze);
    end
  end
    
end
