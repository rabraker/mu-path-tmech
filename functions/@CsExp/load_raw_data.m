function self = load_raw_data(self, cs_paths, opts)
% Loads the raw cs data from the paths defined in cs_paths. This is
% a helper function for the main constructor.
% 
%Arguments
%---------
%  cs_paths :  struct of paths
%              Generate with the function get_cs_paths(). This 
%              function expects that cs-paths will minimally contain 
% 
%  cs_paths.meta_path : (string) 
%                        Path to the matfile saved by laview and contains meta
%                        data generated in the experiment, like
%                        control parameters and settings.
% 
%  cs_paths.parent_meta_path : (string)
%                              This is the .mat file containing the parents
%                              meta data, saved, e.g., by get_cs_paths.m.
%   
%  cs_paths.data_path :        This is the csv file of the actual raw data, saved
%                              by labview.
% 
%  opts :  struct
%          an options structure (typically from name-value pairs
%          passed to the constructor). We exect the following
%          properties.
%  opts.channel_map integer vector
%                   Construct with ChannelMap
%  opts.Ts : double
%            If empty, will use the value in AFM.Ts
%  opts.feature_height: double
%                       Expected feature height of the sample. Used
%                       for visualizing when verbose=true in
%                       process_cs_data().
%  opts.gg :  lti model
%             If this is non-empty, the uz control data will be
%             filtered through this transfer function. Useful for
%             removing drift etc.
% 
% 
%Returns
%-------
%  self : CsExp
%         The current instance with data loaded into the proper
%         properties.    
% 
%See Also
%--------
%  get_cs_paths.m, CsExp.load_mat()
% 
  
  channel_map = opts.channel_map;
  self.channel_map = opts.channel_map;
  self.cs_paths = cs_paths;
  
  dat_meas = csvread(cs_paths.data_path);
  tmp = loadjson(cs_paths.meta_path);  % Provides ExpMetaData
  self.meta_exp = tmp; %.ExpMetaData;
  tmp = loadjson(cs_paths.parent_meta_path); % Provides CsExpMetaIn
  self.meta_in = tmp; %.CsExpMetaIn;

  if size(dat_meas, 2) > 5
      have_friction = true;
  else
      have_friction = false;
  end
  
  if isfield('self.meta_in', 'Ts')
    self.Ts = self.meta_in.Ts;
  elseif ~isnan(opts.Ts)
    self.Ts = opts.Ts;
  else
    self.Ts=AFM.Ts;
  end
  
  self.feature_height = AFM.nm2volts_z * opts.feature_height;
  self.npix = self.meta_in.npix;
  self.width = self.meta_in.width;
   
  
  % Get indices for each state.
  self.met_ind = dat_meas(:, channel_map.met);
  % convert the meta cs-measurment index to -4.
  self.met_ind(self.met_ind > 0) = -4;      
  self.idx_state_s = CsExp.divide_by_state(self.met_ind);
  
  self.t = (0:length(self.x)-1)'*AFM.Ts;
  self.x = dat_meas(:, channel_map.x);
  self.x = self.x; %- min(self.x); % move to positive orthant.
  self.y = dat_meas(:, channel_map.y);
  self.y = self.y; %- min(self.y); % move to positive orthant.
  
  self.ze = dat_meas(:, channel_map.ze);
  
  if have_friction
      self.z_friction = dat_meas(:, channel_map.friction);
  end
  
  if ~isempty(opts.gg) && isa(opts.gg, 'lti')
    fprintf('Performing Dynamic detrend...');
    self.gg = opts.gg;
    self.uz = lsim(opts.gg, dat_meas(:, channel_map.uz), self.t);
    fprintf('done\n');
  elseif  ~isempty(opts.gg) && isa(opts.gg, 'function_handle')
    self.uz = opts.gg(dat_meas(:, channel_map.uz), self.idx_state_s);
  else
    self.uz = dat_meas(:, channel_map.uz);
  end
  
  self.Img_raw = zeros(self.npix, self.npix);
  self.Img_smp1d = zeros(self.npix, self.npix);
  self.Img_bp = zeros(self.npix, self.npix);
  self.pix_mask = zeros(self.npix, self.npix);
  self.Gz = zpk([], [], 1, self.Ts);
  
  state_ticks = self.meta_exp.state_counts;
  self.state_times = state_ticks*self.Ts;
  self.time_total = sum(self.state_times);

  
  
end

% if ~isa(channel_map, 'ChannelMap')
%   error(['channel_map must be of class ChannelMap, but is a ' ...
%          '%s'], class(channel_map));
% end
