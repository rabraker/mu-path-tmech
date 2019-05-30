function self = load_raw_data(self, raster_paths, opts)
  % Loads the raw raster data from the paths defined in raster_paths. This is
  % a helper function for the main constructor.
  %
  %Arguments
  %---------
  %  raster_paths :  struct of paths
  %                  Generate with the function get_raster_paths(). This
  %                  function expects that cs-paths will minimally contain
  %
  %  raster_paths.meta_path : (string)
  %                           Path to the matfile saved by laview and contains meta
  %                           data generated in the experiment, like
  %                           control parameters and settings.
  %
  %  raster_paths.parent_meta_path : (string)
  %                                  This is the .mat file containing the parents
  %                                  meta data, saved, e.g., by get_raster_paths.m.
  %
  %  raster_paths.data_path : This is the csv file of the actual raw data, saved
  %                           by labview.
  %  opts :  struct
  %          an options structure (typically from name-value pairs
  %          passed to the constructor). We exect the following
  %          properties.
  %  opts.channel_map ChannelMap object.
  %                   Defines the indeces for x,y,ze,uz data
  %                   channells in the csv file.
  %  opts.Ts : double
  %            If empty, will use the value in AFM.Ts
  %  opts.feature_height: double
  %                       Expected feature height of the sample. Used
  %                       for visualizing when verbose=true in
  %                       process_raster_data().
  %  opts.gg :  lti model
  %             If this is non-empty, the uz control data will be
  %             filtered through this transfer function. Useful for
  %             removing drift etc.
  %
  %Returns
  %-------
  %  self : RasterExp
  %         The current instance with data loaded into the proper
  %         properties.
  %
  %See Also
  %--------
  %  get_raster_paths.m, RasterExp.load_mat()
  %
  
  channel_map = opts.channel_map;
  self.channel_map = opts.channel_map;
  fprintf('Loading Meta file...\n%s\n', raster_paths.meta_path)
  meta_data = loadjson(raster_paths.meta_path);
  %   self.meta_data = meta_data.scan_meta;
  self.meta_in = meta_data.scan_meta;
  self.meta_exp = meta_data;
  self.Ts = AFM.Ts;
  
  % Load parent data.
  parent_dat = loadjson(raster_paths.parent_path);
  xyref = reshape(parent_dat.fpga_input', 2, [])';
  self.xref = xyref(:,1);
  self.yref = xyref(:,2);
  
  fprintf('Loading Data file...\n%s\n', raster_paths.data_path)
  datmat = csvread(raster_paths.data_path);
  if size(datmat, 2) > 5
      have_friction = true;
  else
      have_friction = false;
  end
  
  
  npix = meta_data.scan_meta.npix;
  width = meta_data.scan_meta.width;
  % The experiment doesnt care about npix_x. Until we are told otherwise, assume
  % we have as many x and y pixels.
  self.npix_y = npix;
  self.npix_x = npix;
  self.width = width;
  
  
  % Do a sanity check on width:
  xmax = max(datmat(:, channel_map.x));
  xmin = min(datmat(:, channel_map.x));
  % xref is volts. Convert to microns.
  width_exp = round( (xmax - xmin) * AFM.volts2mic_xy, 0);
  
  if abs(width_exp - width) > .2
    warning(['Experiment width is %f, but claimed width from meta data is %f'....
      'Are you loading the right files?'], width_exp, width);
  end
  
  % Compute unit conversions
  micron2pix_x = self.npix_x/width;
  micron2pix_y = self.npix_y/width;
  volts2pix_x = AFM.volts2mic_xy  * micron2pix_x;
  volts2pix_y = AFM.volts2mic_xy  * micron2pix_y;
  self.volts2pix_x = volts2pix_x;
  self.volts2pix_y = volts2pix_y;
  self.micron2pix_x = micron2pix_x;
  self.micron2pix_x = micron2pix_y;
  
  self.samps_per_period = self.meta_in.points_per_period;
  self.samps_per_line = self.samps_per_period/2;
  if floor(self.samps_per_line)~= self.samps_per_line
    warning('Non integer number of samples per line = %f.', ...
      self.samps_per_line)
  end
  
  % Pull out data. First, drop anything extra that got collected.
  
  % [1:npix*self.samps_per_period]
  met_ind = datmat(:, channel_map.met);
  
  start_idx = find(met_ind == 1, 1, 'first');
  end_idx = find(met_ind == 1, 1, 'last');
  slice = start_idx:end_idx;
  
  % We have an off by one error: I think this comes from not taking the last
  % sample when we lift the tip. Try to preserve the total number of specified
  % samples, so we get the full image of pixels.
  slice_ = start_idx:npix*self.samps_per_period+start_idx-1;
  if length(slice_) ~= length(slice)
    warning('Should have %d scanning samples, but only have %d',...
      length(slice_), length(slice));
  end
  if size(datmat, 1) >= slice_(end)
    warning("Taking %d samples from state %d.\n", length(slice_)-length(slice),...
      met_ind(slice_(end)+1));
    
    slice = slice_;
  end
  
  self.x = datmat(slice, channel_map.x);
  self.y = datmat(slice, channel_map.y);
  self.ze = datmat(slice, channel_map.ze);
  self.uz = datmat(slice, channel_map.uz);
  self.met_ind = datmat(slice, channel_map.met);
  if have_friction
      self.z_friction = datmat(slice, channel_map.friction);
  end
  
end
