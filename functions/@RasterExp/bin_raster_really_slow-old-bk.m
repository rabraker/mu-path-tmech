% [ self] = bin_raster_really_slow(self, line_detrender)
% 
% Bins the raw raster data based on x,y, and uz. The data is binned
% and averaged according to the actual values x, converted to
% pixels. This means that it is possible, due to noise, for the
% data to be binned non-contiguously.
%
% This method will populate the properties of a RasterExp
% instance of self.pix_mat and self.pix_mask.
% The pix_mask is a mask, similar to the CS mask, which will have
% ones where we actually have data and zeros elsewhere. This may be
% usefull for detrending the image.
%
% If you wish to detrend each line, provide a handle to function to
% do so. The prototype for such a function is 
% [x] = line_detrender_local(x)

function [ self] = bin_raster_really_slow(self, line_detrender, use_error)
  if isempty(self.x) || isempty(self.y) || isempty(self.ze) || isempty(self.uz)
    warning(['Raw dat properties are empty. To process raw data, load',...
      'the raw data with the load_full flag. Skipping'])
    return
  end
  % nperiods, samps_per_period, volt2pix,
  % if size(xyu_dat, 1) ~= nperiods*samps_per_period
  %   s = sprintf('Expected length(xyu_dat) == nperiods*samps_per_period')
  %   s = sprintf('%s \n but have length(xyu_dat)=%d, nperiods*samps_per_period = %d',...
  %     s, length(xyu_dat), nperiods*samps_per_period)
  %   error(s)
  % end
  if ~exist('line_detrender', 'var') || isempty(line_detrender)
    line_detrender = @(x) line_detrender_local(x);
  end
  
  if ~exist('use_error', 'var')
      use_error = false;
  end
  xpix = self.npix; 
  ypix = self.npix; % TODO: make this work with rectangular image.
    
  % Get the indeces corresponding to trace data only.
  trace_inds = self.get_trace_indeces();
    
  xdat_trace = self.x(trace_inds);
  ydat_trace = self.y(trace_inds);

  if use_error
      udat_trace = self.ze(trace_inds);
  else
      udat_trace = self.uz(trace_inds);
  end
  
  % Make the image be at (0,0, --) and convert to pixel coordinates.
  xdat_trace = (xdat_trace - min(xdat_trace))*self.volts2pix;
  ydat_trace = (ydat_trace - min(ydat_trace))*self.volts2pix;
  % make
  
  self.pix_mask = zeros(xpix, ypix);
  self.pix_mat = zeros(xpix,ypix);
  
  samps_per_line = self.samps_per_line;
  for j_row = 0:ypix-1
    % Guard against the possibility that samps_per_line is not an integer.
    % Need to fix this in the raster trajectory generator.
    indy_start = ceil(j_row*(samps_per_line)+1);
    indy_end = floor((j_row+1)*(samps_per_line));
    ind_y = indy_start:indy_end;
    try
    x_dat_j = xdat_trace(ind_y)';
    catch
      warning('Exited early from processing raster data (j_row=%d)', jrow);
      break
    end
    U_dat_j_init = udat_trace(ind_y)';
    
    [U_dat_j] = line_detrender(U_dat_j_init);

    for i_col = 0:xpix-1
      ind_x = find(x_dat_j >= i_col & x_dat_j < i_col+1);
      u_mean_ij = mean(U_dat_j(ind_x));
      if ~isnan(u_mean_ij)
        self.pix_mat(j_row+1, i_col+1) = u_mean_ij;
        self.pix_mask(j_row+1, i_col+1) = 1;
      else
        %             keyboard
      end
      
    end
    
  end

end

function [x] = line_detrender_local(x)
    return
end
