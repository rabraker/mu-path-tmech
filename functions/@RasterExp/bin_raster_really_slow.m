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

function [ self] = bin_raster_really_slow(self, line_detrender, use_error, npix_x)
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
  
  if ~exist('use_error', 'var') || isempty(use_error)
      use_error = false;
  end
  if ~exist('npix_x', 'var') || isempty(npix_x)
      npix_x = self.npix_x;
  else
      micron2pix_x = npix_x/self.width;
      volts2pix_x = AFM.volts2mic_xy  * micron2pix_x;
      
      % Update the instance properties
      self.npix_x = npix_x;
      self.volts2pix_x = volts2pix_x;
      self.micron2pix_x = micron2pix_x;
  end
  
  xpix = self.npix_x; 
  ypix = self.npix_y_og; 
  
    
  if use_error
      udat = self.ze; 
  else
      udat = self.uz; 
  end
  
  % Make the image be at (0,0, --) and convert to pixel coordinates.
  %xdat_trace = (xdat_trace - min(xdat_trace))*self.volts2pix_x;
  %ydat_trace = (ydat_trace - min(ydat_trace))*self.volts2pix_y;

  x_dat_pix = self.x * self.volts2pix_x;
  
  pix_mask = zeros(ypix, xpix);
  pix_mat = zeros(ypix, xpix);
  
  samps_per_line = self.samps_per_line;
  samps_per_period = self.samps_per_period;
  
  % The raster scans are offset from the reference by approximately 
  % ess = abs(freqresp(Her*Int_z, 1) samples, where Int_z is an integrator.
  offset = self.find_lag();
  
  for j_row = 0:ypix-1
    % Guard against the possibility that samps_per_line is not an integer.
    % Need to fix this in the raster trajectory generator.
    
    indy_start = ceil(j_row*(samps_per_period)+1);
%     indy_end = floor((j_row+1)*(samps_per_period));
    indy_end = indy_start + samps_per_line - 1;
    ind_y = (indy_start:indy_end) + offset;
    assert(length(ind_y) == samps_per_line);
    try
      x_dat_j = x_dat_pix(ind_y);
    catch
      warning('Exited early from processing raster data (j_row=%d)', j_row);
      break
    end
    U_dat_j_init = udat(ind_y)';
    
    [U_dat_j] = line_detrender(U_dat_j_init);

    for i_col = 0:xpix-1
      ind_x = find(x_dat_j >= i_col & x_dat_j < i_col+1);
      u_mean_ij = mean(U_dat_j(ind_x));
      if ~isnan(u_mean_ij)
        pix_mat(j_row+1, i_col+1) = u_mean_ij;
        pix_mask(j_row+1, i_col+1) = 1;
      else
        %             keyboard
      end
      
    end
    
  end
  
  if use_error
      self.pix_mat_ze = pix_mat;
      self.pix_mask_ze = pix_mask;
  else
      self.pix_mat = pix_mat;
      self.pix_mask = pix_mask;
  end

end


function [x] = line_detrender_local(x)
    return
end
