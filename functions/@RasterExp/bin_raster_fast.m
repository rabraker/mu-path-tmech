function [ pixmat ] = bin_raster_fast(udat_vec, nperiods, samps_per_period)
udat_vec = udat_vec(:);  % make sure we have a column

if length(udat_vec) ~= nperiods*samps_per_period
   s = sprintf('Expected length(udat_vec) == nperiods*samps_per_period')
   s = sprintf('%s \n but have length(udat_vec)=%d, nperiods*samps_per_period = %d',...
       s, length(udat_vec), nperiods*samps_per_period)
  
    error(s)
end

npix = nperiods;


% Get the indeces corresponding to trace data only.
trace_inds = get_trace_indeces(nperiods, samps_per_period);

% And take a slice of those indeces.
udat_trace = udat_vec(trace_inds);

% make each trace a row in a matrix. 
udat_mat = reshape(udat_trace, [], nperiods);


pixmat = zeros(nperiods, nperiods);
for i_row=1:256
   img_row = udat_mat(:,i_row)'; 
    
   u_pix_row = pixelize(img_row, npix);
   
   pixmat(i_row, :) = u_pix_row(1:npix)';
    
end


end

