% [sig_psd, freqs, k] = psd_from_intervals(self, signal,state, starts, ends)
%
% Given a vector of start times, starts, and end times ends,
% computes the PSD of 'signal' by averaging the periodograms of
% each CS cycle for 'state' between all start and end times.
%
% Inputs
% -----
%   signal : string = ('x' | 'y' | 'uz'| 'ze' | )
%   state  : string = ('scan' | 'move' | 'tdown' | 'tsettle' |'tup')
%   
%   starts : vector containing all the start times. Must have the
%            same length as ends.
%   ends   : vector containing all the end times.
%   no_detrend: boolean, optional. If set to true, the data in each
%               segment will not be first detrended before taking the periodogram.
% 
% Outputs
% ------
%   sig_psd : the PSD estimate of signal based on averaging the
%             periodograms signal during state 'state' between the
%             intervals defined by starts and ends. The estimates
%             are compuated based on the shortest data set.
%   freqs  :  The frequencies associated with the PSD estimate.
%   k      :  The total number of averages taken.
%
%
%Example
%-------
%  This example demonstrates the use case this method was designed
%  for. Suppose we wish to compute the PSD of the control signal during
%  the intervals 1.2 to 2.5 seconds and 7.1 to 9.4 seconds. We wish
%  to compute the PSD during the scanning state. The idea here is
%  to see what the frequency content looks like. If we selected the
%  time intervals to correspond to flat portions of the sample,
%  then we know that anything extra is dynamics from poor control.
% 
%  To accomplish this, we set:
%  signal = 'uz';
%  state = 'scan';
%  starts = [1.2, 7.1];
%  ends  = [2.5, 9.4];
% 
% [sig_psd, freqs] = cs_exp.psd_from_intervals(signal, state, starts, ends);
% figure;
% semilogx(freqs, 10*log10(sig_psd));
%
% See Also : get_idx_by_state_in_time_range
function [sig_psd, freqs, k] = psd_from_intervals(self, signal, ...
                                                  state, starts, ...
                                                  ends, no_detrend)
  if ~exist('no_detrend', 'var')
    no_detrend = false;
  end
  
  % ensure inputs are sane
  if length(starts) ~= length(ends)
    error('Must have the same number of starts and ends')
  end
  if ~any(strcmp(signal,  {'x', 'y', 'uz', 'ze'}))
    error(['Expected signal = (''x'' | ''y'' | ''uz''| ''ze'' | ), ' ...
           'but recieved %s'], signal)
  end
  if ~any(strcmp(state, {'scan', 'move', 'tdown', 'tsettle', ...
                        'tup'}))
    error(['Expted state = (''scan'' } ''move'' | ''tdown'', ' ...
           '''tsettle''| ''tup'') but recieved %s '], state)
  end
  
  
  idx_cell = {};
  
  N = Inf;
  for k=1:length(starts)
    [idx_s_k, Nk] = self.get_idx_by_state_in_time_range(state, starts(k), ends(k));
    % Concatenate the cell array we just got onto the running set.
    idx_cell(end+1:end+length(idx_s_k)) = idx_s_k(:);
    % [uz_mat_k, nk] = extract_uz_in_time_range(self, starts(k), ends(k));
    N = min(N, Nk);
    
  end
  
    
  sig_psd = 0;
  win = window(@hann, N);
  for k=1:length(idx_cell)
    sig_k = self.(signal)(idx_cell{k});
    if ~no_detrend  
      sig_k = detrend(sig_k);
    end
    
    % Compute single realization power spectrum. Truncate to the 
    % shortest of all the realizations.
    try 
    [sig_psd_k, freqs] = power_spectrum(sig_k(1:N), self.Ts, win);
    catch
      keyboard
    end
    
    sig_psd = sig_psd + sig_psd_k;
    
  end
  sig_psd = sig_psd/k;


end
