%[idx_cell, N_min]= get_idx_by_state_in_time_range(self, state_name, t_start, t_end)
% Given a start time, t_start and end time t_end, this function
% will extract all indeces of state_name within that time range.
%
% For example: suppose t_start = 1.2, and t_end = 5 
% and state_name = 'scan'. Then the function returns a matrix of
% indeces idx_mat such that 
% self.uz(idx_mat{k} returns the control signal for the kth
% scan in the range t_start to t_end.
%
% Inputs
% ------
% state_name:  string = ('scan' | 'move' | 'tdown' | 'tsettle' |'tup')
%              Note that these correspond to the field names in idx_state_s
% 
% t_start   :  Time at which the search begins.
% 
% t_end     :  Time at which the search ends.
%
% Outputs
% ------
% idx_cell  : A cell array. Each cell contains a set of indeces
%             such indexing the raw data on the kth index set returns
%             the data for the kth scan in range t_start to t_end.
% 
% N_min    :  Integer. The length of the smallest index set.
%
% See Also : psd_from_intervals
function [idx_cell, N_min]= get_idx_by_state_in_time_range(self, state_name, t_start, t_end)

  CS_idx1 = self.find_cycle_idx(t_start);
  CS_idx2 = self.find_cycle_idx(t_end);
  
  N_min = Inf;
  idx_cell = {};
  j = 1;
  for k = CS_idx1:CS_idx2
    state_idx = self.idx_state_s.(state_name){k};
    % uz_scan = self.uz(state_idx);
    % uz_scan = self.ze(state_idx);
    % uz_scan = uz_scan - mean(uz_scan);
    N_min = min(N_min, length(state_idx));
    
    idx_cell{j} = state_idx;
    j = j+1;
  end

end
