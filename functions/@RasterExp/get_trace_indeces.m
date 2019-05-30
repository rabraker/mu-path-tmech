% [ trace_ind ] = get_trace_indeces(nperiods, samps_per_period)

% Returns a set of indeces corresponding only to the trace data, so that we
% can discard the retrace.
%
% use like 
% [ trace_ind ] = get_trace_indeces(nperiods, samps_per_period)
% xdat_trace = xdat(trace_ind)
function [ trace_ind ] = get_trace_indeces(self) 
%nperiods, samps_per_period)

trace_ind = [];
for k=0:self.npix_y-1
    ind = [k*self.samps_per_period+1:(k*self.samps_per_period + floor(self.samps_per_line))]';
    trace_ind = [trace_ind, ind(:)];
end

end

