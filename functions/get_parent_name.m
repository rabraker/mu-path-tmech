function [parent_name ] = get_parent_name(child_name, child_handle, varargin)

  if length(varargin) == 1
    parent_ext = varargin{1};
  else
    parent_ext = '.csv';
  end

  k = regexp(child_name, child_handle);
  parent_name = sprintf('%s%s', child_name(1:k-1), parent_ext);



end

