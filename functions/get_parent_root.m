% parent_root = get_parent_root(data_root)
%  
% Given a path data_root, e.g.,
% /media/labserver/afm-cs/imaging/cs-imaging/5microns/11-13-2018
%
% produces the location of the parent data
%
% /media/labserver/afm-cs/imaging/cs-imaging/5microns/parents

function parent_root = get_parent_root(data_root, parent_folder)
  if nargin < 2
      parent_folder = 'parents';
  end    
  % remove trailing slash
  data_root = strip(data_root, 'right', filesep);
  
  % break path apart and drop the final term
  fparts = strsplit(data_root, filesep);
  
  fparts = {fparts{1:end-1}, parent_folder};
  parent_root = strjoin(fparts, filesep);
  
end