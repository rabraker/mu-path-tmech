% [raster_paths] = get_raster_paths(data_root, data_name,  parent_folder)
%
% Given data_root and data_name, returns a struct raster_paths with the 
% following fields:
% 
% raster_paths.data_path
% raster_paths.meta_path
% raster_paths.parent_path
% raster_paths.data_path_mat

function [raster_paths] = get_raster_paths(data_root, data_name, parent_folder)
  if nargin < 3
      parent_folder = 'parents';
  end
  parent_root = get_parent_root(data_root, parent_folder);
  parent_name = get_parent_name(data_name, '_out_', '.json');
  
  raster_exp_meta_name = strrep(data_name, '.csv', '-meta.json');
  % Location of the raw, experimental data.
  raster_paths.data_path = fullfile(data_root, data_name);
  % Location of the post-experiment meta-data.
  raster_paths.meta_path = fullfile(data_root, raster_exp_meta_name);
  % Where the input file is
  raster_paths.parent_path = fullfile(parent_root, parent_name);
  % Where to save the post-processed data.
  raster_paths.data_path_mat = strrep(raster_paths.data_path, '.csv', '.mat');
end