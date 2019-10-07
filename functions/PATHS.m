classdef PATHS
  %PATHS Summary of this class goes here
  %   Detailed explanation goes here
  
  properties
  end
  
  methods (Static)
    function [ PATH ] = exp()
      % PATH constant to where all experimental data is stored for the
      if ispc
        PATH = 'Z:\afm-cs';
      else
        PATH = '/media/labserver/afm-cs';
      end
    end
    function path = tmech_fig()
       path = '/home/arnold/gradschool/publications/tmech-afm-cs/latex/figures/'; 
    end
    function path = tmech_table()
      path = '/home/arnold/gradschool/publications/tmech-afm-cs/latex/tables';
    end
    
    function path = cs_final_fig()
      path = '/home/arnold/gradschool/thesis/plots-afm-cs-final/figures';
    end
    function path = cs_final_table()
      path = '/home/arnold/gradschool/thesis/plots-afm-cs-final/tables';
    end
    function path = thesis_root()
      path = '/home/arnold/gradschool/thesis';
    end
    function path = defense_fig()
      path = '/home/arnold/gradschool/thesis/defense/presentation/figures';
    end
    function [ PATH ] = step_exp()
      % PATH constant to where all experimental data for step experiments 
      % (for the x-y plane) is stored for the.
        PATH = fullfile(PATHS.exp(), 'step-exps');
    end
    
    function PATH = cs_image_data(size, date_dir)
      PATH = fullfile(PATHS.exp(), 'imaging', 'cs-imaging', size, date_dir);
    end
    
    function PATH = raster_image_data(size, date_dir)
      PATH = fullfile(PATHS.exp(), 'imaging', 'raster', size, date_dir);
    end
    
    function [ PATH ] = CS_root()
    % PATH constant to where all experimental data is stored for the
    % MPC journal paper.
      if ~ispc()
        PATH = fullfile(getMatPath, 'afm-cs');
      else
        PATH = 'C:\Users\arnold\Documents\afm-cs';
      end
    end
    
    function [ PATH ] = labview()
    % PATH constant to where all experimental data is stored for the
    % MPC journal paper.
      
      PATH = fullfile(PATHS.CS_root, 'labview');
    end
    
    function path = reconstruction_BP()
      path = fullfile(PATHS.CS_root, 'reconstruction', 'BP');
    end
    function path = reconstruction_SMP1D()
      path = fullfile(PATHS.CS_root, 'reconstruction', 'SMP_1D');
    end
    
    function [ PATH ] = sysid()
      % PATH constant to where the system ID data is stored
      
      PATH = fullfile(PATHS.exp, 'sysID');
    end
    
    function PATH = tuesday_fig_path(date_folder, fname)
    % PATH = tuesday_fig_path(date_folder, fname)
    % Constructs a figure path to the tuesday figs subfolder, date_folder,
    % for figure file fname.
    % If date_folder does not exist, creates it.
      folder_path = fullfile(PATHS.CS_root, 'matlab-code', 'tuesday-figs',...
        date_folder);
      if ~exist(folder_path, 'dir')
        mkdir(folder_path);
      end
      PATH = fullfile(folder_path, fname);
    end
    
    function PATH = thesis_fig_final
      PATH = fullfile(PATHS.thesis_root, 'plots-afm-cs-final', 'figures');
    end
    
    function PATH = note_fig(fname)
      PATH = fullfile(PATHS.CS_root, 'matlab-code', 'notes', 'figures');
      if nargin >0
        PATH = fullfile(PATH, fname);
      end

    end
    
  end
  
end

