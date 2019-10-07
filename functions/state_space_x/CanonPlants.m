classdef CanonPlants
  
  properties
    path = fullfile(PATHS.sysid, 'FRF_data_current_stage2.mat');
    plants;
  end
  
  
  methods (Static)
    function [plants, frf_data, modelFit] = plants_ns14(nd, version)
    % [plants, frf_data] = plants_ns14(nd)
      
    %       load(fullfile(PATHS.sysid, ['hysteresis/steps_hyst_model.mat']));
    %       hyst.rp = rp;
    %       hyst.wp = wp;
    %       hyst.r = r;
    %       hyst.w = w;
    %
    %       modelFit_file = fullfile(cs_root,'xy-stage', 'x-axis_sines_infoFourierCoef_9-11-2018-01.mat');
    if ~exist('version', 'var')
      version = 1;
    end
      
      if version == 1
        modelFit_file = fullfile(PATHS.sysid,'xy-stage', 'x-axis_sines_infoFourierCoef_9-11-2018-01.mat');
      elseif version == '5micron'
        modelFit_file = fullfile(PATHS.sysid,'xy-stage', 'x-axis_sines_infoFourierCoef_9-11-2018-01_5michyst.mat');
      else
        error('version number $d not recocnized');
      end
      
      load(modelFit_file, 'modelFit')
      if exist('nd', 'var') && ~isempty(nd)
        modelFit.models.Gvib.IOdelay = 0;
        modelFit.models.G_uz2stage.IOdelay = 0;
        modelFit.models.Gvib.InputDelay = nd;
        modelFit.models.G_uz2stage.InputDelay = nd;
      else
        modelFit.models.Gvib.IOdelay = 0;
        modelFit.models.G_uz2stage.IOdelay = 0;
        modelFit.models.Gvib.InputDelay = 9;
        modelFit.models.G_uz2stage.InputDelay = 9;
      end
      plants = modelFit.models;
      SYS = ss(modelFit.models.Gvib);
      
      %plants.gdrift = modelFit.models.gdrift;
      plants.gdrift_inv = 1/plants.gdrift;
      
      SYS = balreal(SYS);
      Nx = SSTools.getNxNu(SYS);
      T = diag(1./Nx)/10;
      SYS = ss2ss(SYS, T);
      PLANT = SYS;
      
      Nd = SYS.InputDelay;

      plants.sys_nodelay = SYS;
      
      SYS = absorbDelay(SYS);
      PLANT = absorbDelay(PLANT);
      
      plants.PLANT = PLANT;
      plants.SYS = SYS;
%       plants.hyst = hyst;
      plants.hyst = modelFit.models.hyst;
      plants.hyst_sat = modelFit.models.hyst_sat;
      plants.sys_recyc=SSTools.deltaUkSys(SYS);
      plants.sys_recyc_nodelay=SSTools.deltaUkSys(plants.sys_nodelay);
      plants.Nd = Nd;      
      
      
      frf_data = modelFit.frf;
    end

    
  end
  
  
  
end
