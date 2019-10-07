classdef AFM
  % Static methods that provide constants associated with the nPoint stage and
  % z-axis parameters.
  
  methods (Static)
    function ts = Ts()
    % Ts = Ts() 
    % Standard sample time.
          
      ts = 40e-6;
    end
    
    function mic = volts2mic_xy()
    % mic = volts2mic_xy()
    % 
    % Unit conversion of volts 2 microns for xy-stage
      mic = 5;
    end
    
    function volt = mic2volt_xy()
    % volt = mic2volt_xy()
    % 
    % Unit conversion of microns to volts for xy-stage
        
      volt = 1/AFM.volts2mic_xy();
      
    end

    function mic = volts2mic_z()
    % mic = volts2mic_z()
    % 
    % Unit conversion of volts 2 microns for z-axis
        
      mic = 7/20; % 7 micron range in +-20 volts
    end

    function volts = mic2volts_z()
    % volts = mic2volts_z()
    % 
    % Unit conversion of volts 2 microns for z-axis
        
      volts = 1/AFM.volts2mic_z(); 
    end

    function nm = volts2nm_z()
    % nm = volts2nm_z()
    % 
    % Unit conversion of volts 2 nanometers for z-axis
        
      nm= AFM.volts2mic_z()*1000; 
      
    end

    function volts = nm2volts_z()
    % volts = nm2volts_z()
    % 
    % Unit conversion of volts 2 microns for z-axis
        
      volts = 1/AFM.volts2nm_z(); 
      
    end
    
    function map = nametrans()
    % map = nametrans()
    % 
    % A translation map between x, y, z corridinates to the ADC/DAC
    % number they are hooked up to.
        
      keySet ={ 'AI-0', 'AI-1', 'AI-2', 'AI-3',...
        'Input-0', 'Input-1', 'Input-2', 'Input-3',...
        'DAC control', 'Reference'};
      valSet = {'y_X', 'y_Y', 'y_Z', 'null',...
        'X', 'Y', 'Z', 'null', 'u_', 'ref_'};
      map = containers.Map(keySet, valSet);
    end
    
  end % Static Methods
  
end