classdef ImshowDataView

  properties
    
  end
  
  methods (Static)
    function setup(Fig)
    % Call this method with a figure handle Fig to initialize. 
      dcm = datacursormode(Fig);
      fun = @(x, event_obj) ImshowDataView.dataView_callback_function(x, event_obj);
      set(dcm, 'UpdateFcn', fun);
               
    end
    
    function h = imshow(IM, thresh, ax1, ax2, fun)
      if nargin < 5
        fun = @(event_obj) ImshowDataView.plotter(event_obj, ax1, ax2);
      end
      
      if ~isempty(thresh)
        h_im = imagesc(ax1, IM, thresh);
      else
        h_im = imagesc(ax1, IM);
      end
      colormap('gray');
%       ax1 = h_im.Parent;
      
      % create a structure which contains the callback function
      % and an extra  field for the auxilary line handle.
      
      h_im.UserData =struct('fun', fun,...        
                            'line', []);
                          
      if nargout >0
        h = h_im;
      end
      
    end
    
    function output_txt = dataView_callback_function(~,event_obj)
    % Display the position of the data cursor
    % obj          Currently not used (empty)
    % event_obj    Handle to event object

      pos = get(event_obj,'Position');
      ind = get(event_obj, 'DataIndex');
      if isa(event_obj.Target.UserData.fun, 'function_handle')
        % check to see if the target line has the user data set. If so,
        % assume it will return the text to update the data-tip.
        output_txt = event_obj.Target.UserData.fun(event_obj);
      else % Fallback to default.
        yidx = pos(2);
        xidx = pos(1);
        height = event_obj.Target.CData(yidx, xidx);
        
        output_txt = {['x-pix: ',num2str(pos(1),4)],...
                      ['y-pix: ',num2str(pos(2),4)],...
                      ['Index: ', num2str(ind)],...
                      ['value:  ', num2str(height)]};
        
        
      end
    end    

    function output_txt = plotter(event_obj, ax1, ax2)
    % This function should get stored in UserData.
      pos = get(event_obj,'Position');

      ind = get(event_obj, 'DataIndex');
      
      yidx = pos(2);
      xidx = pos(1);
      
      dat = event_obj.Target.CData(yidx, :);
      height = dat(xidx);
      
      zmin = min(event_obj.Target.CData(:));
      zmax = max(event_obj.Target.CData(:));
      cla(ax2);
      plot(ax2, dat);
      hold(ax2, 'on')
      plot(ax2, xidx, dat(xidx), 'xk')
      
      grid(ax2, 'on');
      ylim(ax2, [zmin, zmax]);
      hold(ax1, 'on')
      
      dims = size(event_obj.Target.CData);
      h_line = plot(ax1, [1, dims(1)], [yidx, yidx], 'r');
      if ~isempty(event_obj.Target.UserData.line) ...
          && isa(event_obj.Target.UserData.line, 'matlab.graphics.chart.primitive.Line')
        delete(event_obj.Target.UserData.line)
      end
      event_obj.Target.UserData.line = h_line;
      
      output_txt = {['x-pix: ',num2str(pos(1),4)],...
                    ['y-pix: ',num2str(pos(2),4)],...
                    ['Index: ', num2str(ind)],...
                    ['value:  ', num2str(height)]};
    end    
  end
  
  
  
end




