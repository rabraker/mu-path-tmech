function imshow_dataview(IM, thresh, ax1, ax2)
  axes(ax1)

  if length(thresh) ==2
    h_im = imagesc(ax1, IM, thresh);
  else
    h_im = imagesc(ax1, IM);
  end
  colormap('gray');
  ax1 = h_im.Parent;
  
  if exist('ax2', 'var')
    dcm = datacursormode(gcf);
    set(dcm, 'UpdateFcn', @(x, event_obj) local_callback_function(x, event_obj, ax1, ax2));
  end
  
end

function output_txt = local_callback_function(~,event_obj, ax1, ax2)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object

    if isa(event_obj.Target.UserData, 'function_handle')
        % check to see if the target line has the user data set. If so,
        % assume it will return the text to update the data-tip.
        output_txt = event_obj.Target.UserData();
    else % Fallback to default.
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
        if ~isempty(event_obj.Target.UserData) ...
            && isa(event_obj.Target.UserData, 'matlab.graphics.chart.primitive.Line')
          delete(event_obj.Target.UserData)
        end
        event_obj.Target.UserData = h_line;
        
        output_txt = {['x-pix: ',num2str(pos(1),4)],...
            ['y-pix: ',num2str(pos(2),4)],...
            ['Index: ', num2str(ind)],...
            ['value:  ', num2str(height)]};


    end
end

function datatip_txt = make_tip_txt(sz_loc, freq, Wn, zet, mag)
if abs(Wn - freq) < 1e-13
    unit = '[rad/s]';
else
    unit = '[Hz]';
end
    datatip_txt = {['Freq: ',num2str(freq,4), ' ',unit],...
        ['Mag: ',num2str(mag,4), ' [dB]'],...
        ['Wn: ',num2str(Wn,4), ' rad/s'],...
        ['damping: ',num2str(zet,4)],...
        ['Re: ', num2str(real(sz_loc))],...
        ['Im: +-', num2str(abs(imag(sz_loc)))]};
end