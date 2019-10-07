
function output_txt = dataview_callback(self, event_obj, ax1, ax2)
% output_txt = dataview_callback(self, event_obj, ax1, ax2)
%
% use as: fun = @(event_obj)(dataview_callback(self, event_obj, ax1, ax2)
% where you have provided ax1 and ax2.
% 
% This should be used as the callback function when Using
% ImshowDataView.imshow().
%
% It will plot, in dashed red, the original, unconstructed CS data.
%
%
  pos = get(event_obj,'Position');
  
  ind = get(event_obj, 'DataIndex');
  
  yidx = pos(2);
  xidx = pos(1);
  
  dat = event_obj.Target.CData(yidx, :);
  height = dat(xidx);
  
  npix = size(self.pix_mask_uz,1);
  dat_og = self.pix_mat_uz(yidx, :);
  mask_row = self.pix_mask_uz(yidx, :);
  dat_og(mask_row~=1) = NaN; % matlab wont plot this.
  
  zmin = min(event_obj.Target.CData(:));
  zmax = max(event_obj.Target.CData(:));
  cla(ax2);
  plot(ax2, dat);
  hold(ax2, 'on')
  plot(ax2, xidx, dat(xidx), 'xk')
  plot(ax2, dat_og, '--r')
  grid(ax2, 'on');
  ylim(ax2, [zmin, zmax]);
  hold(ax1, 'on')
  
  dims = size(event_obj.Target.CData);
  h_line = plot(ax1, [1, dims(1)], [yidx, yidx], '-r');
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