function remove_ticks(ax)
   % Remove the x and y axis ticks labels.
   for k=1:length(ax)
       set(ax(k), 'XTickLabel', [])
       set(ax(k), 'yTickLabel', [])
   end
end