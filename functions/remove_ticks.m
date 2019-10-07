function remove_ticks(ax)
   % Remove the x and y axis ticks labels.
   for k=1:size(ax, 1)
       for j=1:size(ax, 2)
           set(ax(k, j), 'XTickLabel', [])
           set(ax(k, j), 'yTickLabel', [])
       end
   end
end