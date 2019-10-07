classdef ChannelMap
  properties
    names;
    x;
    y;
    ze;
    uz;
    met;
    friction;
  end
    methods
      function self = ChannelMap(idx_s)
        self.names = {'x', 'y', 'ze', 'uz', 'met', 'friction'};
        for k=1:min(length(self.names), length(idx_s))
          self.(self.names{k}) = idx_s(k);
        end
      end
    end
      
  
end
  
  