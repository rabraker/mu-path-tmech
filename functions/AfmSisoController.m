classdef AfmSisoController
% This class is just a container for the siso axis controllers, to enforce
% consistent naming, rather than a struct.
% 
%            R-prime
%              ^  (Mss is sort of fictitious with state-space, and this is what
%              |   we usually want from r-prime)
%              |
%    R-->[ M ]--> [M_ss ]--->O---->[Dki*D ]-->[ G ]----+---> y
%                            |                      |
%                            +-----[ Dss ]<---------+
% 
% -- For transfer function designs, Dss=Mss=1. Dki is an integrator and D is likely
%    a notch etc.
% -- For state space designs, D=1.
% 

   properties
    G;  % System plant
    Sens; % Sensitivity function
    Hyr; 
    Hy_rprime;
    Loop;
    plants; % All the plants
    Ki;     % Integrator gain.
    Dki;   % Integrator.
    D;
    M;
    Dss;
    Mss;
   end
   
   methods
       function self = AfmSisoController(G, varargin)
          p = inputParser();
          g_static = zpk([], [], 1, AFM.Ts);
          
          p.addParameter('plants', []);
          p.addParameter('Ki', 0);
          p.addParameter('Dki', g_static);
          p.addParameter('D', g_static);
          p.addParameter('M', g_static);
          p.addParameter('Dss', g_static);
          p.addParameter('Mss', g_static);
          
          p.parse(varargin{:});
          self.plants = p.Results.plants;
          self.Ki = p.Results.Ki;
          self.Dki = p.Results.Dki;
          self.D = p.Results.D;
          self.M = p.Results.M;
          self.Dss = p.Results.Dss;
          self.Mss = p.Results.Mss;

          self.G = absorbDelay(ss(G));
          
          self.Loop = self.Dki * self.D * self.G *self.Dss;
%           self.Sens = balreal(minreal(1/(1+self.Loop)));
          self.Sens = feedback(1, self.Loop);
%           self.Hy_rprime = minreal(self.Mss * (self.Dki * self.D * self.G) * self.Sens);
          self.Hy_rprime = minreal(self.Mss * feedback(self.Dki*self.D*self.G, self.Dss));
          
          self.Hyr = minreal(self.M * self.Hy_rprime);
          
       end
   end
end