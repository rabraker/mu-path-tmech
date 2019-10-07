classdef AfmSisoSSController < AfmSisoController
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
%     G;  % System plant
%     Sens; % Sensitivity function
%     Hyr; 
%     Hy_rprime;
%     Loop;
%     plants; % All the plants
%     Ki;     % Integrator gain.
%     Dki;   % Integrator.
%     D;
%     M;
%     Dss;
%     Mss;
      K_lqr;
      L_dist;
      Nx;
      Nu;
      sys_recyc;
      sys_obsDist;
      
   end
   
   methods
       function self = AfmSisoSSController(G, G_recyc, G_obsdist, L_dist, K_lqr, varargin)
           self = self@AfmSisoController(G, varargin{:});
          
          [Sens, ~, Hy_rprime, ~, Loop, Dss, Mss] = ss_loops_delta_dist(G, G_recyc,...
              G_obsdist, K_lqr, L_dist);
          
          self.G = balreal(absorbDelay(ss(G)));
          self.sys_recyc = G_recyc;
          self.sys_obsDist = G_obsdist;
          self.Sens = Sens;
          self.Hy_rprime = Hy_rprime;
          
          self.Dss = Dss;
          self.Mss = Mss;
          self.K_lqr = K_lqr;
          self.L_dist = L_dist;
          
          [Nx, Nu] = SSTools.getNxNu(G_recyc);
          self.Nx = Nx;
          self.Nu = Nu;
          
          g_static = zpk([], [], 1, AFM.Ts);
          self.Ki = 0;
          self.Dki = g_static;
          self.D = g_static;
          
          self.Hyr = minreal(self.M * self.Hy_rprime);
          
       end
   end
end