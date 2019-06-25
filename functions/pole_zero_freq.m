function [pl_ws, zr_ws, pl_zet, zr_zet] = pole_zero_freq(G)
   pls = pole(G);
   zrs = tzero(G);
   
   
   if ~(G.Ts > 0)
       error('This function only works with discrete time models');
   end
   
   Ts = G.Ts;
   sloc_pl = -log(pls)/Ts;
   sloc_zr = -log(zrs)/Ts;
   pl_ws = abs(sloc_pl)/2/pi;
   zr_ws = abs(sloc_zr)/2/pi;
   
   pl_zet = real(pls)./pl_ws;
   zr_zet = real(zrs)./zr_ws;
end