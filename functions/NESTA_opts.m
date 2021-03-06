%               This is a structure that contains additional options,
%               some of which are optional.
%               The fieldnames are case insensitive.  Below
%               are the possible fieldnames:
%               
%           .U and opts.Ut - Analysis/Synthesis operators as function
%               handles.
%           .normU - if opts.U is provided, this should be norm(U)
% 
%           .alpha_v - weight for vertical TV (default 0)
%           .alpha_h - weight for horizontal TV (default 0)
% 
%           .MaxIntIter - number of continuation steps.
%                 default is 5
%           .maxiter - max number of iterations in an inner loop.
%                 default is 10,000
%           .muf - The desired value of mu at the last continuation step.
%                 A smaller mu leads to higher accuracy.
%           .sigma - l2 error bound.  This enforces how close the variable
%               must fit the observations b, i.e. || y - Ax ||_2 <= delta
%               If delta = 0, enforces y = Ax
%               Common heuristic: delta = sqrt(m + 2*sqrt(2*m))*sigma;
%               where sigma=std(noise).
%               opts.TolVar - tolerance for the stopping criteria, i.e.,
%                             when the relative change in the objective 
%                             function is less than TolVar.
%           .opts.Verbose - if this is 0 or false, then very
%                   little output is displayed.  If this is 1 or true,
%                   then output every iteration is displayed.
%                   If this is a number p greater than 1, then
%                   output is displayed every pth iteration.
%            .opts.errFcn - if this is a function handle,
%                   then the program will evaluate opts.errFcn(xk)
%                   at every iteration and display the result.
%                   ex.  opts.errFcn = @(x) norm( x - x_true )

function opts = NESTA_opts(varargin)
    
   p = inputParser();
   U_default_fun = @(x) x;
   p.addParameter('Verbose', true);
  
  
   p.addParameter('MaxIntIter',5);
   p.addParameter('TypeMin','L1');
   p.addParameter('TolVar',1e-5);
   
   p.addParameter('xplug',[]);
   p.addParameter('normU',[]);  % so we can tell if it's been set
   
   p.addParameter('errFcn',[]);

   p.addParameter('U', U_default_fun);
   p.addParameter('Ut', U_default_fun);
   
   p.addParameter('alpha_v', 0);
   p.addParameter('alpha_h', 0);
   p.addParameter('mu', 1e-5);
   p.addParameter('sigma', 0.01);
   
   p.parse(varargin{:});
   
   opts = p.Results;
   
   if ~isequal(opts.U, U_default_fun)
       opts.U_userSet = true;
   else
       opts.U_userSet = false;
       opts.normU = 1;
   end
end
