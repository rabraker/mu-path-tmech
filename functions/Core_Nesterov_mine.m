function [xk,niter] = Core_Nesterov_mine(...
            A,At,b,mu,delta,opts)
% [xk,niter,residuals,outputData,opts] =Core_Nesterov(A,At,b,mu,delta,opts)
%
% Solves a L1 minimization problem under a quadratic constraint using the
% Nesterov algorithm, without continuation:
%
%     min_x || U x ||_1 s.t. ||y - Ax||_2 <= delta
% 
% If continuation is desired, see the function NESTA.m
%
% The primal prox-function is also adapted by accounting for a first guess
% xplug that also tends towards x_muf 
%
% The observation matrix A is a projector
%
% Inputs:   A and At - measurement matrix and adjoint (either a matrix, in which
%               case At is unused, or function handles).  m x n dimensions.
%           b   - Observed data, a m x 1 array
%           muf - The desired value of mu at the last continuation step.
%               A smaller mu leads to higher accuracy.
%           delta - l2 error bound.  This enforces how close the variable
%               must fit the observations b, i.e. || y - Ax ||_2 <= delta
%               If delta = 0, enforces y = Ax
%               Common heuristic: delta = sqrt(m + 2*sqrt(2*m))*sigma;
%               where sigma=std(noise).
%           opts - see NESTA_opts.
%  Outputs:
%           xk  - estimate of the solution x
%           niter - number of iterations
%           residuals - first column is the residual at every step,
%               second column is the value of f_mu at every step
%           opts - the structure containing the options that were used      
%
% Written by: Jerome Bobin, Caltech
% Email: bobin@acm.caltech.edu
% Created: February 2009
% Modified: May 2009, Jerome Bobin and Stephen Becker, Caltech
% Modified: Nov 2009, Stephen Becker
%
% NESTA Version 1.1
%   See also NESTA

    maxiter = 10000; %setOpts('maxiter',10000,0);

    if delta <= 0 
      error('delta must be greater than zero'); 
    end

    Atb = At(b);

    if isempty(opts.xplug)
        error('Must provide initial guess in opts.xplug.')
    end
    
    %---- Initialization
    N = length(opts.xplug);
    wk = zeros(N, 1); 
    xk = opts.xplug;


    %---- Init Variables

    % Evidently, ||D||_2 ~= sqrt(8), and ||Dv||_2 ~= 2, ||Dh||_2 ~=2
    % In general, we take L=||W||_2^2, so when we multiply alp*D,
    % we need alp^2.
    % Lmu = (opts.normU + alp_v^2*4 + alp_h^2*4)/mu; %Lmu;
    Lmu = opts.normU/mu;
    n = floor(sqrt(N));
    
    if opts.alpha_v > 0
        alp_v = opts.alpha_v;
        Lmu = Lmu + (alp_v^2 * 4)/mu;

        Dv = alp_v*spdiags([reshape([-ones(n-1,n); zeros(1,n)],N,1) ...
            reshape([zeros(1,n); ones(n-1,n)],N,1)], [0 1], N, N);
    else
        Dv = sparse([]);
    end
    
    if opts.alpha_h > 0
        alp_h = opts.alpha_h;%.125;
        Lmu = Lmu + (alp_h^2 * 4)/mu;
        Dh = alp_h*spdiags([reshape([-ones(n,n-1) zeros(n,1)],N,1) ...
            reshape([zeros(n,1) ones(n,n-1)],N,1)], [0 n], N, N);
    else
        Dh = sparse([]);
    end
    
    residuals = zeros(maxiter, 2);
    if opts.Verbose > 0
        st = sprintf('Iter |     fmu     |  Rel. Vartn fmu |  Residual  |');
        fprintf('%s\n%s\n', st, repmat('-', 1, length(st)));
    end
    % ---------------------------- MAIN ITERATION -------------------------- %
    max_means = 10;
    fmu_s = zeros(max_means+1, 1);
    for k = 0:maxiter-1

       %---- Dual problem
       [df, fx] = Perform_L1TVanis_Constraint(xk, mu, opts.U, opts.Ut, Dv, Dh);

       %---- Primal Problem

       %--------------- Update yk --------------------------------
       yk = nesta_project(xk, df, b, Lmu, A, At, Atb, delta);

       %-------------- Update zk ----------------------------------
       apk = 0.5*(k+1);
       tauk = 2/(k+3); 

       wk =  apk*df + wk;

       zk = nesta_project(opts.xplug, wk, b, Lmu, A, At, Atb, delta);

       %----------------- Update xk --------------------------
       xk = tauk*zk + (1 - tauk)*yk;
       
       %------------------- Check for exit ------------------------
       % but wait to do so, so we can print.
       n_mean = min(max_means, k)+1;
       % (option 1) look at the relative change in function value
       % eq. 3.13. Always infinite here for k==0. Unclear what is intended in
       % the first iteration from 3.13, since sum is over x_{k - \ell}.
       rel_del_fmu = abs(fx - mean(fmu_s(1:n_mean)))/mean(fmu_s(1:n_mean));

       fmu_s(n_mean) = fx;
       if n_mean > max_means
           fmu_s(1:end-1) = fmu_s(2:end);
       end
       
       %--- display progress if desired
       if ~mod(k+1, opts.Verbose )
         % fprintf('Iter: %3d  ~ fmu: %.3e ~ Rel. Variation of fmu: %.2e ~
         % Residual: %.2e',...
         fprintf('%3d     %.3e       %.2e        %.2e\n',...
                 k+1,fx,rel_del_fmu,residuals(k+1,1) ); 
       end
       if abs(fx)>1e20 || abs(residuals(k+1,1)) >1e20 || isnan(fx)
         error('Nesta: possible divergence or NaN.  Bad estimate of ||A''A||?');
       end

       if rel_del_fmu <= opts.TolVar
           break;
       end
      
       
    end

    niter = k+1; 

end %% end of main Core_Nester  ov routine

function [vk] = nesta_project(xx, g, b, Lmu, A, At, Atb, delta)
% Performs the projections for either yk or zk.
% 
% For yk:
% % yk = Argmin_x Lmu/2 ||x - xk||_l2^2 + <df,x-xk> s.t. ||b-Ax||_l2 < delta
% Let xp be sqrt(Lmu) (x-xk), dfp be df/sqrt(Lmu), 
% bp be sqrt(Lmu)(b- Axk) and deltap be sqrt(Lmu)delta
% Then
% 
% yk =  xk + 1/sqrt(Lmu) Argmin_xp 1/2 || xp ||_2^2 + <dfp,xp> s.t. || bp - Axp ||_2 < deltap
% 
% For zk:
% zk = Argmin_x Lmu/2 ||b - Ax||_l2^2 + Lmu/2||x - xplug ||_2^2 + <wk,x-xk> 
%   s.t. ||b-Ax||_l2 < delta
%

  q = xx - 1/Lmu*g; % eq. 3.7 or 3.12

  Aq = A( q );
  AtAq = At( Aq );

  lambda = max(0, Lmu*(norm(b - Aq)/delta - 1));
  gamma = lambda/(lambda + Lmu);
  vk = lambda/Lmu*(1 - gamma)*Atb + q - gamma*AtAq;
  
end

%%%%%%%%%%%% PERFORM THE L1+TV CONSTRAINT %%%%%%%%%%%%%%%%%%
function [df,fx] = Perform_L1TVanis_Constraint(xk,mu,U,Ut, Dv, Dh)
    df = 0;
    if ~isempty(Dh)
        Dhx = Dh*xk;
        w_h = max(mu, abs(Dhx));
        uh = Dhx ./ w_h;
        df = df + Dh'*uh;
    else
       Dhx = [];
       uh = [];
    end
    
    if ~isempty(Dv)
        Dvx = Dv*xk;
        w_v = max(mu, abs(Dvx));
        uv = Dvx ./ w_v;
        df = df + Dv'*uv;
    else
       Dvx = [];
       uv = [];
    end
    
    Ux = U(xk);
    w_U = max(mu, abs(Ux));
    uU = Ux./ w_U;
    
    df = df + Ut(uU);
    u = [uh;uv; uU];
  
    fx = u' * [Dhx; Dvx; Ux] - (mu/2) * norm(u)^2;
  
%   fx_l1 = uk_U'*Uxk - mu/2*norm(uk_U)^2;  
%   df = Dh'*uh + Dv'*uv + Ut(uU);

end
  
%%%%%%%%%%%% PERFORM THE L1+TV CONSTRAINT %%%%%%%%%%%%%%%%%%
function [df,fx] = Perform_L1TV_Constraint(xk,mu,U,Ut, Dv, Dh, D)
    Dhx = Dh*xk;
    Dvx = Dv*xk;
    
    % tvx = sum(sqrt(abs(Dhx).^2+abs(Dvx).^2));
    w = max(mu, sqrt(abs(Dhx).^2 + abs(Dvx).^2));
    uh = Dhx ./ w;
    uv = Dvx ./ w;
    if sum(abs(uh)) == 0
        len_u = numel(uv);
    elseif sum(abs(uv)) == 0
        len_u = numel(uh);
    else
        len_u = numel(uv) + numel(uh);
    end
%     len_u = 1;
    u = [uh;uv];
    fx_tv = u'*D*xk - (mu/2) * (1/len_u)*(u'*u);
    df_tv = D'*u;

  uk = U(xk);
  fx = uk;
    
  uk = uk./max(mu, abs(uk));
  fx_l1 = real(uk'*fx - mu/2*norm(uk)^2);
    
  df_l1 = Ut(uk);
  
  fx = fx_l1 + fx_tv;
  df = df_tv + df_l1;
end

function [df,fx,val,uk] = Perform_L1_Constraint(xk,mu,U,Ut)

  uk = U(xk);
  fx = uk;
    
  uk = uk./max(mu,abs(uk));
  val = real(uk'*fx);
  fx = real(uk'*fx - mu/2*norm(uk)^2);
    
  df = Ut(uk);
end


    