function [Ir] = bregman_split(I,E, mu, gam,maxiter,tol)
% [Ir] = bregman_split(I,E, mu, gam,maxiter,tol)
% Solve the CS Anylsis problem
% 
% min |\Psi x|_1  s.t. ||\Phi x - y||_2 < \sigma
%
% through Bregman Splitting.
% 
% Arguments
% ---------
%   I - sampled image with size n*m
%   E - sample matrix with size n*m, 1 if sampled, 0 otherwise
%   mu,gam - tuning vaiables - they are turned out to be very
%   maxiter - max # of iteration
%   tol - allowed error
% 
% Outputs
% ------
%   Ir - reconstrcuted image
% 


  
  [n, m] = size(I);
  
  I_vec = CsTools.pixmat2vec(I);
  E_vec = CsTools.pixmat2vec(E);
  pix_idx = find(E_vec==1);
  I = speye(n^2);
  E_mat = I(pix_idx, :);
  
  u = I_vec.*E_vec;

  y_full_sparse = zeros(n^2, 1);
  w = zeros(n^2,1);
  b_w = zeros(n^2,1);

  y = I_vec(E_vec>0.5);
  
  
  % Build the matrix
  % A = mu *Phi^T*Phi + lambda * Del^T*Del + gamma*I.

   A = mu*(E_mat'*E_mat) +  gam * I;

  yk = y;
  
  max_outer_iter = 1000;

  j = 1;
  %fprintf('| outer-iter |    || Phi*x-y||_2     |\n');
  fprintf(' | out-iter  | Residual |\n')
  fprintf('----------------------------\n')
  for outer_iter = 1:max_outer_iter
    
    %fprintf('loop time   | Lin Solver time |   dct time |   remainder    |\n')
    for i = 1:maxiter
      Ety = y_full_sparse;
      Ety(pix_idx) = yk;
      
      rhs = mu*Ety + .1*gam*idct(w-b_w);

      %uk_1 = A\rhs; % Inv(A)*rhs
      % use conjugate gradient with pre-conditioner instead. WAY faster.
      [uk_1, ~, ~] = pcg(A, rhs, 1e-8, 100); %, L1, L1T);
     
      dct_uk1 = 10*dct(uk_1);

      w = shrink(dct_uk1 + b_w, 1/gam);
      b_w = b_w + dct_uk1 - w;
      
    end
 
    % -------- update -------
    yk_est = uk_1(pix_idx);
    yk = yk + y - yk_est; %Phi*x
    
    Nrm = norm(yk_est-y);
    if mod(j-1, 100) == 0
      fprintf('     %05d     %.4e\n',  outer_iter, Nrm);
    end

    if 1
      Ir = CsTools.pixvec2mat(uk_1, n);
      figure(5)
      imagesc(Ir)
      colormap('gray')
      drawnow();
    end
    if Nrm  < tol
      break
    end
    
  end
  
  % only compute this once we're done.
  Ir = CsTools.pixvec2mat(uk_1, n);
  
end


function [xs] = shrink(x, lambda)
  
%   nrm_x = norm(x, 1);
  xs = (x./abs(x)).* max(abs(x) - lambda, 0);
end