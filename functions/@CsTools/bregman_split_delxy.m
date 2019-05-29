function [Ir] = bregman_split_delxy(I,E,epsy, epsx, mu, lam,gam,maxiter,tol)
% [Ir] = bregman_split_delxy(I,E,weight,mu,lamy,lamx,gam,maxiter,tol)
% Solve the CS Anylsis problem
% 
% min |\Psi x|_1  +epsx|Delx x|_1 + epsy|Dely x|_1 s.t. ||\Phi x - y||_2 < \sigma
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


  
  if isempty(epsx) || epsx==0
    do_x = false;
  else
    do_x = true;
  end
  
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
  d_y = zeros(n^2,1);
  d_x = zeros(n^2,1);
  
  b_y = zeros(n^2,1)+1;
  b_x = zeros(n^2,1)+1;
  
  
  y = I_vec(E_vec>0.5);
  
  
  % Build the matrix
  % A = mu *Phi^T*Phi + lambda * Del^T*Del + gamma*I.
  % We should only use the matrices Delx and Dely to form A. It is faster to
  % compute Dely*y and Delx*x based on indeces.
  %Dely = spalloc(n^2,n^2, 2*(n^2-n) );
  Delx = spalloc(n^2,n^2, 2*(n^2-1));
  
  j1 = [1:n^2]';
  j2 = mod([n:n^2+n-1]', n^2)+1;
  i = [1:n^2]';
  one = ones(n^2,1);
  Dely = sparse([i;i], [j1;j2], [-one; one], n^2, n^2) ;
%   for i = 1:n^2-n
%     Dely(i,i) = -1;
%     Dely(i,i+n) = 1;
% %   Dely(i,i) = -1;
% %   Dely(i, mod(i-1 + n, n^2)+1 ) = 1;
%   end
  fprintf('Finished building Dely\n')
  if do_x
    j1 = [1:n^2-1]';
    j2 = [2:n^2]';
    i1 = [1:n^2-1]';
    one = ones(n^2-1,1);
    Delx = sparse([i1; i1; n^2], [j1; j2; n^2], [-one; one; 0]);
    % for i = 1:n^2-1
    %  Delx(i, i)   = -1;
    %  Delx(i, i+1) = 1;
    % end
    fprintf('Finished building Delx\n')
    Dely = Dely*epsy;
    Delx = Delx*epsx;
    A = mu*(E_mat'*E_mat) + 1 * lam*(Dely'*Dely) +...
      1 * lam*(Delx'*Delx)      +      gam * I;
    

  else
    Dely = Dely*epsy;
    A = mu*(E_mat'*E_mat) + 1 * lam*(Dely'*Dely) + gam * I;
    
  end
%    A = mu*(E_mat'*E_mat) +  gam * I;

  
  L1 = ichol(A); % preconditioner
  L1T = L1';

  
  uk = u;
  uk_1 = uk;
  yk = y;
  
  max_outer_iter = 1000;
  tloop = zeros(maxiter*max_outer_iter,1);
  tdct = zeros(maxiter*max_outer_iter,1);
  tsolve = zeros(maxiter*max_outer_iter,1);
  tshrink = zeros(maxiter*max_outer_iter,1);
  trem   = zeros(maxiter*max_outer_iter,1);
  
  j = 1;
  %fprintf('| outer-iter |    || Phi*x-y||_2     |\n');
  fprintf('loop time   | Lin Solver time |    dct time  |  shrink  |  remainder    | total iter |  out-iter  | Residual |\n')
  fprintf('-------------------------------------------------------------------------------------------------\n')
  for outer_iter = 1:max_outer_iter
    uk = uk_1;
    %fprintf('loop time   | Lin Solver time |   dct time |   remainder    |\n')
    for i = 1:maxiter
      t_loop_start = tic();
      Ety = y_full_sparse;
      Ety(pix_idx) = yk;
      if do_x
        rhs = mu*Ety + lam*Delx'*(d_x - b_x) +...
          lam*Dely'*(d_y-b_y)+gam*idct(w-b_w);
      else
        rhs = mu*Ety + lam*Dely'*(d_y-b_y)+gam*idct(w-b_w);
      end
      

      t_solve_start = tic();
      %uk_1 = A\rhs; % Inv(A)*rhs
      % use conjugate gradient with pre-conditioner instead. WAY faster.
      [uk_1, ~, ~] = pcg(A, rhs, 1e-6, 10); %, L1, L1T);
      t_solve = toc(t_solve_start);
      
      t_dct_start = tic();
      dct_uk1 = dct(uk_1);
      t_dct = toc(t_dct_start);
      
      t_shrink_start=tic();
      if do_x
        Dely_uk = Dely*uk_1;
        Delx_uk = Delx*uk_1;
        % How is this the same shrink(Del*u_k + b_y, lambda) ?????
        s = sqrt( (Dely_uk + b_y).^2 + (Delx_uk + b_x).^2 );
        d_x = max(s-1/lam, 0).*(Delx_uk + b_x)./s;
        d_y = max(s-1/lam, 0).*(Dely_uk + b_y)./s;
      else
        Dely_uk = Dely*uk_1;
        % How is this the  shrink(Del*u_k + b_y, lambda) ?????
        s = abs(Dely_uk + b_y);
        d_y = max(s-1/lam, 0).*(Dely_uk + b_y)./s;
        
      end
      t_shrink = toc(t_shrink_start);

      w = shrink(dct_uk1 + b_w, 1/gam);
      if do_x
       b_x = b_x + Delx*uk_1 - d_x;
      end
      b_y = b_y + Dely*uk_1 - d_y;
      
      b_w = b_w + dct_uk1 - w;
      
      % ---------- Bench mark ----------
      tloop(j) = toc(t_loop_start);
      
      tsolve(j) = 100*t_solve/tloop(j);
      tdct(j) = 100*t_dct/tloop(j);
      tshrink(j) = 100*t_shrink/tloop(j);
      trem(j) = 100*(tloop(j) - t_dct - t_solve - t_shrink)/tloop(j);
      %fprintf('%.6e  %.6e      %.6e    %.6e \n', tloop(i), tsolve(i), tdct(i), trem(i));
      j = j + 1;
    end
 
    % -------- update -------
    yk_est = uk_1(pix_idx);
    yk = yk + y - yk_est; %Phi*x
    
    Nrm = norm(yk_est-y);
    if mod(j-1, 100) == 0
      fprintf('%.4e  (%.2f %%,  %03.2f [s])   %.2f %%       %.2f %%      %.2f %%       %05d        %05d     %.4e\n',...
        sum(tloop), sum(tsolve)/j, 0.01*sum(tsolve),...
        sum(tdct)/j, sum(tshrink)/j, sum(trem)/j, j, outer_iter, Nrm);
      %fprintf('| outer-iter |    || Phi*x-y||_2     |\n');
    end
    if mod(j-1, 5000) == 0
      fprintf('loop time   | Lin Solver time |    dct time  |  shrink  | remainder    | total iter |  out-iter  | Residual |\n')
      fprintf('-------------------------------------------------------------------------------------------------\n')
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
  % N.b. The notation in their paper is WRONG. On page 1 they define |x| as the
  % l1 norm of x. That is NOT what is meant by |x| in equation 3.10 etc, which
  % you can only figure out by inspecting their code.
  xs = (x./abs(x)).* max(abs(x) - lambda, 0);
end