classdef CsTools
  
  
  methods (Static)
    xp = l1qc_logbarrier(x0, A, At, b, epsilon, lbtol, mu, ...
                              cgtol, cgmaxiter, verbose);
    [xp, up, niter, cgtot_iter] = l1qc_newton(x0, u0, A, At, b, epsilon, tau, newtontol,...
      newtonmaxiter, cgtol, cgmaxiter, Tii, verbose)
    [x, res, iter] = cgsolve(A, b, tol, maxiter, verbose, x0);
    
    
    vpix = pixmat2vec(Mat);
    
    Mpix = pixvec2mat(vpix, nrows, ncol);
    
    [pix_mask, pix_idx, npaths] = mu_path_mask(mupath_len, n, m, samplingRatio, repeat_sampling);
    
    [Ir] = bregman_split_delxy(I,E,epsy, epsx, mu, lam,gam,maxiter,tol)
    
    [Ir] = bregman_split(I,E, mu, gam,maxiter,tol);
    

    function [ b ] = Afun_dct(eta, pix_idx)
    % Given the CS equation
    % b = E * M * eta
    % where E is the subsampling matrix and M is the idct.
      b = idct(eta);
      b = b(pix_idx);
    end
    function [ y ] = E_fun1(x, pix_idx)
    % Given the CS equation
    % y = E * x
    % where E is the subsampling matrix and M is the idct.
      y = x(pix_idx);
    end
    function [ x ] = Et_fun1(y, pix_idx, N, M)
    % Computes the adjoint of the CS equation
    % x = E^T y
    % where E is the subsampling matrix and M is the idct. That is
    %
    % eta = M'*E'*b
    %
    % N2 is the length eta.
      x = zeros(N*M, 1);
      x(pix_idx) = y;
    end

    function [ eta ] = Atfun_dct(b, pix_idx, N, M)
    % Computes the adjoint of the CS equation
    % b = E * M * eta
    % where E is the subsampling matrix and M is the idct. That is
    %
    % eta = M'*E'*b
    %
    % N2 is the length eta.
      
      eta = zeros(N*M, 1);
      eta(pix_idx) = b;
      eta = dct(eta);
      
    end

    function [ eta ] = Atfun_dct2(b, pix_idx, N, M)
    % Computes the adjoint of the CS equation
    % b = E * M * eta
    % where E is the subsampling matrix and M is the 2D-idct. That is
    %
    % eta = M'*E'*b
    %
    % N2 is the length eta.
    
      z_sparse = zeros(N*M,1);
      z_sparse(pix_idx) = b;
      Z_sparse_mat = CsTools.pixvec2mat(z_sparse, N);
      
      Eta_mat = dct2(Z_sparse_mat);
      
      eta = CsTools.pixmat2vec(Eta_mat);
    end
    
    function [b] = Afun_dct2(eta, pix_idx, N, M)
    % Computes the adjoint of the CS equation
    % b = E * M * eta
    % where E is the subsampling matrix and M is the 2D-idct. That is
    %
    
      Eta_mat = CsTools.pixvec2mat(eta, N);
      B_mat = idct2(Eta_mat);
      b = CsTools.pixmat2vec(B_mat);
      b = b(pix_idx);

    end
    function [ x ] = Utfun_dct2(z, N)
      Z_mat = CsTools.pixvec2mat(z, N);
      X_mat = dct2(Z_mat);
      x = CsTools.pixmat2vec(X_mat);
    end
    
    function [z] = Ufun_dct2(x, N)
      X_mat = CsTools.pixvec2mat(x, N);
      Z_mat = idct2(X_mat);
      z = CsTools.pixmat2vec(Z_mat);

    end
  end
  
end