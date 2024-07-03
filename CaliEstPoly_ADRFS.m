function [gEst,gammahat]= CaliEstPoly_ADRFS(t,T,PS,weight,Y,p)
% PURPOSE: compute the estimator of the ADRFS in the model: E{Y*(t)|T \in S}
% = g(t,gamma_S), where g(t,gamma_S) =
% gamma_0+gamma_1*t+...+gamma_{p-1}*t^{p-1}.
%--------------------------------------------------------------------------
% USAGE: [gEst,gammahat,PS]= CaliEstPoly_ADRFS(t,T,X,IS,weight,betahat,Y,p,K2)
% where: t is a vector of values of treatment where the prediction of the
% potential outcome are wanted (n x 1)
%        T is the vector of the treatment data (N x 1)
%        X is the confounding variable (N x r)
%        IS is the subpopulation condition for T; I(T in S) (N x 1)
%        Y is the vector of the observed outcome (N x 1)
%        weight is the estimated pi_0(T,X) (N x 1)
%        betahat is the estimated beta (r x 1)
%        p is the order of the polynomials in the model
%        K2 is the number of sieve basis used to estimate PS
%--------------------------------------------------------------------------
% RETURNS: gEst (n x 1) is the vector of regression estimator of
% g(t0,gammahat_S).    
%          gammahat (p x 1) is the vector of estimated gamma.
%          PS (N x 1) is the estimated P_S(X)
%--------------------------------------------------------------------------
% SEE ALSO:
% -------------------------------------------------------------------------
% References: 
%    
% -------------------------------------------------------------------------
% Written by:
%    Wei Huang
%    Lecturer
%    School of Mathematics and Statistics, The University of Melbourne
%--------------------------------------------------------------------------
% Last updated:
%    Feb 19, 2024.
% -------------------------------------------------------------------------

N = length(Y);
T = reshape(T,N,1);
Y = reshape(Y,N,1);
weight = reshape(weight,N,1);

%% Estimate gamma
T_poly_mat = repmat(T, 1, p).^repmat(0:(p-1), N, 1);     % N x p
w = weight.*PS;
W = diag(w);
TT = T_poly_mat'*W*T_poly_mat;

gammahat = TT\T_poly_mat'*W*Y;

g = @(theta,x_mat) x_mat*theta;

%ini = zeros(p,1);
%gT = @(theta) g(theta,T_poly_mat);

%fun = @(ini) Obj(ini,weight.*PS,Y,gT);

%gs = GlobalSearch('Display','off');
%options=optimset('Display','off');
%problem = createOptimProblem('fmincon','x0',ini,'objective',fun,'options',options);
%gammahat = run(gs,problem);

n = length(t);
t_mat = repmat(t,1,p).^repmat(0:(p-1),n,1);
gEst = g(gammahat,t_mat);

end

%function f = Obj(theta,weight,Y,gT)

%gval = gT(theta);
%diff = Y - gval;
%f = sum(weight.*diff.^2);

%end

