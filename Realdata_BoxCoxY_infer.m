clear all;
clc;

addpath '/Users/zyl/Library/CloudStorage/OneDrive-Personal/Non_CTE20220226 _2/Matlab Code/28-April-2024'
addpath '/Users/zyl/Library/CloudStorage/OneDrive-Personal/Non_CTE20220226 _2/Matlab Code/real-data'

tic;                                           % start timing
start_date = now;                              % start date
disp(['Job-Realdata started: ', datestr(start_date)]);  % display start date

%%

% load('Realdata_BoxCox_CV_p1u3p2u3k2u4_estu-02-Jun-2024.mat');
%filename=sprintf('Realdata_CV_infer_p1u3p2u3k2u4-%s.mat',datetime('today'));

 % load('Realdata_BoxCox_CV_tc0_p1u3p2u3k2u3_estu-02-Jun-2024.mat');
 % filename=sprintf('Realdata_CV_infer_tc0_p1u3p2u3k2u3-%s.mat',datetime('today'));
 
 load('Realdata_BoxCox_CV_tc20less_p1u3p2u3k2u4_estu-02-Jun-2024.mat');
 filename=sprintf('Realdata_CV_infer_tc20less_p1u3p2u3k2u4-%s.mat',datetime('today'));


%% Step 1: Set-up
%B = 500;              % # of bootstrap samples  
p = 3;                 % link function is quadratic if p = 3.
alpha = 0.05;          % nominal size

% options_fminsearch = optimset('Display','off');              % options for fminsearch
% options_fsolve = optimoptions('fsolve','Display','off');     % options for fsolve

% fix seed so that we always get the same result
defaultStream = RandStream.getGlobalStream;
defaultStream.reset(3);

%% Step 2: Read Data
% load zip_ads_v4.txt
% zip = zip_ads_v4(:,1); % zip code
% ads = zip_ads_v4(:,2); % # of political advertisements aired in each zipcode (treatment of interest)
% pop = zip_ads_v4(:,3); % population
% pct_over65 = zip_ads_v4(:,4); % percent over age 65
% median_income = zip_ads_v4(:,5); % median household income
% pct_black = zip_ads_v4(:,6); % percent black
% pct_hispanic = zip_ads_v4(:,7); % percent hispanic
% pop_density = zip_ads_v4(:,8); % population density
% pct_colgrads = zip_ads_v4(:,9)/100; % percent college graduates
% can_commute = zip_ads_v4(:,10); % binary indicator of whether it is possible to commute to the zip code from a competitive state
% contribution = zip_ads_v4(:,11); % total amount of campaign contribution 
% log_pop = log(pop); % log of population
% log_median_income = log(median_income + 1); % some elements of median_income are 0
% log_pop_density = log(pop_density + 1); % log of population density
% clear median_income pop_density
% Y = contribution;
% T = ads; % treatment
% X = [log_pop, log_pop_density, log_median_income, pct_over65,...
% pct_hispanic, pct_black, pct_colgrads,can_commute];
% clear log_pop log_pop_density log_median_income pct_over65 pct_hispanic pct_black pct_colgrads can_commute
% N = length(Y); % sample size
% data = [Y, T, X];
% dimX = size(X,2);
% clear contribution
% 
% lambda = FindLamBoxCox(Y); % find the best Box-Cox transformation parameter for outcome
% Y = BoxCox(Y,lambda); % outcome 
% Y = Y-min(Y);
% T = log(T+1);
% Nt=100;
% t0 = min(T):(max(T)-min(T))/(Nt-1):max(T); 
% 
% IS = double(T>= log(tc+1));



%% Step 2: Empirical analysis (point estimation)
% Y = log(contribution + 1);                     % outcome 
% T = log(ads + 1);                              % treatment
% t = (min(T) : 0.05 : max(T))';
% X = [log_pop, log_pop_density, log_median_income, pct_over65, pct_hispanic, pct_black, pct_colgrads, can_commute];    % covariates
% dimX = size(X,2);                              % dimension of covariates
% N = size(zip,1);                               % sample size

%% Obtain Estimator
    % %[p1(j),p2(j),K2(j)] = GCV_K1K2v2(T,X,Y,IS,p); 
    % [p1,p2,K2] = CV_K1K2(T,X,Y,IS,p);
    % [weight, betahat, ~] =  get_weight_semiclosedv3(T,X,T,X,N,p1,p2);
    % 
    % PShat= PS_estv2(X,betahat,IS,K2);
    % [~,gammahat]=CaliEstPoly_ADRFS(T,T,PShat,weight,Y,p);

    
    %% Step 3: Empirical analysis (interval estimation)
    
    % Obtain undersmoothing parameters
    % p1u = p1+1;
    % p2u = p2+1;
    % K2u = K2+1;
    % 
    % [weightu, betau, lambdahatu] =  get_weight_semiclosedv3(T,X,T,X,N,p1u,p2u);
    % PSu= PS_estv2(X,betau,IS,K2u);
    % [gTu,gammahatu]=CaliEstPoly_ADRFS(T,T,PSu,weightu,Y,p);

    %p1u = p1(j);
    %p2u = p2(j);
    %K2u = K2(j);
    %weightu = weight(:,j);
    %hatbeta = betahat(:,j);
    %PSu = PS(:,j);
    %gammahatu = gammahat(:,j);
    %hatgamma = gammahatu;
%%
    Z = X*betau;
[uii,par2_uii,par22_uii]...
            = BasisFUNv8(T,Z,p1u,p2u);

    %[uii,uij,par2_uii,par2_uij,par22_uii,par22_uij]...
     %       = BasisFUNv4(T,Z,N,p1u,p2u);
%%
    %orthogonomalize uK1 basis matrix
    Ku = size(uii,2);
    %exp_uu = (1/N) * (uii' * uii);               % K x K
    %exp_uu_half = chol(exp_uu, 'lower');
    %exp_uu_half_inv = exp_uu_half \ eye(Ku);
    %uii_std = uii * exp_uu_half_inv'; %N x K

    %uij_std = uij * exp_uu_half_inv'; %N^2 x K
    %par2_uii_std = par2_uii * exp_uu_half_inv'; %N^2 x K
    %par2_uij_std = par2_uij * exp_uu_half_inv'; %N^2 x K
    %par22_uii_std = par22_uii * exp_uu_half_inv'; %N^2 x K
    %par22_uij_std = par22_uij * exp_uu_half_inv'; %N^2 x K

    lambdahatu = reshape(lambdahatu,Ku,1); %Ku x 1

    %orthogonomalize vK2 basis matrix
    vZ = Z.^(0:K2u);
    dvZ = [zeros(N,1),Z.^(0:K2u-1).*repmat((1:K2u),[N,1])];
    
    Kv = size(vZ,2);
    exp_vv = vZ' * vZ;               % K x K
    %exp_vv_half = chol(exp_vv, 'lower');
    %exp_vv_half_inv = exp_vv_half \ eye(Kv);
    %vZ_std = vZ * exp_vv_half_inv'; %N x K

    %dvZ_std = dvZ * exp_vv_half_inv'; %N x K

    hatmu = exp_vv\vZ'*IS;
    hatmu = reshape(hatmu,Kv,1);

    T_poly_mat = repmat(T, 1, p).^repmat(0:(p-1), N, 1);     % N x p

    [indT,indZ] = ndgrid((1:N)', (1:N)');
    indT = indT(:);
    indZ = indZ(:);
%%
    % More auxiliary quantities 
    ulambda1ii = uii*lambdahatu-1;
    %ulambda1ij =  ulambda1ii(indT,:);
    par2ulambda = par2_uii*lambdahatu;
    %par2ulambdaij = par2ulambda(indT,:);
    par22ulambda = par22_uii*lambdahatu;
    ulambda1 = uii*lambdahatu-1;
    vmu = vZ*hatmu;
    %vmuij = vmu(indT,:);
    YgTT = ( Y - gTu).*T_poly_mat;

    %% hatQ2
    par_beta_g4p  = (par2ulambda+ulambda1).*vmu.*YgTT;
    par_lambda_g4p = vmu.*YgTT;
    par_mu_g4p = ulambda1.*YgTT;

    hatQ2 = [par_beta_g4p'*X, par_lambda_g4p'*uii, par_mu_g4p'*vZ]./N;

    %% Q3
    par_gamma_g4p =-ulambda1.*vmu.*T_poly_mat;
    hatQ3 = par_gamma_g4p'*T_poly_mat./N;
    
%% elements ij

 % Estimate E(gK*gK')



 g1r_ii = (ulambda1+1).*par2ulambda.* X;
 %g1r_ij = g1r_ii(indT,:); 

 g2K1_ii = ulambda1.* uii + uii;
 
 g3K2_ii = (vmu - IS ).*vZ;
 %g3K2_ij = g3K2_ii(indT,:);

 g4p_ii = ulambda1.*vmu.*YgTT;
 %g4p_ij = g4p_ii(indT,:);

 gK_ii = [g1r_ii, g2K1_ii, g3K2_ii, g4p_ii];

%%
Nvset = N*1000;
Kfold = floor(N^2/Nvset);
rem = N - Nvset*Kfold;

ngk = dimX + Ku + Kv + p; 
hatEgKgK_f1 =zeros(ngk);
for f = 1:Kfold
        id = (f-1)*Nvset+1:f*Nvset;
        indTf=indT(id);
        indZf=indZ(id);
        g1r_ijf = g1r_ii(indTf,:);
        g3K2_ijf = g3K2_ii(indTf,:);
        g4p_ijf = g4p_ii(indTf,:);
        
        Tf=T(indTf);
        Zf=Z(indZf);
        uijf = BasisFUNv8(Tf,Zf,p1u,p2u);
        g2K1_ijf = ulambda1ii(indTf,:).* uii(indTf,:) + uijf;

     % hatEgKgK
    gK_ijf = [g1r_ijf, g2K1_ijf, g3K2_ijf, g4p_ijf];
    
    hatEgKgK_f1 = hatEgKgK_f1 + gK_ijf'*gK_ijf./N ;
end
if rem >0
       id = Kfold*Nvset+1:N^2;
       indTf=indT(id);
        indZf=indZ(id);
        g1r_ijf = g1r_ii(indTf,:);
        g3K2_ijf = g3K2_ii(indTf,:);
        g4p_ijf = g4p_ii(indTf,:);
        
        Tf=T(indTf);
        Zf=Z(indZf);
        uijf = BasisFUNv8(Tf,Zf,p1u,p2u);
        g2K1_ijf = ulambda1ii(indTf,:).* uii(indTf,:) + uijf;

     % hatEgKgK
    gK_ijf = [g1r_ijf, g2K1_ijf, g3K2_ijf, g4p_ijf];
    
    hatEgKgK_f1 = hatEgKgK_f1 + gK_ijf'*gK_ijf./N ;
end
    hatEgKgK =  hatEgKgK_f1./(N-1) - gK_ii'*gK_ii./N./(N-1);

 %% inv_Q1
    hatQ13 = vZ'*vZ./N;
    hatQ13 = (hatQ13 + hatQ13')/2;


    par_beta_g3K2 = (dvZ*hatmu).*vZ+(vmu-IS).*dvZ;
    hatQ12 = [par_beta_g3K2'*X./N , zeros(Kv,Ku) ];

    %%

    hatQ11_1_1 = (par2ulambda.^2.*X)'* X./N;
    hatQ11_1_3 = (ulambda1.*par22ulambda.*X)'* X./N  - (par22ulambda.*X)'* X./N ;
    hatQ11_2_1 =  (par2ulambda.*X)'* uii./N;
    hatQ11_2_3 = ( ulambda1.*X)'* par2_uii./N - X'* par2_uii./N;
%%
nQ11_1_2=[dimX,dimX];
nQ11_2_2 = [dimX ,Ku ];


hatQ11_1_2 = zeros(nQ11_1_2);
hatQ11_2_2 = zeros(nQ11_2_2);
for f = 1:Kfold
        id = (f-1)*Nvset+1:f*Nvset;
        indTf=indT(id);
        indZf=indZ(id);
         Tf=T(indTf);
        Zf=Z(indZf);
        Xf = X(indZf,:);

        [ ~, par2_uijf,par22_uijf ]= BasisFUNv8(Tf,Zf,p1u,p2u);
  hatQ11_1_2 = hatQ11_1_2 +((par22_uijf*lambdahatu).*Xf)'*Xf./N;
   hatQ11_2_2 =  hatQ11_2_2 + Xf'* par2_uijf./N;
end
if rem >0
       id = Kfold*Nvset+1:N^2;
       indTf=indT(id);
        indZf=indZ(id);
         Tf=T(indTf);
        Zf=Z(indZf);
        Xf = X(indZf,:);

        [ ~, par2_uijf,par22_uijf ]= BasisFUNv8(Tf,Zf,p1u,p2u);
      hatQ11_1_2 = hatQ11_1_2 +((par22_uijf*lambdahatu).*Xf)'*Xf./N;
   hatQ11_2_2 =  hatQ11_2_2 + Xf'* par2_uijf./N;
end
 hatQ11_1_2 =  hatQ11_1_2./(N-1) - (par22ulambda.*X)'*X./N./(N-1);
 hatQ11_2_2 = hatQ11_2_2./(N-1) - X'* par2_uii./N./(N-1);

    % hatQ11_1_2 =((par22_uij*lambdahatu).*X(indZ,:))'*X(indZ,:);
    % hatQ11_2_2 = X(indZ,:)'* par2_uij./N./(N-1);
    % 

%%



    hatQ11_1 = hatQ11_1_1 + hatQ11_1_2 + hatQ11_1_3;
    hatQ11_2 = hatQ11_2_1 + hatQ11_2_2 + hatQ11_2_3;
    hatQ11_3 = hatQ11_2';


    hatQ11_4 = uii'*uii./N;
    hatQ11_4 = (hatQ11_4 + hatQ11_4')/2;

    hatQ11 = [hatQ11_1,hatQ11_2;hatQ11_3,hatQ11_4];
    hatQ11 = (hatQ11+hatQ11')/2;
    inv_hatQ11 = inv(hatQ11);

    l_Q11 = size(hatQ11,1);
    l_Q13 = size(hatQ13,1);
    %%
    inv_hatQ1 = [inv_hatQ11,zeros(l_Q11,l_Q13);...
        hatQ13\hatQ12/hatQ11,inv(hatQ13)];
    inv_hatQ3=inv(hatQ3);

    VQ = [-inv_hatQ3*hatQ2*inv_hatQ1,inv_hatQ3];
    var= VQ*hatEgKgK*VQ'/N;
    varhat = (1/2)*(var +var');
    %[V,D] = eig(varhat);
    %half_varhat = V*sqrt(D)*V';
    stdhat = diag(sqrt(diag(varhat)));
   

    CI_up= gammahatu + stdhat*norminv([1-alpha/2 ;1-alpha/2; 1-alpha/2] );
    CI_low= gammahatu - stdhat*norminv([1-alpha/2 ;1-alpha/2 ; 1-alpha/2]);



save(filename)

%% Step 4: Report computational time 
time = toc;        % finish timing
end_date = now;    % end date
disp('*****************************************');
disp(['Job-Realdata started: ', datestr(start_date)]);
disp(['Job-Realdata finished: ', datestr(end_date)]);
disp(['Computational time: ', num2str(time), ' seconds.']);
disp(['Computational time: ', num2str(time / 60), ' minutes.']);
disp(['Computational time: ', num2str(time / 3600), ' hours.']);
disp('*****************************************');
disp(' ');
