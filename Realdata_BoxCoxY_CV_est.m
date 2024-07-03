clear all;
close all;
clc;

addpath '/Users/zyl/Library/CloudStorage/OneDrive-Personal/Non_CTE20220226 _2/Matlab Code/28-April-2024'


tic;                                           % start timing
start_date = now;                              % start date
disp(['Job-Realdata started: ', datestr(start_date)]);  % display start date

%% Step 1: Set-up
 
p = 3;                 % link function is quadratic if p = 3.
alpha = 0.05;          % nominal size

% fix seed so that we always get the same result
defaultStream = RandStream.getGlobalStream;
defaultStream.reset(3);

%% Step 2: Read Data
load zip_ads_v4.txt
zip = zip_ads_v4(:,1); % zip code
ads = zip_ads_v4(:,2); % # of political advertisements aired in each zipcode (treatment of interest)
pop = zip_ads_v4(:,3); % population
pct_over65 = zip_ads_v4(:,4); % percent over age 65
median_income = zip_ads_v4(:,5); % median household income
pct_black = zip_ads_v4(:,6); % percent black
pct_hispanic = zip_ads_v4(:,7); % percent hispanic
pop_density = zip_ads_v4(:,8); % population density
pct_colgrads = zip_ads_v4(:,9)/100; % percent college graduates
can_commute = zip_ads_v4(:,10); % binary indicator of whether it is possible to commute to the zip code from a competitive state
contribution = zip_ads_v4(:,11); % total amount of campaign contribution 
log_pop = log(pop); % log of population
log_median_income = log(median_income + 1); % some elements of median_income are 0
log_pop_density = log(pop_density + 1); % log of population density
clear median_income pop_density
Y = contribution;
T = ads; % treatment
X = [log_pop, log_pop_density, log_median_income, pct_over65,...
pct_hispanic, pct_black, pct_colgrads,can_commute];
clear log_pop log_pop_density log_median_income pct_over65 pct_hispanic pct_black pct_colgrads can_commute
N = length(Y); % sample size
data = [Y, T, X];
dimX = size(X,2);
clear contribution

lambda = FindLamBoxCox(Y); % find the best Box-Cox transformation parameter for outcome
Y = BoxCox(Y,lambda); % outcome 
Y = Y-min(Y);
T = log(T+1);
%Nt=100;
%t0 = min(T):(max(T)-min(T))/(Nt-1):max(T); 
 tc=0;
IS = double(T>=log(tc+1));



%% Step 2: Empirical analysis (point estimation)
% Y = log(contribution + 1);                     % outcome 
% T = log(ads + 1);                              % treatment
% t = (min(T) : 0.05 : max(T))';
% X = [log_pop, log_pop_density, log_median_income, pct_over65, pct_hispanic, pct_black, pct_colgrads, can_commute];    % covariates
% dimX = size(X,2);                              % dimension of covariates
% N = size(zip,1);                               % sample size

%% Obtain Estimator
    %[p1(j),p2(j),K2(j)] = GCV_K1K2v2(T,X,Y,IS,p); 
    %[p1,p2,K2] = CV_K1K2(T,X,Y,IS,p);
    p1 = 3;
    p2 =2;
    K2=2;
    [weight, betahat, ~] =  get_weight_semiclosedv_largeN(T,X,T,X,N,p1,p2);

    PShat= PS_estv2(X,betahat,IS,K2);
    [~,gammahat]=CaliEstPoly_ADRFS(T,T,PShat,weight,Y,p);

filename=sprintf('Realdata_BoxCox_CV_tc0_est-%s.mat',DGP,N,datetime('today'));
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
