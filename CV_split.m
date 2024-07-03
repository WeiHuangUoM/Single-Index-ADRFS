function CV_split(p1,p2,K2)

% clear all;
% close all;
% clc;

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

tc=20;
IS = double(T> log(tc+1));


%% Obtain Estimator
    %[p1(j),p2(j),K2(j)] = GCV_K1K2v2(T,X,Y,IS,p); 
%ndgrid of K1, K2

% pg = 1:3;
% Kg2 = 1:2;
% 
% [pg1,pg2,Kg2] = ndgrid(pg,pg,Kg2);
% pg1 = pg1(:);
% pg2 = pg2(:);
% Kg2 = Kg2(:);
% lK = length(pg1);

Kfold = 5;
Nvset = floor(N/Kfold);
rem = N - Nvset*Kfold;
% CV = zeros(1,lK);

% for j = 1:lK
%     p1 = pg1(j);
%     p2 = pg2(j);
%     K2 = Kg2(j);
    D = zeros(N,1);
    for i = 1:Kfold
        i
        idx = (i-1)*Nvset+1:i*Nvset;
        D(idx) = diffi(idx,p1,p2,K2,Y,T,X,IS,p);
    end
    if rem >0
       idx = Kfold*Nvset+1:N;
       D(idx) = diffi(idx,p1,p2,K2,Y,T,X,IS,p);
    end
        
    CV = mean(D);
% end
  
% index = find(CV == min(CV),1);
% p1 = pg1(index);
% p2 = pg2(index);
% K2 = Kg2(index);


filename=sprintf('CV-tc%d-%d-%d-%d-%s.mat', tc, p1,p2,K2,datetime('today'));
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

end


function D = diffi(i,p1,p2,K2,Y,T,X,IS,p)

Tv = T(i);
Yv = Y(i);

Xt = X;
Tt = T;
ISt = IS;
Xt(i,:)=[];
Tt(i)=[];
Y(i)=[];
ISt(i) = [];
Nt = length(Tt);

[weight, betahat] =  get_weight_semiclosed_largeN(T,X,Tt,Xt,Nt,p1,p2);
wv = weight(i);
weight(i) = [];
PS= PS_estCV(X,Xt,betahat,ISt,K2);
PSv = PS(i);
PS(i) = [];

gT=CaliEstPoly_ADRFS(Tv,Tt,PS,weight,Y,p);

D = wv.*PSv.*(Yv - gT).^2;

end
