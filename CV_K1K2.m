function [p1,p2,K2] = CV_K1K2(T,X,Y,IS,p)

N = length(Y);

%ndgrid of K1, K2

pg = 1:3;
Kg2 = 1:2;

[pg1,pg2,Kg2] = ndgrid(pg,pg,Kg2);
pg1 = pg1(:);
pg2 = pg2(:);
Kg2 = Kg2(:);
lK = length(pg1);

Kfold = 5;
Nvset = floor(N/Kfold);
rem = N - Nvset*Kfold;
CV = zeros(1,lK);

for j = 1:lK
    p1 = pg1(j);
    p2 = pg2(j);
    K2 = Kg2(j);
    D = zeros(N,1);
    for i = 1:Kfold
        idx = (i-1)*Nvset+1:i*Nvset;
        D(idx) = diffi(idx,p1,p2,K2,Y,T,X,IS,p);
    end
    if rem >0
       idx = Kfold*Nvset+1:N;
       D(idx) = diffi(idx,p1,p2,K2,Y,T,X,IS,p);
    end
        
    CV(j) = mean(D);
end
  
index = find(CV == min(CV),1);
p1 = pg1(index);
p2 = pg2(index);
K2 = Kg2(index);
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

[weight, betahat] =  get_weight_semiclosedv3(T,X,Tt,Xt,Nt,p1,p2);
wv = weight(i);
weight(i) = [];
PS= PS_estCV(X,Xt,betahat,ISt,K2);
PSv = PS(i);
PS(i) = [];

gT=CaliEstPoly_ADRFS(Tv,Tt,PS,weight,Y,p);

D = wv.*PSv.*(Yv - gT).^2;

end

