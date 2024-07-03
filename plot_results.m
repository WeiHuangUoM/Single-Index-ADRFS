gammau20plus = cell2mat(struct2cell( load('Realdata_BoxCox_CV_p1u3p2u3k2u4_estu-02-Jun-2024.mat','gammahatu')));
gammau0plus = cell2mat(struct2cell( load('Realdata_BoxCox_CV_tc0_p1u3p2u3k2u3_estu-02-Jun-2024.mat','gammahatu')));
gammau20less = cell2mat(struct2cell( load('Realdata_BoxCox_CV_tc20less_p1u3p2u3k2u4_estu-02-Jun-2024.mat','gammahatu')));





gamma20plus = cell2mat(struct2cell( load('Realdata_BoxCox_CV_p13p22k23_est-02-Jun-2024.mat','gammahat')));
gamma0plus = cell2mat(struct2cell( load('Realdata_BoxCox_CV_tc0_p13p22k22_est-02-Jun-2024.mat','gammahat')));
gamma20less = cell2mat(struct2cell( load('Realdata_BoxCox_CV_tc20less_p13p22k23_est-02-Jun-2024.mat','gammahat')));



t0= 1:100;
t0=reshape(t0,100,1);


g = @(theta,x_mat) x_mat*theta;
n = length(t0);
t_mat = repmat(log(t0+1),1,3).^repmat(0:(3-1),n,1);
guEst20plus = g(gammau20plus,t_mat);
guEst0plus = g(gammau0plus,t_mat);
guEst20less = g(gammau20less,t_mat);

figure(1)
plot(t0,guEst20plus)
hold on 
plot(t0,guEst0plus)
hold on 
plot(t0,guEst20less)
legend('S=\{t>20\}','S=\{t>=0\}','S=\{t<=20\}')
title('undersmoothing:p1=p1cv   p2=p2cv+1   k2=k2cv+1')


gEst20plus = g(gamma20plus,t_mat);
gEst0plus = g(gamma0plus,t_mat);
gEst20less = g(gamma20less,t_mat);
figure(2)
plot(t0,gEst20plus)
hold on 
plot(t0,gEst0plus)
hold on 
plot(t0,gEst20less)
legend('S=\{t>20\}','S=\{t>=0\}','S=\{t<=20\}')
title('p1=p1cv  p2=p2cv   k2=k2cv')




