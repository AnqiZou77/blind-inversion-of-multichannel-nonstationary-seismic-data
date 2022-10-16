%% code for impadence inversion from nonstationary seismic data
clc; clear; close all;
%%  读取数据（load Synthetic Data: Acoustic Impedance, Seismic data, wavelet）
load('.\SyntheticData\AI_1000_1361.mat');
AI = double(AI);
load('.\SyntheticData\Dnons.mat');
Dnons = double(Dnons);
load('.\SyntheticData\Wa.mat');
load('.\SyntheticData\Ws.mat');
load('.\SyntheticData\rw.mat');
Z = AI;
[m,n] = size(AI);

Dt = zeros(m,m);
for i = 1:m-1
    Dt(i,i:i+1) = [-1,1]; % Differential matrix along t axis
end

D_true = Dnons;
M_true = 0.5*log(AI);
R = Dt*M_true;
S = D_true;
Wtrue = Wa;

%% parameter 
para.mu = 0.001; % l2 0.001
para.beta = 0.01;  % 0.01
para.lambda = 0; % l1
para.gama = 0.05; % DtDx 0.005
para.sigma = 0.0005;  % Dt 0.0005
para.yita = 0.0005; % Dx 0.005
para.maxIter = 20;
para.tol = 1e-6;
para.method = 1; % 1 : ATV ,2 : TV
dt = 0.002;
Q = 800;

%% initialization
Lw = 63;
W0 = eye(m);
Ra = NonstatReflectivity(R, dt,Q);

%% 反演(ADMM algorithm)
tic
[WQ,err] = ADMM4WaveletATV(Ra,S,Ws,W0,Lw,para);
L = 129;
WQa = NonstatWave( WQ,L,m,dt,Q);
save('.\output4WEst\WQa_ADMM_TV.mat','WQa');
save('.\output4WEst\WQ_ADMM_TV.mat','WQ');
t = toc
norm(WQa-Wa)/norm(Wa)

%% plot
figure(1)
imagesc(WQa)
colormap('jet')
title('True nonstationary wavelet matrix')
figure(2)
imagesc(Wa)
colormap('jet')
title('Estimated nonstationary wavelet matrix')
figure(3)
rwq = zeros(1,129);
for i = 65:m-64
    rwq = rwq+WQ(i,i-64:i+64);
end
rwq = rwq/(m-L+1);
plot(rwq,'r')
hold on
plot(rw,'b')
title('Comparison of true and estimated wavelet in time domain')

Da = Wa*Dt*M_true;
Ds = Ws*Ra;
norm(Da-Ds)
figure(4)
plot(Ds(:,650))
hold on
plot(Dnons(:,650))
title('Comparison of true and sythetic sismic data')

