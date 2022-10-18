%% code for nonstatianary seismic inversion
clc;
clear all;
close all;
%%  读取数据（load Synthetic Data: Acoustic Impedance, Seismic data, wavelet）
load('SyntheticData\AI_1000_1361.mat');
AI = double(AI);
load('SyntheticData\Dnons.mat');
Dnons = double(Dnons);
load('SyntheticData\Wa.mat');
load('SyntheticData\Ws.mat');
load('SyntheticData\W0.mat');
Z = AI;
[m,n] = size(Z);

D_obs = Dnons;
U_true = 0.5*log(Z);
Wtrue = Wa;
Dt = zeros(m,m);
for i = 1:m-1
    Dt(i,i:i+1) = [-1,1]; % Differential matrix along t axis
end

%% low frequency filter, get the initial AI model
G1 = dftmtx(m);
G2 = dftmtx(n);
u=0:(m-1);v=0:(n-1);
idx = find(u > m/2);u(idx) = u(idx) - m;
idy = find(v > n/2);v(idy) = v(idy) - n;
[V,U]=meshgrid(v,u);
D=sqrt(U.^2+V.^2);
D0=0.005*m;
FL=exp(-(D.^2)./(2*(D0)^2));
Ulf = FL.*(G1*U_true*G2);
ulf = real(G1'*Ulf*G2'/m/n);
Zlf = exp(2*ulf);

%% inversion parameter for WITV
% paraW.lambda = 0; % l1 0
% paraW.beta = 0.01;  % 0.01 
% paraW.mu = 0.001; % l2  0.001
% paraW.gama = 0.05; % DtDx 0.05
% paraW.sigma = 0.0005;  % Dt 0.0005
% paraW.yita = 0.0005; % Dx 0.0005
% paraW.maxIter = 20; 
% paraW.tol = 1e-3;
% paraW.method = 1; % 1 : ATV ,2 : TV
% Lw = 63;
% paraZ.tol = 1e-5; % 1-6
% paraZ.lambda = 2; % 2
% paraZ.mu = 12; % 15
% paraZ.beta = 0.05; % 0.05 parameter for Low frequency
% paraZ.gama = 0.05; % 0.05 parameter for low frequency 
% paraZ.tao = 2; % parameter for wtv Ht 0.9
% paraZ.muw = 0.0018; % parameter for wtv 0.0018
% paraZ.maxIter = 20; % 200
% outfilepre = '.\output4BDWITV\WITV_';
% kmax = 8;

%% inversion parameter for TV
paraW.lambda = 0; % l1
paraW.beta = 0.01;  % l1 0.01
paraW.mu = 0.001; % l2  0.005
paraW.gama = 0.05; % DtDx 0.05
paraW.sigma = 0.0005;  % Dt 0.0005
paraW.yita = 0.0005; % Dx 0.0005
paraW.maxIter = 20;
paraW.tol = 1e-5;
paraW.method = 1; % 1 : ATV ,2 : TV
Lw = 63;

paraZ.tol = 1e-6; % 1-6
paraZ.lambda = 2; %2
paraZ.mu = 300; % 150
paraZ.beta = 0.05; % low frequency 3.5
paraZ.gama = 0.05; % low frequency 3.5
paraZ.maxIter = 20; % 200
paraZ.method = 1; % 1 for TV; 2 for ATV
outfilepre = 'output4BDTV\TV_';
kmax = 20;

%% 交替迭代反演 (alternating minimization algorithm)
k = 0;

%% initialization
Unew = ulf;
WQa = W0;
WQ = W0;
dt = 0.002;
Q = 800;
L = 129;
errZ = zeros(paraZ.maxIter,kmax);
errW = zeros(paraW.maxIter,kmax);
k = 0;
tic
while k<kmax
    % AI inversion
    k = k+1
    G = WQa*Dt;
%     [Z_inversion,Unew,err,~] = ZinversionWITV( Z,G,D_obs,Unew,Ulf,FL,paraZ);
    [Z_inversion,Unew,err,~] = ZinversionTV( Z,G,D_obs,Unew,Ulf,FL,paraZ);
    errZ(1:length(err),k) = err;
    Zfilename = [outfilepre 'Z_inversion_' num2str(k) '.mat'];
    save(Zfilename,'Z_inversion');
    Rs= Dt*Unew;
    Ra = NonstatReflectivity(Rs, dt,Q);
    % wavelet estimation
    paraW.maxIter = paraW.maxIter-1; 
    [WQ,err] = ADMM4WaveletATV(Ra,D_obs,Ws,WQ,Lw,paraW);
    WQa = NonstatWave( WQ,L,m,dt,Q);
    errW(1:length(err),k) = err;
    Wfilename = [outfilepre 'W_inversion_' num2str(k) '.mat'];
    save(Wfilename,'WQa');
   %     prints = 'W inversion done'
    eZ(k) = norm(Z_inversion-Z)/norm(Z)
    eW(k) = norm(WQa-Wtrue)/norm(Wtrue)
    
end
time = toc
Zerrfilename = [outfilepre 'errorZ.mat'];
Werrfilename = [outfilepre 'errorW.mat'];
save(Zerrfilename,'eZ')
save(Werrfilename,'eW')


figure(1)
imagesc(Z); colormap(1-gray);
title('Ture AI')
tempC = caxis;
ylabel('Time(s)')
xlabel('Distance(km)')
ax = gca;
ax.XAxisLocation = 'origin';
figure(2)
imagesc(Z_inversion); colormap(1-gray);
title('Estimated AI')
caxis(tempC)
ylabel('Time(s)')
xlabel('Distance(km)')
ax = gca;
ax.XAxisLocation = 'origin';
figure(3)
imagesc(Zlf); colormap(1-gray);
title('Initial AI')
caxis(tempC)
ylabel('Time(s)')
xlabel('Distance(km)')
ax = gca;
ax.XAxisLocation = 'origin';

trace = 650;
figure(4)
plot(Z(:,trace),'k');
axis tight;
xlabel('t/ms');
ylabel('v/ms^{-1}');
hold on;
plot(Zlf(:,trace),'b');
axis tight;
xlabel('t/ms');
ylabel('v/ms^{-1}');
hold on
plot(Z_inversion(:,trace),'r');
axis tight;
xlabel('t/ms');
ylabel('v/ms^{-1}');
title('Comparison AI of single channe')

figure(5)
imagesc(WQa)
colormap('jet')
title('Estimated nonstationary wavelet matrix')
figure(6)
imagesc(Wtrue)
colormap('jet')
title('True nonstationary wavelet matrix')

