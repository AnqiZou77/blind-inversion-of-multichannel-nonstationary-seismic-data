%% code for statianary impadence inversion
clc; 
clear all;
close all;

%%  ¶ÁÈ¡Êý¾Ý£¨load Synthetic Data: Acoustic Impedance, Seismic data, wavelet£©
load('SyntheticData\AI_1000_1361.mat');
AI = double(AI);
load('SyntheticData\Dnons.mat');
Dnons = double(Dnons);
load('SyntheticData\Wa.mat');
load('SyntheticData\Ws.mat');

Z = AI;
[m,n] = size(AI);

D_obs = Dnons;
U_true = 0.5*log(AI);
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

%% inversion parameter for TV done
% para.tol = 1e-5; % 1-6
% para.lambda = 2; % 1.5
% para.mu = 150; % 150
% para.beta = 0.05; % 3.5
% para.gama = 0.05; %3.5
% para.maxIter = 20; % 200
% tic
% G = Wa*Dt;
% U0 = ulf;
% [Z_inversion,U,err,errIter] = ZinversionTV( Z,G,D_obs,U0,Ulf,FL,para);
% Zfile = '.\output4TV\Z_nonSinversion_TV.mat';
% errfile = '.\output4TV\error_TV.mat';
% errIterfile = '.\output4TV\error_iter_TV.mat';

%% inversion parameter for ATV
% para.tol = 1e-5; % 1-6
% % para.lambda = 13; %13
% para.lambda1 = 15; %4.5
% para.lambda2 = 3.5; %2.5
% para.mu = 150; % 10
% para.beta = 0.05; % 3.5
% para.gama = 0.05; %3.5
% para.maxIter = 20; % 200
% 
% tic
% G = Wa*Dt;
% U0 = ulf;
% [Z_inversion,U,err,errIter] = ZinversionATV( Z,G,D_obs,U0,Ulf,FL,para);
% toc
% Zfile = '.\output4ATV\Z_nonSinversion_TV.mat';
% errfile = '.\output4ATV\error_TV.mat';
% errIterfile = '.\output4ATV\error_iter_TV.mat';
   

%% inversion parameter for WITV done
para.tol = 1e-5; % 1-6
para.lambda = 2; %1.5
para.mu = 15; % 15
para.beta = 0.05; % 0.05 parameter for Low frequency
para.gama = 0.05; % 0.05 parameter for low frequency 
para.tao = 2; % 1.5 parameter for wtv Ht 
para.muw = 0.0018; % parameter for wtv 0.0018
para.maxIter = 20; % 20
G = Wa*Dt;
U0 = ulf;
tic
[Z_inversion,U,err,errIter] = ZinversionWITV( Z,G,D_obs,U0,Ulf,FL,para);
t = toc

Zfile = '.\output4WITV\Z_nonSinversion_WITV.mat';
errfile = '.\output4WITV\error_WITV.mat';
errIterfile = '.\output4WITV\error_iter_WITV.mat';
%% inversion parameter for WITV with stationary model
% para.tol = 1e-5; % 1-6
% para.lambda = 1; %9.3
% para.mu = 6.3; % 200
% para.beta = 0.05; % 0.35 parameter for Low frequency
% para.gama = 0.05; % 0.35 parameter for low frequency 
% para.tao = 0.9; % parameter for wtv Ht 9.2
% para.muw = 0.0018; % parameter for wtv 0.0018
% para.maxIter = 20; % 20
% % para.method = 3; % 1 for TV; 2 for ATV; 3 for WTV
% % load('.\output4BD\W_inversion_1.mat')
% G = Ws*Dt;
% U0 = ulf;
% tic
% [Z_inversion,U,err,errIter] = ZinversionWITV( Z,G,D_obs,U0,Ulf,FL,para);
% t = toc
% 
% Zfile = '.\output4WITVS\Z_nonSinversion_WITV_S.mat';
% errfile = '.\output4WITVS\error_WITV_S.mat';
% errIterfile = '.\output4WITVS\error_iter_WITV_S.mat';
%% inversion parameter for WITV without lowfrequency constrain
% para.tol = 1e-5; % 1-6
% para.lambda = 1; %1
% para.mu = 6.3; % 6.3
% para.beta = 0; % 0.05 parameter for Low frequency
% para.gama = 0; % 0.05 parameter for low frequency 
% para.tao = 0.9; % parameter for wtv Ht 0.9
% para.muw = 0.0018; % parameter for wtv 0.0018
% para.maxIter = 20; % 20
% G = Wa*Dt;
% U0 = ulf;
% tic
% [Z_inversion,U,err,errIter] = ZinversionWITV( Z,G,D_obs,U0,Ulf,FL,para);
% t = toc
% Z_inversion = Z_inversion/norm(Z_inversion(1,742))*norm(Z(1,742));
% Zfile = '.\output4WITVnoLF\Z_nonSinversion_WITV_noLF.mat';
% errfile = '.\output4WITVnoLF\error_WITV_noLF.mat';
% errIterfile = '.\output4WITVnoLF\error_iter_WITV_noLF.mat'
%% inversion parameter for WATV done
% para.tol = 1e-5; % 1-6
% para.lambda1 = 9.3; %9.3
% para.lambda2 = 15.3; %9.3
% para.mu = 200; % 200
% para.beta = 0.35; % 0.35 parameter for Low frequency
% para.gama = 0.35; % 0.35 parameter for low frequency 
% para.tao1 = 9.3; % parameter for wtv Ht 9.2
% para.tao2 = 15.3; % parameter for wtv Ht 9.2
% para.muw = 0.018; % parameter for wtv 0.0018
% para.maxIter = 20; % 20
% G = Wa*Dt;
% U0 = ulf;
% tic
% [Z_inversion,U,err,errIter] = ZinversionWATV( Z,G,D_obs,U0,Ulf,FL,para);
% t = toc
% 
% Zfile = '.\output4WATV\Z_nonSinversion_WATV.mat';
% errfile = '.\output4WATV\error_WATV.mat';
% errIterfile = '.\output4WATV\error_iter_WATV.mat';

%% post-process
P = 50;
Z  = Z(P+1:m-P,P+1:n-P);
% save('.\output\Ztrue.mat','Z');
Z_inversion = Z_inversion(P+1:m-P,P+1:n-P);
norm(Z_inversion-Z)/norm(Z)
save(Zfile,'Z_inversion');
Zlf = Zlf(P+1:m-P,P+1:n-P);

[m,n] = size(Z);

x = (1:m)*2/1000;
y = (1:n)*12.5/1000;
%% plot
trace = 650;
figure
imagesc(y,x,Dnons); colormap(1-gray);
title('Nonstationary Seismic Data')
colorbar 
figure(1)
imagesc(y,x,Z); colormap(1-gray);
title('Ture AI')
tempC = caxis;
ylabel('Time(s)')
xlabel('Trace Number')
ax = gca;
ax.XAxisLocation = 'origin';
figure(2)
imagesc(y,x,Z_inversion); colormap(1-gray);
title('Estimated AI')
caxis(tempC)
ylabel('Time(s)')
xlabel('Trace Number')
ax = gca;
ax.XAxisLocation = 'origin';
figure(3)
imagesc(y,x,Zlf); colormap(1-gray);
title('Initial AI')
caxis(tempC)
ylabel('Time(s)')
xlabel('Trace Number')
ax = gca;
ax.XAxisLocation = 'origin';
figure(4)
plot(x,Z(:,trace),'k');
axis tight;
hold on;
plot(x,Zlf(:,trace),'b');
axis tight;
xlabel('t/ms');
ylabel('v/ms^{-1}');
hold on
plot(x,Z_inversion(:,trace),'r');
axis tight;
xlabel('t/ms');
ylabel('v/ms^{-1}');
title('Comparison AI of single channe')

figure(5)
plot(err)

figure(6)
plot(errIter(1:end))


