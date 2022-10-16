function [Z_inversion,Mz,err,errIter] = ZinversionWATV(Z,W,S,Z0,Zlf,FL,para )
% S = WZ
[m,n] = size(S);
%% prepare
Dt= zeros(m,m);
Dx= zeros(n,n);
for i = 1:m-1
    Dt(i,i:i+1) = [-1,1]; % Differential matrix along t axis
end
for j = 1:n-1
    Dx(j:j+1,j) = [-1,1].'; % Differential matrix along x axis
end
% G1 = dftmtx(m);
% G2 = dftmtx(n);
%% inversion R : ER+RH=F
E = para.mu*(W.'*W)+para.lambda1*(Dt.'*Dt)+para.gama*eye(m);
[U,X] = eig(E); 
for k = 1:m
    U(:,k) = U(:,k)/norm(U(:,k));
end 
V = zeros(n,n);
Y = diag(2-2*cos(pi*(0:n-1)/n));
for k = 1:n
    V(:,k) = cos(pi*(k-1)*(2*(1:n)-1)/(2*n)).'/norm(cos(pi*(k-1)*(2*(1:n)-1)/(2*n)).');
end   
tmp1 = zeros(m,n);
pinv_tmp1 = zeros(m,n);
for i = 1:m
    for j = 1:n
        tmp1(i,j) = X(i,i)+para.lambda2*Y(j,j);
        if tmp1(i,j) == 0
            pinv_tmp1(i,j) = 0;
        else
            pinv_tmp1(i,j) = 1/tmp1(i,j);
        end
    end
end
Mz = Z0;
Mtemp = zeros(m,n);
Bt = Dt*Z0;
Bx = Z0*Dx;
Bttemp = zeros(m,n);
Bxtemp = zeros(m,n);
C = Z0;
Ctemp = zeros(m,n);

Wt = 1/(para.muw*log(2))*1./(1+exp(abs(Dt*Mz)/para.muw));
Wx = 1/(para.muw*log(2))*1./(1+exp(abs(Mz*Dx)/para.muw));
Ht = Wt.*(Dt*Mz);
Hx = Wx.*(Mz*Dx);
Httemp = zeros(m,n);
Hxtemp = zeros(m,n);

k = 0;
while (norm(Mz-Mtemp)/norm(Mz)>para.tol&&k<para.maxIter) % the stop criterion can be improved
    k = k+1
    %% subproblem R
    Mtemp = Mz;
    F = para.mu*W.'*S+para.lambda1*Dt.'*(Bt-Bttemp)+para.lambda2*(Bx-Bxtemp)*Dx.'+para.gama*(C-Ctemp);
    tmp2 = U.'*F*V;
    Mz = U*(pinv_tmp1.*tmp2)*V.';
    
    %% B subproblem TV
    
    Wt = 1/(para.muw*log(2))*1./(1+exp(abs(Dt*Mz)/para.muw));
    Wx = 1/(para.muw*log(2))*1./(1+exp(abs(Mz*Dx)/para.muw));
    %% H subproblem in B
    
    for iter = 1:10

        Bt = para.lambda1*(Dt*Mz+Bttemp)+para.tao1*Wt.*(Ht-Httemp);
        Bt = Bt./(para.lambda1*ones(m,n)+para.tao1*Wt.*Wt);
        Ht = max(1-1/para.tao1./abs(Wt.*Bt+Httemp),0).*(Wt.*Bt+Httemp);
        Httemp = Httemp+Wt.*Bt-Ht;
        Bx = para.lambda2*(Mz*Dx+Bxtemp)+para.tao2*Wx.*(Hx-Hxtemp);
        Bx = Bx./(para.lambda2*ones(m,n)+para.tao2*Wx.*Wx);
        Ht = max(1-1/para.tao2./abs(Wx.*Bx+Hxtemp),0).*(Wx.*Bx+Hxtemp);
        Hxtemp = Hxtemp+Wx.*Bx-Hx;

    end
    Bttemp = Bttemp+Dt*Mz-Bt;
    Bxtemp = Bxtemp+Mz*Dx-Bx;
    %% C subproblem
    C = para.beta*m*n*(FL.*Zlf)+para.gama*fft2(Mz+Ctemp);
    C = C./(para.beta*m*n*FL.*FL+para.gama*ones(m,n));
    C = max(C,0);
    C = real(ifft2(C));
    Ctemp = Ctemp+Mz-C;

    Z_inversion = exp(2*Mz);
    err(k) = norm(Z_inversion-Z)/norm(Z);
    errIter(k) = norm(Mz-Mtemp)/norm(Mz);
    norm(Z_inversion-Z)/norm(Z);

end
    
end

