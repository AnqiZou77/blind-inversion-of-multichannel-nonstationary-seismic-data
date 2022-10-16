function [Z_inversion,Mz,err,errIter] = ZinversionWITV(Z,W,S,Z0,Zlf,FL,para )
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

%% inversion R : ER+RH=F
E = para.mu*(W.'*W)+para.lambda*(Dt.'*Dt)+para.gama*eye(m);
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
        tmp1(i,j) = X(i,i)+para.lambda*Y(j,j);
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

[Wt,Wx] = WtWx( Dt,Dx,Mz,para );
Ht = Wt.*(Dt*Mz);
Hx = Wx.*(Mz*Dx);
Httemp = zeros(m,n);
Hxtemp = zeros(m,n);
k = 0;
while (norm(Mz-Mtemp)/norm(Mz)>para.tol&&k<para.maxIter) % the stop criterion can be improved
    k = k+1;
    %% subproblem R
    Mtemp = Mz;
    F = para.mu*W.'*S+para.lambda*Dt.'*(Bt-Bttemp)+para.lambda*(Bx-Bxtemp)*Dx.'+para.gama*(C-Ctemp);
    tmp2 = U.'*F*V;
    Mz = U*(pinv_tmp1.*tmp2)*V.';
    %% subproblem B
    [Wt,Wx] = WtWx( Dt,Dx,Mz,para );
    for iter = 1:10
        sk = sqrt((Wt.*Bt+Httemp).^2+(Wx.*Bx+Hxtemp).^2);
%         sk(find(sk<1e-8))=1e-6;
        sk(find(sk<10^-6.5))=1e-5;
        Bt = para.lambda*(Dt*Mz+Bttemp)+para.tao*Wt.*(Ht-Httemp);
        Bt = Bt./(para.lambda*ones(m,n)+para.tao*Wt.*Wt);
        Ht = max(1-1/para.tao./sk,0).*(Wt.*Bt+Httemp);
        Httemp = Httemp+Wt.*Bt-Ht;
        Bx = para.lambda*(Mz*Dx+Bxtemp)+para.tao*Wx.*(Hx-Hxtemp);
        Bx = Bx./(para.lambda*ones(m,n)+para.tao*Wx.*Wx);
        Hx = max(1-1/para.tao./sk,0).*(Wx.*Bx+Hxtemp);
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
    %% Z computation
    Z_inversion = exp(2*Mz);
    err(k) = norm(Z_inversion-Z)/norm(Z);
    errIter(k) = norm(Mz-Mtemp)/norm(Mz);
    norm(Z_inversion-Z)/norm(Z);

end
    



end

