function [WQ,err] = ADMM4WaveletATV(R,S,W,W0,L,para)
% WR = S
[m,n] = size(R);
Dt = zeros(m,m);
Dx = zeros(m,m);
for i = 1:m-1
    Dt(i,i:i+1) = [-1,1]; % 纵向差分算子
end
for j = 1:m-1
    Dx(j:j+1,j) = [-1,1].'; % 横向差分算子
end
%% l1
mask0 = zeros(L+m-1, m);
onew = ones(L,1);
onew(1:10) = linspace(0,1,10);
onew(end-9:end) = linspace(1,0,10);
for k = 1:m
    mask0(k:L+k-1, k) = onew;
end
mask = mask0((L-1)/2+1:(L-1)/2+m,:);
%% Dx*Dx' = UYU'  Dt'*Dt = UYU'
% A = para.sigma*(Dt'*Dt)=UYU';
U = zeros(m,m);
Y = diag(2-2*cos(pi*(0:m-1)/m));
for k = 1:m
    U(:,k) = cos(pi*(k-1)*(2*(1:m)-1)/(2*m)).'/norm(cos(pi*(k-1)*(2*(1:m)-1)/(2*m)).');
end
% B = VXV'
B = para.mu*eye(m)+para.beta*eye(m)+R*R'+para.yita*(Dx*Dx');
% [U,Y] = eig(A);
[V,X] = eig(B);
for k = 1:m
    V(:,k) = V(:,k)/norm(V(:,k));
end
invAB = zeros(m);
for i = 1:m
    for j = 1:m
        temp = X(j,j)+para.sigma*Y(i,i);
        if(temp == 0)
            invAB(i,j) = 0;
        else
            invAB(i,j) = 1/temp;
        end
    end
end

%% 

SRT = S*R';
k = 0;
Q = W0;
Qtemp = zeros(m);
Ht = Dt*W0;
Httemp = zeros(m);
Hx = W0*Dx;
Hxtemp = zeros(m);
WQ = W0;
Wtemp = zeros(m);
while norm(WQ-Wtemp)/norm(WQ)>para.tol&&k<para.maxIter
    k = k+1;
    Wtemp = WQ;
    C = SRT+para.beta*(Q-Qtemp)+para.sigma*Dt'*(Ht-Httemp)+para.yita*(Hx-Hxtemp)*Dx';
    WQ = (U'*C*V).*invAB;
    WQ = U*WQ*V';
    %  Q sub-problem
    Q= mask.*WQ;
    Qtemp = Qtemp+WQ-Q;
    if para.method == 1 % ATV
        %  Ht sub-problem
        Ht = max(1-para.gama./(para.sigma*abs(Dt*WQ+Httemp)),0).*(Dt*WQ+Httemp);
        Httemp = Httemp+Dt*WQ-Ht;
        %  Hx sub-problem
        Hx = max(1-para.gama./(para.sigma*abs(WQ*Dx+Hxtemp)),0).*(WQ*Dx+Hxtemp);
        Hxtemp = Hxtemp+WQ*Dx-Hx;
    else % TV
        s = sqrt((Dt*WQ+Httemp).^2+(WQ*Dx+Hxtemp).^2);
        s(find(s<1e-8))=1e-6;
        Ht = max(s-para.gama/para.sigma,0).*(Dt*WQ+Httemp)./s;
        Httemp = Httemp+Dt*WQ-Ht;
        Hx = max(s-para.gama/para.sigma,0).*(WQ*Dx+Hxtemp)./s;
        Hxtemp = Hxtemp+WQ*Dx-Hx;
%             sk = sqrt((DRow*u_new+bRow_old).^2+(u_new*DColumn+bColumn_old).^2);
%     sk(find(sk==0)) = 1e-6;
%     dRow_old = max(sk-1/lambda,0).*(DRow*u_new+bRow_old)./sk;
%     dColumn_old = max(sk-1/lambda,0).*(u_new*DColumn+bColumn_old)./sk;
%     
%     bRow_old = bRow_old+DRow*u_new-dRow_old;
%     bColumn_old = bColumn_old+u_new*DColumn-dColumn_old;
    end
    err(k) = norm(WQ-W)/norm(W);
end


