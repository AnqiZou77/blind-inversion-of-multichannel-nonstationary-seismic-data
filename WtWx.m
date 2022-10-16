function [Wt,Wx] = WtWx( Dt,Dx,Mz,para )
%% WITV 1
Wt = 1/(para.muw*log(2))*1./(1+exp(abs(Dt*Mz)/para.muw));
Wx = 1/(para.muw*log(2))*1./(1+exp(abs(Mz*Dx)/para.muw));

%% WITV 2
% [m,n] = size(Mz);
% Dt= zeros(m,m);
% Dx= zeros(n,n);
% for i = 2:m-1
%     Dt(i,i-1:i+1) = [-1,0,1]; % 纵向差分算子
% end
% for j = 2:n-1
%     Dx(j-1:j+1,j) = [-1,0,1].'; % 横向差分算子
% end
% Wt = 1/(para.muw*log(2))*1./(1+exp(abs(Dt*Mz)/para.muw));
% Wx = 1/(para.muw*log(2))*1./(1+exp(abs(Mz*Dx)/para.muw));

%% WITV 2-2
% [m,n] = size(Mz);
% Wt = ones(m,n);
% Wx = ones(m,n);
% for i = 2:m-1
%     for j = 2:n-1
%         Wt_i = Mz(i,j-1:j+1);
%         Wt_i_1 = Mz(i+1,j-1:j+1);
%         
%         Wx_j = Mz(i-1:i+1,j);
%         Wx_j_1 = Mz(i-1:i+1,j+1);
%         
%         if norm(Wt_i)==0||norm(Wt_i_1)==0
%             d = 1;
%             Wt(i,j) = 1/(para.muw*log(2))*1./(1+exp(abs(d)/para.muw));
%         else
%             d = 1-(Wt_i*Wt_i_1')/(norm(Wt_i)*norm(Wt_i_1));
% %             Wt(i,j) = exp(-d^2/para.muw);
%             Wt(i,j) = 1/(para.muw*log(2))*1./(1+exp(abs(d)/para.muw));
%         end
%         if  norm(Wx_j)==0||norm(Wx_j_1)==0
%             d = 1;
%             Wx(i,j) = 1/(para.muw*log(2))*1./(1+exp(abs(d)/para.muw));
%         else
%             d = 1-(Wx_j'*Wx_j_1)/(norm(Wx_j)*norm(Wx_j_1));
%             Wx(i,j) = 1/(para.muw*log(2))*1./(1+exp(abs(d)/para.muw));
%         end
%     end
% end
%% WITV3
% dt = 0.002;
% Q = 800;
% Wt = 1/(para.muw*log(2))*1./(1+exp(abs(NonstatReflectivity(Dt*Mz, dt,Q))/para.muw));
% Wx = 1/(para.muw*log(2))*1./(1+exp(abs(Mz*Dx)/para.muw));


end

