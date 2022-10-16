function Ra = NonstatReflectivity(Rs, dt,Q)
% Lm = floor(L/2);
[m,n] = size(Rs);
t = (0:m-1)*dt; % 1*1000
f = (0:m-1)*1/(m*dt); % 1*1000
tf = t'*f;
Aq = exp(-pi*tf/Q);
clear i
Etf = exp(-2*pi*i*tf);
Ra = zeros(m,n);
for j = 1:n
    Rtemp1 = Rs(:,j);
    Rtemp2 = Rtemp1*ones(1,m);
    Rfa = sum(Aq.*Rtemp2.*Etf,1);
    Ra(:,j) = real(ifft(Rfa));

end

end

