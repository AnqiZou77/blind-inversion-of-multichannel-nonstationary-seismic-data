function Wa = NonstatWave( Ws,Lw,m,dt,Q)
% Lw = length(rw);
Lm = floor(Lw/2);
rw = zeros(Lw,1);
for i = Lm+1:m-Lm
    rw = rw+Ws(i,i-Lm:i+Lm).';
end
rw = rw/(m-2*Lm);
Fw = fft(rw,m);
t = (0:m-1)*dt;
f = (0:m-1)*1/(m*dt);
% Q = 800;
Aq = t'*f;
Aq = exp(-pi*Aq/Q);
tu = t'*ones(1,m);
tu = tu-tu';
Wa = zeros(m,m);
for j = 1:m
Eu = exp(2*pi*1i*f'*tu(j,:));
Wa(:,j) = Aq.*(Eu.')*Fw;
end
Wa = real(Wa)/m;
Wa = circshift(Wa,[0,-Lm]);


end

