function yn = constPhaseRotation( xn, phi )
% -------------------------------------------------------------------------
% Make a constant phase rotation for input data.
% 
% Syntax    : yn = constPhaseRotation( xn, phi )
% 
% Input
%       xn  : input signal(column vector)
%       phi : constant phase to be rotated in degrees
% Output
%       yn  : signal after phase rotation
% 
% Written by Zhang Ming, 2009-10.
% -------------------------------------------------------------------------

phi = phi * pi/180;
Ls = length(xn);
nf = floor(Ls/2);

if( rem(Ls,2) )
    Phase = exp( -i * phi*[0; ones(nf,1); -ones(nf,1)] );
else
    Phase = exp( -i * phi*[0; ones(nf,1); -ones(nf-1,1)] );
end

rotXk = fft(xn) .* Phase;
yn = real( ifft(rotXk) );
