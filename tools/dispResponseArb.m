function [resp_a varargout] = dispResponseArb(f_vec,k_a,fs,sign,x)
%
% function resp_a = dispResponseArb(f_vec,k_a,fs,sign,x)
% 
% optional:
% function [resp_a t] = dispResponseArb(f_vec,k_a,fs,sign,x)
% 
% function to generate a dispersive response based on input data: 
% f_vec, k_a - dispersion curves in form of f, k pairs 
% fs      - sampling frequency
% sign    - arbitrary signal subject to processing 
% x       - propagation distance
% 
%   09.2019 L. Ambrozinski ambrozin@agh.edu.pl

S = fft(sign);
n_sampl = length(S); 
df = fs/ n_sampl;
f = 0:df:fs-df;
f=f(1:end/2+1);

% interpolate dispersion curve to have values where FT samples are
k_int_a=interp1(f_vec,k_a,f);       k_int_a(end)=0;

% Create structure's transfer function     
G_a=[ exp(1).^(-i.*k_int_a.*x*1e-3) exp(1).^(i.*fliplr(k_int_a(2:end-1)).*x*1e-3)];

% convolve and inverse Fourier transform 
StimesGa=S.*G_a; 

resp_a=ifft(StimesGa);  
% varargout{1} = t; 

