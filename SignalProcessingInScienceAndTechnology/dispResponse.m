function [resp_a varargout] = dispResponse(f_vec,k_a,fs,n_sampl,f_exc,n_cycl,x)
%
% function resp_a = dispResponse(f_vec,k_a,fs,n_sampl,f_exc,n_cycl,x)
% 
% function to generate a dispersive response based on input data: 
% f_vec, k_a - dispersion curves in form of f, k pairs 
% fs      - sampling frequency
% n_sampl - number of samples
% f_exct  - excitation frequency 
% n_cycl  - number of periods in the excitation signal
% x       - propagation distance
% 
%   09.2019 L. Ambrozinski ambrozin@agh.edu.pl



t = [0:(n_sampl-1)]./fs;           % time vector

N_win = round(n_cycl/f_exc*fs);
w = window(@hann,N_win);
win = zeros(1,length(t));
win(1:length(w)) = w;
sign = sin(2*pi*f_exc*t).*win ;


S = fft(sign);
df = fs/ n_sampl;
f = 0:df:fs-df;
f=f(1:end/2+1);

% interpolate dispersion curve to have values where FT samples are
k_int_a=interp1(f_vec,k_a,f);       k_int_a(end)=0;

% Create structure's transfer function     
G_a=[ exp(1).^(-i.*k_int_a.*x*1e-3) exp(1).^(i.*fliplr(k_int_a(2:end-1)).*x*1e-3)];

% convolve and inverse Fourier transform 
StimesGa=S.*G_a; 


resp_a=real(ifft(StimesGa));  
varargout{1} = t; 
varargout{2} = real(sign);
end
