%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%       This script illustrates disperison compensation based on method
%       proposed by Wilcox.  
% 
%       The example uses simulated response.
%             
%           
%                   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../tools')

clear all 
% close all
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%                  Load the dispersion curves data
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('../dispersion_data/disp_curve')  % disp curves for 2 mm alu plates 

a_cg = diff(2.*pi.*f_vec)./diff(k_a);       a_cg(end+1) = a_cg(end);
% figure(1),  plot(f_vec,k_a,f_vec,k_s)
% xlabel('f [Hz]'), ylabel('k [rad/m]')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%                           signal parameters 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

exct.fs = 5.12e6;               % sampling f 
exct.n_sampl = 4096;            % no of samples
exct.f_exc = 100e3;             % excitation f
exct.n_cycl = 2;                % no of cycles

[w indF] = min(abs(exct.f_exc - f_vec));
wave.vph = 2*pi*f_vec(indF)./k_a(indF);     % phase velocity in m/s for excitation freq
wave.lambda = wave.vph / f_vec(indF);       % wavelength in m for excitation freq

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%       Simulate  responses using structure's transfer function
%                          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d = 300;

[resp_a t sig] = dispResponse(f_vec,k_a,exct.fs,exct.n_sampl,exct.f_exc,exct.n_cycl,d);
figure(1), plot(t,sig,t,resp_a)


%% dispersion compensation steps:
        %% 1. padding the signal with zeros to be at least 8 times longer
        num_of_zeros = 2^(nextpow2(8*length(resp_a))) - length(resp_a);
        sig_padded = padarray(resp_a(round(exct.fs/exct.f_exc*exct.n_cycl/2):end),[0 num_of_zeros], 'post'); %padded g(t)

        %% 2. performing an FFT on the zero-padded time-domain signal
        fft_sig = fft(sig_padded); % G(omega)
        f_v = [0:length(fft_sig)-1]*exct.fs/length(fft_sig);
        fft_sig(end/2:end) = 0;

        %% 3. calculating:
        % the size of the distance step (less or equals 1/(2*k_Nyq))
        f_Nyq = exct.fs/2;
        k_Nyq = k_a(find(f_vec == f_Nyq));
        dx = 2*pi/(2*k_Nyq);
        % the size of the wavenumber step (less than the reciprical of length of original signal
        % multiplied by max. group velocity of a guided wave)
        m = length(resp_a);
        dt = 1/exct.fs;
        vg_max = max(a_cg);
        dk = 1./(m*dt*vg_max);
        dk = fix(dk * 100)/100;
        % the number of points (bigger than 2 times k_Nyq over dk)
        n_pts = ceil(2*k_Nyq/(2*pi)/dk);

        %% omega(k)

        %k_vec = [0:(n_pts-1)]*dk;
        k_vec = linspace(0,max(k_a),n_pts-1);
        omega_k = interp1(k_a,f_vec,k_vec);
        % figure(323)
        % plot(freq*1e3, k_a, omega_k,k_vec)
        %% 4. Interpolating G(omega) to find its value, G(k), at points equally spaced in k for the desired mode

        Gk = interp1(f_v,fft_sig,omega_k);
        Gk(end/2+1:end) = 0;
        % plot(f_v,abs(fft_sig), omega_k,abs(Gk))
        %% 5. Calculating the group velocity of the guided wave mode at the same wavenumber points, vgr(k).

        vgr_k = interp1(f_vec*1e3,a_cg*1e3,omega_k);
        % plot(omega_k,vgr_k)
        %% 6. Computing H(k) = G(k)vgr(k).
        Hk = Gk.*vgr_k;
        %% 7. Applying inverse FFT to H(k) to obtain the dispersion compensated distance-trace h(x).
        h_x = ifft(Hk,'nonsymmetric');

        h_x = h_x(1:2048);
        h_x = h_x ./ max(h_x); 
        h_dist = [0:length(real(h_x))-1].*2*pi./k_vec(end);
        plot(h_dist*1e3,abs(h_x), h_dist*1e3, real(h_x)*30)
        
        xlabel('x [mm]')
