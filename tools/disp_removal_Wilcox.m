function [out h_dist] = disp_removal_Wilcox(resp_a,exct, disp)

%% read dispersion data from the structure 
k_a = disp.k_a; 
f_vec = disp.f_vec; 

% calculate group velocity 
a_cg = diff(2.*pi.*f_vec)./diff(k_a);       a_cg(end+1) = a_cg(end);


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
        h_dist = [0:length(real(h_x))-1].*2*pi./k_vec(end);
%         h_x = h_x(1:2048);
%         h_x = h_x ./ max(h_x); 
        
        
        
%         plot(h_dist*1e3,abs(h_x), h_dist*1e3, real(h_x)*30)
%         xlabel('x [mm]')
        
        out = h_x; 
