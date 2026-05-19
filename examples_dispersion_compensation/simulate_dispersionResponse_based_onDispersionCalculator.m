%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    based on
%           simulated square array response.
%             
%           
%                   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('..\tools')

clear all 
close all
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%                  Load the dispersion curves data
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try 
    load('..\dispersion_data\alu_4mm_DispersionCurves')  % disp curves for 2 mm alu plates 
catch % mac book
    load('../dispersion_data/alu_4mm_DispersionCurves')  % disp curves for 2 mm alu plates
end
% figure(1),  plot(f_vec,k_a,f_vec,k_s)
% xlabel('f [Hz]'), ylabel('k [rad/m]')
%%


f_A0 = A{:,1};      % frequency of A0 mode
v_A0 = A{:,2};      % phase velcity of A0 mode
vg_A0 = A{:,3};      % phase velcity of A0 mode
k_A0 = A{:,10}
f_A1 = A{:,12};      % frequency of A1 mode
v_A1 = A{:,13};      % phase velcity of A1 mode
f_S0 = S{:,1};      % frequency of S0 mode
v_S0 = S{:,2};      % phase velcity of S0 mode
k_S0 = S{:,10}
vg_S0 = S{:,3};      % phase velcity of A0 mode

plot(f_A0,v_A0, f_A1, v_A1, f_S0, v_S0)
xlabel('frequency (kHz)')
ylabel('Phase velocity (km/s)')

 % load('..\dispersion_data\disp_curve'); 
 % dispC.f_v = f_vec; 
 % dispC.k_a = k_a; 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%                           signal parameters 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

exct.fs = 8e6;               % sampling f 
exct.n_sampl = 4096;            % no of samples
exct.f_exc = 1000e3;             % excitation f
exct.n_cycl = 1;                % no of cycles
% 
% [w indF] = min(abs(exct.f_exc - f_vec));
% wave.vph = 2*pi*f_vec(indF)./k_a(indF);     % phase velocity in m/s for excitation freq
% wave.lambda = wave.vph / f_vec(indF);       % wavelength in m for excitation freq


% create interpolated disp. curves to span exct.fs = 5.12e6;  exct.n_sampl = 4096; 
% create interpolated disp. curves to span exct.fs = 5.12e6;  exct.n_sampl = 4096; 
% Build frequency vector for simulation (Hz)
df = exct.fs / exct.n_sampl;
f_sim = (0:exct.n_sampl-1)' * df;        % column vector, 0..fs-df

% Original dispersion data use f_vec (Hz) and k_a, k_s (rad/m) from loaded file.
% Interpolate k vectors onto f_sim. Use NaN for out-of-range values then fill with nearest.
k_a_sim = interp1(f_A0*1e3, k_A0*1e3, f_sim, 'pchip', NaN);
k_s_sim = interp1(f_S0*1e3, k_S0*1e3, f_sim, 'pchip', NaN);

% Fill NaNs at ends by nearest valid value to avoid gaps in further processing
if any(isnan(k_a_sim))
    valid = find(~isnan(k_a_sim));
    if ~isempty(valid)
        k_a_sim(1:valid(1)-1) = k_a_sim(valid(1));
        k_a_sim(valid(end)+1:end) = k_a_sim(valid(end));
    end
end
if any(isnan(k_s_sim))
    valid = find(~isnan(k_s_sim));
    if ~isempty(valid)
        k_s_sim(1:valid(1)-1) = k_s_sim(valid(1));
        k_s_sim(valid(end)+1:end) = k_s_sim(valid(end));
    end
end

% Save interpolated vectors into variables expected later (f_vec, k_a, k_s)
% Overwrite f_vec to the simulation frequency grid and ensure column orientation
f_vec = f_sim;
k_a = k_a_sim;
k_s = k_s_sim;


% plot(f_vec,k_a, dispC.f_v, dispC.k_a)
% return
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%       Simulate  responses using structure's transfer function
%                          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d = 500;

[resp_a t sig] = dispResponse(f_vec,k_a,exct.fs,exct.n_sampl,exct.f_exc,exct.n_cycl,d);
[resp_s t sig] = dispResponse(f_vec,k_s,exct.fs,exct.n_sampl,exct.f_exc,exct.n_cycl,d);
resp_a = resp_a+ resp_s;
figure(1), plot(t,sig,t,resp_a), 
xlabel('time [s]')
legend('excitation','response')

%% perform STFT of the resulting signal 


%%
figure(323)
spectrogram(real(resp_a), 128,120,2048, exct.fs, 'onesided','yaxis')
hold on 
% plot group velocity
plot(d./vg_A0, f_A0/1e3)
plot(d./vg_S0, f_S0/1e3)
xlim([0 500])
ylim([0 5])
% plot(d/A0_vg, f_A0)
hold off

%%
% figure(3)
% plot(d./vg_A0/1e3, f_A0*1e3)
% xlim([0 800e-3])