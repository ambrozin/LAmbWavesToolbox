clear all 
close all 
%%

load("dispersion_data_4_mm_498mm.mat")

%%

figure(1)
plot(t,y)
y(1,:) = y(1,:) - mean(y(1,:)); 
y(2,:) = y(2,:) - mean(y(2,:)); 
[R lags] = xcorr(y(2,:),y(1,:)); 
R = R(end/2:end); 
lags = lags(end/2:end);
figure(2)
plot(R)
%
 load("alu_4mm_DispersionCurves.mat")
% figure(1),  plot(f_vec,k_a,f_vec,k_s)
% xlabel('f [Hz]'), ylabel('k [rad/m]')
%

%%
f_A0 = A{:,1};      % frequency of A0 mode
v_A0 = A{:,2};      % phase velcity of A0 mode
vg_A0 = A{:,3};      % phase velcity of A0 mode
k_A0 = A{:,10};
f_A1 = A{:,12};      % frequency of A1 mode
v_A1 = A{:,13};      % phase velcity of A1 mode
vg_A1 = A{:,14};      % phase velcity of A1 mode
f_S0 = S{:,1};      % frequency of S0 mode
v_S0 = S{:,2};      % phase velcity of S0 mode
k_S0 = S{:,10};
vg_S0 = S{:,3};      % phase velcity of A0 mode

plot(f_A0,v_A0, f_A1, v_A1, f_S0, v_S0)
xlabel('frequency (kHz)')
ylabel('Phase velocity (km/s)')

%%
figure(3)
spectrogram(R(1:4000), 128, 127, 4096,Fs, "onesided", "yaxis")
caxis([-60 0]+max(caxis))
d = 510;

hold on 
% plot group velocity
plot(d./vg_A0, f_A0/1e3/4)
plot(d./vg_S0, f_S0/1e3/4)
plot(d./vg_A1, f_A1/1e3/4)
xlim([0 500])
ylim([0 1.5])

hold off

%% apply reassigned spectrogram
d = 515; 
figure(4)
spectrogram(R(1:4000), 256, 250, 4096,Fs, "onesided", "yaxis", "reassigned")
caxis([-60 0]+max(caxis))
hold on 
% plot group velocity
plot(d./vg_A0, f_A0/1e3/4)
plot(d./vg_S0, f_S0/1e3/4)
plot(d./vg_A1, f_A1/1e3/4)
xlim([0 500])
ylim([0 1.5])

hold off

%% Hilbert-Huang Transform using built-in MATLAB functions
d = 515;
signal = double(R(1:4000));
dt = 1/Fs;
time_vec = (0:length(signal)-1) * dt;

% EMD decomposition (requires Signal Processing Toolbox)
[imf, residual] = emd(signal, 'MaxNumIMF', 8 );

% HHT spectrum (figure 5)
figure(5)
hht(imf, Fs)
title('Hilbert-Huang Transform Spectrum')

% EMD components (figure 6)
figure(6)
num_imfs = size(imf, 2);
for i = 1:num_imfs
    subplot(num_imfs+1, 1, i)
    plot(time_vec*1000, imf(:,i))
    ylabel(['IMF ' num2str(i)])
    grid on
end
subplot(num_imfs+1, 1, num_imfs+1)
plot(time_vec*1000, residual)
ylabel('Residual')
xlabel('Time (ms)')
grid on
sgtitle('EMD Components')

%% Wavelet Synchrosqueezed Transform
figure(7)
% wsst(double(R(1:4000)), Fs)

twin = hamming(511);
fwin = hamming(511);
[d f tt] = wvd(R(1:4000),Fs,"smoothedPseudo", twin, fwin);
imagesc(tt,f,mag2db(d))
axis xy
xlabel('Time [s]')
ylabel('Frequency [Hz]')
title('Wigner-Ville Distribution')
colorbar
caxis([- 100 0]+max(caxis))
