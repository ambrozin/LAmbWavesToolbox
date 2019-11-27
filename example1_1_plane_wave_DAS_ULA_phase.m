%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%           this script illustrates time-domain DAS beamforming based on
%           simulated ULA response.
%             
%           DAS is implemented in two ways: 
%           1  - shift of digital signal         
%           2  - phase shift of analitical signal
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('tools')

clear all 
close all
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%                  Load the dispersion curves data
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('disp_curve')  % disp curves for 2 mm alu plates 

% figure(1),  plot(f_vec,k_a,f_vec,k_s)
% xlabel('f [Hz]'), ylabel('k [rad/m]')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%                          Generate array topology
% 
%                                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%    ULA   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
array.dx = 6.5;               % pitch in mm
x = (1:10).*array.dx;         x = x - mean(x);
y = zeros(size(x));

array.X = x;                  array.Y = y; 
array.name = 'ULA';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot array elements 
figure(2),  plot(array.X,array.Y,'s'), 
xlabel('x [mm]'), ylabel('y [mm]')
title(array.name)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%                           signal parameters 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

exct.fs = 2.5e6;               % sampling f 
exct.n_sampl = 2048;            % no of samples
exct.f_exc = 100e3;             % excitation f
exct.n_cycl = 5;                % no of cycles

[mCphs indF] = min(abs(exct.f_exc - f_vec));
wave.vph = 2*pi*f_vec(indF)./k_a(indF);     % phase velocity in m/s for excitation freq
wave.lambda = wave.vph / f_vec(indF);       % wavelength in m for excitation freq

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%                           Target to be imaged
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

target.r = 300;
target.theta = 90;      % azimuth
target.x = target.r*cosd(target.theta); 
target.y = target.r*sind(target.theta); 

figure(2),  plot(array.X,array.Y,'s', target.x,target.y,'*'),  xlabel('x [mm]')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%       Simulate array responses using structure's transfer function
%                          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arrayRes = simulate_array_response(array.X,array.Y,target, exct,f_vec,k_a, k_s);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%                       DAS beamforming
%                          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theta = 0:1:360;                      % azimuth in deg
r = (1:exct.n_sampl)';

im = plane_wave_DAS_ULA_phs(array, arrayRes, theta, r, wave, exct) ;




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%                      Plot the results 
%                          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axisSize = [-500, 500, -500, 500];
figure(10), subplot(2,2,1)
pcolor(im.x, im.y, im.Cphs'), shading flat, axis equal tight, axis(axisSize)
title('Phase shift DAS')
colorbar 

figure(10), subplot(2,2,2)
pcolor(im.x, im.y, im.C'), shading flat, axis equal tight, axis(axisSize)
title('Discrete-time DAS')
colorbar 

% beam pattern 
[mC indC] = max(max(im.C)); 
[mCphs indCphs] = max(max(im.Cphs)); 

figure(10), subplot(2,2,3:4)
plot(theta, 20*log10(im.C(:,indC)./mC), theta, 20*log10(im.Cphs(:,indCphs)./mCphs))
xlabel('angle [\circ]'), xlim([0 360])
ylabel('BP [dB]'), ylim([-40 0])
hlgd = legend('Discrete-time DAS','Phase shift DAS'); 
set(hlgd, 'location', 'northoutside', 'orientation', 'horizontal', 'box', 'off')

%%
figure(11), subplot(2,2,1)
pcolor(im.x, im.y, im.C'), shading flat, axis equal tight, axis(axisSize)
title('Discrete-time DAS')
colorbar 

figure(11), subplot(2,2,2)
pcolor(im.x, im.y, im.coh_C'), shading flat, axis equal tight, axis(axisSize)
title('Coherence of discr. time DAS')
colorbar 

coh_C_image = im.coh_C.* im.C; 
[m_coh_C_image indcoh_C_image] = max(max(coh_C_image)); 
BP_coh_C_image = coh_C_image(:,indcoh_C_image); 

figure(11), subplot(2,2,3)
pcolor(im.x, im.y, coh_C_image'), shading flat, axis equal tight, axis(axisSize)
title('Coherence of discr. time DAS')
colorbar 

figure(11), subplot(2,2,4)
plot(theta, 20*log10(im.C(:,indC)./mC),theta, 20*log10(BP_coh_C_image./max(BP_coh_C_image)))
xlabel('angle [\circ]'), xlim([0 360])
ylabel('BP [dB]'), ylim([-40 0])
hlgd = legend('Discrete-time DAS','Phase coherence weighted'); 
set(hlgd, 'location', 'northoutside', 'orientation', 'horizontal', 'box', 'off')

%%
figure(12), subplot(2,2,1)
pcolor(im.x, im.y, im.Cphs'), shading flat, axis equal tight, axis(axisSize)
title('Phase-shift DAS')
colorbar 

figure(12), subplot(2,2,2)
pcolor(im.x, im.y, im.coh_Cphs'), shading flat, axis equal tight, axis(axisSize)
title('Coherence of phase-shift DAS')
colorbar 

coh_Cphs_image = im.coh_Cphs.* im.Cphs; 
[m_coh_Cphs_image ind_coh_Cphs_image] = max(max(coh_Cphs_image)); 
BP_coh_Cphs_image = coh_C_image(:,ind_coh_Cphs_image)./m_coh_Cphs_image; 

figure(12), subplot(2,2,3)
pcolor(im.x, im.y, coh_Cphs_image'), shading flat, axis equal tight, axis(axisSize)
title('Coherence of discr. time DAS')
colorbar 

figure(12), subplot(2,2,4)
plot(theta, 20*log10(im.Cphs(:,indCphs)./mCphs),theta, 20*log10(BP_coh_Cphs_image./max(BP_coh_Cphs_image)))
xlabel('angle [\circ]'), xlim([0 360])
ylabel('BP [dB]'), ylim([-40 0])
hlgd = legend('Discrete-time DAS','Phase coherence weighted'); 
set(hlgd, 'location', 'northoutside', 'orientation', 'horizontal', 'box', 'off')