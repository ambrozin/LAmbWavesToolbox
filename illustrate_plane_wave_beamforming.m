%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%           this script illustrates time-domain DAS beamforming based on
%           simulated ULA response.
%             
%           Additionally, it performes phase coherence imaging anaysing
%           standard deviation of phase shift 
%                   
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
%                          Generate arrau topology
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ULA  
dx = 6.5;               % pitch in mm
x = (1:10).*dx;         x = x - mean(x);
y = zeros(size(x));

X = x;                  Y = y; 

figure(2),  plot(X,Y,'s'), title('ULA'), xlabel('x [mm]')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%                           signal parameters 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

exct.fs = 5.12e6;               % sampling f 
exct.n_sampl = 4096;            % no of samples
exct.f_exc = 100e3;             % excitation f
exct.n_cycl = 5;                % no of cycles

[w indF] = min(abs(exct.f_exc - f_vec));
wave.vph = 2*pi*f_vec(indF)./k_a(indF);     % phase velocity in m/s for excitation freq
wave.lambda = wave.vph / f_vec(indF);       % wavelength in m for excitation freq

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%                           Target to be imaged
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

target.r = 300;
target.theta = 50;
target.x = target.r*cosd(target.theta); 
target.y = target.r*sind(target.theta); 

figure(2),  plot(X,Y,'s', target.x,target.y,'*'),  xlabel('x [mm]')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%       Simulate array responses using structure's transfer function
%                          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arrayRes = simulate_array_response(X,Y,target, exct,f_vec,k_a, k_s);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%                       DAS beamforming
%                          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theta = 0:360;                      % azimuth in deg
r = (1:exct.n_sampl)';

im = plane_wave_beamforming_ULA(X, arrayRes, theta, r, wave, exct) ;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%                      Plot the results 
%                          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(10)
pcolor(im.x, im.y, im.C'), shading flat, axis equal tight, axis([-1500, 1500, -1500,1500])
title('DAS result')
colorbar 

figure(11)
pcolor(im.x, im.y, im.Cph'), shading flat, axis equal tight, axis([-1500, 1500, -1500,1500])
title('weights based on Phase coherence ')
colorbar 

figure(12)
pcolor(im.x, im.y, im.Cph'.*im.C'), shading flat, axis equal tight, axis([-1500, 1500, -1500,1500])
title('Weighted image DAS and phase weights ')
colorbar 

figure(13)
plot(theta,max(im.C')./max(max(im.C')),  theta, max(im.Cph'.*im.C')./max(max(im.Cph'.*im.C')))
title('Beam patterns comparision')
legend('DAS','Phase Cohrerence')
colorbar