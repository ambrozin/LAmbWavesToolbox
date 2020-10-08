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
clear all, close all

try, addpath('..\tools');   catch, addpath('../tools') , end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%                  Load the dispersion curves data
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try 
    load('..\dispersion_data\disp_curve.mat')  % disp curves for 2 mm alu plates 
catch 
    load('../dispersion_data/disp_curve.mat')  % disp curves for 2 mm alu plates
end

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

exct.fs = 5.12e6;               % sampling f 
exct.n_sampl = 4096;            % no of samples
exct.f_exc = 100e3;             % excitation f
exct.n_cycl = 2;                % no of cycles

[w indF] = min(abs(exct.f_exc - f_vec));
wave.vph = 2*pi*f_vec(indF)./k_a(indF);     % phase velocity in m/s for excitation freq
wave.lambda = wave.vph / f_vec(indF);       % wavelength in m for excitation freq

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%                           Targets to be imaged
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

target.r = [300 200];
target.theta = [90 60];
target.reflectivity = [1 0];
target.x = target.r.*cosd(target.theta); 
target.y = target.r.*sind(target.theta); 

figure(2),  plot(array.X,array.Y,'s', target.x,target.y,'*'),  xlabel('x [mm]')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%       Simulate array responses using structure's transfer function
%                          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arrayRes = simulate_array_response(array.X,array.Y,target, exct,f_vec,k_a, k_s,'a0');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%                       DAS beamforming
%                          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theta = 0:360;                      % azimuth in deg
r = (1:exct.n_sampl)';

sig  = squeeze(arrayRes(1,:)); 
disp.k_a = k_a; 
disp.k_s = k_s; 
disp.f_vec = f_vec; 

% [out h_dist] = disp_removal_Wilcox(sig,exct, disp);

im = plane_wave_DAS_ULA_disp_compens(array, arrayRes, theta, r, wave, exct, disp) ;
% return
% im = plane_wave_DAS_ULA(array, arrayRes, theta, r, wave, exct) ;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%                      Plot the results 
%                          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scale_axis = [-500 500 -500 500]*1e-3; 


figure(10)
pcolor(im.x, im.y, im.C'), shading flat, axis equal tight, axis(scale_axis)
title('DAS result')
xlabel('x [m]')
ylabel('y [m]')
colorbar 

% return
figure(11)
pcolor(im.x, im.y, im.coh_C'), shading flat, axis equal tight, axis(scale_axis)
title('weights based on Phase coherence ')
colorbar 
xlabel('x [m]')
ylabel('y [m]')

figure(12)
pcolor(im.x, im.y, im.coh_C'.*im.C'), shading flat, axis equal tight, axis(scale_axis)
title('Weighted image DAS and phase weights ')
colorbar
xlabel('x [m]')
ylabel('y [m]')

figure(13)
plot(theta,20*log10(max(im.C')./max(max(im.C'))),  theta, 20*log10(max(im.coh_C'.*im.C')./max(max(im.coh_C'.*im.C'))))
% title('Beam patterns comparision')
box off
hlgd = legend('DAS','Phase Cohrerence'); 
set(hlgd, 'location', 'northoutside', 'orientation', 'horizontal', 'box', 'off')
xlabel('angle [\circ]'), xlim([0 360])
ylabel('[dB]')