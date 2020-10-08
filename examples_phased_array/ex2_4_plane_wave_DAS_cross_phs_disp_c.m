%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%           this script illustrates time-domain DAS beamforming based on
%           simulated square array response.
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

disp_data.k_a = k_a; 
disp_data.k_s = k_s; 
disp_data.f_vec = f_vec; 
% figure(1),  plot(f_vec,k_a,f_vec,k_s)
% xlabel('f [Hz]'), ylabel('k [rad/m]')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%                          Generate array topology
% 
%                                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%    cross   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
array.dx = 6.5;               % pitch in mm
x = (1:10).*array.dx;         x = x - mean(x);
y = zeros(size(x));

% [x,y] = meshgrid(x,x); 
array.X = [x y];
array.Y = [y x];
array.u = [diag(ones(1,10)) zeros(10,10)];
array.v = [ zeros(10,10)     ones(10,10)];
% array.Y = reshape (y,1,size(y,1)*size(y,2));
array.name = 'cross';
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
exct.n_cycl = 5;                % no of cycles

[w indF] = min(abs(exct.f_exc - f_vec));
wave.vph = 2*pi*f_vec(indF)./k_a(indF);     % phase velocity in m/s for excitation freq
wave.lambda = wave.vph / f_vec(indF);       % wavelength in m for excitation freq

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%                           Target to be imaged
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

target.r = [300 300];
target.theta = [-25 25];
target.reflectivity = [1 0.9]; 
target.x = target.r.*cosd(target.theta); 
target.y = target.r.*sind(target.theta); 

figure(2),  plot(array.X,array.Y,'s', target.x,target.y,'*'),  xlabel('x [mm]')
 axis equal
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%       Simulate array responses using structure's transfer function
%                          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arrayRes = simulate_sparse_array_response(array,target, exct,f_vec,k_a, k_s, 'a0');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%                       DAS beamforming
%                          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theta = 0:360;                      % azimuth in deg
r = (1:exct.n_sampl)';

% im = plane_wave_DAS_sparse(array, arrayRes, theta, r, wave, exct) ;
im = plane_wave_DAS_sparse_dispC(array, arrayRes, theta, r, wave, exct,disp_data) ;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%                      Plot the results 
%                          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(10)
pcolor(im.x, im.y, abs(im.C)'), shading flat, axis equal tight, axis([-1500, 1500, -1500,1500])
title('DAS result')
colorbar 

figure(11)
pcolor(im.x, im.y, im.coh_C'), shading flat, axis equal tight, axis([-1500, 1500, -1500,1500])
title('weights based on Phase coherence ')
colorbar 


figure(12)
pcolor(im.x, im.y, im.coh_C'.*abs(im.C)'), shading flat, axis equal tight, axis([-1500, 1500, -1500,1500])
title('Weighted image DAS and phase weights ')
colorbar 

BP1 = imagePolar2BP(abs(im.C), 'norm');
BP2 = imagePolar2BP((im.coh_C'.*abs(im.C)')','norm'); 
figure(13)
plot(theta,BP1, theta, BP2)


return

figure(13)
plot(theta,20*log10(max(im.C')./max(max(im.C'))),  theta, 20*log10(max(im.Cph'.*im.C')./max(max(im.Cph'.*im.C'))))
title('Beam patterns comparision')
legend('DAS','Phase Cohrerence')
colorbar