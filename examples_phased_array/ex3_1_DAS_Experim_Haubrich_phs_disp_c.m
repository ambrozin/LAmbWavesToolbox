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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

error('This script does not work yet')

load('/Users/LukaszA/Dropbox/Lukasz (1)/M_Malarz/Clean/experimental data/Haubrich1.svd.mat')
pos = VibData.XYZ;
x = pos(:,1)';
y = pos(:,2)';
array.X = x; 
array.Y = y; 

array.u = [diag(ones(1,10)) ];
array.v = [ ones(10,10)];
load('/Users/LukaszA/Dropbox/Lukasz (1)/M_Malarz/Clean/experimental data/HaubrichFMC.mat')
array.name = 'Haubrich';

% plot array elements 
figure(2),  plot(array.X,array.Y,'s'), 
xlabel('x [mm]'), ylabel('y [mm]')
title(array.name)

return
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%                           signal parameters 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

exct.fs = 1/(VibData.t(2) -  VibData.t(1)) ;            % sampling f 
exct.n_sampl = 745;            % no of samples
exct.f_exc = 100e3;             % excitation f
exct.n_cycl = 5;                % no of cycles

[w indF] = min(abs(exct.f_exc - f_vec));
wave.vph = 2*pi*f_vec(indF)./k_a(indF);     % phase velocity in m/s for excitation freq
wave.lambda = wave.vph / f_vec(indF);       % wavelength in m for excitation freq


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%                       DAS beamforming
%                          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theta = 0:360;                      % azimuth in deg
r = (1:exct.n_sampl)';

% im = plane_wave_DAS_sparse(array, arrayRes, theta, r, wave, exct) ;
im = plane_wave_DAS_sparse_dispC(array, Response, theta, r, wave, exct,disp_data) ;

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