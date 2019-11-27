%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           this script illustrates time-domain DAS beamforming based on
%           simulated ULA response.
%
%           The script permits using selected elements as transmitter /
%           receivers
%
%           Additionally, it performes phase coherence imaging anaysing
%           standard deviation of phase shift
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% error('This code is not finished yet!')
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
x = (1:9).*array.dx;          x = x - mean(x);
y = zeros(size(x));

array.X = x;                  array.Y = y;
array.name = 'ULA';

% select transmitting / receiving elements - comment suitable lines

%____________________ case 1 - two firings___________________________
array.u = zeros(2,numel(x));  % two firings
array.u(1) = 1; array.u(end) = 1;   % only outermost transmitters
array.v = ones(2,numel(x));   % all receivers


% w = hanning(17);
% w(9) = 0.5;
% array.v(2,:) = w(9:end);
% array.v(1,:) = w(1:9);
% array.v(2,:) = [4.5 8:-1:1];
% array.v(1,:) = [1: 8 4.5 ];

%___________________ case 2 - TFM
%
% array.u = zeros(numel(x),numel(x));  % two firings
% array.u = eye(numel(x)) ;  % only outermost transmitters
% array.v = ones(numel(x),numel(x));   % all receivers


%___________________ case 3    STMR
%         (single transmit multiple receivers)
%
% array.u = zeros(1,numel(x));  % one firing
% array.u(5) = 1;
% array.v = ones(size(array.u ));   % all receivers

%___________________ case 4   all transmit 1 receive
% array.u = eye(numel(x)) ;
% array.v = zeros(size(array.u ));
% array.v(:,5) = 1;
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

target.r = 300;
target.theta = 45;
target.x = target.r*cosd(target.theta);
target.y = target.r*sind(target.theta);

figure(2),  plot(array.X,array.Y,'s', target.x,target.y,'*'),  xlabel('x [mm]')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       Simulate array responses using structure's transfer function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

arrayRes = simulate_sparse_array_response(array,target, exct,f_vec,k_a, k_s);


% returnr
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                       DAS beamforming
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


theta = 0:360;                      % azimuth in deg
r = (1:exct.n_sampl)';

im = plane_wave_DAS_sparse(array, arrayRes , theta, r, wave, exct) ;

figure(10),
pcolor(im.x, im.y, abs(im.C)'), shading flat, axis equal tight, axis([-1500, 1500, -1500,1500])
title('DAS result')
colorbar

figure(11),
BP = imagePolar2BP(im.C,'norm') ;
plot(theta, 20*log10(BP))



hold off
return
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                      Plot the results
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% figure(9)
% pcolor(im.x, im.y, im.Cphs'), shading flat, axis equal tight, axis([-1500, 1500, -1500,1500])
% title('DAS result')
% colorbar

figure(10), subplot(1,2,1)
pcolor(im.x, im.y, im.C'), shading flat, axis equal tight, axis([-1500, 1500, -1500,1500])
title('DAS result')
colorbar


figure(10), subplot(1,2,2)
BP = imagePolar2BP(im.C,'norm') ;
plot(theta, 20*log10(BP))

return

figure(11)
pcolor(im.x, im.y, im.Cph'), shading flat, axis equal tight, axis([-1500, 1500, -1500,1500])
title('weights based on Phase coherence ')
colorbar

figure(12)
pcolor(im.x, im.y, im.Cph'.*im.C'), shading flat, axis equal tight, axis([-1500, 1500, -1500,1500])
title('Weighted image DAS and phase weights ')
colorbar

figure(13)
plot(theta,20*log10(max(im.C')./max(max(im.C'))),  theta, 20*log10(max(im.Cph'.*im.C')./max(max(im.Cph'.*im.C'))))
% title('Beam patterns comparision')
box off
hlgd = legend('DAS','Phase Cohrerence');
set(hlgd, 'location', 'northoutside', 'orientation', 'horizontal', 'box', 'off')
xlabel('angle [\circ]'), xlim([0 360])