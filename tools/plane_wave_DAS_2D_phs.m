function im = plane_wave_DAS_2D_phs(array, Res, theta, r, wave, exct) 
% function im = plane_wave_DAS_2D_phs(array, Res, theta, r, wave, exct)
%  
% function to perform beamforming based on plane wave assumption
% 
%  im = plane_wave_DAS_2D_phs(array, Res, theta, r, wave, exct)
% 
% where: 
% im - structure with 
%                        
%         DAS result: 
%         im.C         % discrete-time C
%         im.Cphs      % phase-shift C
% 
%         COHERENCE  - standard deviation of phase at subsequent points 
%         im.coh_C        % coherence calculated in discrete time 
%         im.coh_Cphs     % coherence calculated using phase shift
%         
%         GEOMETRY
%         im.x, im.y  - x,y vectors to scale the images 
%
% array - structure of array parameters 
%         array.X, array.Y - coordinates of array elements
%         Res         - inpust signals from array elements 
%         theta       - vector angles to evaluate beamforming 
%         r           - vector of radial distances 
%         wave        - structure of wave parameters 
%         wave.vph    - phase velocity
%         exct        - excitation parameters 
%         exct.fs     - sampling frequency 
%
%
%        2019 ambrozin@agh.edu.pl


%% plane wave beamforming for a 2D shaped array

% theta = 0:360;     % azimuth in deg
% r = (1:exct.n_sampl)';

im.x = r*cosd(theta);
im.y = r*sind(theta); 
im.C = zeros(size(im.x,2), size(im.x,1));           % discrete-time C
im.Cphs = zeros(size(im.x,2), size(im.x,1));        % phase-shift C
im.coh_C = zeros(size(im.x,2), size(im.x,1));       % coherence calculated in discrete time 
im.coh_Cphs = zeros(size(im.x,2), size(im.x,1));    % coherence calculated using phase shift

for i = 1:length (theta) 
   % delay 
   delay = round( (array.X .* cosd(theta(i)) + array.Y .* sind(theta(i))) * 1e-3 / wave.vph * exct.fs);
   Res2 = zeros(size(Res));
   for j = 1:length(delay) 
        Res2(j,:) = circshift(Res(j,:),delay(j)); 
   end
%    and sum
   resp = sum(Res2);
   C_coh = std(Res2);
   
    %% phase  shift approach 
   % calculate delay in radians 
   delayRad = 2*pi*( (array.X .* cosd(theta(i)) + array.Y .* sind(theta(i))) * 1e-3 / wave.vph * exct.f_exc);
%    delayRad = 2*pi*( array.X .* cosd(theta(i)) * 1e-3 / wave.vph * exct.f_exc);
   % apply the phase shift 
   Res3 = diag(exp(-1i.*delayRad))*Res;
   Cphs_coh = std(Res3);
%    figure(3), plot(real(Res3)'), xlim([500 1200]), drawnow
   
   im.C(i,:) = abs(sum(Res2));
   im.coh_C(i,:) = 1-(C_coh);
   
   im.Cphs(i,:) = abs(sum(Res3));
   im.coh_Cphs(i,:) = 1-(Cphs_coh);
   
   
   clc, disp(['calculating DAS: ' num2str(i/length (theta) *100), '%'])
end