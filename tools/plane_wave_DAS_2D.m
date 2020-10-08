function im = plane_wave_DAS_2D(array, Res, theta, r, wave, exct) 
% function im = plane_wave_DAS_2D(array, Res, theta, r, wave, exct)
%  
% function to perform beamforming based on plane wave assumption
% 
%  im = plane_wave_DAS_2D(array, Res, theta, r, wave, exct)
% 
% where: 
% im - structure with 
%         im.C        - DAS result 
%         im.Cph      - standard deviation of phase at subsequent points 
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
im.C = zeros(size(im.x,2), size(im.x,1));
im.Cph = zeros(size(im.x,2), size(im.x,1));

for i = 1:length (theta) 
   % delay 
   delay = round( (array.X .* cosd(theta(i)) + array.Y .* sind(theta(i))) * 1e-3 / wave.vph * exct.fs);
   Res2 = zeros(size(Res));
   for j = 1:length(delay) 
        Res2(j,:) = circshift(Res(j,:),delay(j)); 
   end
%    and sum
   resp = sum(Res2);
   Cphs_coh =std(angle( Res2.*exp(-1i*angle(mean(Res2)))));
   Cphs_coh = 1 - Cphs_coh; 
   Cphs_coh(Cphs_coh<=0) = 1e-12;
   
   
%    figure(3), plot(real(Res2)'), drawnow
   
   im.C(i,:) = abs(resp);
   im.Cph(i,:) = Cphs_coh;
   clc, disp(['calculating DAS: ' num2str(i/length (theta) *100), '%'])
end