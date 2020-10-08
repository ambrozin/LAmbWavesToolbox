function im = plane_wave_DAS_ULA_disp_compens(array, Res, theta, r, wave, exct, disp_data) 
% function im = plane_wave_DAS_ULA(array, Res, theta, r, wave, exct)
%  
% function to perform beamforming based on plane wave assumption
% 
%  im = plane_wave_DAS_ULA(array, Res, theta, r, wave, exct)
% 
% where: 
% im - structure with 
%         im.C        - DAS result 
%         im.Cph      - standard deviation of phase at subsequent points 
%         im.x, im.y  - x,y vectors to scale the images 
%
% array - structure of array parameters 
%         array.X,    - coordinates of array elements
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
%
%% plane wave beamforming for ULA

% theta = 0:360;     % azimuth in deg
% r = (1:exct.n_sampl)';

% compensate for dispersion 
for i = 1: size(Res,1) 

    [out(i,:) h_dist] = disp_removal_Wilcox(Res(i,:),exct, disp_data) ; 

end

Res = out;
dx = (h_dist(2)  - h_dist(1) ) ; 

im.x = r*cosd(theta).*dx;
im.y = r*sind(theta).*dx; 
im.C = zeros(size(im.x,2), size(im.x,1));
im.Cph = zeros(size(im.x,2), size(im.x,1));
im.r = r.*dx; 

disp('calculating DAS.... ')

for i = 1:length (theta) 
   % delay 
%    delay = round( array.X .* cosd(theta(i)) * 1e-3 / wave.vph * exct.fs);
   delay = round( array.X .* cosd(theta(i))*1e-3 /dx );
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
   
   im.C(i,:) = abs(resp(:,1:size(im.x,1)));
%    im.Cph(i,:) = 1-(respPhs);
   im.coh_C  (i,:) = Cphs_coh(:,1:size(im.x,1));
   
%  figure(4), pcolor(im.x, im.y, im.C'), shading interp, drawnow
    clc, disp(['calculating DAS: ' num2str(i/length (theta) *100), '%'])
end
disp('done')