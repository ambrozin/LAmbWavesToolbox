function im = plane_wave_beamforming_ULA(X, Res, theta, r, wave, exct) 
% function im = plane_wave_beamforming(X, theta, r, wave, exct) 
%  
% funtion to perform beamforming based on plane wave assumptions 


%% plane wave beamforming for ULA

% theta = 0:360;     % azimuth in deg
% r = (1:exct.n_sampl)';

im.x = r*cosd(theta);
im.y = r*sind(theta); 
im.C = zeros(size(im.x,2), size(im.x,1));
im.Cph = zeros(size(im.x,2), size(im.x,1));

for i = 1:length (theta) 
    
   delay = round( X .* cosd(theta(i)) * 1e-3 / wave.vph * exct.fs);
   Res2 = zeros(size(Res));
   for j = 1:length(delay) 
        Res2(j,:) = circshift(Res(j,:),delay(j)); 
   end
   resp = sum(Res2);
   respPhs = std(Res2);
   
%    figure(3), plot(real(Res2)'), drawnow
   
   im.C(i,:) = abs(resp);
   im.Cph(i,:) = 1-(respPhs);
end