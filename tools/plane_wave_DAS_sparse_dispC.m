function im = plane_wave_DAS_sparse(array, Res, theta, r, wave, exct,disp_data)
% function im = plane_wave_DAS_sparse(array, Res, theta, r, wave, exct)
%
% function to perform beamforming based on plane wave assumption
%
%  im = plane_wave_DAS_2D(array, Res, theta, r, wave, exct)
%
% where:
%  DAS result: 
%         im.C         % discrete-time C
%         
% 
%         COHERENCE  - standard deviation of phase at subsequent points 
%         im.coh_C        % coherence calculated in discrete time 
%      
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

% dispersion removal 
for i = 1: size(Res,1)
    for j = 1 : size(Res,2)
        [out(i,j,:) h_dist] = disp_removal_Wilcox(squeeze(Res(i,j,:))',exct, disp_data) ;
    end
end

Res = out;
dx = h_dist(2) - h_dist(1); 
im.x = r*cosd(theta);
im.y = r*sind(theta);
im.C = zeros(size(im.x,2), size(im.x,1));
im.Cph = zeros(size(im.x,2), size(im.x,1));


for i = 1:length (theta)
    
    % DELAY 
    Res2 = zeros(size(array.u,1)*size(array.u,2),size(Res,3));
    for j = 1: size(array.u,1)  % for subsequent transmitters
        % delay for receivers 
%         delayR = round( (array.X .* cosd(theta(i)) + array.Y .* sind(theta(i)))...
%             * 1e-3 / wave.vph * exct.fs);

%        shift in spatial domain for receivers in [mm]
        delayR = round( (array.X .* cosd(theta(i))*1e-3 + array.Y .* sind(theta(i)) * 1e-3  )/dx);
%        shift in spatial domain for transmitters in [mm]        
        delayT = round( (array.X(array.u(j,:)>0) .* cosd(theta(i)) + array.Y(array.u(j,:)>0) .* sind(theta(i)))* 1e-3 / dx);
        
        
        apodiz = array.u(j,array.u(j,:)>0) .* array.v(j,:); 
        
        for k = 1:length(delayR)
            Res2(k + (j-1).*size(Res,2),:) = circshift(Res(j,k,:).*apodiz(k),delayR(k) + delayT );
        end
        %    and sum
%        figure(232), plot(apodiz)
    end
    % AND SUM 
    resp = sum(Res2);
    im.C(i,:) = (resp(1:size(im.x,1)));
    
    
    % phase coherence 
    % rotate phase to minimalize error resulting from phase wrapping 
    % calculate std 
    C_coh =std(angle( Res2.*exp(-1i*angle(mean(Res2)))));  
    
    C_coh = 1 - C_coh(1:size(im.x,1)); 
  	C_coh(C_coh<=0) = 1e-12;
    im.coh_C(i,:) = C_coh(1:size(im.x,1));
        
    clc, disp(['calculating DAS: ' num2str(i/length (theta) *100), '%'])
end
