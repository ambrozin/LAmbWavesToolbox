function im = plane_wave_DAS_sparse(array, Res, theta, r, wave, exct)
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

im.x = r*cosd(theta);
im.y = r*sind(theta);
im.C = zeros(size(im.x,2), size(im.x,1));
im.Cph = zeros(size(im.x,2), size(im.x,1));



for i = 1:length (theta)
    
    Res2 = zeros(size(array.u,1)*size(array.u,2),size(Res,3));
    for j = 1: size(array.u,1)
        % delay for receivers
        delayR = round( (array.X .* cosd(theta(i)) + array.Y .* sind(theta(i)))...
            * 1e-3 / wave.vph * exct.fs);
         delayT = round( (array.X(array.u(j,:)>0) .* cosd(theta(i)) + array.Y(array.u(j,:)>0) .* sind(theta(i)))...
            * 1e-3 / wave.vph * exct.fs);
        
        apodiz = array.u(j,array.u(j,:)>0) .* array.v(j,:); 
        
        for k = 1:length(delayR)
            Res2(k + (j-1).*size(Res,2),:) = circshift(Res(j,k,:).*apodiz(k),delayR(k) + delayT );
        end
        %    and sum
%        figure(232), plot(apodiz)
    end
    
     resp = sum(Res2);

    C_coh =std(angle( Res2.*exp(-1i*angle(mean(Res2)))));
    C_coh = 1 - C_coh; 
  	C_coh(C_coh<=0) = 1e-12;
        
%            figure(3), plot(real(Res2([4:6],:))'), drawnow
%            figure(3), plot(real(Res2)'), drawnow
        
        im.C(i,:) = (resp);
        im.coh_C(i,:) = C_coh;
        clc, disp(['calculating DAS: ' num2str(i/length (theta) *100), '%'])
end
return
%%

phs2 = Res2(:,1482); 

phs3 = phs2*exp(-1i*angle(mean(phs2)));
std(angle(phs2))


phs4 = Res2.*exp(-1i*angle(mean(Res2)));

figure(323) 
plot(std(angle(phs3)), '.')




