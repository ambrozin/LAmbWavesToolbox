 addpath('disp_curves_tool')
%%
clear all 
% close all

load('disp_curve')

% figure(1)
% plot(f_vec,k_a,f_vec,k_s)
% xlabel('f [Hz]')
% ylabel('k [rad/m]')
%%  Topologia
% ISSUE: replace this  this with a more general function 


% kwadrat 
dx = 6.5; 
x = (1:10).*dx;         x = x - mean(x);
y = (1:10).*dx;         y = y - mean(y);
y = zeros(size(x));

X = x; 
Y = y; 

[X Y]  = meshgrid(x,y);
X= reshape(X,1, length(x)*length(y));
Y= reshape(Y,1,length(x)*length(y));


sel_topol = 8; 
predefined_coarrays;

% X = x; Y = y; 

figure(2),
plot(X,Y,'.'), axis equal tight
%% signal parameters 
exct.fs = 5.12e6;
exct.n_sampl = 4096; 

exct.f_exc = 100e3; 
exct.n_cycl = 3; 

[w indF] = min(abs(exct.f_exc - f_vec));
wave.vph = 2*pi*f_vec(indF)./k_a(indF);     % phase velocity in m/s
wave.lambda = wave.vph / f_vec(indF);

%% target 

target.r = 500;
target.theta = 30;
target.x = target.r*cosd(target.theta); 
target.y = target.r*sind(target.theta); 

target2.r = 500;
target2.theta = 200;
target2.x = target2.r*cosd(target2.theta); 
target2.y = target2.r*sind(target2.theta); 

figure(2),
plot(X,Y,'.', target.x,target.y,'*')

%% simulate dispersive responses 


Res = zeros(exct.n_sampl,size(X,2),size(X,1)) ;
for i = 1:size(X,2)
    for j = 1:size(X,1)
        d = sqrt((X(i)-target.x).^2+(Y(i)-target.y).^2); 
        d2 = sqrt((X(i)-target2.x).^2+(Y(i)-target2.y).^2);
        resp_a = dispResponse(f_vec,k_a,exct.fs,exct.n_sampl,exct.f_exc,exct.n_cycl,d);
        resp_a2 = dispResponse(f_vec,k_a,exct.fs,exct.n_sampl,exct.f_exc,exct.n_cycl,d2);
%         resp_s = dispResponse(f_vec,k_s,exct.fs,exct.n_sampl,exct.f_exc,exct.n_cycl,d);
        Res(:,i,j) =  hilbert(resp_a + resp_a2  );
    end
    disp(i/size(X,2)*100)
end
Res = reshape(Res, exct.n_sampl, size(X,2)* size(X,1)); 
%% plane wave beamforming for ULA


theta = 0:360;     % azimuth in deg
r = (1:exct.n_sampl)';

im.x = r*cosd(theta);
im.y = r*sind(theta); 
im.C = zeros(size(im.x,2), size(im.x,1));
im.Cph = zeros(size(im.x,2), size(im.x,1));
im.respPhs = zeros(size(im.x,2), size(im.x,1));
for i = 1:length (theta) 
%    dist = abs(-tand(90-theta(i)).*X + Y) ./ sqrt( (-tand(90-theta(i))).^2 + (1).^2); 
   
%    delay = round( X .* cosd(theta(i)) * 1e-3 / wave.vph * exct.fs);
   
   delay = round( (X .* cosd(theta(i)) + Y .* sind(theta(i))) * 1e-3 / wave.vph * exct.fs);
   
   
   Res2 = zeros(size(Res));
   for j = 1:length(delay) 
        Res2(:,j) = circshift(Res(:,j),delay(j)); 
   end
   resp    = sum(Res2');
   respPhs = std((angle(Res2)'));
   
   respPhs = 1-sqrt( var(cos(angle(Res2)')) + var(sin(angle(Res2)')) );
%    respPhs(respPhs < 0 ) = 0; 
%    figure(3), plot(respPhs), 
%    xlim([500 1500])
   
%    drawnow
   i
   im.C(i,:) = abs(resp);
   im.respPhs(i,:) = respPhs;
end
%%

im.Cph = im.respPhs;
% im.Cph(im.Cph<0) = 0; 
figure(10)
pcolor(im.x, im.y, 20*log10(im.C'./max(max(im.C)))), shading flat, axis equal tight, axis([-1500, 1500, -1500,1500])
colorbar
caxis([-50 0])
figure(11)
pcolor(im.x, im.y, im.respPhs'./max(max(im.respPhs'))), shading flat, axis equal tight, axis([-1500, 1500, -1500,1500])
colorbar

resIm = im.Cph'.*im.C' ; 
resIm = resIm ./ max(max(resIm)); 
figure(12)
pcolor(im.x, im.y, 20*log10(abs(resIm))), shading flat, axis equal tight, axis([-1500, 1500, -1500,1500])
colorbar
caxis([-50 0])

figure(13)
plot(theta,20*log10(max(im.C')./max(max(im.C'))),  theta, 20*log10(max(im.Cph'.*im.C')./max(max(im.Cph'.*im.C'))))
xlabel('[\circ]')
ylabel('[dB]')