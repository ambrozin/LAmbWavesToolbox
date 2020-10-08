 addpath('disp_curves_tool')
%%
clear all 
close all

load('disp_curve')

% figure(1)
% plot(f_vec,k_a,f_vec,k_s)
% xlabel('f [Hz]')
% ylabel('k [rad/m]')
%%  Topologia
% ISSUE: replace this  this with a more general function 


% kwadrat 
dx = 5; 
x = (1:10).*dx;         x = x - mean(x);
y = (1:10).*dx;         y = y - mean(x);
y = zeros(size(x));

X = x; 
Y = y; 

% [X Y]  = meshgrid(x,y);
% X= reshape(X,1, length(x)*length(y));
% Y= reshape(Y,1,length(x)*length(y));

figure(2),
plot(X,Y,'.')
%% signal parameters 
exct.fs = 2.56e6;
exct.n_sampl = 2048; 

exct.f_exc = 100e3; 
exct.n_cycl = 3; 

[w indF] = min(abs(exct.f_exc - f_vec));
wave.vph = 2*pi*f_vec(indF)./k_a(indF);     % phase velocity in m/s
%% simulate dispersive responses 
target.x = 500; 
target.y = 500; 

Res = zeros(size(X,1), size(X,2), exct.n_sampl);
for i = 1:size(X,1)
 
    for j = 1:size(X,2)
        
        d = sqrt((X(i,j)-target.x).^2+(Y(i,j)-target.y).^2); 
        resp_a = dispResponse(f_vec,k_a,exct.fs,exct.n_sampl,exct.f_exc,exct.n_cycl,d);
        Res(i,j,:) =  resp_a;
    end
    disp(i/size(X,1)*100)
end

figure(3) 
plot(squeeze(Res)')
%% plane wave beamforming for ULA


theta = 0:360;     % azimuth in deg
r = 0:1000; 

[im.x im.y] = meshgrid(theta,r); 
im.C = zeros(size(im.x,1), size(im.x,2));

for i = 45%:length (theta) 
    
  
   delay = round( X .* sind(theta(i)) * 1e-3 / wave.vph * exct.fs);
    
    
   return
end




