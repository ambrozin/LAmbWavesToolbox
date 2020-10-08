function [Y,t]=remove_dispersion(signal,fs,k,cg,f_k, f0)
% function [Y,t]=remove_dispersion(signal,fs,k,f_k,f0)
% 
%
% k   - [1/m]  wavenumber
% cg  - [m/s]  group velocity ;
% f_k - [Hz]    frequency where group velocity and wavenumber is specified
% f0 -  [Hz]      central frequency
% 

%make signal even numbers of samples

    parzysta=mod(length(signal)-1,2) ;          
    if ~parzysta
        signal(end+1)=0;
    end
    
%interpolate wavenumber and group velocity for frequencies specified for
%spectrum of input signal.

    f= fs*linspace(0,1,(length(signal)));
    f=f(1:end/2+1);
    
    k_int=interp1(f_k,k,f);
    k_int(end)=0;
  
    
    cg_int=interp1(f_k,cg,f);
  
%wavenumber linearisation for f0 

     [w i_f0]=min(abs(f-f0));
     k_lin=k_int(i_f0)+2*pi*(f-f0)./cg_int(i_f0);
%      cg_int(i_f0)

%operations on signal spectrum

     G=fft(signal) ;
     G1=G(1:end/2) ; 
     G1_klin=interp1(k_int(1:end-1),G1,k_lin(1:end-1)); 

     G2=fliplr(G(end/2+1:end));
     G2_klin=interp1(k_int(2:end),G2,k_lin(2:end)); 

     G_klin=[G1_klin fliplr(G2_klin)];
     G_klin(isnan(G_klin))=0;
     G_klin(1)=0;
 
 %output signals
 
     Y=ifft(G_klin);
     t=0:1/fs:(length(Y)-1)*1/fs;
    
%  plot(ifft(G_klin))
    