%%
%
%       this script probaly does not work.
%%
%% input data
    meas=(pas_measurement(3).measurement(1600:2800));
    meas=[meas zeros(1,8048-length(meas))];
%     meas=(paq_measurement(2).measurement(5706:7000));
%     meas=[meas zeros(1,4048-length(meas))];
    
%     meas=(paq_measurement(1).measurement(1:3000));
%     meas=[meas zeros(1,4048-length(meas))];
    fs=2.5e6;

%% remove dispersion
%     meas=wavelet_filter(meas,100e3,fs);
    [Y,X]=remove_dispersion(meas,fs,k,a_cg*1e3,freq*1e3,340e3);
%     return
%     Y=wavelet_filter(Y,100e3,fs);
%      meas=wavelet_filter(meas,100e3,fs);
     figure(99)
    plot(X,Y/max(Y),'r','linewidth',2)
    hold on
    plot(X,meas/max(meas),'b--')
    hold off
    