clear all
close all


%% 1.2 Unsupervised artifact removal using CCA
%% 1.2.1
%%1.2.1.1 Visualize EEG
load eegdata_artifacts.mat
figure('Name','Raw EEG'),
eegplot_simple(eegdata,fs)

%%1.2.1.3 Apply CCA
[channel , time] = size(eegdata);
eegdata_delay    = [zeros(channel,1),eegdata(:,1:end-1)];
% eegdata_delay = delayseq(eegdata',ones(channel,1))';
[A,B,r,U,V]  = canoncorr(eegdata', eegdata_delay');
figure('Name','Estimated Sources'),
eegplot_simple(U',fs) % muscle artifacts at 39:end

%%1.2.1.4 Reconstruct EEG
%%A and B correspond to W_x and W_y, U and V correspond to x and y, U=XA
selected_channels = 38;
removed_channels  = channel-selected_channels;
U_clear = [U(:,1:selected_channels)' ; zeros(removed_channels,time)];
eegdata_clear = U_clear'*inv(A);
figure('Name','Reconstructed EEG with muscle artifacts removal'),
eegplot_simple(eegdata_clear',fs)



%% 1.3 Supervised artifact removal using MWF
%% 1.3.1
%%1.3.1.1 creat eye blink mask in fist 30 Sec, MWF process with delay=0
% [mask_eyeblink] = mwf_getmask(eegdata,fs);
% save('mask_eyeblink.mat','mask_eyeblink')
load mask_eyeblink.mat
[n, d, W, SER, ARR, p] = mwf_process(eegdata, mask_eyeblink,0, []);
figure('Name','MWF filtered EEG with eye blink mask(delay=0)'),
eegplot_simple(n,fs)

%%1.3.1.2
ch = 1;
T = 0:1/fs:30;
figure('Name','Eye blink artifact removal using MWF with delay=0'),
hold on
plot(T,eegdata(ch,1:numel(T)))
plot(T,n(ch,1:numel(T)))
plot(T,d(ch,1:numel(T)))
xlabel('Time (sec)')
title(sprintf('SER = %.2f , ARR= %.2f ',SER,ARR))
legend('Raw EEG','MWF filtered EEG (delay=0)','Estimated artifact')

%%1.3.1.3
%%use same mask as before, MWF process with delay=3
[n, d, W, SER, ARR, p] = mwf_process(eegdata, mask_eyeblink, 3,[]);
figure('Name','MWF filtered EEG with eye blink mask(delay=3)'),
eegplot_simple(n,fs)

ch = 1;
T = 0:1/fs:30;
figure('Name','Eye blink artifact removal using MWF with delay=3'),
hold on
plot(T,eegdata(ch,1:numel(T)))
plot(T,n(ch,1:numel(T)))
plot(T,d(ch,1:numel(T)))
xlabel('Time (sec)')
title(sprintf('SER = %.2f [dB], ARR= %.2f [dB]',SER,ARR))
legend('Raw EEG','MWF filtered EEG (delay=3)','Estimated artifact')

%%SER:Signal to Error Ratio, measures clean EEG distortion
%%ARR:Artifact to Residue Ratio, measures artifact estimation
%%For delay 3, the ARR is increased by 2.3 which means that we have better 
%%artifact estimation.But, SER is decreased by 30 which means that we have 
%%more distortion on clean EEG. 


%%1.3.1.5 retaining only the first few generalized eigenvalues 
delay = 0;
SER_all = zeros(1,64);
ARR_all = zeros(1,64);
for j=1:64
    p = mwf_params('rank', 'first', 'rankopt', j,...
                    'delay', delay);
    [n, d, W, SER, ARR, p] = mwf_process(eegdata, mask_eyeblink, delay,p);
    SER_all(j)= SER;
    ARR_all(j)=ARR;
end

[~,index]= max(ARR_all);
figure('Name','Find ideal number of retained eigenvalues in MWF'),
hold on
plot(SER_all)
plot(ARR_all)
xline(index,'m')
legend('SER','ARR',sprintf('ideal number = %d',index))
xlabel('Number of retained eigenvalues')
ylabel('[dB]')
title('Compare SER and ARR of different number of retained eigenvalues')


%%1.3.1.6 creat muscle mask, MWF process with delay=0
% [mask_muscle] = mwf_getmask(eegdata,fs);
% save('mask_muscle.mat','mask_muscle')
load mask_muscle.mat
[n, d, W, SER, ARR, p] = mwf_process(eegdata, mask_muscle,0, []);
figure('Name','MWF filtered EEG with muscle mask(delay=0)'),
eegplot_simple(n,fs)

ch = 6;
T = 60:1/fs:100;
figure('Name','Muscle artifact removal using MWF, delay=0'),
hold on
plot(T,eegdata(ch,T(1)*fs:T(end)*fs))
plot(T,n(ch,T(1)*fs:T(end)*fs))
plot(T,d(ch,T(1)*fs:T(end)*fs))
xlabel('Time (sec)')
title(sprintf('SER = %.2f [dB] , ARR= %.2f [dB] ',SER,ARR))
legend('Raw EEG','MWF filtered EEG (delay=0)','Estimated artifact')

%%muscle artifact removal using MWF process with delay=3
[n, d, W, SER, ARR, p] = mwf_process(eegdata, mask_muscle,3, []);
figure('Name','MWF filtered EEG with muscle mask(delay=3)'),
eegplot_simple(n,fs)

ch = 6;
T = 60:1/fs:100;
figure('Name','Muscle artifact removal using MWF, delay=3'),
hold on
plot(T,eegdata(ch,T(1)*fs:T(end)*fs))
plot(T,n(ch,T(1)*fs:T(end)*fs))
plot(T,d(ch,T(1)*fs:T(end)*fs))
xlabel('Time (sec)')
title(sprintf('SER = %.2f [dB] , ARR= %.2f [dB] ',SER,ARR))
legend('Raw EEG','MWF filtered EEG (delay=3)','Estimated artifact')
%%we have better performance(higher SER and ARR) by increasing delay muscle
%%artifact removal using MWF

%%1.3.1.7 Compute SER and ARR of the muscle artifact removal using CCA
d_CAA = eegdata - eegdata_clear';
[SER_CAA, ARR_CAA] = mwf_performance(eegdata, d_CAA, mask_muscle);
fprintf('SER CAA: %.2f , SER MWF: %.2f',SER_CAA,SER)
fprintf('\nARR CAA: %.2f , ARR MWF: %.2f \n',ARR_CAA,ARR)


%%1.3.1.8 Union of eye blink and muscle mask, 
%%create new mask for eye blink without 30 sec limitation 
% [mask_eyeblink_full] = mwf_getmask(eegdata,fs);
% save('mask_eyeblink_full.mat','mask_eyeblink_full')
load mask_eyeblink_full.mat
% mask_union = mask_eyeblink + mask_muscle;
mask_union = mask_eyeblink_full + mask_muscle;
mask_union(mask_union==2) = 1;

[n, d, W, SER, ARR, p] = mwf_process(eegdata, mask_union, 3,[]);
figure('Name','MWF filtered EEG with union mask(delay=3)'),
eegplot_simple(n,fs)

ch = 35;
T = 0:1/fs:(time/fs) - (1/fs);
figure('Name','Eye blink and muscle artifact removal using MWF with delay=3'),
hold on
plot(T,eegdata(ch,:))
plot(T,n(ch,:))
plot(T,d(ch,:))
xlabel('Time (sec)')
title(sprintf('SER = %.2f [dB], ARR= %.2f [dB]',SER,ARR))
legend('Raw EEG','MWF filtered EEG (delay=3)','Estimated artifact')


