function [slope] = freq_drop_off(RR,ID)

fs = 2;
time = cumsum(RR/1000); % Time in seconds
time_new = 0:1/fs:time(end)-time(1);  
RR_resamp = spline(time-time(1),RR,time_new); 

n = length(RR_resamp);
RR_fft = fft(RR_resamp);
power = (abs(RR_fft).^2)/n;

freq = (0:n-1)*(fs/n);
indx = find(freq<=10^(-2) & freq>=10^(-4));

freq_log  = log(freq(indx));
power_log = log(power(indx));

p = polyfit(freq_log,power_log,1);
power_log_fit = polyval(p,freq_log);
slope = p(1);

figure,
hold on,
plot(freq_log,power_log,'b')
xlabel('log(freq)')
ylabel('log(power)')
plot(freq_log,power_log_fit,'r')
legend('power','regression line')
title(sprintf('slope = %.2f , ID = %s',slope,ID))




end