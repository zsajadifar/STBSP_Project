function [SD1 , SD2] = poincare_plot(RR,ID)
xp=RR;
xp(end)=[];
xm=RR;
xm(1)=[];

%SD1
SD1 = std(xp-xm)/sqrt(2);

%SD2
SD2 = std(xp+xm)/sqrt(2);

figure,
plot(xm,xp,'.')
title(sprintf('poincare plot , SD1 = %.2f , SD2 = %.2f , ID = %s',SD1,SD2,ID))

end