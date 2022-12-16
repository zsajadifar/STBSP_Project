close all
clear all

low_anx_ID = [507, 514, 571, 619, 621];
high_anx_ID = [518, 542, 547, 562, 576];

for i=1:length(low_anx_ID)
    filename = "RR" + num2str(low_anx_ID(i)) + ".mat";
    load(filename)
    ID = num2str(low_anx_ID(i));

    [FD,SE,LE,CD]=nonlin_measures(RR);
    [alpha1 , alpha2]= DFA(RR,ID);
    slope = freq_drop_off(RR,ID);
    [SD1 , SD2] = poincare_plot(RR,ID);
    
    FD_all_low_anx(i)=FD;
    SE_all_low_anx(i)=SE;
    LE_all_low_anx(i)=LE;
    CD_all_low_anx(i)=CD;
    alpha1_all_low_anx(i)=alpha1;
    alpha2_all_low_anx(i)=alpha2;
    freq_dropoff_all_low_anx(i)=slope;
    SD1_all_low_anx(i) = SD1;
    SD2_all_low_anx(i) = SD2;
end

for i=1:length(high_anx_ID)
    filename = "RR" + num2str(high_anx_ID(i)) + ".mat";
    load(filename)
    ID = num2str(high_anx_ID(i));

    [FD,SE,LE,CD]=nonlin_measures(RR);
    [alpha1 , alpha2]= DFA(RR,ID);
    slope = freq_drop_off(RR,ID);
    [SD1 , SD2] = poincare_plot(RR,ID);

    FD_all_high_anx(i)=FD;
    SE_all_high_anx(i)=SE;
    LE_all_high_anx(i)=LE;
    CD_all_high_anx(i)=CD;
    alpha1_all_high_anx(i)=alpha1;
    alpha2_all_high_anx(i)=alpha2;
    freq_dropoff_all_high_anx(i)=slope;
    SD1_all_high_anx(i) = SD1;
    SD2_all_high_anx(i) = SD2;

end

% SD
figure
subplot(1,2,1)
bar([SD1_all_high_anx;SD1_all_low_anx])
xlabel('1 = High Anxiety  |   2 = Low Anxiety ')
title('SD1')

subplot(1,2,2)
bar([SD2_all_high_anx;SD2_all_low_anx])
xlabel('1 = High Anxiety  |   2 = Low Anxiety ')
title('SD2')

% alpha
figure
subplot(1,2,1)
bar([alpha1_all_high_anx;alpha1_all_low_anx])
xlabel('1 = High Anxiety  |   2 = Low Anxiety ')
title('alpha1')

subplot(1,2,2)
bar([alpha2_all_high_anx;alpha2_all_low_anx])
xlabel('1 = High Anxiety  |   2 = Low Anxiety ')
title('alpha2')

% drop-off
figure,
bar([freq_dropoff_all_high_anx;freq_dropoff_all_low_anx])
xlabel('1 = High Anxiety  |   2 = Low Anxiety ')
title('frequency drop-off')

% FD
figure,
bar([FD_all_high_anx;FD_all_low_anx])
xlabel('1 = High Anxiety  |   2 = Low Anxiety ')
title('FD')

% SE
figure,
bar([SE_all_high_anx;SE_all_low_anx])
xlabel('1 = High Anxiety  |   2 = Low Anxiety ')
title('SE')

% LE
figure,
bar([LE_all_high_anx;LE_all_low_anx])
xlabel('1 = High Anxiety  |   2 = Low Anxiety ')
title('LE')

% CD
figure,
bar([CD_all_high_anx;CD_all_low_anx])
xlabel('1 = High Anxiety  |   2 = Low Anxiety ')
title('CD')
