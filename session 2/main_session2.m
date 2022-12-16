close all
clear all

%% 2.1
load ex1data.mat
T1_rankest = mlrankest(T1);
[U,S,sv] = mlsvd(T1); % 6 9 7 
% lmlra_hooi
% lmlra_nls

%% 2.2
A = 1/sqrt(5).*[1,-2;-2,1];
SNR = 0:5:50;
iter = 200;
SIR_ICA=[];
SIR_PCA=[];

for i=SNR
    sir_ica_tot=0;
    sir_pca_tot=0;
    for j = 1:iter
        s1 = rand(1,800);
        s2 = rand(1,800);
        S = [s1;s2];
        x = A*S;
        [x,n] = noisy(x,i);

        % ICA
        [F,delta] = aci(x); 
        [sir_ica,P,D] = sir(A,F); %Signal-to-interference ratio
        sir_ica_tot = sir_ica_tot + sir_ica;

        % PCA
        x_corr = cov(x');
        Q = eye(2);
        [V,D] = eig(x_corr,Q);
        [sir_pca,P,D] = sir(A,V');
        sir_pca_tot = sir_pca_tot + sir_pca;
    end
    SIR_ICA = [SIR_ICA sir_ica_tot/iter];
    SIR_PCA = [SIR_PCA sir_pca_tot/iter];
end

figure
hold on
plot(SNR, SIR_ICA)
plot(SNR, SIR_PCA)
xlabel('SNR')
ylabel('SIR')
legend('SIR (ICA)','SIR (PCA)')

%% 2.3
load ex3data.mat

% Plotting
figure,
subplot(3,1,1)
hold on,
plot(A(:,1))
plot(A(:,2))
plot(A(:,3))
plot(A(:,4))
title('Components of matrix A')
legend('Component1','Component2','Component3','Component4')

subplot(3,1,2)
hold on,
plot(B(:,1))
plot(B(:,2))
plot(B(:,3))
plot(B(:,4))
title('Components of matrix B')
legend('Component1','Component2','Component3','Component4')

subplot(3,1,3)
hold on,
plot(C(:,1))
plot(C(:,2))
plot(C(:,3))
plot(C(:,4))
title('Components of matrix C')
legend('Component1','Component2','Component3','Component4')

%% 2.3.1
T = squeeze(T3(1,:,:));
[U,S,V] = svd(T);

projC = (U(:,1:4)'*T)';
projB = T*V(:,1:4);

% Plotting, SVD components
figure,
subplot(2,1,1)
hold on,
plot(U(:,1))
plot(U(:,2))
plot(U(:,3))
plot(U(:,4))
title('Components of matrix B SVD')
legend('Component1','Component2','Component3','Component4')

subplot(2,1,2)
hold on,
plot(V(:,1))
plot(V(:,2))
plot(V(:,3))
plot(V(:,4))
title('Components of matrix C SVD')
legend('Component1','Component2','Component3','Component4')


% Plotting, projected C and B
figure,
subplot(2,1,1)
hold on,
plot(projB(:,1))
plot(projB(:,2))
plot(projB(:,3))
plot(projB(:,4))
title('Components of matrix B PCA')
legend('Component1','Component2','Component3','Component4')

subplot(2,1,2)
hold on,
plot(projC(:,1))
plot(projC(:,2))
plot(projC(:,3))
plot(projC(:,4))
title('Components of matrix C PCA')
legend('Component1','Component2','Component3','Component4')

%% 2.3.2

T3_noise = noisy(T3,25);
[U_noise, S_noise, sv_noise] = mlsvd(T3_noise);

for R = 1:5 %rank-1 term
    
    % NLS
    options = struct;
    options.Compression = @mlsvd_rsi;
    options.Initialization = @cpd_gevd;
    options.Algorithm = @cpd_nls;
    options.Refinement = @cpd_nls;
    options.ExploitStructure = true;
    options.AlgorithmOptions.MaxIter = 100;      % Default 200
    options.AlgorithmOptions.TolFun = eps^2;     % Default 1e-12
    options.AlgorithmOptions.TolX = eps;         % Default 1e-6
    [U_NLS,output_NLS] = cpd(T3_noise,R,options);
    
    % ALS
    options.Algorithm = @cpd_als;
    options.Refinement = @cpd_als;
    options.ExploitStructure = true;
    options.AlgorithmOptions.MaxIter = 100;      % Default 200
    options.AlgorithmOptions.TolFun = eps^2;     % Default 1e-12
    options.AlgorithmOptions.TolX = eps;         % Default 1e-6
    [U_ALS,output_ALS] = cpd(T3_noise,R,options);
    
    
    figure,
    hold on
    semilogy(output_NLS.Algorithm.fval)
    semilogy(output_ALS.Algorithm.fval);
    ylabel('Objective function')
    xlabel('Iteration');
    title(sprintf('Convergence plot, R = %d',R))
    legend('cpd (NLS)','cpd (ALS)')
end

%% 2.4
load ex4data.mat

%% 2.5
fs = 256;
load('demosignal3_963.mat');
eeg = demosignal3_963;
eegplot_simple(eeg,fs);

% Normalise and wavelet-transform data near t=52
[data_3D,m,s] = normalise_3D(eeg,51,53);

% Decompose the tensor in two rank one terms.
R = 2;
U = cpd(data_3D,R);
A = U{1}; 
B = U{2}; 
C = U{3};

result = transform_back(A,s);

% Frequency signatures.
figure;
subplot(2,1,1)
plot(B);
title('Frequency signatures');
xlabel('Frequency')
ylabel('Normalized amplitude')
% Temporal signatures.
subplot(2,1,2)
plot(C);
title('Temporal signatures');
xlabel('Time index');
ylabel('Normalized amplitude')

