function lagged_noise_covariance = noise_covariance(noiseSegments,trainingData,spike_window)
lag = spike_window;
cov=zeros(size(trainingData,2));
for i=1:size(noiseSegments,1)-lag
    segm = trainingData(noiseSegments(i):noiseSegments(i+lag),:);
    cov = cov + segm'*segm;
end

lagged_noise_covariance = cov./(size(noiseSegments,1)-lag);





