function [ output ] = applyMultiChannelFilter( data, filter )
%APPLYMULTICHANNELFILTER Calculate multi-channel FIR filter output on given data.
%
% inputs:
%   data : multi-channel data (numberSamples, numberChannels)
%   filter : multi-channel FIR filter (numberTaps, numberChannels)
%
% output:
%   output : filter output
%

% extract the temporal window length from the give matrix filter
% coefficients
L = size(filter, 1);
% initialize the filter output
output = zeros(size(data, 1), 1);

% calculate output
for idx=L:length(output)
    
    % TODO: Calculate the filter output(idx). At time idx consider only idx
    % and its past samples (i.e., we implement a causal filter).

   output(idx) = sum(sum(data(idx-L+1:idx,:).*filter));
end

end

