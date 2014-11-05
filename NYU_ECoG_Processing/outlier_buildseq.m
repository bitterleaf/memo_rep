% Create binarized and state sequence-transformed HGP peak outlier matrices
%
% Created by Xi Jiang, Oct. 1st, 2014
%
% Dependencies: N/A
%
% Outputs:
%  hgp_bins_orig:    binarized channel-by-time_bin matrix
%  hgp_bins_numeric: hgp_bins_orig converted to TraMineR-ready state
%                    sequences, with zero-only time bins removed
%  ind:              index of time bins in hgp_bins_numeric, denoting their 
%                    positions within the 12-hour period covered by the 
%                    .edf file
%        
% Inputs:
%  good_peaks:  one of the outputs of outlier_chopperAll.m, containing 
%               HGP peak outliers that passed the artifact rejection stage
%  data_raw:    raw data matrix (channel-by-sampling point)
%  methods:     a structure that contains at least the following fields:
%                -time_bin: integer representing size of time bins, in seconds
%                -sfreq: sampling frequency (in Hz)
%                -synthresh: threshold for too synchronous events (percentage of
%                       "coactive" channels), can take on any value in [0,1]

function [hgp_bins_orig,hgp_bins_numeric,ind] = ... 
                        outlier_buildseq(good_peaks,data_raw,methods)

    %outlier_chopperAll
    step = methods.sfreq*methods.time_bin;
    num_bins = length(data_raw)/step;
    hgp_bins = zeros(size(data_raw,1),num_bins);
    for i = 1:size(data_raw,1)
        temp = good_peaks{i};
        for j = 1:num_bins
            if ~isempty(intersect(temp,step*(j-1)+1:step*j));
                hgp_bins(i,j) = 1;
            end
        end
        disp(['Channel ' num2str(i) ' processed!'])
    end

    % remove too synchronous events (machine blip?)
    for j = 1:num_bins
        if sum(hgp_bins(:,j)) > round(size(data_raw,1)*methods.synthresh);
            hgp_bins(:,j) = 0;
        end
    end
    hgp_bins_orig = hgp_bins;

    ind = 1:num_bins;
    ind = ind(any(hgp_bins));

    % make it so that traminer considers subsequences properly
    for i = 1:size(data_raw,1)
        for j = 1:size(hgp_bins,2)
            if hgp_bins(i,j) ~= 0
                hgp_bins(i,j) = i;
            end
        end
        disp(['Channel ' num2str(i) ' processed!'])
    end
    hgp_bins_numeric = hgp_bins(:,any(hgp_bins));
    clear hgp_bins
    
end