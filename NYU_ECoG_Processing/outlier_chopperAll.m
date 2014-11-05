% chop the raw and preprocessed data into R-digestable chunks of cleaned 
% HGP peak outlier-containing time bins (default 2 seconds). While this
% function is designed to create inputs for the dynamic time warping
% algorithm in R, it also screens the previously obtained HGP peak outlier 
% indices for artifacts in the process.
%
% Created by Xi Jiang, Sept. 23th, 2014
% Last edited by Xi Jiang, Nov. 5th, 2014
%
% Dependencies: outlier_wavedecompC.m
%
% Outputs:
%  good_bins:  logical indices of (allegedly) artifact-free time bins, with
%              respect to the elements in all_bins_seq
%  bad_bins:   logical indices of (allegedly) artifact-containing time bins, 
%              with respect to the elements in all_bins_seq
%  data_LFP:   61-by-1 cell containing 61 matrices, each row of which 
%              represents (by default) a 2-second time bin containing at 
%              least one HGP peak outlier that may serve as the bin's center
%  data_HGP:   similar to the above, except the data has been
%              bandpass-filtered, hilbert-transformed, turned into z-score,
%              and smoothed with a 0.1s Gaussian kernel.
%  good_peaks: indices of HGP peak outliers in each channel that happen to
%              fall into one of the time bins
% Do note that both data_LFP and data_HGP contains all HGP peak 
%  outlier-containing bins, regardless of artifact presence. To select only
%  the non-artifact bins, use the indices in good_bins.
%        
% Inputs:
%  data_out: one of the outputs of outlier_cleaning.m
%  data_raw: raw data in matrix form
%  methods: a structure that contains at least the following fields:
%           -time_bin:  integer representing size of time bins, in seconds
%           -sfreq:     sampling frequency (in Hz)
%           -wavthresh: threshold (in z-score) for artifact detection
%           -synthresh: threshold for too synchronous events (percentage of
%                       "coactive" channels), can take on any value in [0,1]



function [good_bins,bad_bins,data_LFP,data_HGP,good_peaks] = ...
            outlier_chopperAll(data_out,data_raw,methods)


    good_peaks = data_out.newPOI;    
        
    % preassign space for outputs
    good_bins = cell(size(data_out.newPOI,1),1);
    bad_bins = cell(size(data_out.newPOI,1),1);
    data_LFP = cell(size(data_out.newPOI,1),1);
    data_HGP = cell(size(data_out.newPOI,1),1);
    
    if ~isfield(methods,'time_bin')
        methods.time_bin = 2;
    end

    for i = 1:size(data_out.newPOI,1)

    % find needed indices for the given channel and pre-assign space for
    % time bins   
    
        step = methods.time_bin * methods.sfreq;  % time bin size in sampling points
        num_bins = length(data_raw)/step;
        bins_LFP = cell(1,num_bins);
        bins_HGP = cell(1,num_bins);
        step = methods.time_bin * methods.sfreq;  % time bin size in sampling points
    
    % obtain LFP and convert data to matrix form
    
        for j = 1:num_bins
            bins_LFP{j} = data_raw(i,step*(j-1)+1:step*j);
            bins_HGP{j} = data_out.smoothed(i,step*(j-1)+1:step*j);
        end

        all_bins_LFP = cell2mat(bins_LFP);
        all_bins_HGP = cell2mat(bins_HGP);
    % new matrix: time bins by sampling points (note the transpose)
        all_bins_LFP = reshape(all_bins_LFP,step,num_bins)';
        all_bins_HGP = reshape(all_bins_HGP,step,num_bins)';
    
    % reject artifacts based on wavelet decomposition
        [good_bin_ind,bad_bin_ind] = ...
            outlier_wavedecompC(all_bins_LFP,methods);
        
        good_bins{i} = good_bin_ind;
        bad_bins{i} = bad_bin_ind;
        data_LFP{i} = all_bins_LFP;
        data_HGP{i} = all_bins_HGP;
        for j = 1:size(all_bins_LFP,1)
            if good_bin_ind(j) ~= 1
                good_peaks{i}(good_peaks{i} >= step*(j-1)+1 & good_peaks{i} <= step*j) = 0;
            end
        end
        
    % make sure centering does not re-introduce bad bins
    
        % prepare for removing machine artifact
        LFP_elevation_thresh = median(data_raw(i,:)) + 5*std(data_raw(i,:));
        
        for j = 1:size(all_bins_LFP,1)
%             if sum(abs(diff(data_LFP{i}(j,:))) < methods.tol) > methods.contthresh
%                 % remove flatline artifacts: BAD! DO NOT TOUCH!
%                 good_peaks{i}(good_peaks{i} >= step*(j-1)+1 & good_peaks{i} <= step*j) = 0;
%                 good_bins{i}(j) = 0;
%                 bad_bins{i}(j) = 1;
                
            if max(data_HGP{i}(j,:)) > methods.peakthresh
                % remove "too high" HGP peaks
                good_peaks{i}(good_peaks{i} >= step*(j-1)+1 & good_peaks{i} <= step*j) = 0;
                good_bins{i}(j) = 0;
                bad_bins{i}(j) = 1;
                
            elseif any(data_LFP{i}(j,:) > LFP_elevation_thresh)
                % remove extreme LFP elevations (electronic in origin? A common problem with NY394)
                good_peaks{i}(good_peaks{i} >= step*(j-1)+1 & good_peaks{i} <= step*j) = 0;
                good_bins{i}(j) = 0;
                bad_bins{i}(j) = 1;
                
            end
        end
        
        good_peaks{i}(good_peaks{i} == 0) = [];
      
        
        disp([num2str(i) ' channel(s) processed!'])
    
    end
    
    % reject too synchronous time bins
    good_bins = cell2mat(good_bins);
    for i = 1:length(good_bins)
        if sum(good_bins(:,i)) > round(size(data_raw,1)*methods.synthresh);
            good_bins(:,i) = 0;
        end
    end
    
end