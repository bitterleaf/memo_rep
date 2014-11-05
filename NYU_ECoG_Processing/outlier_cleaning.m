% Obtain outliers and remove artifact periods to create ROI
%
% Created by Xi Jiang, Aug. 22nd, 2014
% Last edited by Xi Jiang, Nov. 5th, 2014
%
% Dependency: N/A
%
% Outputs:
%   newPeakInd:  by-channel indices (in sampling points) of HGP peak 
%                outliers that do not fall into one of the 
%                artifact-containing time bins
%   binaryInd:   by-channel binary indices of whether a given time point 
%                corresponds to a HGP peak outlier (1 for yes)
%   all_bins:    by-channel segmented HGP envelope, with artifact-containing
%                time bins removed
%   all_bin_seq: by-channel indices of the time bins in all_bins, i.e. 
%                markers for where they are in the original data
%   data_out: output structure of outlier_finder
%
% Inputs:
%   data_out: output structure of outlier_finder, containing at least the
%             following field:
%             -smoothed: HGP envelope data matrix (channel-by-sampling points)
%   data_raw: raw data matrix (channel-by-sampling points)
%   methods:  structure (potentially) containing the following fields:
%              -time_sep:   minimum separation between HGP peaks (in seconds). 
%                           Default: 0.1s
%              -sfreq:      sampling frequency in Hz. Default: 512Hz
%              -time_bin:   size of time bins (in seconds) into which HGP 
%                           peaks are assigned. Default: 10s
%              -peakthresh: threshold (in z-score) for HGP peak outliers 
%                           that are likely to be artefactual
%   bad_bins: time bins marked by hand as containing artifacts, epileptic or
%             otherwise. If no such index is available, set it to 0.
%

function [newPeakInd,binaryInd,all_bins,all_bin_seq,data_out] = ...
outlier_cleaning(data_out,data_raw,methods,bad_bins)

    if nargin < 4
        bad_bins = 0;
    end

    if ~isfield(methods,'time_sep')
        methods.time_sep = 0.1;
    end 
    if ~isfield(methods,'sfreq')
        sfreq = 512;
    else
        sfreq = methods.sfreq;
    end 
    
    if ~isfield(methods,'time_bin')
        time_bin = 10;
    else
        time_bin = methods.time_bin;
    end
    if ~isfield(methods,'tol')
        methods.tol = 1e-12;
    end
    
    % preallocate space for output
    newPeakInd = cell(size(data_out.smoothed,1),1);
    binaryInd = zeros(size(data_out.smoothed));
    all_bins = cell(size(data_out.smoothed,1),1);
    all_bin_seq = cell(size(data_out.smoothed,1),1);
    
    for i = 1:size(data_out.smoothed,1)
        % locate and binarize HGP peak index   
        Ind = sort(data_out.POI{i});
        inter_peak_distance = [0 diff(Ind)];
        Ind = Ind(inter_peak_distance > methods.time_sep*sfreq);
    
        binary_ind = zeros(1,length(data_out.smoothed));
        binary_ind(Ind) = 1;

    
        segment_size = length(data_out.smoothed)/(time_bin*sfreq);
        bins = cell(1,segment_size);
        LFP_elevation_thresh = median(data_raw(i,:)) + 5*std(data_raw(i,:));
        
        % find time bins that have non-artifact HGP peak(s)
        for j = 1:segment_size
            bins{j} = data_out.smoothed(i,(j-1)*time_bin*sfreq+1:time_bin*sfreq*j);
            if ~any(binary_ind((j-1)*time_bin*sfreq+1:time_bin*sfreq*j) == 1)
                bins{j} = {};
            elseif any(bad_bins == j)
                % remove hand-marked artifacts
                bins{j} = {};
                binary_ind((j-1)*time_bin*sfreq+1:time_bin*sfreq*j) = 0;
                [~,IA,~] = intersect(Ind,((j-1)*time_bin*sfreq+1:time_bin*sfreq*j));
                Ind(IA) = 0;
                Ind(Ind == 0) = [];
%             elseif sum(abs(diff(data_raw(i,(j-1)*time_bin*sfreq+1:time_bin*sfreq*j))) < methods.tol) > methods.contthresh
%                 % remove flatline artifacts: DEFUNCT, DO NOT TOUCH
%                 bins{j} = {};
%                 binary_ind((j-1)*time_bin*sfreq+1:time_bin*sfreq*j) = 0;
%                 [~,IA,~] = intersect(Ind,((j-1)*time_bin*sfreq+1:time_bin*sfreq*j));
%                 Ind(IA) = 0;
%                 Ind(Ind == 0) = [];
            elseif max(bins{j}) > methods.peakthresh
                % remove "too high" HGP peaks
                bins{j} = {};
                binary_ind((j-1)*time_bin*sfreq+1:time_bin*sfreq*j) = 0;
                [~,IA,~] = intersect(Ind,((j-1)*time_bin*sfreq+1:time_bin*sfreq*j));
                Ind(IA) = 0;
                Ind(Ind == 0) = [];
            elseif any(data_raw(i,(j-1)*time_bin*sfreq+1:time_bin*sfreq*j) > LFP_elevation_thresh)
                % remove extreme LFP elevations (electronic in origin? A common problem with NY394)
                bins{j} = {};
                binary_ind((j-1)*time_bin*sfreq+1:time_bin*sfreq*j) = 0;
                [~,IA,~] = intersect(Ind,((j-1)*time_bin*sfreq+1:time_bin*sfreq*j));
                Ind(IA) = 0;
                Ind(Ind == 0) = [];
            end
        end
    
    
        % obtain index of non-empty bins and remove empty bins
        full_bin_ind = ~cellfun('isempty',bins);
        full_bin_seq = 1:segment_size;
        full_bin_seq = full_bin_seq(full_bin_ind);
        bins(cellfun('isempty',bins)) = [];
        
        newPeakInd{i} = Ind;
        binaryInd(i,:) = binary_ind;
        all_bins{i} = bins;
        all_bin_seq{i} = full_bin_seq;
        disp([num2str(i) ' channel(s) processed!'])
    end
    data_out.newPOI = newPeakInd;
    
end