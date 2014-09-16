% chop the raw and preprocessed data into R-digestable chunks of cleaned 
% HGP peak outlier-containing time bins (default 2 seconds)
%
% Created by Xi Jiang, Sept. 5th, 2014
% Last edited by Xi Jiang, Sept.10th, 2014
%
% Dependencies: outlier_wavedecompC.m
%
% Outputs:
%  good_bins: logical indices of (allegedly) artifact-free time bins, with
%             respect to the elements in all_bins_seq
%  bad_bins: logical indices of (allegedly) artifact-containing time bins, 
%            with respect to the elements in all_bins_seq
%  data_LFP: # of channels-by-1 cell containing one matrix per channel, 
%            each row of which representing (by default) a 2-second time 
%            bin containing at least one HGP peak outlier that may serve as
%            the bin's center
%  data_HGP: similar to the above, except the data has been
%            bandpass-filtered, hilbert-transformed, turned into z-score,
%  new_bins: indices of time bin start/end points
%
% Do note that both data_LFP and data_HGP contains all HGP peak 
%  outlier-containing bins, regardless of artifact presence. To select only
%  the non-artifact bins, use the indices in good_bins.
%        
% Inputs:
%  data_out: one of the outputs of outlier_cleaning.m, with required field
%            "newPOI", representing numerical indices of all HGP peak 
%            outliers
%  data_raw: raw data in matrix form
%  methods: a structure that contains at least the following fields:
%           -time_bin: integer representing size of time bins, in seconds
%           -sfreq: sampling frequency (in Hz)
%           -wavthresh: threshold (in z-score) for artifact detection
%           -contthresh: threshold (in sampling points) for flat-lining artifact 
%                        (i.e. how long the signal needs to stay at the exact same 
%                        value before the time bin can be considered contaminated)
%           -peakthresh: threshold (in z-score) for HGP peak outliers that are
%                        likely to be artefactual
%           -tol: floating number difference tolerance, for removal of flatline 
%                        artifacts (consecutively occuring exact same values). 
%


function [good_bins,bad_bins,data_LFP,data_HGP,new_bins] = ...
            outlier_chopperC_nonincl(data_out,data_raw,methods)
        
    % preassign space for outputs
    good_bins = cell(size(data_out.newPOI,1),1);
    bad_bins = cell(size(data_out.newPOI,1),1);
    data_LFP = cell(size(data_out.newPOI,1),1);
    data_HGP = cell(size(data_out.newPOI,1),1);
    new_bins.start = cell(size(data_out.newPOI,1),1);
    new_bins.end = cell(size(data_out.newPOI,1),1);
    
    if ~isfield(methods,'time_bin')
        methods.time_bin = 2;
    end

    for i = 1:size(data_out.newPOI,1)

    % find needed indices for the given channel and pre-assign space for
    % time bins   
        ind_all = data_out.newPOI{i};
        
        num_bins = length(ind_all);       % number of bins
        bins_LFP = cell(1,num_bins);
        bins_HGP = cell(1,num_bins);
        step = methods.time_bin * methods.sfreq;  % time bin size in sampling points
        
        new_bins.start{i} = [];
        new_bins.end{i} = [];
    
    % obtain LFP and convert data to matrix form
        for j = 1:num_bins
            if ind_all(j) <= methods.sfreq-1
                bins_LFP{j} = data_raw(i,1:step);
                bins_HGP{j} = data_out.smoothed(i,1:step);
                new_bins.start{i} = horzcat(new_bins.start{i},1);
                new_bins.end{i} = horzcat(new_bins.end{i},step);
            elseif ind_all(j) > length(data_out.smoothed)-methods.sfreq
                bins_LFP{j} = data_raw(i,...
                    length(data_raw)-(step-1):length(data_raw));
                bins_HGP{j} = data_out.smoothed(i,...
                    length(data_out.smoothed)-(step-1):length(data_out.smoothed));
                new_bins.start{i} = horzcat(new_bins.start{i},length(data_raw)-(step-1));
                new_bins.end{i} = horzcat(new_bins.end{i},length(data_raw));
            else
                bins_LFP{j} = data_raw(i,...
                    (ind_all(j)-(round(step/2)-1)):ind_all(j)+round(step/2));
                bins_HGP{j} = data_out.smoothed(i,...
                    (ind_all(j)-(round(step/2)-1)):ind_all(j)+round(step/2));
                new_bins.start{i} = horzcat(new_bins.start{i},ind_all(j)-(round(step/2)-1));
                new_bins.end{i} = horzcat(new_bins.end{i},ind_all(j)+round(step/2));
            end
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
        
    % make sure centering does not re-introduce bad bins
    
        % prepare for removing machine artifact
        LFP_elevation_thresh = median(data_raw(i,:)) + 5*std(data_raw(i,:));
        
        for j = 1:size(all_bins_LFP,1)
            if sum(abs(diff(data_LFP{i}(j,:))) < methods.tol) > methods.contthresh
                % remove flatline artifacts
                good_bins{i}(good_bins{i} == j) = 0;
                good_bins{i}(good_bins{i} == 0) = [];
                bad_bins{i} = unique(sort([bad_bins{i} j]));
                
            elseif max(data_HGP{i}(j,:)) > methods.peakthresh
                % remove "too high" HGP peaks
                good_bins{i}(good_bins{i} == j) = 0;
                good_bins{i}(good_bins{i} == 0) = [];
                bad_bins{i} = unique(sort([bad_bins{i} j]));
                
            elseif any(data_raw(i,(j-1)*time_bin*sfreq+1:time_bin*sfreq*j) > LFP_elevation_thresh)
                % remove extreme LFP elevations (electronic in origin? A common problem with NY394)
                good_bins{i}(good_bins{i} == j) = 0;
                good_bins{i}(good_bins{i} == 0) = [];
                bad_bins{i} = unique(sort([bad_bins{i} j]));
                
            end
        end
        
        disp([num2str(i) ' channel(s) processed!'])
    
    end
    
end