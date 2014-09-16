% chop the raw and preprocessed data into R-digestable chunks of cleaned 
% HGP peak outlier-containing time bins (default 2 seconds)
%
% Created by Xi Jiang, Sept. 9th, 2014
%
% Dependencies: outlier_wavedecompC.m
%
% Outputs:
%  data_LFP: # of channels-by-1 cell containing one matrix per channel, 
%            each row of which representing (by default) a 2-second time 
%            bin containing at least one HGP peak outlier that may serve as
%            the bin's center
%  data_HGP: similar to the above, except the data has been
%            bandpass-filtered, hilbert-transformed, turned into z-score,
%  good_peaks: indices of HGP peak outliers in each channel that happen to
%              fall into one of the new time bins (cleaned with extreme
%              prejudice)
%  new_starts: one of the outputs from outlier_unifier.m, denoting the
%              starting indices of clean time bins (cleaned)
%  new_ends: similar to the variable above, except containing ends of time
%            bins (cleaned)
%        
% Inputs:
%  data_out: one of the outputs of outlier_cleaning.m, with required field
%            "newPOI", representing numerical indices of all HGP peak 
%            outliers
%  data_raw: raw data in matrix form
%  new_starts: one of the outputs from outlier_unifier.m, denoting the
%              starting indices of clean time bins
%  new_ends: similar to the variable above, except containing ends of time
%            bins
%  methods: a structure that contains at least the following fields:
%           -time_bin: integer representing size of time bins, in seconds
%           -sfreq: sampling frequency (in Hz)
%           -contthresh: threshold (in sampling points) for flat-lining artifact 
%                        (i.e. how long the signal needs to stay at the exact same 
%                        value before the time bin can be considered contaminated)
%           -peakthresh: threshold (in z-score) for HGP peak outliers that are
%                        likely to be artefactual
%           -tol: floating number difference tolerance, for removal of flatline 
%                        artifacts (consecutively occuring exact same values). 
%           Optionally:
%           -savedir: string path of directory for saving outputs AND 
%                     filename, e.g. '/home/test/output.mat'
%  good_peaks: indices of HGP peak outliers in each channel that happen to
%              fall into one of the new time bins 


function [data_LFP,data_HGP,new_starts,new_ends,good_peaks,bad_bins] = ...
           outlier_chopperC(data_out,data_raw,new_starts,new_ends,good_peaks,methods)
%% prepare reusable variables, allocate spaces
    
    step = methods.sfreq * methods.time_bin;

    chan_num = size(data_raw,1);
    bin_num = length(new_starts);
    
    data_LFP = cell(chan_num,1);
    data_HGP = cell(chan_num,1);
    
    
%% assign data chunks to bins

    for i = 1:chan_num
        all_bins_LFP = cell(1,bin_num);
        all_bins_HGP = cell(1,bin_num);
        
        for j = 1:bin_num
            all_bins_LFP{j} = data_raw(i,new_starts(j):new_ends(j));
            all_bins_HGP{j} = data_out.smoothed(i,new_starts(j):new_ends(j));
        end
        
        all_bins_LFP = cell2mat(all_bins_LFP);
        all_bins_HGP = cell2mat(all_bins_HGP);
        
            % new matrix: time bins by sampling points (note the transpose)
        all_bins_LFP = reshape(all_bins_LFP,step,bin_num)';
        all_bins_HGP = reshape(all_bins_HGP,step,bin_num)';
        
        data_LFP{i} = all_bins_LFP;
        data_HGP{i} = all_bins_HGP;
        
        
    % make sure centering does not re-introduce bad bins
    
        % prepare for removing machine artifact
        LFP_elevation_thresh = median(data_raw(i,:)) + 5*std(data_raw(i,:));
        bad_bins = [];
        
        for j = 1:size(all_bins_LFP,1)
            if sum(abs(diff(data_LFP{i}(j,:))) < methods.tol) > methods.contthresh
                % remove flatline artifacts
                good_peaks{i}(good_peaks{i} >= new_starts(j) && good_peaks{i} <= new_ends(j)) = 0;
                bad_bins = sort([bad_bins j]);
                
            elseif max(data_HGP{i}(j,:)) > methods.peakthresh
                % remove "too high" HGP peaks
                good_peaks{i}(good_peaks{i} >= new_starts(j) && good_peaks{i} <= new_ends(j)) = 0;
                bad_bins = sort([bad_bins j]);
                
            elseif any(data_raw(i,(j-1)*time_bin*sfreq+1:time_bin*sfreq*j) > LFP_elevation_thresh)
                % remove extreme LFP elevations (electronic in origin? A common problem with NY394)
                good_peaks{i}(good_peaks{i} >= new_starts(j) && good_peaks{i} <= new_ends(j)) = 0;
                bad_bins = sort([bad_bins j]);
                
            end
            good_peaks{i}(good_peaks{i} == 0) = [];
            data_LFP{i}(bad_bins,:) = [];
            data_HGP{i}(bad_bins,:) = [];
        end
        
        disp([num2str(i) ' channel(s) processed!'])
    end
    
    new_starts(bad_bins) = [];
    new_ends(bad_bins) = [];
    
    if isfield(methods,'savedir')
        save(methods.savedir,'data_LFP','data_HGP',...
                            'new_starts','new_ends','good_peaks')
    end

end