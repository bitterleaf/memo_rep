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


function [data_LFP,data_HGP] = ...
           outlier_chopperC(data_out,data_raw,new_starts,new_ends,methods)
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
    
        disp([num2str(i) ' channel(s) processed!'])
    end

end