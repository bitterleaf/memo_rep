% Read fuzzy_cluster_wrap.R output text files as cells of matrices
%
% Created by Xi Jiang, Aug. 1st, 2014
%
% Last edited by Xi Jiang, Aug. 5th, 2014

function [cluster,membership] = fuzz_reader
    
% cluster: hard-cluster results for all time bins (N cells of n by 2
%   matrices, N being the number of clustering runs, n being the number of 
%   time bins)
% membership: percentage likelihood of the bins being in any given cluster
%   (N cells of n by m matrices, N being the number of clustering runs, n
%   being the number of time bins, and m being the number of clusters used
%   in a given run
% ALSO REQUIRED: have files named 'fuzz_cluster' and 'fuzz_membership'
%   under the working directory

%% Opening files and loading data

    C = textscan(fopen('fuzz_cluster'),'%q %f');
    
    C{1} = str2double(C{1});
    C = cell2mat(C);
    nan_ind = find(isnan(C(:,1)));  % index of NaNs that separate the fuzzy 
                                    % clustering runs. Length of this
                                    % vector is equal to the number of runs
                                    
    % obtain segment size for loading fuzz_membership in blocks
    if length(nan_ind) > 1
        step = nan_ind(2) - nan_ind(1);
    else
        step = length(C);
    end
    
    % load the irregular text file as rows of characters
    M = textread('fuzz_membership', '%s','delimiter', '\n');
    
%% Re-organizing data into cells of matrices
    
    % converting hard cluster results
    cluster = cell(1,length(nan_ind));  % pre-allocate space
                                    
    for i = 1:(length(nan_ind)-1)
        cluster{i} = C([nan_ind(i)+1:nan_ind(i+1)-1],:);
    end
    cluster{length(nan_ind)} = C([nan_ind(end)+1:end],:);
    
    membership = cell(1,length(nan_ind));   % pre-allocate space
    
    % figure out the number of clusters used 
    n_vec = zeros(1,length(nan_ind));
    for k = 1:length(nan_ind)
        n_cluster_line = M{nan_ind(k)};
        n_cluster_line = double(n_cluster_line);
        n_cluster = n_cluster_line(end-1);
        n_cluster = str2num(char(n_cluster));
        n_vec(k) = n_cluster;
    end
    
    % converting fuzzy (degree-of-belonging) results
    for i = 1:(length(nan_ind))
        if i ~= length(nan_ind)
            M_segment = [M{nan_ind(i)+1:nan_ind(i+1)-1}];
        else
            M_segment = [M{nan_ind(i)+1:end}];
        end
        M_segment = double(M_segment);
        
        % locating quotation marks and "remove" the numbers within
        % quotation marks, i.e. the numerical labels for observations
        quote_loc = find(M_segment == 34);
        for k = 1:length(quote_loc)
            if M_segment(quote_loc(k)+1) == 32
                M_segment(quote_loc(k-1):quote_loc(k)) = 48;  % char(48)=0
            end
        end
        
        space_loc = find(M_segment == 32);

        for j = 1:length(space_loc)
            if j ~= length(space_loc)
                M_slice = M_segment(space_loc(j)+1:space_loc(j+1)-1);
            else
                M_slice = M_segment(space_loc(j)+1:end);
            end
            M_slice = char(M_slice);
            M_slice = str2num(M_slice);
            membership{i}(j) = M_slice;
        end
        
        membership{i} = reshape(membership{i},n_vec(i),length(membership{i})/n_vec(i))';
    end
    
    
end