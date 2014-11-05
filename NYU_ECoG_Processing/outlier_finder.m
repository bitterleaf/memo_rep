% Locates high gamma power peaks
% Fits inverse Gaussian distributions over each channel
% Locates "outlying" populations
% Thresholds those populations according to given rules
% Binarizes the binned data 
%
% Created by Xi Jiang, Aug.5th, 2014
% Last edited by Xi Jiang, Aug 20th, 2014
%
% Dependencies: ft_preproc_bandpassfilter.m
%               extrema.m
%
% Output:
%   data_out: output structure with the following fields:
%       -POI: indices of "outlying" HGP-hilbert-Gaussian-smoothed peaks
%       -hgpZscore/hgp: hilbert envelope of high gamma power
%       -smoothed: hgpZscore/hgp smoothed with Gaussian kernel
%       -allpeakInd: all local maxima in data_out.smoothed (UNSORTED!)
%       -distributions: ProbDist objects from each channel's distribution fit
%       -outlyingthresh: outlier thresholds for each channel based on 
%                       the fitted distributions' interquartile ranges
%
% Inputs:
%   data: pre-processed continuous iEEG data (LFP) in fieldtrip format
%   methods: structure of parameters that determine how HGP is defined and
%            how thresholds are set for HGP peaks, with the following fields:
%       -freq_range: range (in Hz) of HG frequency band, e.g. [70 190]
%       -filter_order: order of Butterworth filter (needs to be even integer)
%       -dist: type of distribution supported by MATLAB's "fithist", 
%              e.g. 'gamma', 'inversegaussian' 
%       -zscore: numerical "0" or "1". To convert pre-hgp data to z-score, 
%                use "1". Set to "0" otherwise.
%       -olc: outlier criterion, can be 'IQR' or 'ESD' ("three-sigma"
%           rule). default is 'IQR'.
%       -smooth: smoothing dur. of Gaussian kernel in seconds. Default: 0.1
%   params: structure of data labeling-related parameters 
%       -savedir: directory to which the output file is saved
%       -ID: subject ID
%       -day: dayfile number
%       -period: morning or evening ('morn' or 'eve')
%       -clinsys: clinical system from which the data is obtained, can be
%                 'clin1' or 'clin2'

function [data_out] = outlier_finder(data,methods,params)
%% Xiaojing's linenoise selection
    linenoise_freq = [60 120 180]; % Hz, for dft filter
    linenoise_bsfreq = [58 62; 118 122; 178 182]; % Hz, has to be an Nx2 matrix!!
    
    % Get rid of line noise using DFT filter
    cfg = [];
    cfg.dftfilter = 'yes';
    cfg.dftfreq = linenoise_freq;
    data = ft_preprocessing(cfg, data);
       
    % Get rid of line noise using BS filter
    cfg = [];
    cfg.bsfilter = 'yes';
    cfg.bsfreq = linenoise_bsfreq;
    data = ft_preprocessing(cfg, data);

%% Finding HGP peaks

    sfreq = data.fsample;
    freq_range = methods.freq_range;
    highGammaRange = [freq_range(2) freq_range(1)];
    
    hgp = ft_preproc_bandpassfilter(data.trial{1},sfreq,highGammaRange,...
        methods.filter_order/2,'but');  % 'filter_order/2' due to fieldtrip 
                                       % idiosyncracies. 'but': Butterworth
                                       
    for i = 1:size(hgp,1)
        hgp(i,:) = hilbert(hgp(i,:));
    end
    hgp = abs(hgp);             % obtain envelope
    if methods.zscore == 1   % converting data points to z-score
        [hgp, ~, ~] = zscore(hgp,0,2);  
    else
        methods.zscore = 0;
    end
    
    disp('Bandpassed! Smoothing...')
    
    % smoothing the hilbert
    %these lines of code were borrowed from ekaestner
    % Gaussian
        HGP_smooth_msec = 0.1;         % Extra Smooth, 100ms
        if methods.smooth ~= HGP_smooth_msec
                HGP_smooth_msec = methods.smooth;
        end
        w = sfreq*HGP_smooth_msec;

        gauss_x = -w:1:w;
        gauss_y = normpdf(gauss_x,0,w/4);
        gauss = gauss_y/sum(gauss_y);

    % create the smooth matrix, and convolve
        a_mov_gauss_gamma = zeros(size(hgp));
        
        for j = 1:size(hgp,1)
            a_mov_gauss_gamma(j,:) = conv(hgp(j,:),gauss,'same');
        end
        
    disp('Smoothed! Finding peaks...')
    
    % peak-finding with extrema.m
%     Ind = cell(size(hgp,1),1);
%     for k = 1:size(hgp,1)
%         [~,Ind{k},~,~] = extrema(a_mov_gauss_gamma(k,:));
%     end
%     disp('Peaks found. ')

%     clear hgp

%% Fitting distributions

    dim = size(a_mov_gauss_gamma,1);
    PD_list = cell(dim,1);          % preassign space for distributions and
    peak_Ind_list = cell(dim,1);    % peak indices
                                    
    % find peaks and fit distributions over the empirical dist.
    for i = 1:dim
        disp(['Processing channel ',num2str(i),':'])
        [~,Ind,~,~] = extrema(a_mov_gauss_gamma(i,:));
        peak_Ind_list{i} = Ind;
        PD_list{i} = fitdist(a_mov_gauss_gamma(i,Ind)',methods.dist);   
        disp(['Channel ',num2str(i),' processed.'])
    end
    
%% Setting thresholds based on interquartile range and find outliers
    
    if ~isfield(methods,'olc')
        methods.olc = 'IQR';
    end

    if strcmp(methods.olc,'ESD')
        % make alternate thresholds (3 times the std above mean of the 
        % probability distribution-generated data)
        POI_list = cell(dim,1);
        outlying_thresh = zeros(dim,1);
        for j = 1:dim
            temp = a_mov_gauss_gamma(j,:); % placeholder for a given channel
            temp_mean = zeros(1,methods.clt_iterations); % placeholder for sample means (calculated for CLT application)
            for k = 1:methods.clt_iterations
                Y = random(PD_list{j},1,length(peak_Ind_list{j}));
                temp_mean(k) = mean(Y);
            end
            outlying_thresh(j) = mean(temp_mean) + 3*std(temp_mean)*sqrt(length(Y));
            POI_list{j} = sort(peak_Ind_list{j}(temp(peak_Ind_list{j}) ...
                >= outlying_thresh(j)));
        end
    else        
        POI_list = cell(dim,1);
        outlying_thresh = zeros(dim,1);
        for j = 1:dim
            IQ_range = iqr(PD_list{j});     % No, not that IQ. "InterQuartile".
            IQ_bounds = icdf(PD_list{j},[0.25,0.75]);   % peak height values
            outlying_thresh(j) = IQ_bounds(2) + IQ_range*1.5;
            temp = a_mov_gauss_gamma(j,:); % placeholder for a given channel
        
        % find indices of outliers
            temp_Ind = sort(peak_Ind_list{j});
        
        % POI: points of interest (indices of outliers for all channels)
            POI_list{j} = temp_Ind(temp(temp_Ind) >= outlying_thresh(j));
        end
    end
    

    
%% Organizing output and saving

    data_out.POI = POI_list;
    if methods.zscore == 1
        data_out.hgpZscore = hgp;
    else
        data_out.hgp = hgp;
    end
    data_out.smoothed = a_mov_gauss_gamma;
    data_out.allPeakInd = peak_Ind_list;
    data_out.distributions = PD_list;
    data_out.outlyingthresh = outlying_thresh;
    save(sprintf('%s/%s/%s_D%s_%s_%s_Zxsmooth2.mat',params.savedir,params.ID,...
        params.ID,num2str(params.day),params.period,params.clinsys),'data_out','-v7.3')

end
