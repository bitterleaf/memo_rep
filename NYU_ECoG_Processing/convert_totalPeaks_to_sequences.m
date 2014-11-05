% put the peaks in temporal order over all the channels as well as bin these orders into time bins
%
% Dependency: N/A
%
% Created by Isaac Shamie
% Last edited by Xi Jiang, Nov. 5th
%
% Output:
%   sequences:           numerical (channel-representing, e.g. "21" means 
%                        "MT03") sequences of temporal order in which HGP 
%                        peak outliers occur within a time bin
%   sequences_plus_time: same as above, with each time bin represented by
%                        two columns instead, one containing the channel 
%                        sequence, the other the sampling point indices 
%                        denoting where the peaks were found
%   order_peaks:         sequences of all time bins concatenated
%   bins_with_same_peak_times: time bins where some/all HGP peaks occur at 
%                              the exact same sampling point
%
% Input:
%   good_peaks: an N cell array where N = num chans and in each cell is a matrix with the HGP peak_times (in edf sample points)  
%   methods:    structure that contains at least the following fields:
%                 -time_bin: size of the time bin in seconds
%                 -sfreq: sampling frequency
%   params:     structure that contains the following fields:
%                 -savedir: directory to which the output file is saved
%                 -ID:      subject ID
%                 -day:     dayfile number
%                 -period:  morning or evening ('morn' or 'eve')
%                 -clinsys: clinical system from which the data is obtained, 
%                           can be 'clin1' or 'clin2'
function [sequences,sequences_plus_time, order_peaks, bins_with_same_peak_times] = convert_totalPeaks_to_sequences(good_peaks,methods,params)
%% params
fsample = methods.sfreq; %512;
seconds = methods.time_bin; %2;

dir = params.savedir; % '/space/mdeh4/1/halgdev/projects/ishamie/24_7/NY394/394_D3_eve';
ID = params.ID; %'394';
day = params.day; %'3';
period = params.period; %'eve';
clinsys = params.clinsys;  %'clin1';


%% put all cells in a huge matrix 
totalOrder = zeros(1,2);%zeros(sum(cellfun(@length,good_peaks)),2);
current = 1; %the current index in the totalOrder
for i = 1:length(good_peaks)
  chan_peaks = good_peaks{i};
  num_peaks = length(chan_peaks);
  temp = [chan_peaks;i * ones(1,num_peaks)]'; 
  totalOrder(current:current+num_peaks-1,:) = temp;
  current = current + num_peaks;  
end

%% sort the peaks, along with its respective channel
[Y,I] = sort(totalOrder(:,1));
B = totalOrder(I,2);
order_peaks = [Y,B];

%% construct the sequence for every time bin
edfLength = fsample*seconds;
current = 1; %this is to save time looking through the totalOrder matrix
num_sequences = ceil(order_peaks(end,1)/edfLength);
sequences_plus_time = cell(1,num_sequences);
sequences = cell(1,num_sequences);
bins_with_same_peak_times = {};
for bin = 1:edfLength:num_sequences*edfLength
  sequences_plus_time{(bin-1)/edfLength + 1} = order_peaks(order_peaks(:,1) >= bin & order_peaks(:,1) < bin+edfLength,:);
  sequences{(bin-1)/edfLength + 1} = order_peaks(order_peaks(:,1) >= bin & order_peaks(:,1) < bin+edfLength,2);
  if size(unique(sequences_plus_time{(bin-1)/edfLength+1}(:,1)),1) ~= size(sequences_plus_time{(bin-1)/edfLength+1}(:,1),1)
     bins_with_same_peak_times{end+1} = (bin-1)/edfLength + 1; 
  end
end

fname =  sprintf('%s/%s_D%s_%s_%s_temporal_sequence_2s_bins.mat',dir, ID,day,period,clinsys);
save(fname, 'sequences', 'sequences_plus_time', 'order_peaks', 'bins_with_same_peak_times','-v7.3')