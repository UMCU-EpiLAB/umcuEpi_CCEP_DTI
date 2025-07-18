function [dataBase] = detect_n1peak_ECoG_ccep(dataBase, cfg)

% Function for detecting the the N1 peaks of CCEPs

% INPUTS:
% - database
%   Structure with metadata and averaged epoched data split into 
%   [channels x stimulations x time].
% - amplitude_thresh
%   Threshold factor that needs to be exceeded to detect a peak as
%   significant peak. This factor is multiplied by the minimal 
%   pre-stimulation standard deviation of 50 uV. 
%       * A standard deviation of 50uV, amplitude_thresh of 3.4 are 
%       recommended for a conservative algorithm
%       * A standard deviation of 50uV, amplitude_thresh of 2.6 are 
%       recommended for a more sensitive algorithm.
% - n1_peak_range
%   End point (in ms) of the range in which the algorithm will search for the 
%   peak of the N1. Algorithm will search between 10ms and the end point.

% OUTPUTS:
% - database
%   To the database structure the following matrices will be added:
% - n1_peak_sample 
%   matrix[channels x stimulations] containing latency of the significant
%   n1 peaks that are detected.
% - n1_peak_amplitude
%   matrix[channels x stimulations] containing the amplitude of the
%   significant n1 peaks that are detected.

% Peakfinder function is used for peak detection, add to path

% This function is heavily based on the function and algorithm developed by
% Dorien van Blooijs during her master's thesis 'Improving the SPES protocol
% by automating ER and DR detection and evaluation of the spatial relation
% between ERs and DRs. Algorithm is validated during the thesis and used in
% publication (van Blooijs et al., 2018). 

% Some slight differences:
% - While algorithm of van Blooijs detects both positive as negative peaks,
%   this algorithm only selects negative peaks
% - After validation the algorithm of van Blooijs used 125uV as amplitude threshold
%   and 100ms as n1 peak range, this algorithm performs best with 140 uV
%   and 80 ms. A reason for this possibly is the other ages of patients which are
%   used during validation. The first validation is done in only younger
%   (9yrs - 12 yrs.) patients, while later validation is also done in older
%   patients. Age seems to effect the characteristics of CCEPs and
%   therefore also the detection. 

% This version of the algorithm is validated in three patients 
% (age 9, 21 and 50) by S Jelsma. Validation done by comparing the 
% performance of the code and visual assesment. Algorithm optimized by 
% setting different parameters and comparing their performances. 

% FIXED PARAMETERS (that are validated by both van Blooijs & van der Aar):
% - sel = 20, which is how many samples around peak not considered as another peak
% - minSD = 50, minimum standard deviation (in uV) of prestimulus baseline

% original author: Dorien van Blooijs, UMC Utrecht, January 2018
% modified by: Jaap van der Aar, Dora Hermes, Dorien van Blooijs, Giulio Castegnaro, UMC Utrecht, 2019

%%
amplitude_thresh = cfg.amplitude_thresh;
n1_peak_range = cfg.n1_peak_range;
epoch_prestim = cfg.epoch_prestim;
epoch_length = cfg.epoch_length;
minSD = cfg.minSD;
sel = cfg.sel;

% iterate over all subjects in database
for nSubj = 1:length(dataBase)
    % iterate over all their runs
    for nRun = 1:size(dataBase(nSubj).metadata,2)
        
        bad_ch = find(strcmp(dataBase(nSubj).metadata(nRun).tb_channels.status,'bad')); % this is also when in silicon is yes (so on top of each other)
        
        signal = dataBase(nSubj).metadata(nRun).cc_epoch_sorted_reref_avg;
        
        % pre-allocation: output in [channels X stimulations X latency amplitude]
        n1_peak = NaN(size(signal,1), ...
            size(signal,2),2);

        % for every averaged stimulation
        for nStimp = 1:size(dataBase(nSubj).metadata(nRun).cc_epoch_sorted_avg,2)
            % for every channel
            for nChan = 1:size(dataBase(nSubj).metadata(nRun).cc_epoch_sorted_avg,1)

                % create time struct 
                tt = (1:epoch_length*dataBase(nSubj).metadata(nRun).ccep_header.Fs) / ...
                    dataBase(nSubj).metadata(nRun).ccep_header.Fs - epoch_prestim;
               
                % baseline subtraction: take median of part of the averaged signal for
                % this stimulation pair before stimulation, which is the half of the
                % epoch                
                baseline_tt = tt>-2 & tt<-.1;
                signal_median = median(signal(nChan,nStimp,baseline_tt),3,'omitnan');

                % subtract median baseline from signal
                new_signal = squeeze(signal(nChan,nStimp,:)) - signal_median;
                % testplot new signal: plot(tt,squeeze(new_signal))

                % save the corrected signal again
                signal(nChan,nStimp,:) = new_signal;

                % take area before the stimulation of the new signal and calculate its SD
                pre_stim_sd = std(new_signal(baseline_tt));
                
                % if the pre_stim_sd is smaller that the minimally needed SD, 
                % which is validated as 50 uV, use this the minSD as pre_stim_sd
                if pre_stim_sd < minSD
                    pre_stim_sd = minSD;
                end

                % when the electrode is stimulated
                if nChan == dataBase(nSubj).metadata(nRun).cc_stimsets(nStimp,1) || ...
                        nChan == dataBase(nSubj).metadata(nRun).cc_stimsets(nStimp,2)
                    n1_peak_sample = NaN;
                    n1_peak_amplitude = NaN;

                elseif ismember(nChan , bad_ch)
                    % do nothing because noisy/bad electrode
                    n1_peak_sample = NaN;
                    n1_peak_amplitude = NaN;
                    
                else % in other electrode

                    % use peakfinder to find all negative peaks and their
                    % amplitude.
                    % sel = 20 , which is how many samples around a peak not considered as another peak
                    [all_sampneg, all_amplneg] = peakfinder(new_signal(find(tt>0,1):find(tt>0.5,1)),sel,[],-1);

                    % If the first selected sample is a peak, this is not a real peak,
                    % so delete
                    all_amplneg(all_sampneg==1) = [];
                    all_sampneg(all_sampneg==1) = [];

                    % convert back timepoints based on tt, substract 1 because before
                    % the first sample after stimulation is taken
                    all_sampneg = all_sampneg + find(tt>0,1) - 1;

                    % set the starting range in which the algorithm looks
                    % for a peak. At least 18 samples are necessary because
                    % of the micromed amplifier does not record the
                    % stimulated electrode before this. Peak detection
                    % start 9 ms after stimulation, which is 19 samples
                    % after stimulation (with sample frequency = 2048 Hz)
                    n1_samples_start = find(tt>(19/dataBase(nSubj).metadata(nRun).ccep_header.Fs),1);
                    
                    % find first sample that corresponds with the given n1
                    % peak range
                    n1_samples_end = find(tt>(n1_peak_range/1000),1);
                    
                    % for N1, first select the range in which the N1 could appear, and
                    % select the peaks found in this range
                    temp_n1_peaks_samp = all_sampneg((n1_samples_start <= all_sampneg) & (all_sampneg <= n1_samples_end));
                    temp_n1_peaks_ampl = all_amplneg((n1_samples_start <= all_sampneg) & (all_sampneg <= n1_samples_end));

                    % if peak(s) found, select biggest peak
                    if ~isempty(temp_n1_peaks_samp)
                        max_n1_ampl = find(abs(temp_n1_peaks_ampl) == max(abs(temp_n1_peaks_ampl)));
                        n1_peak_sample = temp_n1_peaks_samp(max_n1_ampl(1));
                        n1_peak_amplitude = temp_n1_peaks_ampl(max_n1_ampl(1));
                        % otherwise give the amplitude the value NaN
                    elseif isempty(temp_n1_peaks_samp)
                        n1_peak_sample = NaN;
                        n1_peak_amplitude = NaN;
                    end

                    % if N1 exceeds positive threshold, it is deleted
                    if temp_n1_peaks_ampl > 0
                       n1_peak_sample = NaN;
                       n1_peak_amplitude = NaN;   
                    end

                    % when peak amplitude is saturated, it is deleted
                    if abs(n1_peak_amplitude) > 3000
                        n1_peak_sample = NaN;
                        n1_peak_amplitude = NaN;
                    end

                    % if the peak is not big enough to consider as a peak, assign NaN
                    if abs(n1_peak_amplitude) < amplitude_thresh* abs(pre_stim_sd)
                        n1_peak_sample = NaN;
                        n1_peak_amplitude = NaN;
                    end
                end
              
                % add properties to output frame 
                n1_peak(nChan,nStimp,1) = n1_peak_sample;
                n1_peak(nChan,nStimp,2) = n1_peak_amplitude;

            end
        end
        
        % write n1_peak (sample and amplitude) to database
        dataBase(nSubj).metadata(nRun).ccep.n1_peak_sample = n1_peak(:,:,1);
        dataBase(nSubj).metadata(nRun).ccep.n1_peak_amplitude = n1_peak(:,:,2);
        
        % baseline corrected signal
        dataBase(nSubj).metadata(nRun).cc_epoch_sorted_reref_avg = signal;

        dataBase(nSubj).metadata(nRun).ccep.amplitude_thresh = cfg.amplitude_thresh;
        dataBase(nSubj).metadata(nRun).ccep.n1_peak_range = cfg.n1_peak_range;
        dataBase(nSubj).metadata(nRun).ccep.minSD = cfg.minSD;
        dataBase(nSubj).metadata(nRun).ccep.sel = cfg.sel;

    end
end