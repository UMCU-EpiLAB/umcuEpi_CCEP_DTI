% function to split data into stimulation trials and average
 
% author: Dorien van Blooijs
% date: June 2019

% Epoch stimuli of each SPES trial in time windows defined in
% cfg.epoch_length and cfg.epoch_prestim, time-locked to the stimulus artifact. 
% Average these epochs for each trial per electrode contact.

function dataBase = preprocess_ECoG(dataBase,cfg)
epoch_length = cfg.epoch_length;
epoch_prestim = cfg.epoch_prestim;

% minimal number of stimulations needed to be included in further analysis
minstim = 5;

for nSubj = 1:size(dataBase,2)
    
    for nRun = 1:size(dataBase(nSubj).metadata,2)
   
        %% unique stimulation pairs
        stimpair = dataBase(nSubj).metadata(nRun).tb_events.electrical_stimulation_site(contains(dataBase(nSubj).metadata(nRun).tb_events.sub_type,'SPES') & ...
            ~contains(dataBase(nSubj).metadata(nRun).tb_events.electrical_stimulation_site,'n/a')) ;
        
        stimnum = NaN(size(stimpair,1),2);
        for stimp = 1:size(stimpair,1)
            stimchans = strsplit(stimpair{stimp},'-');
            for chan = 1:2
                stimnum(stimp,chan) = find(strcmp(stimchans{chan},dataBase(nSubj).metadata(nRun).ch)==1);
            end
        end
        
        stimelek = sort(stimnum,2);
        [cc_stimsets,~,IC] = unique(stimelek,'rows');
        
        nTrial = histcounts(IC,'BinMethod','integers');
        
        if any(diff(nTrial) ~= 0)
            indivstimremove = find(nTrial<minstim); % remove al stimulation pairs that are stimulated less than 5 times
            
            stimremove = sum(IC==indivstimremove,2);
            stimelek(stimremove>0,:) = [];
            
            % if there are still pairs with more or less stimuli than
            % median, give warning
            [cc_stimsets,~,IC] = unique(stimelek,'rows');
            nTrial = histcounts(IC,'BinMethod','integers');
            if any(diff(nTrial) ~= 0)
                unequalstim = find(nTrial ~= median(nTrial));
                for i=1:size(unequalstim,2)
                    
                    if nTrial(unequalstim(i)) < median(nTrial)
                        
                        warning('%s %s: %s-%s is stimulated less (%dx) than all others (%dx) \n',...
                            dataBase(nSubj).sub_label,dataBase(nSubj).metadata(nRun).run_label,...
                            dataBase(nSubj).metadata(nRun).ch{cc_stimsets(unequalstim(i),1)},...
                            dataBase(nSubj).metadata(nRun).ch{cc_stimsets(unequalstim(i),2)},...
                            nTrial(unequalstim(i)),...
                            median(nTrial))
                    elseif nTrial(unequalstim(i)) > median(nTrial)
                        warning('%s %s: %s-%s is stimulated more (%dx) than all others (%dx) \n',...
                            dataBase(nSubj).sub_label,dataBase(nSubj).metadata(nRun).run_label,...
                            dataBase(nSubj).metadata(nRun).ch{cc_stimsets(unequalstim(i),1)},...
                            dataBase(nSubj).metadata(nRun).ch{cc_stimsets(unequalstim(i),2)},...
                            nTrial(unequalstim(i)),...
                            median(nTrial))
                    end
                end
            end
        end
        
        % determine stimulus channels (instead of numbers)
        % pre-allocation
        cc_stimchans = cell(size(cc_stimsets,1),2);
        
        for stimp = 1:size(cc_stimsets,1)
            for chan =1:2
                cc_stimchans{stimp,chan} = dataBase(nSubj).metadata(nRun).ch{cc_stimsets(stimp,chan)};
            end
        end
        
        max_stim = median(nTrial);
        
        % write stimsets etc to dataBase struct
        dataBase(nSubj).metadata(nRun).cc_stimsets = cc_stimsets;
        dataBase(nSubj).metadata(nRun).cc_stimchans = cc_stimchans;
        dataBase(nSubj).metadata(nRun).max_stim = max_stim;
        
        %% select epochs
        t = round(epoch_length*dataBase(nSubj).metadata(nRun).ccep_header.Fs);
        
        % allocation
        cc_epoch_sorted = NaN(size(dataBase(nSubj).metadata(nRun).data,1),dataBase(nSubj).metadata(nRun).max_stim,size(dataBase(nSubj).metadata(nRun).cc_stimsets,1),t); % [channels, trials, stimpairs, samples]
        tt_epoch_sorted = NaN(dataBase(nSubj).metadata(nRun).max_stim,size(dataBase(nSubj).metadata(nRun).cc_stimsets,1),t); % samplenumbers for each epoch
        cc_epoch_sorted_avg = NaN(size(dataBase(nSubj).metadata(nRun).data,1),size(dataBase(nSubj).metadata(nRun).cc_stimsets,1),t); % [channels, stimpairs, samples]
        
        for nElec = 1:size(dataBase(nSubj).metadata(nRun).data,1) % for all channels
            for nStimp = 1:size(dataBase(nSubj).metadata(nRun).cc_stimsets,1) % for all epochs with >4 stimuli
                                
                eventnum1 = find(strcmp(dataBase(nSubj).metadata(nRun).tb_events.electrical_stimulation_site,...
                    [dataBase(nSubj).metadata(nRun).cc_stimchans{nStimp,1}, '-',dataBase(nSubj).metadata(nRun).cc_stimchans{nStimp,2}]));
                eventnum2 = find(strcmp(dataBase(nSubj).metadata(nRun).tb_events.electrical_stimulation_site,...
                    [dataBase(nSubj).metadata(nRun).cc_stimchans{nStimp,2}, '-',dataBase(nSubj).metadata(nRun).cc_stimchans{nStimp,1}]));
                eventnum = [eventnum1;eventnum2];

                if size(eventnum,1) > dataBase(nSubj).metadata(nRun).max_stim
                    events = dataBase(nSubj).metadata(nRun).max_stim;
                else
                    events = size(eventnum,1);
                end
                
                for nTrial = 1:events
                    
                    if dataBase(nSubj).metadata(nRun).tb_events.sample_start(eventnum(nTrial))-round(epoch_prestim*dataBase(nSubj).metadata(nRun).ccep_header.Fs)+1< 0
                        % do nothing, the start of the selected epoch is before
                        % the start of the data-file (stimulus is less than 2s
                        % (epoch_prestim) after start of the data recording)
                    elseif dataBase(nSubj).metadata(nRun).tb_events.sample_start(eventnum(nTrial))+round((epoch_length-epoch_prestim)*dataBase(nSubj).metadata(nRun).ccep_header.Fs)+1 > size(dataBase(nSubj).metadata(nRun).data,2)
                        % do nothing, the end of the selected epoch is after
                        % the end of the data-file (stimulus is less than 3s
                        % (epoch_length - epoch_prestim) before end of the data
                        % recording)
                    else
                        
                        singleTrial = dataBase(nSubj).metadata(nRun).data(nElec,dataBase(nSubj).metadata(nRun).tb_events.sample_start(eventnum(nTrial))-round(epoch_prestim*dataBase(nSubj).metadata(nRun).ccep_header.Fs)+1:...
                            dataBase(nSubj).metadata(nRun).tb_events.sample_start(eventnum(nTrial))+round((epoch_length-epoch_prestim)*dataBase(nSubj).metadata(nRun).ccep_header.Fs));
                        
                        % create time struct
                        tt = (1:epoch_length*dataBase(nSubj).metadata(nRun).ccep_header.Fs) / ...
                            dataBase(nSubj).metadata(nRun).ccep_header.Fs - epoch_prestim;
                        
                        % baseline subtraction: take median of part of the averaged signal for
                        % this stimulation pair before stimulation, which is the half of the
                        % epoch
                        baseline_tt = tt>-2 & tt<-.1;
                        
                        cc_epoch_sorted(nElec,nTrial,nStimp,:) = singleTrial - median(singleTrial(baseline_tt));
                        tt_epoch_sorted(nTrial,nStimp,:) = dataBase(nSubj).metadata(nRun).tb_events.sample_start(eventnum(nTrial))-round(epoch_prestim*dataBase(nSubj).metadata(nRun).ccep_header.Fs)+1:...
                            dataBase(nSubj).metadata(nRun).tb_events.sample_start(eventnum(nTrial))+round((epoch_length-epoch_prestim)*dataBase(nSubj).metadata(nRun).ccep_header.Fs);
                    end
                end
            end
        end
        
        cc_epoch_sorted_avg(:,1:nStimp,:) = squeeze(mean(cc_epoch_sorted,2,'omitnan'));
        
        dataBase(nSubj).metadata(nRun).cc_epoch_sorted = cc_epoch_sorted;
        dataBase(nSubj).metadata(nRun).tt_epoch_sorted = tt_epoch_sorted;
        dataBase(nSubj).metadata(nRun).cc_epoch_sorted_avg = cc_epoch_sorted_avg;
        
        fprintf('...%s %s has been epoched and averaged... \n',...
            dataBase(nSubj).sub_label,dataBase(nSubj).metadata(nRun).run_label)
    end
end