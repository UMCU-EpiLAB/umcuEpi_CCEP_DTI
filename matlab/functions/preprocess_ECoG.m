function dataBase = preprocess_ECoG(dataBase,cfg)
epoch_length = cfg.epoch_length;
epoch_prestim = cfg.epoch_prestim;

% minimal number of stimulations needed to be included in further analysis
minstim = 5;

for subj = 1:size(dataBase,2)
    
    for run = 1:size(dataBase(subj).metadata,2)
        %% define artefact period
        ev_artefact_start = dataBase(subj).metadata(run).tb_events.sample_start(strcmp(dataBase(subj).metadata(run).tb_events.trial_type,'artefact'));
        ev_artefact_stop = dataBase(subj).metadata(run).tb_events.sample_end(strcmp(dataBase(subj).metadata(run).tb_events.trial_type,'artefact'));
        
        if iscell(ev_artefact_start)
            ev_artefact_start = str2double(ev_artefact_start);
        end
        if iscell(ev_artefact_stop)
            ev_artefact_stop = str2double(ev_artefact_stop);
        end
        
        ev_artefact = [];
        for i=1:size(ev_artefact_start,1)
            ev_artefact = [ev_artefact, ev_artefact_start(i):ev_artefact_stop(i)]; %#ok<AGROW>
        end
        
        %% unique stimulation pairs
        stimpair = dataBase(subj).metadata(run).tb_events.electrical_stimulation_site(contains(dataBase(subj).metadata(run).tb_events.sub_type,'SPES') & ...
            ~contains(dataBase(subj).metadata(run).tb_events.electrical_stimulation_site,'n/a')) ;
        
        stimnum = NaN(size(stimpair,1),2);
        for stimp = 1:size(stimpair,1)
            stimchans = strsplit(stimpair{stimp},'-');
            for chan = 1:2
                stimnum(stimp,chan) = find(strcmp(stimchans{chan},dataBase(subj).metadata(run).ch)==1);
            end
        end
        
        stimcur = str2double(dataBase(subj).metadata(run).tb_events.electrical_stimulation_current(contains(dataBase(subj).metadata(run).tb_events.sub_type,'SPES') & ...
            ~contains(dataBase(subj).metadata(run).tb_events.electrical_stimulation_site,'n/a')));
        
        if strcmp(cfg.dir,'yes') && strcmp(cfg.amp,'yes')
            stimelek = [stimnum stimcur];
        elseif strcmp(cfg.dir,'yes') && strcmp(cfg.amp,'no')
            stimelek = stimnum;
        elseif strcmp(cfg.dir,'no') && strcmp(cfg.amp,'yes')
            stimelek = [sort(stimnum,2) stimcur];
        elseif strcmp(cfg.dir,'no') && strcmp(cfg.amp,'no')
            stimelek = sort(stimnum,2);
        end
        
        [cc_stimsets,~,IC] = unique(stimelek,'rows');
        
        n = histcounts(IC,'BinMethod','integers');
        
        if any(diff(n) ~= 0)
            indivstimremove = find(n<minstim); % remove al stimulation pairs that are stimulated less than 5 times
            
            stimremove = sum(IC==indivstimremove,2);
            stimelek(stimremove>0,:) = [];
            
            % if there are still pairs with more or less stimuli than
            % median, give warning
            [cc_stimsets,~,IC] = unique(stimelek,'rows');
            n = histcounts(IC,'BinMethod','integers');
            if any(diff(n) ~= 0)
                unequalstim = find(n ~= median(n));
                for i=1:size(unequalstim,2)
                    
                    if n(unequalstim(i)) < median(n)
                        
                        warning('%s %s: %s-%s is stimulated less (%dx) than all others (%dx) \n',...
                            dataBase(subj).sub_label,dataBase(subj).metadata(run).run_label,...
                            dataBase(subj).metadata(run).ch{cc_stimsets(unequalstim(i),1)},...
                            dataBase(subj).metadata(run).ch{cc_stimsets(unequalstim(i),2)},...
                            n(unequalstim(i)),...
                            median(n))
                    elseif n(unequalstim(i)) > median(n)
                        warning('%s %s: %s-%s is stimulated more (%dx) than all others (%dx) \n',...
                            dataBase(subj).sub_label,dataBase(subj).metadata(run).run_label,...
                            dataBase(subj).metadata(run).ch{cc_stimsets(unequalstim(i),1)},...
                            dataBase(subj).metadata(run).ch{cc_stimsets(unequalstim(i),2)},...
                            n(unequalstim(i)),...
                            median(n))
                    end
                end
            end
        end
        
        % determine stimulus channels (instead of numbers)
        % pre-allocation
        cc_stimchans = cell(size(cc_stimsets,1),2);
        
        for stimp = 1:size(cc_stimsets,1)
            for chan =1:2
                cc_stimchans{stimp,chan} = dataBase(subj).metadata(run).ch{cc_stimsets(stimp,chan)};
            end
            
        end
        
        max_stim = median(n);
        
        % write stimsets etc to dataBase struct
        dataBase(subj).metadata(run).cc_stimsets = cc_stimsets;
        dataBase(subj).metadata(run).cc_stimchans = cc_stimchans;
        dataBase(subj).metadata(run).max_stim = max_stim;
        
        %% select epochs
        t = round(epoch_length*dataBase(subj).metadata(run).ccep_header.Fs);
        
        % allocation
        cc_epoch_sorted = NaN(size(dataBase(subj).metadata(run).data,1),dataBase(subj).metadata(run).max_stim,size(dataBase(subj).metadata(run).cc_stimsets,1),t); % [channels, trials, stimpairs, samples]
        tt_epoch_sorted = NaN(dataBase(subj).metadata(run).max_stim,size(dataBase(subj).metadata(run).cc_stimsets,1),t); % samplenumbers for each epoch
        cc_epoch_sorted_avg = NaN(size(dataBase(subj).metadata(run).data,1),size(dataBase(subj).metadata(run).cc_stimsets,1),t); % [channels, stimpairs, samples]
        
        for elec = 1:size(dataBase(subj).metadata(run).data,1) % for all channels
            for ll = 1:size(dataBase(subj).metadata(run).cc_stimsets,1) % for all epochs with >4 stimuli
                
                % use either all negative and positive stimuli as one trial, or
                % split these stimuli (C01-C02 and C02-C01 for example)
                if strcmp(cfg.dir,'no')
                    eventnum1 = find(strcmp(dataBase(subj).metadata(run).tb_events.electrical_stimulation_site,...
                        [dataBase(subj).metadata(run).cc_stimchans{ll,1}, '-',dataBase(subj).metadata(run).cc_stimchans{ll,2}]));
                    eventnum2 = find(strcmp(dataBase(subj).metadata(run).tb_events.electrical_stimulation_site,...
                        [dataBase(subj).metadata(run).cc_stimchans{ll,2}, '-',dataBase(subj).metadata(run).cc_stimchans{ll,1}]));
                    eventnum = [eventnum1;eventnum2];
                    
                elseif strcmp(cfg.dir,'yes')
                    eventnum = find(strcmp(dataBase(subj).metadata(run).tb_events.electrical_stimulation_site,...
                        [dataBase(subj).metadata(run).cc_stimchans{ll,1}, '-',dataBase(subj).metadata(run).cc_stimchans{ll,2}]));
                end
                
                if size(eventnum,1) > dataBase(subj).metadata(run).max_stim
                    events = dataBase(subj).metadata(run).max_stim;
                else
                    events = size(eventnum,1);
                end
                
                for n = 1:events
                    
                    if dataBase(subj).metadata(run).tb_events.sample_start(eventnum(n))-round(epoch_prestim*dataBase(subj).metadata(run).ccep_header.Fs)+1< 0
                        % do nothing, the start of the selected epoch is before
                        % the start of the data-file (stimulus is less than 2s
                        % (epoch_prestim) after start of the data recording)
                    elseif dataBase(subj).metadata(run).tb_events.sample_start(eventnum(n))+round((epoch_length-epoch_prestim)*dataBase(subj).metadata(run).ccep_header.Fs)+1 > size(dataBase(subj).metadata(run).data,2)
                        % do nothing, the end of the selected epoch is after
                        % the end of the data-file (stimulus is less than 3s
                        % (epoch_length - epoch_prestim) before end of the data
                        % recording)
%                     elseif ismember(dataBase(subj).metadata(run).tb_events.sample_start(eventnum(n)),ev_artefact)
                        % do nothing, because part of artefact
                    else
                        
                        singleTrial = dataBase(subj).metadata(run).data(elec,dataBase(subj).metadata(run).tb_events.sample_start(eventnum(n))-round(epoch_prestim*dataBase(subj).metadata(run).ccep_header.Fs)+1:...
                            dataBase(subj).metadata(run).tb_events.sample_start(eventnum(n))+round((epoch_length-epoch_prestim)*dataBase(subj).metadata(run).ccep_header.Fs));
                        
                        % create time struct
                        tt = (1:epoch_length*dataBase(subj).metadata(run).ccep_header.Fs) / ...
                            dataBase(subj).metadata(run).ccep_header.Fs - epoch_prestim;
                        
                        % baseline subtraction: take median of part of the averaged signal for
                        % this stimulation pair before stimulation, which is the half of the
                        % epoch
                        baseline_tt = tt>-2 & tt<-.1;
                        
                        cc_epoch_sorted(elec,n,ll,:) = singleTrial - median(singleTrial(baseline_tt));
                        tt_epoch_sorted(n,ll,:) = dataBase(subj).metadata(run).tb_events.sample_start(eventnum(n))-round(epoch_prestim*dataBase(subj).metadata(run).ccep_header.Fs)+1:...
                            dataBase(subj).metadata(run).tb_events.sample_start(eventnum(n))+round((epoch_length-epoch_prestim)*dataBase(subj).metadata(run).ccep_header.Fs);
                    end
                end
            end
        end
        
        cc_epoch_sorted_avg(:,1:ll,:) = squeeze(mean(cc_epoch_sorted,2,'omitnan'));
        
        dataBase(subj).metadata(run).cc_epoch_sorted = cc_epoch_sorted;
        dataBase(subj).metadata(run).tt_epoch_sorted = tt_epoch_sorted;
        dataBase(subj).metadata(run).cc_epoch_sorted_avg = cc_epoch_sorted_avg;
        
        fprintf('...%s %s has been epoched and averaged... \n',...
            dataBase(subj).sub_label,dataBase(subj).metadata(run).run_label)
    end
end