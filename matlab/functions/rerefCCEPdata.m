function dataBase = rerefCCEPdata(dataBase,subj,run,cfg)

fs = dataBase(subj).metadata(run).ccep_header.Fs;
tt = (1:cfg.epoch_length*fs)/fs - cfg.epoch_prestim;

for stimp = 1:size(dataBase(subj).metadata(run).cc_epoch_sorted,3) % for each stimulus pair

    these_epochs_data = squeeze(dataBase(subj).metadata(run).cc_epoch_sorted(:,:,stimp,:));
    these_epochs_data_reref = NaN(size(these_epochs_data));

    stimChan = dataBase(subj).metadata(run).cc_stimsets(stimp,:);
    these_epochs_data(stimChan,:,:) = NaN;

    for numstim = 1:size(these_epochs_data,2) % for each of the trials of one stimulus pair

        % calculate common avarage reference
        CAR = squeeze(mean(these_epochs_data(:,numstim,:),1,'omitnan'));

        % calculate the variance in 400ms prior to stimulus and
        % during 100ms post stimulus
        varCARprio = var(CAR(tt>-0.5&tt<-0.1));
        varChanprio = var(squeeze(these_epochs_data(:,numstim,tt>-0.5&tt<-0.1)),[],2,'omitnan');
        varChanccep = var(squeeze(these_epochs_data(:,numstim,find(tt>0.01,1):find(tt<0.1,1,'last'))),[],2,'omitnan');

        % determine which channels have a variance lower than
        % the CAR
        idx_usedCAR = varChanccep < varCARprio & varChanprio < varCARprio;

        % I want minimal 5% of the electrodes included in the CARlowvar. If
        % the above does not result in at least 10% of the electrodes, I
        % will include the electrodes with lowest variance in the period
        % for ccep (0.01:0.1s), but first I remove the channels that are
        % already used in CARlowvar in the previous step
        [~,sort_varChan] = sort(varChanccep,'ascend');
        sort_varChan_excl_usedCAR = sort_varChan;
        sort_varChan_excl_usedCAR(ismember(sort_varChan,find(idx_usedCAR==1))) = [];

        % determine which channels are included in CAR, with a minimum of
        % 10% of the total amount of channels
        totChan = sum(~isnan(these_epochs_data(:,1,1)))+2; %total good channels, including the two stimulus channels
        usedCAR = [find(idx_usedCAR ==1); sort_varChan_excl_usedCAR(1:(ceil(0.1*totChan)-sum(idx_usedCAR)))];
        CARlowvar = squeeze(median(these_epochs_data(usedCAR,numstim,:),1,'omitnan'));

        these_epochs_data_reref(:,numstim,:) = squeeze(these_epochs_data(:,numstim,:)) - repmat(transpose(CARlowvar),[size(these_epochs_data,1),1]);

        dataBase(subj).metadata(run).ref(:,numstim,stimp,:) = CARlowvar;

    end

    % re-reference the individual trials
    dataBase(subj).metadata(run).cc_epoch_sorted_reref(:,:,stimp,:) = these_epochs_data_reref;

    fprintf('--- Stimpair %s-%s ---\n', ...
        dataBase(subj).metadata(run).cc_stimchans{stimp,1},...
        dataBase(subj).metadata(run).cc_stimchans{stimp,2})

end

end