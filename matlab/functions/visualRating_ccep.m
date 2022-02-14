
function dataBase = visualRating_ccep(dataBase,myDataPath,subj,run,cfg,endstimp)

%%
filefolder = fullfile(myDataPath.CCEPpath,dataBase(subj).sub_label,dataBase(subj).ses_label,dataBase(subj).metadata(run).run_label);
if ~exist(filefolder,'dir')
    mkdir(filefolder)
end

if any(strcmp(fieldnames(dataBase(subj).metadata(run)),'ccep'))
    ccep = dataBase(subj).metadata(run).ccep;
end
fs = dataBase(subj).metadata(run).ccep_header.Fs;
tt = -cfg.epoch_prestim+1/fs:1/fs:cfg.epoch_length-cfg.epoch_prestim;

cfg.n1Detected = 'y'; % --> if n1 peaks are detected
% if n1 peaks are not detected yet: pre-allocation
if sum(strcmp(fieldnames(dataBase(subj).metadata(run)),'ccep'))==0
    ccep.n1_peak_sample = ones(size(dataBase(subj).metadata(run).cc_epoch_sorted_avg,1),size(dataBase(subj).metadata(run).cc_epoch_sorted_avg,2));
    ccep.n1_peak_amplitude = zeros(size(dataBase(subj).metadata(run).cc_epoch_sorted_avg,1),size(dataBase(subj).metadata(run).cc_epoch_sorted_avg,2));
    cfg.n1Detected = 'n';
elseif  sum(strcmp(fieldnames(dataBase(subj).metadata(run).ccep),'n1_peak_amplitude'))==0
    ccep.n1_peak_sample = ones(size(dataBase(subj).metadata(run).cc_epoch_sorted_avg,1),size(dataBase(subj).metadata(run).cc_epoch_sorted_avg,2));
    ccep.n1_peak_amplitude = zeros(size(dataBase(subj).metadata(run).cc_epoch_sorted_avg,1),size(dataBase(subj).metadata(run).cc_epoch_sorted_avg,2));
    cfg.n1Detected = 'n';
end

% if no n1 peaks are checked yet: pre-allocation
if sum(strcmp(fieldnames(dataBase(subj).metadata(run)),'ccep'))==0
    ccep.checked = zeros(size(dataBase(subj).metadata(run).cc_epoch_sorted_avg,1),size(dataBase(subj).metadata(run).cc_epoch_sorted_avg,2));
elseif  sum(strcmp(fieldnames(dataBase(subj).metadata(run).ccep),'checked'))==0
    ccep.checked = zeros(size(dataBase(subj).metadata(run).cc_epoch_sorted_avg,1),size(dataBase(subj).metadata(run).cc_epoch_sorted_avg,2));
end

% if seeg, use only gray matter channels, hippocampus, amygdala, lesion,
% gliosis
if any(contains(fieldnames(dataBase(subj).tb_electrodes),'graymatter'))
    idx_screw = strcmpi(dataBase(subj).tb_electrodes.screw,'yes');
    idx_csf = strcmpi(dataBase(subj).tb_electrodes.csf,'yes');
    idx_whitematter = strcmpi(dataBase(subj).tb_electrodes.whitematter,'yes') & strcmpi(dataBase(subj).tb_electrodes.graymatter,'no');
    idx_all = sum([idx_screw, idx_csf, idx_whitematter],2);

    elec_include = cell(size(dataBase(subj).tb_electrodes,1),1);
    [elec_include{idx_all>0}] = deal('no');
    [elec_include{idx_all==0}] = deal('yes');

else
    elec_include = cell(size(dataBase(subj).tb_electrodes,1),1);
    [elec_include{1:end}] = deal('yes');
end

n=numel(ccep.checked(:,1:endstimp))+1;
for stimp = endstimp+1:size(dataBase(subj).metadata(run).cc_epoch_sorted_avg,2)
    for chan = 1:size(dataBase(subj).metadata(run).cc_epoch_sorted_avg,1)

        if strcmpi(elec_include{chan},'yes') && ~ismember(chan,dataBase(subj).metadata(run).cc_stimsets(stimp,:))
            if ~isnan(ccep.n1_peak_sample(chan,stimp)) % NaN is als hij niet gedetecteerd is, dus kijkt alleen naar de gedetecteerde als er hiervoor al gedetecteerd is
                H = plot_ccep(dataBase,subj,run,cfg,tt,chan,stimp);

                perc = n / size(ccep.n1_peak_amplitude(:),1) *100;

                x = input(sprintf('%2.1f %% --- stimpair = %s-%s chan = %s --- Is this a CCEP (y/n) and correct N1 detection (y/yn)? [y/yn/n]: ',...
                    perc,dataBase(subj).metadata(run).cc_stimchans{stimp,:},dataBase(subj).metadata(run).ch{chan}),'s');

                if strcmp(x,'yn') || strcmp(x,'cyn') % so it is an ER, but N1 is not correctly detected
                    currkey = 0;
                    while ~strcmp(currkey,'c') % currkey == 0 %
                        w = waitforbuttonpress;

                        if w == 0     %0 for mouse click, 1 for button press
                            cp = get(gca,'CurrentPoint');
                            % find sample number closest to the selected point
                            [~,sampnum] = min(abs(tt-cp(1,1)));

                            % find nearby peak
                            [~,locs] = findpeaks(-1*squeeze(dataBase(subj).metadata(run).cc_epoch_sorted_reref_avg(chan,stimp,...
                                sampnum-round(0.01*fs):sampnum+round(0.01*fs))),'NPeaks',1,'SortStr','descend');

                            % find x-position of nearby peak
                            locsamp = sampnum-round(0.01*fs)+locs-1;


                            hold on
                            plot(tt(locsamp),dataBase(subj).metadata(run).cc_epoch_sorted_reref_avg(chan,stimp,locsamp),'o','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',6); drawnow;
                            hold off

                            fprintf('sample: %1.4f, amplitude: %2.1f, Press "c" when correct \n',tt(locsamp),dataBase(subj).metadata(run).cc_epoch_sorted_avg(chan,stimp,locsamp))

                        elseif w==1
                            currkey = get(gcf,'CurrentCharacter');
                        end
                    end

                    ccep.n1_peak_amplitude(chan,stimp) = dataBase(subj).metadata(run).cc_epoch_sorted_avg(chan,stimp,locsamp);
                    ccep.n1_peak_sample(chan,stimp) = locsamp;
                    ccep.checked(chan,stimp) = 1;

                elseif strcmp(x,'y') || strcmp(x,'cy')
                    ccep.checked(chan,stimp) = 1 ;
                else
                    ccep.checked(chan,stimp) = 0 ;
                end 
            close(H)
            end
            %show the scoring of another scorer
%             ccep_scorer2= load(fullfile(myDataPath.CCEPpath2, ...
%             dataBase(subj).sub_label,dataBase(subj).ses_label, dataBase(subj).metadata(run).run_label,...
%             [dataBase(subj).sub_label, '_', dataBase(subj).ses_label,'_',...
%             dataBase(subj).metadata(run).task_label,'_',...
%             dataBase(subj).metadata(run).run_label,'_N1sChecked.mat']));
%             checked_scorer2=ccep_scorer2.checked(chan,stimp);
%             checked_scorer2_str=strrep(strrep(sprintf('%d ', checked_scorer2), '1', 'Yes'), '0', 'No');
%             uiwait(msgbox({'CCEP according to other scorer?:' ,'',sprintf('%s',checked_scorer2_str)},'','help'));
            
        else
            %             ccep.checked(chan,stimp) = NaN;
        end
        n=n+1;
    end
    % save also till which stimpair visual N1s are checked.
    ccep.checkUntilStimp = stimp;

    filename = [dataBase(subj).sub_label,'_',dataBase(subj).ses_label,'_',dataBase(subj).metadata(run).task_label,'_',dataBase(subj).metadata(run).run_label,'_N1sChecked.mat'];

    % save file during scoring in case of error
    save(fullfile(filefolder,filename),'-struct','ccep');

end

dataBase(subj).metadata(run).ccep = ccep;
fprintf('Visual rating of %s completed\n',dataBase(subj).sub_label)
end
