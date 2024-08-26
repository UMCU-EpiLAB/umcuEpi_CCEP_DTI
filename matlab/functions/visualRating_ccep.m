
function dataBase = visualRating_ccep(dataBase,myDataPath,subj,run,cfg,endstimp)

filefolder = fullfile(myDataPath.input_dev,dataBase(subj).sub_label);
if ~exist(filefolder,'dir')
    mkdir(filefolder)
end

ccep = dataBase(subj).metadata(run).ccep;

fs = dataBase(subj).metadata(run).ccep_header.Fs; 
tt = -cfg.epoch_prestim+1/fs:1/fs:cfg.epoch_length-cfg.epoch_prestim;

% if no n1 peaks are checked yet: pre-allocation
if  sum(strcmp(fieldnames(dataBase(subj).metadata(run).ccep),'checked'))==0
    ccep.checked = zeros(size(dataBase(subj).metadata(run).ccep.ch,1),size(dataBase(subj).metadata(run).ccep.cc_stimchans,1));
end

n = numel(ccep.checked(:,1:endstimp))+1;
for nStimp = endstimp+1:size(dataBase(subj).metadata(run).ccep.cc_stimchans,1)
    for nChan = 1:size(dataBase(subj).metadata(run).ccep.ch,1)

            if ~isnan(ccep.n1_peak_sample(nChan,nStimp)) 
                H = plot_ccep(dataBase,subj,run,tt,nChan,nStimp);

                perc = n / size(ccep.n1_peak_amplitude(:),1) *100;

                x = input(sprintf('%2.1f %% --- stimpair = %s-%s chan = %s --- Is this a CCEP (y/n) and correct N1 detection (y/yn)? [y/yn/n]: ',...
                    perc,dataBase(subj).metadata(run).ccep.cc_stimchans{nStimp,:},dataBase(subj).metadata(run).ccep.ch{nChan}),'s');

                if strcmp(x,'yn') || strcmp(x,'cyn') % so it is an ER, but N1 is not correctly detected
                    currkey = 0;
                    while ~strcmp(currkey,'c') % currkey == 0 %
                        w = waitforbuttonpress;

                        if w == 0     %0 for mouse click, 1 for button press
                            cp = get(gca,'CurrentPoint');
                            % find sample number closest to the selected point
                            [~,sampnum] = min(abs(tt-cp(1,1)));

                            % find nearby peak
                            [~,locs] = findpeaks(-1*squeeze(dataBase(subj).metadata(run).cc_epoch_sorted_reref_avg(nChan,nStimp,...
                                sampnum-round(0.01*fs):sampnum+round(0.01*fs))),'NPeaks',1,'SortStr','descend');

                            % find x-position of nearby peak
                            locsamp = sampnum-round(0.01*fs)+locs-1;


                            hold on
                            plot(tt(locsamp),dataBase(subj).metadata(run).cc_epoch_sorted_reref_avg(nChan,nStimp,locsamp),'o','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',6); drawnow;
                            hold off

                            fprintf('sample: %1.4f, amplitude: %2.1f, Press "c" when correct \n',tt(locsamp),dataBase(subj).metadata(run).cc_epoch_sorted_avg(nChan,nStimp,locsamp))

                        elseif w==1
                            currkey = get(gcf,'CurrentCharacter');
                        end
                    end

                    ccep.n1_peak_amplitude(nChan,nStimp) = dataBase(subj).metadata(run).cc_epoch_sorted_avg(nChan,nStimp,locsamp);
                    ccep.n1_peak_sample(nChan,nStimp) = locsamp;
                    ccep.checked(nChan,nStimp) = 1;

                elseif strcmp(x,'y') || strcmp(x,'cy')
                    ccep.checked(nChan,nStimp) = 1 ;
                else
                    ccep.checked(nChan,nStimp) = 0 ;
                end 
            close(H)
            end
            
        n=n+1;
    end
    % save also till which stimpair visual N1s are checked.
    ccep.checkUntilStimp = nStimp;

    filename = [dataBase(subj).sub_label,'_',dataBase(subj).ses_label,'_',dataBase(subj).metadata(run).task_label,'_',dataBase(subj).metadata(run).run_label,'_N1sChecked.mat'];

    % save file during scoring in case of error
    save(fullfile(filefolder,filename),'-struct','ccep');

end

dataBase(subj).ccep(run).ccep = ccep;
fprintf('Visual rating of %s completed\n',dataBase(subj).sub_label)
end
