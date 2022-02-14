function data=show_ccep(dataBase,F_stimp,F_chan,subj,run,cfg)
data = zeros(1);
for i = 1:length(F_stimp)
stimp = F_stimp(i);
chan = F_chan(i);
fs = dataBase(subj).metadata(run).ccep_header.Fs;
tt = -cfg.epoch_prestim+1/fs:1/fs:cfg.epoch_length-cfg.epoch_prestim;

H = plot_ccep(dataBase,subj,run,cfg,tt,chan,stimp);

perc = i / length(F_stimp) *100;
x = input(sprintf('%2.1f %% --- stimpair = %s-%s chan = %s --- Is this a CCEP? (y/n): ',...
perc,dataBase(subj).metadata(run).cc_stimchans{stimp,:},dataBase(subj).metadata(run).ch{chan}),'s');
% input if it is indeed a FP/FN or not by scoring the CCEPs
% again

if strcmp(x,'y') 
data(chan,stimp) = 1 ;
else
data(chan,stimp) = 0 ;
end
close(H)
clear('x')
           
end
fprintf('Visual rating of %s completed\n',dataBase(subj).sub_label)