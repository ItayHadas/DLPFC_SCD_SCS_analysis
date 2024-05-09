clear all
cd 'D:\OneDrive - University of California, San Diego Health\DATA\Theta-Burst-Dose_tremblay\REST'
load('D:\MATLAB\LAB_MatlabScripts\Chanlocs\Chanlocs_64Ch-EasyCap_for_BrainAmp_AFz_FCz.mat');
chan_fold='D:\OneDrive - University of California, San Diego Health\DATA\Theta-Burst-Dose_tremblay\channel_locations\';

[fileNames, pathName]=Z_getSetsFileNames('vhdr');
for i=158: size(fileNames,1) % PIPELINE LOOP

    fileName=fileNames{i};
    EEG  = pop_loadbv(pathName , fileName );
    disp([ num2str(i) ' out of ' num2str(size(fileNames,1)) '   -   ' fileName ])

    EEG.filename=fileName([1:3,9:end]);
    EEG.subject=EEG.filename(2:3);
    EEG.setname=EEG.filename(1:end-5);
    if contains(fileName,'open')
        EEG.session='open';
    elseif contains(fileName,'close')
        EEG.session='close';
    end
    if contains(fileName,'pre')
        EEG.condition=[EEG.filename(5) '_pre'];
    elseif contains(fileName,'post')
        EEG.condition=[EEG.filename(5) '_post'];
    end
    EEG.data=double(EEG.data);
    EEG = eeg_checkset( EEG );

    chan=readtable([chan_fold EEG.subject(1:2) '_' num2str(EEG.filename(5)) '.xlsx'],'ReadVariableNames',0);
    chan.Properties.VariableNames{1}='labels';
    remove_chan=chan{strcmp(chan{:,2},'-'),1}';

    chanlocsxx= {};
    for c=1:size(chan,1)
        if sum(ismember(lower(chan.labels(c)),lower({chanlocs62.labels})))==1
            chanlocsxx(c,:)= struct2cell(chanlocs62(ismember(lower({chanlocs62.labels}),lower(chan.labels(c)))));
        else
            chanlocsxx(c,1)=chan.labels(c);
        end

    end
    EEG.chanlocs=cell2struct(chanlocsxx,fieldnames(chanlocs62),2);
    EEG.event=[];  EEG.urevent=[];  EEG.eventdescription=[];  EEG.epoch=[];    EEG.epochdescription=[];
    EEG = pop_select( EEG,'nochannel',remove_chan);
    EEG = eeg_checkset( EEG );
    %  pop_eegplot( EEG, 1, 1, 1);
    EEG = pop_saveset( EEG, 'filename',EEG.filename,'filepath',strrep([pathName 'Chanfixed\'],'\','\\'),'check', 'on','savemode','onefile','version','7.3');
end

