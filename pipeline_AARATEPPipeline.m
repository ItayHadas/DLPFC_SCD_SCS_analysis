% clone AARATEPPioeline, Tesa and EEGLAB repos, and add the directories to matlab paths (change the local dirs below)
% https://github.com/chriscline/AARATEPPipeline
% https://github.com/nigelrogasch/TESA


% matlab -nosplash -nodesktop
% workspace 
clear all
if (ispc)
    sep='\';
    not_sep='/';
    rep_space = ' ';
    GITS='D:\GITs\';
    Path = dir('\\ad.ucsd.edu\ahs\apps\INTERPSYC\DATA\Wellcome_Leap_802232\Neurophysiology_Data\**\*SPD_*.cdt');
    outdir='A:\WorkingSet\WellcomeLeap_TEP';
elseif (ismac || isunix)
    sep='/';
    not_sep='\';
    rep_space = '\ ';
    GITS='/media/ipp/DATA/GITs';
    addpath([GITS sep 'TESA']);
    gitdir=[GITS sep 'ThePipeline'];
    Path = dir('/mnt/INTERPSYC/DATA/Wellcome_Leap_802232/Neurophysiology_Data/**/*SPD_*.cdt');
    outdir='/media/ipp/DATA/EEG_DATA/WL_TEP_9RUNa'; 
    mkdir(outdir)
    eeglabdir='/media/ipp/DATA/Documents/MATLAB/eeglab-2023.1/';
    % brainstormdir='/media/ipp/DATA/Documents/MATLAB/brainstorm3'
end
adjustelec=false;
addpath(gitdir, genpath([gitdir sep 'AARATEPPipeline']));
addpath(eeglabdir); %addpath(brainstormdir); 
addpath([gitdir sep 'addons' sep]); eeglab('nogui')
%%
fileNames={Path.name}';
idx=1:35  %:size(fileNames,1);
for i= idx %%%%%%%%%%%%%%%%% PIPELINE LOOP
    
    fileName=fileNames{i};
    pathName=[Path(i).folder sep];
    disp([ 'Loading ' num2str(find(idx==i)) ' out of ' num2str(length(idx)) '   -   ' fileName '   #####' ...
    '  ##################################@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@' ...
        '  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##################################' ])
    if strcmpi(fileName(:,end-2:end),'set')
        EEG = pop_loadset( [pathName fileName]);
    elseif strcmpi(fileName(:,end-2:end),'cnt')
        EEG = pop_loadcnt([pathName fileName] , 'dataformat', 'int32');
    elseif strcmpi(fileName(:,end-2:end),'cdt')
        EEG = loadcurry([pathName fileName], 'CurryLocations', 'True');
    elseif strcmpi(fileName(:,end-3:end),'vhdr')
        EEG  = pop_loadbv(pathName , fileName );
    end

    if size(EEG.data,2)/EEG.srate<15 
        continue; %if less then 15 sec skip dataset
    end 

    EEG.filename=fileName;
    EEG.subject=EEG.filename(1:6);
    if strcmpi(fileName(:,end-3:end),'vhdr')
        EEG.setname=EEG.filename(1:end-5);
    else
        EEG.setname=EEG.filename(1:end-4);
    end
    % pop_eegplot( EEG, 1, 1, 1);
    if contains(EEG.setname,'_S1','IgnoreCase',true)...
        | contains(EEG.setname,'_BL','IgnoreCase',true)...
            | contains(EEG.setname,'_pre','IgnoreCase',true)
        EEG.condition='pre';
        EEG.session='1';
    elseif contains(EEG.setname,'_S3','IgnoreCase',true)...
        | contains(EEG.setname,'_post','IgnoreCase',true)
        EEG.condition='post';
        EEG.session='2';
    end
    try
        EEG.chanlocs(elecName(EEG,{'afp1'})).labels='Fp1'; 
        EEG.chanlocs(elecName(EEG,{'afp2'})).labels='Fp2'; 
    end
    if adjustelec
    EEG = pop_chanedit(EEG, 'lookup',[eeglabdir sep 'plugins' sep 'dipfit' sep ...
    'standard_BEM' sep 'elec' sep 'standard_1005.elc']);
    end
    EEG.data=double(EEG.data);
    EEG = eeg_checkset( EEG );
    disp(['Dataset is loaded ' num2str(i) '/' num2str(size(fileNames,1))])
    disp(['remove un-needed channels  - dataset ' num2str(i) ...
        '/' num2str(size(fileNames,1))]);
    remove_channels= { 'F11' 'F12' 'FT11' 'FT12' 'CB1' 'CB2'...
        'lz' 'VEO' 'HEO' 'VEOG' 'HEOG' 'EKG' 'EMG' 'EOG' 'HL 1' 'HL 2' 'Trigger'};
    %figure; topoplot([],EEG.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo)
    chanlocs62=EEG.chanlocs;
    EEG = pop_select( EEG,'nochannel',remove_channels);
    evname= '128'; newtrigs=1; method=1; DoublePulseINT=0; % paired-pulse interval, Zero (0) in case of single pulse
    [impotant, ~]=elecName(EEG,{'f5' 'f3' 'f1' 'fc5' ...
        'fc3' 'fc1' 'fcz' 'fc2' 'c1' 'cz' 'c2' 'c4' 'c6' 'cp5' 'cp3'});
    if (contains(EEG.setname,'_RS'))
        [impotant, ~]=elecName(EEG,{'f6' 'f4' 'f2' 'fc6'...
            'fc4' 'fc1' 'fcz' 'fc2' 'c1' 'cz' 'c2' 'c4' 'c6' 'cp5' 'cp3'});
    end
    % pop_eegplot( EEG, 1, 1, 1);
    EEG = IDpulse(EEG,method,DoublePulseINT,impotant,newtrigs,evname);
    if size(EEG.urevent,2)<8 
        continue; %if less then 8 pulses skip dataset
    end 
    EEG = eeg_checkset( EEG );
    EEG_mat = c_TMSEEG_Preprocess_AARATEPPipeline(EEG,...
        'epochTimespan', [-1.5 1.5],...
        'outputDir', [outdir sep],...
        'pulseEvent',evname,...
        'outputFilePrefix',[EEG.setname '_Clean'],...
        'ICAType','fastica',...
        'doDecayRemovalPerTrial',true);   %'ICAType','amica','picard'+ 'notch_freq',[48 52],... %for ANU Australia
    if strcmpi(fileName(:,end-3:end),'vhdr')
        EEG_mat.setname=[fileName(1:end-5) '_AARATEPPipeline'];
    else
        EEG_mat.setname=[fileName(1:end-4) '_AARATEPPipeline'];
    end
    EEG_mat.filename=[EEG_mat.setname '.set'];
    EEG_mat.filepath=[outdir sep];
    EEG_mat.datfile='';
    EEG_mat.report.rank_4={rank(double(EEG_mat.data(:,:)),1e-4)};
    EEG_mat.report.rank_3={rank(double(EEG_mat.data(:,:)),1e-3)};
    EEG_mat = pop_saveset( EEG_mat, 'filename',EEG_mat.filename,...
        'filepath',EEG_mat.filepath,'check', 'on','savemode','onefile','version','7.3');
    if i==min(idx) % EEG_mat.report.orgsamprate=EEG.srate EEG_mat.report.newsamprate=EEG.srate;
        report=EEG_mat.report;
    else
        report=[report ; EEG_mat.report];
    end
    clear EEG_mat EEG
    writetable(report,[outdir sep 'AARATEPPipeline-REPORT.csv'])
end