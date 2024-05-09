%% GRAND Average 
%chan_interp='off';
%load('D:\Google_Drive\MATLAB\LAB_MatlabScripts\chanlocs68.mat');
%load('D:\OneDrive\MATLAB\LAB_MatlabScripts\Chanlocs_64Ch-EasyCap_for_BrainAmp_AFz_FCz.mat'); chanlocs=chanlocs62;
load('D:\OneDrive\DATA\PAS_DATA_temerty\chanlocs62.mat'); %chanlocs=chanlocs62;
chan_interp='on'
Subloc=[4:6];
Z2_grand_average('Pre_0min_PAS_temerty',Subloc,chan_interp,Chanlocs62,[0]) 

figure; plot(EEG.times,squeeze(EEG.data(elecName(EEG,{'f3'}),:,:))); ylim([-40 40]); xlim([-300 300])


[a loc]=max(EEG.data(17,5,:))
EEG.data(:,:,8)=[]


%% SCD SCS comp
%load('D:\WORKINGSET_D\Trial_compare_SCD_SCS_SGC_hipp\conditions_run_loop.mat')
%load('D:\OneDrive\UCSD - Daskalakis\PROJECTS\TremblaySara_thetaBurst\localized_Theta_dose\condition_Theta_dose.mat')
% [fileNames, pathName]=Z_getSetsFileNames('set'); D:\OneDrive\DATA\Theta-Burst-Dose_tremblay\SP_preprocessed
% D:\OneDrive\DATA\Theta-Burst-Dose_tremblay\SP_preprocessed
% conditions{3,6}=fileNames

brainsto_data=['D:\Brainstorm_db\'];
%load('D:\DATA\Trial_compare_SCD_SCS_SGC_hipp\brainstormchannel2.mat')
%load('E:\Google_Drive\MATLAB\LAB_MatlabScripts\chanlocs62.mat');
%load('D:\OneDrive\MATLAB\LAB_MatlabScripts\chanlocations_for_brainstorm.mat');
%load('D:\OneDrive\DATA\Theta-Burst-Dose_tremblay\channel_locations\Chanlocs_64Ch-EasyCap_for_BrainAmp_AFz_FCz.mat'); %ST iTBS dose chanlocs
%chanlocs62=chanlocs;
Subloc=[4:6];
baseline=[-0.95, -0.35];
sigtime=[0.015, 0.615];
AtlasScts={'Destrieux', {'G_front_middle L', 'G_front_middle R', 'G_subcallosal L', 'G_subcallosal R'}};
%   clearvars -except i fileNames pathName Subloc iProtocol showgraphics ProtocolName cond chanlocs62 Channel brainsto_data conditions

if ~brainstorm('status')
    brainstorm nogui
end

%starts brainstrom
for uuu=1:size(conditions,2)
    %uuu=6
    cond=conditions{1,uuu};
    ProtocolName=['TBS_dose_resp_', conditions{1,uuu}];
    pathName='D:\OneDrive\DATA\Theta-Burst-Dose_tremblay\SP_preprocessed\';
    fileNames=conditions{3,uuu};
    
    % Get the protocol index
    iProtocol = bst_get('Protocol', ProtocolName);
    if isempty(iProtocol)
        warning(['Creating new protocol: ' ProtocolName]);
        % Create new protocol
        gui_brainstorm('CreateProtocol', ProtocolName, 1, 0);
    else
        gui_brainstorm('SetCurrentProtocol', iProtocol);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%
    %% Looping trough subjects
    
    
    %[fileNames, pathName]=Z_getSetsFileNames('set');
    % Z2_grand_average('PRE_MST', [1:3])
    % % % EEG = pop_loadset( ['U:\WORKINGset\MST\PRE\PRE_MST_GrandAverage.set']);
    
    for i=1: size(fileNames,1)
        %% Loading the Dataset
        if size(fileNames, 1)==1
            fileName=fileNames{i,1}'
        else
            fileName=fileNames{i,1}
        end
        EEG = pop_loadset( [pathName fileName]);
        EEG = eeg_checkset( EEG );
        %EEG=pop_chanedit(EEG, 'lookup','E:\\Google_Drive\\MATLAB\\EEGLAB\\plugins\\dipfit2.3\\standard_BEM\\elec\\standard_1005.elc');
        remove_channels= {'M1' 'M2' 'VEO' 'HEO' 'EKG' 'EMG' 'HL 1' 'HL 2' 'Trigger'};
        EEG = pop_select( EEG,'nochannel',remove_channels);
        EEG = pop_interp(EEG, chanlocs62, 'spherical');
        EEG = pop_saveset( EEG, [pathName fileName]);
        clear EEG;
        
        %% brainstorm files and folders
        SubjectName =  {fileName(Subloc)}; %{  'NewSubject'};
        SubjectNamestr=fileName(Subloc);
        RawFiles = {[pathName fileName]};
        % '005/MDD_PRE_SP/data_Neuroscan_EEG_data_resampled_pruned_with_ICA_pruned_with_ICA_(_average_171128_1725.mat'};
        brnstmPATH= [brainsto_data ProtocolName '\data\' SubjectNamestr '\' cond '\'];
        
        sFiles =[];
        %% Process: Import MEG/EEG: Existing epochs
        sFiles = bst_process('CallProcess', 'process_import_data_epoch', sFiles, [], ...
            'subjectname',  SubjectName{1}, ...
            'condition',    cond, ...
            'datafile',     {{RawFiles{1}}, 'EEG-EEGLAB'}, ...
            'iepochs',      [], ...
            'eventtypes',   '', ...
            'createcond',   0, ...
            'channelalign', 1, ...
            'usectfcomp',   1, ...
            'usessp',       1, ...
            'freq',         [], ...
            'baseline',     []);
        
        Org_sFiles=sFiles;
        %% Process: Select data files in: Subject01/*/Avg: deviant
        sFiles = bst_process('CallProcess', 'process_select_files_data', sFiles, [], ...
            'subjectname',   SubjectName{1}, ...
            'condition',     cond, ...
            'tag',           '', ...  %Avg: deviant
            'includebad',    0, ...
            'includeintra',  0, ...
            'includecommon', 0);
        %% Process: Set channel file
        sFiles = bst_process('CallProcess', 'process_channel_addloc', sFiles, [], ...
            'channelfile', {'', ''}, ...
            'usedefault',  65, ...  % ICBM152: BrainProducts EasyCap 64
            'fixunits',    1, ...
            'vox2ras',     1);
        
        % 'channelfile',  {'', ''}, ...
        % 'usedefault',   65, ...  % ICBM152: BrainProducts EasyCap 64
        % 'usedefault',   79, ...  % ICBM152: Neuroscan Quik-cap 64
        %save([brainsto_data ProtocolName '\data\' SubjectNamestr '\' cond '\' 'channel.mat'], 'Channel');
        
        %% Process: Extract time: [15ms,615ms]
        
        sFiles= Org_sFiles;
        sFiles = bst_process('CallProcess', 'process_extract_time', sFiles, [], ...
            'timewindow', sigtime, ...
            'overwrite',  0);
        after150_800=  sFiles ;
        
        %% Process: Select data files in: Subject01/*/Avg: deviant
        sFiles = bst_process('CallProcess', 'process_select_files_data', sFiles, [], ...
            'subjectname',   SubjectName{1}, ...
            'condition',     cond, ...
            'tag',           '', ...  %Avg: deviant
            'includebad',    0, ...
            'includeintra',  0, ...
            'includecommon', 0);
        %% Process: Extract time: [-850ms,-215ms]
        sFiles= Org_sFiles;
        sFiles = bst_process('CallProcess', 'process_extract_time', sFiles, [], ...
            'timewindow', baseline, ...
            'overwrite',  0);
        
        before900_115=  sFiles ;
        
        %% Process: Select data files in: Subject01/*/Avg: deviant
        sFiles = bst_process('CallProcess', 'process_select_files_data', sFiles, [], ...
            'subjectname',   SubjectName{1}, ...
            'condition',     '', ...
            'tag',           '', ...  %Avg: deviant
            'includebad',    0, ...
            'includeintra',  0, ...
            'includecommon', 0);
        %% Process: Compute head model
        
        
        sFiles = bst_process('CallProcess', 'process_headmodel', sFiles, [], ...
            'Comment',     '', ...
            'sourcespace', 1, ...  % Cortex surface
            'meg',         1, ...  %
            'eeg',         3, ...  % OpenMEEG BEM
            'ecog',        1, ...  %
            'seeg',        1, ...  %
            'openmeeg',    struct(...
            'BemSelect',    [1, 1, 1], ...
            'BemCond',      [1, 0.0125, 1], ...
            'BemNames',     {{'Scalp', 'Skull', 'Brain'}}, ...
            'BemFiles',     {{}}, ...
            'isAdjoint',    0, ...
            'isAdaptative', 1, ...
            'isSplit',      0, ...
            'SplitLength',  4000), ...
            'channelfile', '');
        
        %% Process: Compute covariance (noise or data) - for baseline
        
        %     a=dir(brnstmPATH);
        %     b='data_';
        %     f={a(find(contains({a.name},b) & ~contains({a.name},'(_average_'))).name};
        %     lnk= repmat({[SubjectNamestr '/' cond '/'  ]},size(f,2),1);
        %     sFiles = join(cat(  2, lnk, f' ),'')';
        
        sFiles = bst_process('CallProcess', 'process_noisecov', {Org_sFiles.FileName}, [], ...
            'baseline',       [-1.9, -0.005], ...
            'datatimewindow', [0.017, 0.8], ...
            'sensortypes',    'EEG', ...
            'target',         1, ...  % Noise covariance     (covariance over baseline time window)
            'dcoffset',       1, ...  % Block by block, to avoid effects of slow shifts in data
            'identity',       0, ...
            'copycond',       0, ...
            'copysubj',       0, ...
            'copymatch',      0, ...
            'replacefile',    1);  % Replace
        
        % Process: Compute covariance (noise or data)
        sFiles = bst_process('CallProcess', 'process_noisecov', {Org_sFiles.FileName}, [], ...
            'baseline',       [-1.9, -0.005], ...
            'datatimewindow', [0.017, 0.8], ...
            'sensortypes',    '', ...
            'target',         2, ...  % Data covariance      (covariance over data time window)
            'dcoffset',       1, ...  % Block by block, to avoid effects of slow shifts in data
            'identity',       0, ...
            'copycond',       0, ...
            'copysubj',       0, ...
            'copymatch',      0, ...
            'replacefile',    1);  % Keep
        %% Process: Compute sources
        %
        
        srcFiles1 = bst_process('CallProcess', 'process_inverse_2018', {Org_sFiles.FileName}, [], ...
            'output',  2, ...  % Kernel only: one per file
            'inverse', struct(...
            'Comment',        'PNAI: EEG', ...
            'InverseMethod',  'lcmv', ...
            'InverseMeasure', 'nai', ...
            'SourceOrient',   {{'fixed'}}, ...
            'Loose',          0.2, ...
            'UseDepth',       1, ...
            'WeightExp',      0.5, ...
            'WeightLimit',    10, ...
            'NoiseMethod',    'median', ...
            'NoiseReg',       0.1, ...
            'SnrMethod',      'rms', ...
            'SnrRms',         1e-06, ...
            'SnrFixed',       3, ...
            'ComputeKernel',  1, ...
            'DataTypes',      {{'EEG'}}));
        
        srcFiles_baseline = bst_process('CallProcess', 'process_inverse_2018', {before900_115.FileName}, [], ...
            'output',  2, ...  % Kernel only: one per file
            'inverse', struct(...
            'Comment',        'PNAI: EEG', ...
            'InverseMethod',  'lcmv', ...
            'InverseMeasure', 'nai', ...
            'SourceOrient',   {{'fixed'}}, ...
            'Loose',          0.2, ...
            'UseDepth',       1, ...
            'WeightExp',      0.5, ...
            'WeightLimit',    10, ...
            'NoiseMethod',    'median', ...
            'NoiseReg',       0.1, ...
            'SnrMethod',      'rms', ...
            'SnrRms',         1e-06, ...
            'SnrFixed',       3, ...
            'ComputeKernel',  1, ...
            'DataTypes',      {{'EEG'}}));
        
        srcFiles_signal = bst_process('CallProcess', 'process_inverse_2018', {after150_800.FileName}, [], ...
            'output',  2, ...  % Kernel only: one per file
            'inverse', struct(...
            'Comment',        'PNAI: EEG', ...
            'InverseMethod',  'lcmv', ...
            'InverseMeasure', 'nai', ...
            'SourceOrient',   {{'fixed'}}, ...
            'Loose',          0.2, ...
            'UseDepth',       1, ...
            'WeightExp',      0.5, ...
            'WeightLimit',    10, ...
            'NoiseMethod',    'median', ...
            'NoiseReg',       0.1, ...
            'SnrMethod',      'rms', ...
            'SnrRms',         1e-06, ...
            'SnrFixed',       3, ...
            'ComputeKernel',  1, ...
            'DataTypes',      {{'EEG'}}));
        
        
        %% Process: Scouts time series: G_front_middle L G_front_middle R G_subcallosal L G_subcallosal R
        actfiles = bst_process('CallProcess', 'process_extract_scout', {srcFiles1.FileName}, [], ...
            'timewindow',     sigtime, ...
            'scouts',         AtlasScts, ...
            'scoutfunc',      5, ...  % All
            'isflip',         1, ...
            'isnorm',         1, ...
            'concatenate',    0, ...
            'save',           1, ...
            'addrowcomment',  1, ...
            'addfilecomment', 1);
        
        % Process: Average: By subject
        ScoutActFile = bst_process('CallProcess', 'process_average', {actfiles.FileName}, [], ...
            'avgtype',       2, ...  % By subject
            'avg_func',      1, ...  % Arithmetic average:  mean(x)
            'weighted',      0, ...
            'keepevents',    0, ...
            'matchrows',     1, ...
            'iszerobad',     1);
        % Process: Extract values: [-2000ms,2000ms] 1 scouts abs
        % actfiles = bst_process('CallProcess', 'process_extract_values', {srcFiles1.FileName}, [], ...
        %     'timewindow', sigtime, ...
        %     'scoutsel',   AtlasScts, ...
        %     'scoutfunc',  5, ...  % All
        %     'isnorm',     1, ...
        %     'avgtime',    0, ...
        %     'dim',        1, ...  % Concatenate signals (dimension 1)
        %     'Comment',    '');
        %
        % % Process: Average: By subject
        % sFiles = bst_process('CallProcess', 'process_average', {actfiles.FileName}, [], ...
        %     'avgtype',       2, ...  % By subject
        %     'avg_func',      1, ...  % Arithmetic average:  mean(x)
        %     'weighted',      0, ...
        %     'keepevents',    0, ...
        %     'matchrows',     1, ...
        %     'iszerobad',     1);
        
        
        
        %% Process: t-test baseline [-2000ms,2000ms]          H0:(X=Baseline), H1:(X<>Baseline)
        
        Basettest = bst_process('CallProcess', 'process_test_baseline', {srcFiles1.FileName}, [], ...
            'baseline',      [-1.5, -0.03], ...
            'timewindow',    [-1.5, 0.615], ...
            'scoutsel',      AtlasScts, ...
            'scoutfunc',     5, ...  % ALL
            'isnorm',        0, ...
            'avgtime',       0, ...
            'test_type',     'ttest_baseline', ...  % Student's t-test vs baseline        X~N(m,v)Y = mean_trials(X)        Y~N(m,v)t = (Y - mean_time(Y(baseline)) / std_time(Y(baseline)))df = Nbaseline - 1 = length(baseline) - 1
            'tail',          'one+');  % 'two'=Two-tailed 'one+'= One-tailed (+)
        %% Permutation test = baseline<>activation
        % {AtlasScts{1}, AtlasScts{2}(:,3:4)}
        perm_2tail_SGC = bst_process('CallProcess', 'process_test_permutation2p', {srcFiles_signal.FileName}, {srcFiles_baseline.FileName}, ...
            'timewindow',     [], ...
            'scoutsel',       {AtlasScts{1}, AtlasScts{2}(:,3:4)}, ... {'Yeo 7 Networks', {'Visual_L', 'Visual_R', 'Somatomotor_L', 'Somatomotor_R', 'Dorsal_Attention_L', 'Dorsal_Attention_R', 'Ventral_Attention_L', 'Ventral_Attention_R', 'Limbic_L', 'Limbic_R', 'FrontoParietal_L', 'FrontoParietal_R', 'DefaultModeNetwork_L', 'DefaultModeNetwork_R'}}, ...
            'scoutfunc',      5, ...  % ALL
            'isnorm',         0, ...
            'avgtime',        0, ...
            'iszerobad',      1, ...
            'Comment',        '', ...
            'test_type',      'ttest_unequal', ...  % Student's t-test   (unequal variance) t = (mean(A)-mean(B)) / sqrt(var(A)/nA + var(B)/nB)
            'randomizations', 1000, ...
            'tail',           'two');  % Two-tailed
        
        perm_1tail_SGC = bst_process('CallProcess', 'process_test_permutation2p', {srcFiles_signal.FileName}, {srcFiles_baseline.FileName}, ...
            'timewindow',     [], ...
            'scoutsel',       {AtlasScts{1}, AtlasScts{2}(:,3:4)}, ... {'Yeo 7 Networks', {'Visual_L', 'Visual_R', 'Somatomotor_L', 'Somatomotor_R', 'Dorsal_Attention_L', 'Dorsal_Attention_R', 'Ventral_Attention_L', 'Ventral_Attention_R', 'Limbic_L', 'Limbic_R', 'FrontoParietal_L', 'FrontoParietal_R', 'DefaultModeNetwork_L', 'DefaultModeNetwork_R'}}, ...
            'scoutfunc',      5, ...  % all
            'isnorm',         0, ...
            'avgtime',        0, ...
            'iszerobad',      1, ...
            'Comment',        '', ...
            'test_type',      'ttest_unequal', ...  % Student's t-test   (unequal variance) t = (mean(A)-mean(B)) / sqrt(var(A)/nA + var(B)/nB)
            'randomizations', 1000, ...
            'tail',           'one+');  % one-tailed
        
        perm_2tail_L_DLPFC = bst_process('CallProcess', 'process_test_permutation2p', {srcFiles_signal.FileName}, {srcFiles_baseline.FileName}, ...
            'timewindow',     [], ...
            'scoutsel',       {AtlasScts{1}, {'G_front_middle L'}}, ... {'Yeo 7 Networks', {'Visual_L', 'Visual_R', 'Somatomotor_L', 'Somatomotor_R', 'Dorsal_Attention_L', 'Dorsal_Attention_R', 'Ventral_Attention_L', 'Ventral_Attention_R', 'Limbic_L', 'Limbic_R', 'FrontoParietal_L', 'FrontoParietal_R', 'DefaultModeNetwork_L', 'DefaultModeNetwork_R'}}, ...
            'scoutfunc',      5, ...  % ALL
            'isnorm',         0, ...
            'avgtime',        0, ...
            'iszerobad',      1, ...
            'Comment',        '', ...
            'test_type',      'ttest_unequal', ...  % Student's t-test   (unequal variance) t = (mean(A)-mean(B)) / sqrt(var(A)/nA + var(B)/nB)
            'randomizations', 1000, ...
            'tail',           'two');  % Two-tailed
        
        perm_1tail_L_DLPFC = bst_process('CallProcess', 'process_test_permutation2p', {srcFiles_signal.FileName}, {srcFiles_baseline.FileName}, ...
            'timewindow',     [], ...
            'scoutsel',       {AtlasScts{1}, {'G_front_middle L'}}, ... {'Yeo 7 Networks', {'Visual_L', 'Visual_R', 'Somatomotor_L', 'Somatomotor_R', 'Dorsal_Attention_L', 'Dorsal_Attention_R', 'Ventral_Attention_L', 'Ventral_Attention_R', 'Limbic_L', 'Limbic_R', 'FrontoParietal_L', 'FrontoParietal_R', 'DefaultModeNetwork_L', 'DefaultModeNetwork_R'}}, ...
            'scoutfunc',      5, ...  % all
            'isnorm',         0, ...
            'avgtime',        0, ...
            'iszerobad',      1, ...
            'Comment',        '', ...
            'test_type',      'ttest_unequal', ...  % Student's t-test   (unequal variance) t = (mean(A)-mean(B)) / sqrt(var(A)/nA + var(B)/nB)
            'randomizations', 1000, ...
            'tail',           'one+');  % one-tailed
        
        perm_2tail_R_DLPFC = bst_process('CallProcess', 'process_test_permutation2p', {srcFiles_signal.FileName}, {srcFiles_baseline.FileName}, ...
            'timewindow',     [], ...
            'scoutsel',       {AtlasScts{1}, {'G_front_middle R'}}, ... {'Yeo 7 Networks', {'Visual_L', 'Visual_R', 'Somatomotor_L', 'Somatomotor_R', 'Dorsal_Attention_L', 'Dorsal_Attention_R', 'Ventral_Attention_L', 'Ventral_Attention_R', 'Limbic_L', 'Limbic_R', 'FrontoParietal_L', 'FrontoParietal_R', 'DefaultModeNetwork_L', 'DefaultModeNetwork_R'}}, ...
            'scoutfunc',      5, ...  % ALL
            'isnorm',         0, ...
            'avgtime',        0, ...
            'iszerobad',      1, ...
            'Comment',        '', ...
            'test_type',      'ttest_unequal', ...  % Student's t-test   (unequal variance) t = (mean(A)-mean(B)) / sqrt(var(A)/nA + var(B)/nB)
            'randomizations', 1000, ...
            'tail',           'two');  % Two-tailed
        
        perm_1tail_R_DLPFC = bst_process('CallProcess', 'process_test_permutation2p', {srcFiles_signal.FileName}, {srcFiles_baseline.FileName}, ...
            'timewindow',     [], ...
            'scoutsel',       {AtlasScts{1}, {'G_front_middle R'}}, ... {'Yeo 7 Networks', {'Visual_L', 'Visual_R', 'Somatomotor_L', 'Somatomotor_R', 'Dorsal_Attention_L', 'Dorsal_Attention_R', 'Ventral_Attention_L', 'Ventral_Attention_R', 'Limbic_L', 'Limbic_R', 'FrontoParietal_L', 'FrontoParietal_R', 'DefaultModeNetwork_L', 'DefaultModeNetwork_R'}}, ...
            'scoutfunc',      5, ...  % all
            'isnorm',         0, ...
            'avgtime',        0, ...
            'iszerobad',      1, ...
            'Comment',        '', ...
            'test_type',      'ttest_unequal', ...  % Student's t-test   (unequal variance) t = (mean(A)-mean(B)) / sqrt(var(A)/nA + var(B)/nB)
            'randomizations', 1000, ...
            'tail',           'one+');  % one-tailed
        %% Loading brainstorm vars for SCD/SCS
        % significance map comparing same sizes vector of baseline to a vector after
        % the stimulation
        % source_stats = load('E:\brainstorm_db\MDD_PRE_SP\data\034\MDD_PRE_SP\presults_no_171206_1243.mat')% significance mask
        %  brnstmPATH='E:\brainstorm_db2\MDD_PRE_MST_SP_DLPFC_DEEP\data\194\MST_re-run_test3\'
        
        
        %      clear a b f c h source_stats;
        %      a=dir(brnstmPATH);
        %      b='results_';
        %      f=a(contains({a.name},b)).name;
        %      source_stats = load( [brnstmPATH f]);
        
        % source_stats =load('E:\brainstorm_db\MDD_PRE_MST_SP_DLPFC_rerun-test4\data\049\MST_re-run_test4\presults_no_190709_1528.mat')
        
        % TEP wave averaged wave - the same size as the significanse map
        % example file:      timedata_avg = load('E:\brainstorm_db\MDD_PRE_SP\data\034\MDD_PRE_SP\data_Merged_datasets_resampled_pruned_with_ICA_pruned_with_ICA_(_average_171129_1247_time.mat')
        %     b= '_average_';   c= '_time';
        %     f1=find(contains({a.name},b)); f2 = find(contains({a.name},c)); f=intersect(f1,f2);
        %     %f1=contains({a.name},b);
        %     %f1=a(f1).name;
        %     timedata_avg = load( [brnstmPATH f1]);
        % timedata_avg =load('E:\brainstorm_db\MDD_PRE_MST_SP_DLPFC_rerun-test4\data\049\MST_re-run_test4\data_ECT049_PRE_TMS_DLPFC_SP_Epo_Dec_interp_ICA1_EXP_NOfreq_pass1_55_ICA2_CLEEG_(_average_190709_1425.mat')
        clear b f loc timedata_avg;
        a=dir(brnstmPATH);
        %b= '_average_';  c= '_time';
        %f1=find(contains({a.name},b)); %f2 = find(contains({a.name},c)); f=intersect(f1,f2);
        %f=a(f1).name;
        
        
        loc=strfind(ScoutActFile.FileName,'/');
        f=ScoutActFile.FileName(loc(end)+1:end);
        timedata_avg = load( [brnstmPATH f]);
        clear b f1 f loc head;
        
        %head=load([brnstmPATH 'headmodel_mix_eeg_3sphereberg.mat']);
        % brainstorm headmodel
        head=load([brnstmPATH 'headmodel_surf_openmeeg.mat']); % headmodel_surf_eeg_3sphereberg
        %head=load([brnstmPATH 'headmodel_mix_openmeeg.mat']);
        
        % head=load('E:\brainstorm_db\MDD_PRE_MST_SP_DLPFC_rerun-test4\data\049\MST_re-run_test4\headmodel_mix_openmeeg.mat')
        
        sensor=load([brnstmPATH 'channel.mat']);
        %sensor=load('E:\brainstorm_db\MDD_PRE_MST_SP_DLPFC_rerun-test4\data\049\MST_re-run_test4\channel.mat')
        
        % brainstorm cortex anatomy
        % [brainsto_data ProtocolName '\anat\' SubjectNamestr '\' cond '\' ]
        cortex=load([brainsto_data ProtocolName '\anat\@default_subject\tess_cortex_pial_low.mat']);
        % brainsto_data=['E:\brainstorm_db2\'];
        % ProtocolName = 'MDD_PRE_MST_SP_DLPFC_DEEP';
        %cortex=load([brainsto_data ProtocolName '\anat\@default_subject\tess_cortex_mixed.mat']);
        % cortex=load('E:\brainstorm_db\MDD_PRE_MST_SP_DLPFC_rerun-test4\anat\@default_subject\tess_cortex_mixed.mat');
        
        % link|009/MST_Pre_DEEP/results_sLORETA_EEG_KERNEL_190910_1318.mat
        %|009/MST_Pre_DEEP/data_Merged_datasets_pruned_with_ICA_pruned_with_ICA_(_average_190910_1317_time.mat
        % brainstorm sLORETA
        
        %     b='results_';%sLORETA_EEG_
        %     f=a(contains({a.name},b)).name;
        %     inversion= load( [brnstmPATH f]);
        
        
        % inversion= load('E:\brainstorm_db\MDD_PRE_MST_SP_DLPFC_rerun-test4\data\049\MST_re-run_test4\results_sLORETA_EEG_190709_1431.mat')
        %% Permutation test variable for saving
        
        ee=strfind(perm_2tail_SGC.FileName,'/');
        f=a(contains({a.name},perm_2tail_SGC.FileName(ee(end)+1:end))).name;
        pmat= load( [brnstmPATH f]);
        Permutation_test.SGC.scouts2=pmat.Description;
        Permutation_test.SGC.tail2.pmap=pmat.pmap;
        Permutation_test.SGC.tail2.time=pmat.Time;
        Permutation_test.SGC.tail2.tmap=pmat.tmap;
        
        ee=strfind(perm_1tail_SGC.FileName,'/');
        f=a(contains({a.name},perm_1tail_SGC.FileName(ee(end)+1:end))).name;
        pmat= load( [brnstmPATH f]);
        Permutation_test.SGC.scouts1=pmat.Description;
        Permutation_test.SGC.tail1.pmap=pmat.pmap;
        Permutation_test.SGC.tail1.time=pmat.Time;
        Permutation_test.SGC.tail1.tmap=pmat.tmap;
        
        ee=strfind(perm_1tail_L_DLPFC.FileName,'/');
        f=a(contains({a.name},perm_1tail_L_DLPFC.FileName(ee(end)+1:end))).name;
        pmat= load( [brnstmPATH f]);
        Permutation_test.L_DLPFC.scouts1=pmat.Description;
        Permutation_test.L_DLPFC.tail1.pmap=pmat.pmap;
        Permutation_test.L_DLPFC.tail1.time=pmat.Time;
        Permutation_test.L_DLPFC.tail1.tmap=pmat.tmap;
        
        ee=strfind(perm_2tail_L_DLPFC.FileName,'/');
        f=a(contains({a.name},perm_2tail_L_DLPFC.FileName(ee(end)+1:end))).name;
        pmat= load( [brnstmPATH f]);
        Permutation_test.L_DLPFC.scouts2=pmat.Description;
        Permutation_test.L_DLPFC.tail2.pmap=pmat.pmap;
        Permutation_test.L_DLPFC.tail2.time=pmat.Time;
        Permutation_test.L_DLPFC.tail2.tmap=pmat.tmap;
        
        ee=strfind(perm_2tail_R_DLPFC.FileName,'/');
        f=a(contains({a.name},perm_2tail_R_DLPFC.FileName(ee(end)+1:end))).name;
        pmat= load( [brnstmPATH f]);
        Permutation_test.R_DLPFC.scouts2=pmat.Description;
        Permutation_test.R_DLPFC.tail2.pmap=pmat.pmap;
        Permutation_test.R_DLPFC.tail2.time=pmat.Time;
        Permutation_test.R_DLPFC.tail2.tmap=pmat.tmap;
        
        ee=strfind(perm_1tail_R_DLPFC.FileName,'/');
        f=a(contains({a.name},perm_1tail_R_DLPFC.FileName(ee(end)+1:end))).name;
        pmat= load( [brnstmPATH f]);
        Permutation_test.R_DLPFC.scouts1=pmat.Description;
        Permutation_test.R_DLPFC.tail1.pmap=pmat.pmap;
        Permutation_test.R_DLPFC.tail1.time=pmat.Time;
        Permutation_test.R_DLPFC.tail1.tmap=pmat.tmap;
        
        ee=strfind(Basettest.FileName,'/');
        f=a(contains({a.name},Basettest.FileName(ee(end)+1:end))).name;
        pmat= load( [brnstmPATH f]);
        Permutation_test.SGC_DLPFC.scoutstt=pmat.Description;
        Permutation_test.SGC_DLPFC.tailtt.pmap=pmat.pmap;
        Permutation_test.SGC_DLPFC.tailtt.time=pmat.Time;
        Permutation_test.SGC_DLPFC.tailtt.tmap=pmat.tmap;
        
        clear pmat f ee
        %% preparing anatomy data and variables
        %
        %     % Calculating distances between stimulation site and sources
        %     % make sure you delete blank channels
        chanlocs = sensor.Channel;
        %
        %
        %sourcelocs = inversion.GridAtlas.Vert2Grid'*head.GridLoc;
        sourcelocs = head.GridLoc;
        %
        %
        total_sources = size(sourcelocs,1);
        %
        %     % chan 7 = F5; chan 26 = C3
        %     % chan_stim = 7; %Stimulation channel
        %     % chan_probe = 7; %electrode of interest
        %
        
        
        chan_stim = find(strcmpi({sensor.Channel.Name},'F3')); %7; %Stimulation channel
        %chan_probe = find(strcmpi({sensor.Channel.Name},'F3')); %7; %electrode of interest
        
        %     for chan = 1:size(chanlocs,2) %extract the coordinates from the channel number
        %         x_chan(chan) = chanlocs(1,chan).Loc(1);
        %         y_chan(chan) = chanlocs(1,chan).Loc(2);
        %         z_chan(chan) = chanlocs(1,chan).Loc(3);
        %     end
        x_stim = chanlocs(1,chan_stim).Loc(1);
        y_stim = chanlocs(1,chan_stim).Loc(2);
        z_stim = chanlocs(1,chan_stim).Loc(3);
        %     x_stim = x_chan(chan_stim);
        %     y_stim = y_chan(chan_stim);
        %     z_stim = z_chan(chan_stim);
        
        %     x_probe = x_chan(chan_probe);
        %     y_probe = y_chan(chan_probe);
        %     z_probe = z_chan(chan_probe);
        %
        %verthemis = ones(total_sources,1);
        %GridAtlas.Grid2Source*GridLoc
        %dist=[];
        dist_stim=[];
        for source_num = 1:total_sources
            x_source = sourcelocs(source_num,1);
            y_source = sourcelocs(source_num,2);
            z_source = sourcelocs(source_num,3);
            
            dist_stim(source_num) = sqrt((x_source-x_stim)^2+(y_source-y_stim)^2+(z_source-z_stim)^2);
            %dist(source_num) = sqrt((x_source-x_probe)^2+(y_source-y_probe)^2+(z_source-z_probe)^2);
            
            %     if source_num > total_sources/2
            %         verthemis(source_num) = 2;
            %     end
        end
        
        %     %dist_stim=dist_stim*inversion.GridAtlas.Vert2Grid;
        %     %dist=dist*inversion.GridAtlas.Vert2Grid;
        %
        %     %inversion.GridLoc(inversion.GridAtlas.Vert2Grid')
        %
        %     %   dist_stim=dist_stim* inversion.GridAtlas.Grid2Source';
        %     %dist=dist* inversion.GridAtlas.Grid2Source';
        %% Adding source values from defined anatomical regions
        % categorize the sources according to atlass ROI
        
        %     verthemis = []; % Left and right hemisphere labelling
        %     inversion.GridAtlas.Vert2Grid
        %     for k = 1:size(head.GridLoc,1)
        %         if find(cortex.Atlas(1, 3).Scouts(1, 1).Vertices == k) % left cortex
        %             verthemis(k) = 1;
        %         else
        %             verthemis(k) = 2;
        %         end
        %     end
        
        % Atlas variables
        % 3 = Destrieux, 148 scouts
        %   = Desikan-Killiany, 68 scouts
        
        % atlas = 3;
        %vertlabel = [];
        
        %     for s = 1:size(cortex.Atlas(1,atlas).Scouts,2) %loop through all the scouts
        %         scoutnames{s} = cortex.Atlas(1, atlas).Scouts(1, s).Label;
        %     end
        %     scoutnames = scoutnames';
        atlas='Destrieux'; atlas=find(contains({cortex.Atlas.Name},atlas));
        scoutatlas={cortex.Atlas(1,atlas ).Scouts(1, :).Label}';
        %atlas2='Structures'; atlas2=find(contains({cortex.Atlas.Name},atlas2));
        %scout_Structures={cortex.Atlas(1, atlas2).Scouts(1, :).Label}';
        %scoutnames=[scoutatlas; scout_Structures];
        
        % K=15627; s=3;
        %     vertlabel=[]; vertlabelnames={};
        %     total_sources=size(source_stats.ImagingKernel,1);
        %     for k = 1:total_sources %loop through all sources
        %         for s = 1:size(scoutatlas,1) %loop through all the scouts
        %             if find(cortex.Atlas(1, atlas).Scouts(1, s).Vertices == k)
        %                 vertlabel(k) = find(strcmp(scoutatlas,cortex.Atlas(1, atlas).Scouts(1, s).Label));
        %                 vertlabelnames{k} = cortex.Atlas(1, atlas).Scouts(1, s).Label;
        %             end
        %         end
        %     end
        % AtlasScts{:,2}
        %     vertlabel=[]; vertlabelnames={}; total_sources=[];
        %     %clear vertlabel vertlabelnames total_sources
        %     total_sources=size(source_stats.ImagingKernel,1);
        %     for k = 1:total_sources %loop through all sources
        %         for s = 1:size(AtlasScts{:,2},2) %loop through all the scouts
        %             if any(cortex.Atlas(1, atlas).Scouts(1, find(contains({cortex.Atlas(1, atlas).Scouts(1, :).Label},AtlasScts{:,2}(s)))).Vertices == k)
        %                 vertlabel(k) = find(strcmp(scoutatlas,AtlasScts{:,2}(s)));
        %                 vertlabelnames{k} =AtlasScts{:,2}{s};
        %             else
        %                 vertlabel(k) =0;
        %                 vertlabelnames{k} ='';
        %             end
        %         end
        %     end
        
        
        total_sources=size(cortex.Vertices,1);
        vertlabel=zeros(1,total_sources); vertlabelnames=cell(1,total_sources);
        for s = 1:size(AtlasScts{:,2},2)
            vertlabel(cortex.Atlas(1, atlas).Scouts(1, find(contains({cortex.Atlas(1, atlas).Scouts(1, :).Label},AtlasScts{:,2}(s)))).Vertices)=find(strcmp(scoutatlas,AtlasScts{:,2}(s)));
            vertlabelnames(cortex.Atlas(1, atlas).Scouts(1, find(contains({cortex.Atlas(1, atlas).Scouts(1, :).Label},AtlasScts{:,2}(s)))).Vertices)={AtlasScts{:,2}{s}};
        end
        %find(strcmp(vertlabelnames,'G_subcallosal L'))
        
        
        %find(contains(vertlabelnames,'G_front')) %{'G_front_middle L'}    {'G_front_middle R'}    {'G_subcallosal L'}    {'G_subcallosal R'}
        %find(strcmp(vertlabelnames,'G_subcallosal L'))
        %      vertlabelnames=vertlabelnames( vertlabel*inversion.GridAtlas.Grid2Source');
        %         vertlabel=vertlabel* inversion.GridAtlas.Grid2Source';
        %     cortex.Atlas(1,3).Scouts
        %     for k = 1:size(head.GridLoc,1) %loop through all sources
        %         for s = 1:size(cortex.Atlas(1,atlas).Scouts,2) %loop through all the scouts
        %             if find(cortex.Atlas(1, atlas).Scouts(1, s).Vertices == k)
        %                 vertlabel(k) = s;
        %                 vertlabelnames{k} = cortex.Atlas(1, atlas).Scouts(1, s).Label;
        %             end
        %         end
        %     end
        %% Using imported statistical values
        % check the number of channel in each file - and delete the extra
        % channels
        % sigtime baseline
        [~, zero] = min(abs(Permutation_test.SGC.tail2.time-(sigtime(1))));
        [~, eigh] = min(abs(Permutation_test.SGC.tail2.time-(sigtime(2))));
        pmap2_SGC=Permutation_test.SGC.tail2.pmap(:,zero:eigh);
        %
        [~, zero] = min(abs(Permutation_test.SGC.tail1.time-(sigtime(1))));
        [~, eigh] = min(abs(Permutation_test.SGC.tail1.time-(sigtime(2))));
        pmap1_SGC=Permutation_test.SGC.tail1.pmap(:,zero:eigh);
        %
        [~, zero] = min(abs(Permutation_test.L_DLPFC.tail2.time-(sigtime(1))));
        [~, eigh] = min(abs(Permutation_test.L_DLPFC.tail2.time-(sigtime(2))));
        pmap2_L_DLPFC=Permutation_test.L_DLPFC.tail2.pmap(:,zero:eigh);
        %
        [~, zero] = min(abs(Permutation_test.L_DLPFC.tail1.time-(sigtime(1))));
        [~, eigh] = min(abs(Permutation_test.L_DLPFC.tail1.time-(sigtime(2))));
        pmap1_L_DLPFC=Permutation_test.L_DLPFC.tail1.pmap(:,zero:eigh);
        %
        [~, zero] = min(abs(Permutation_test.R_DLPFC.tail2.time-(sigtime(1))));
        [~, eigh] = min(abs(Permutation_test.R_DLPFC.tail2.time-(sigtime(2))));
        pmap2_R_DLPFC=Permutation_test.R_DLPFC.tail2.pmap(:,zero:eigh);
        %
        [~, zero] = min(abs(Permutation_test.R_DLPFC.tail1.time-(sigtime(1))));
        [~, eigh] = min(abs(Permutation_test.R_DLPFC.tail1.time-(sigtime(2))));
        pmap1_R_DLPFC=Permutation_test.R_DLPFC.tail1.pmap(:,zero:eigh);
        %
        [~, zero] = min(abs(Permutation_test.SGC_DLPFC.tailtt.time-(sigtime(1))));
        [~, eigh] = min(abs(Permutation_test.SGC_DLPFC.tailtt.time-(sigtime(2))));
        pmaptt_SGC_DLPFC=Permutation_test.SGC_DLPFC.tailtt.pmap(:,zero:eigh);
        
        
        [~, SS2_SGC]=fdr(pmap2_SGC,0.05);
        [~, SS1_SGC]=fdr(pmap1_SGC,0.05);
        [~, SS2_L_DLPFC]=fdr(pmap2_L_DLPFC,0.05);
        [~, SS1_L_DLPFC]=fdr(pmap1_L_DLPFC,0.05);
        [~, SS2_R_DLPFC]=fdr(pmap2_R_DLPFC,0.05);
        [~, SS1_R_DLPFC]=fdr(pmap1_R_DLPFC,0.05);
        [~, SStt_SGC_DLPFC]=fdr(pmaptt_SGC_DLPFC,0.05);
        %spy(sparse(SS2_SGC));
        %pmap2=inversion.GridAtlas.Vert2Grid'*source_stats.pmap;
        
        %     SS = zeros(size(source_stats.pmap)); % significance matrix( mask)
        %     threshold = 0.05/(size(source_stats.pmap,1));
        %     % size(pmap2)
        %     for k1 = 1:size(source_stats.pmap,1)
        %         for k2 = 1:size(source_stats.pmap,2)
        %             if source_stats.pmap(k1,k2)<threshold
        %                 SS(k1,k2) = 1;
        %             end
        %         end
        %     end
        
        %     SS2 = zeros(size(pmap2)); % significance matrix( mask)
        %     threshold = 0.05./(size(pmap2,1).*size(pmap2,2));
        %     %threshold = 0.05/15000 1/15000
        %     % size(pmap2)
        %     for k1 = 1:size(pmap2,1)
        %         for k2 = 1:size(pmap2,2)
        %             if pmap2(k1,k2)<threshold
        %                 SS(k1,k2) = 1;
        %             end
        %         end
        %     end
        %sum(SS,'all')
        %SS= inversion.GridAtlas.Vert2Grid'*SS;
        
        %J = inversion.ImagingKernel*timedata_avg.F; %previous version of cortical s-loreta
        %multiplying the Sloretta kernel with the TEP wave averaged wave
        %J= inversion.ImageGridAmp;
        %sLoreta_act=(inversion.GridAtlas.Vert2Grid'*inversion.ImagingKernel);
        %J = sLoreta_act*timedata_avg.F;
        
        timevec = timedata_avg.Time;
        [~, zero] = min(abs(timevec-(sigtime(1))));
        [~, eigh] = min(abs(timevec-(sigtime(2))));
        %J = timedata_avg.Value(:,[zero:eigh]);
        %J = inversion.ImagingKernel*timedata_avg.F(inversion.GoodChannel,[zero:eigh]);
        %timevec = source_stats.Time;
        
        %timevec = timedata_avg.Time;
        %timevec = source_stats.Time;
        %[~, zero] = min(abs(timevec-(0.015)));
        %[~, eigh] = min(abs(timevec-(0.615)));
        %J = J(:,[zero:eigh]);
        timevec=timevec([zero:eigh]);
        
        %q=[]; J_atlas = [];  %_atlas = [];
        % checksum = 0;
        %         for kk = 1:length(scoutnames)
        %for kk = [63 64]%1:length(scoutnames)
        % kk=151
        %     q=[]; J_atlas = []; %SCD_atlas = [];
        %     dist_vox = []; SCD_atlas = [];
        %     checksum = 0;
        %     %contains(scoutatlas,Permutation_test.YeoNET.scouts2)
        %     % AtlasScts{2}
        %     for kk = 1:length(scoutatlas)% [3:8 31:34 55:58 63 64]%1:length(scoutnames)
        %
        %         q=find(vertlabel==kk);
        %         % J_atlas(kk,:) = sum(abs(J(q,:)*10),1);
        %         %SCD_atlas(kk,:) = sum(SCD(q,:)*10,1);
        %         %SCS_atlas(kk,:) = sum(SCS(q,:),1);
        %         dist_vox(kk,:) = sum(dist_stim(:,q),2);
        %         checksum = checksum + length(q);
        %     end
        
        J_SGC=sum(abs(timedata_avg.Value(find(contains(timedata_avg.Description,'G_subcallosal')),[zero:eigh])),1);
        J_L_DLPFC=sum(abs(timedata_avg.Value(find(contains(timedata_avg.Description,'G_front_middle L')),[zero:eigh])),1);
        J_R_DLPFC=sum(abs(timedata_avg.Value(find(contains(timedata_avg.Description,'G_front_middle R')),[zero:eigh])),1);
        % plot(timevec*1000,J_L_DLPFC), hold on, plot(timevec*1000,J_R_DLPFC) ,hold off
        SCD1_SGC=sum(abs(timedata_avg.Value(find(contains(timedata_avg.Description,'G_subcallosal')),[zero:eigh])).*SS1_SGC,1);
        SCD2_SGC=sum(abs(timedata_avg.Value(find(contains(timedata_avg.Description,'G_subcallosal')),[zero:eigh])).*SS2_SGC,1);
        SCDtt_SGC=sum(abs(timedata_avg.Value(find(contains(timedata_avg.Description,'G_subcallosal')),[zero:eigh])).*...
            SStt_SGC_DLPFC(find(contains(Permutation_test.SGC_DLPFC.scoutstt,'G_subcallosal')),:),1);
        % figure;plot(timevec*1000,SCD1_SGC), hold on, plot(timevec*1000,SCD2_SGC) , plot(timevec*1000,SCDtt_SGC) ,hold off
        SCD1_L_DLPFC=sum(abs(timedata_avg.Value(find(contains(timedata_avg.Description,'G_front_middle L')),[zero:eigh])).*SS1_L_DLPFC,1);
        SCD2_L_DLPFC=sum(abs(timedata_avg.Value(find(contains(timedata_avg.Description,'G_front_middle L')),[zero:eigh])).*SS2_L_DLPFC,1);
        SCDtt_L_DLPFC=sum(abs(timedata_avg.Value(find(contains(timedata_avg.Description,'G_front_middle L')),[zero:eigh])).*...
            SStt_SGC_DLPFC(find(contains(Permutation_test.SGC_DLPFC.scoutstt,'G_front_middle L')),:),1);
        %figure; plot(timevec*1000,SCD1_L_DLPFC), hold on, plot(timevec*1000,SCD2_L_DLPFC) , plot(timevec*1000,SCDtt_L_DLPFC) ,hold off
        SCD1_R_DLPFC=sum(abs(timedata_avg.Value(find(contains(timedata_avg.Description,'G_front_middle R')),[zero:eigh])).*SS1_R_DLPFC,1);
        SCD2_R_DLPFC=sum(abs(timedata_avg.Value(find(contains(timedata_avg.Description,'G_front_middle R')),[zero:eigh])).*SS2_R_DLPFC,1);
        SCDtt_R_DLPFC=sum(abs(timedata_avg.Value(find(contains(timedata_avg.Description,'G_front_middle R')),[zero:eigh])).*...
            SStt_SGC_DLPFC(find(contains(Permutation_test.SGC_DLPFC.scoutstt,'G_front_middle R')),:),1);
        %figure; plot(timevec*1000,SCD1_R_DLPFC), hold on, plot(timevec*1000,SCD2_R_DLPFC) , plot(timevec*1000,SCDtt_R_DLPFC) ,hold off
        SCS1_SGC=sum(repmat(dist_stim([cortex.Atlas(2).Scouts(find(contains({cortex.Atlas(2).Scouts.Label},'G_subcallosal'))).Vertices])',1,size(SS1_SGC,2))...
            .*SS1_SGC,1);
        SCS2_SGC=sum(repmat(dist_stim([cortex.Atlas(2).Scouts(find(contains({cortex.Atlas(2).Scouts.Label},'G_subcallosal'))).Vertices])',1,size(SS2_SGC,2))...
            .*SS2_SGC,1);
        SCStt_SGC=sum(repmat(dist_stim([cortex.Atlas(2).Scouts(find(contains({cortex.Atlas(2).Scouts.Label},'G_subcallosal'))).Vertices])',1,size(SS1_SGC,2))...
            .*SStt_SGC_DLPFC(find(contains(Permutation_test.SGC_DLPFC.scoutstt,'G_subcallosal')),:),1);
        %figure; plot(timevec*1000,SCS1_SGC), hold on, plot(timevec*1000,SCS2_SGC) , plot(timevec*1000,SCStt_SGC) ,hold off
        SCS1_L_DLPFC=sum(repmat(dist_stim([cortex.Atlas(2).Scouts(find(contains({cortex.Atlas(2).Scouts.Label},'G_front_middle L'))).Vertices])',1,size(SS1_L_DLPFC,2))...
            .*SS1_L_DLPFC,1);
        SCS2_L_DLPFC=sum(repmat(dist_stim([cortex.Atlas(2).Scouts(find(contains({cortex.Atlas(2).Scouts.Label},'G_front_middle L'))).Vertices])',1,size(SS1_L_DLPFC,2))...
            .*SS2_L_DLPFC,1);
        SCStt_L_DLPFC=sum(repmat(dist_stim([cortex.Atlas(2).Scouts(find(contains({cortex.Atlas(2).Scouts.Label},'G_front_middle L'))).Vertices])',1,size(SS1_L_DLPFC,2))...
            .*SStt_SGC_DLPFC(find(contains(Permutation_test.SGC_DLPFC.scoutstt,'G_front_middle L')),:),1);
        %figure; plot(timevec*1000,SCS1_L_DLPFC), hold on, plot(timevec*1000,SCS2_L_DLPFC) , plot(timevec*1000,SCStt_L_DLPFC) ,hold off
        SCS1_R_DLPFC=sum(repmat(dist_stim([cortex.Atlas(2).Scouts(find(contains({cortex.Atlas(2).Scouts.Label},'G_front_middle R'))).Vertices])',1,size(SS1_R_DLPFC,2))...
            .*SS1_R_DLPFC,1);
        SCS2_R_DLPFC=sum(repmat(dist_stim([cortex.Atlas(2).Scouts(find(contains({cortex.Atlas(2).Scouts.Label},'G_front_middle R'))).Vertices])',1,size(SS1_R_DLPFC,2))...
            .*SS2_R_DLPFC,1);
        SCStt_R_DLPFC=sum(repmat(dist_stim([cortex.Atlas(2).Scouts(find(contains({cortex.Atlas(2).Scouts.Label},'G_front_middle R'))).Vertices])',1,size(SS1_R_DLPFC,2))...
            .*SStt_SGC_DLPFC(find(contains(Permutation_test.SGC_DLPFC.scoutstt,'G_front_middle R')),:),1);
        %figure; plot(timevec*1000,SCS1_R_DLPFC), hold on, plot(timevec*1000,SCS2_R_DLPFC) , plot(timevec*1000,SCStt_R_DLPFC) ,hold off
        
        %% Calculate SCD and SCS measures
        
        % SCD1_atlas=abs(J_atlas).*SS1*10;% ?A/mm2
        % SCD2_atlas=abs(J_atlas).*SS2*10;% ?A/mm2
        % SCDtt_atlas=abs(J_atlas).*SStt*10;% ?A/mm2
        
        % SCD(inversion.GridAtlas.Grid2Source')
        % SCS=SS.*repmat(dist,[1,size(SS,2)]);% mm
        %GridAtlas.Grid2Source*GridLoc
        %SCS=SS.*repmat(dist_stim',[1,size(SS,2)]);% mm
        %     dist_atlas=repmat(dist_vox,[1 size(SS1,2)]);
        %     SCS1_atlas=dist_atlas.*SS1;
        %     SCS2_atlas=dist_atlas.*SS2;
        %     SCStt_atlas=dist_atlas.*SStt;
        %% Calculating SCD and SCS summed over defined atlas regions
        
        % Atlas 2**************Destrieux
        % label 3, 4 = Rostral Anterior Cingulate
        % label 5, 6 = Dorsal Anterior Cingulate
        % label 7, 8 = Middle Posterior Cingulate
        % label 31,32 = Frontal Middle Gyrus
        % label 33,34 = Superior Frontal Gyrus
        % label 55,56 = Postcentral
        % label 57,58 = Precentral
        % label 63,64 = G_Subcallosal
        
        % Atlas 3**************
        % label 17 = inferiortemporal L
        % label 18 = inferiortemporal R
        % label 49 = precentral L
        % label 55 = rostralmiddlefrontal L
        % label 57 = superiorfrontal L
        % label 53 = rACC L
        % label 54 = rACC R
        
        %         'select the scout areas (e.g. BA regions) of interest'
        %         BA=input(['Scout area: ']);
        %         hemisphere=input(['hemisphere: specify either L or R: '],'s');
        
        %     q=[]; J_atlas = []; SCD_atlas = []; SCS_atlas = [];
        %     checksum = 0;
        %     %         for kk = 1:length(scoutnames)
        %     %for kk = [63 64]%1:length(scoutnames)
        %     % kk=151
        %     for kk = 1:length(scoutnames)% [3:8 31:34 55:58 63 64]%1:length(scoutnames)
        %
        %         q=find(vertlabel==kk);
        %         J_atlas(kk,:) = sum(abs(J(q,:)*10),1);
        %         SCD_atlas(kk,:) = sum(SCD(q,:)*10,1);
        %         SCS_atlas(kk,:) = sum(SCS(q,:),1);
        %
        %         checksum = checksum + length(q);
        %     end
        %% Saves the DATA
        
        % clearvars -except J J_sum J_atlas SCD SCD_atlas SCD_sum SCS SCS_atlas
        % SCS_sum timevec dist_stim vertlabel vertlabelnames SS cortex %timedata_avg
        save (['A:\WorkingSet\ThetaBurst_DOSE_ST\SCD-SCS\' SubjectNamestr '_' cond '_SGC-DLPFC_J_SCD_SCS_' date '.mat'],'J_SGC', 'J_L_DLPFC', 'J_R_DLPFC',...
            'SCD1_SGC', 'SCD2_SGC', 'SCDtt_SGC', 'SCD1_L_DLPFC', 'SCD2_L_DLPFC', 'SCDtt_L_DLPFC', 'SCD1_R_DLPFC',...
            'SCD2_R_DLPFC', 'SCDtt_R_DLPFC','SCS1_SGC', 'SCS2_SGC', 'SCStt_SGC', 'SCS1_L_DLPFC', 'SCS2_L_DLPFC', 'SCStt_L_DLPFC', 'SCS1_R_DLPFC',...
            'SCS2_R_DLPFC', 'SCStt_R_DLPFC', 'SS2_SGC', 'SS1_SGC', 'SS2_L_DLPFC', 'SS1_L_DLPFC', 'SS2_R_DLPFC', 'SS1_R_DLPFC', 'timevec');
        %save ([pathName SubjectNamestr 'J_SCD_SCS_DEEP_atlas_' cond '_' date '.mat'], 'J', 'J_atlas', 'SCD', 'SCD_atlas', 'SCS', 'SCS_atlas', 'timevec', 'sLoreta_act' , 'dist_stim', 'vertlabel', 'vertlabelnames', 'SS', 'pmap2', 'scoutnames', 'Permutation_test');
        %     save ([pathName SubjectNamestr '_' cond 'SGC-DLPFC_J_SCD_SCS_' date '.mat'],  'J_atlas', 'SCD1_atlas',...
        %         'SCD2_atlas', 'SCDtt_atlas', 'SCS1_atlas', 'SCS2_atlas', 'SCStt_atlas' , 'scoutatlas', 'SS1',...
        %         'SS2', 'SStt', 'pmap2', 'pmap1', 'pmaptt');
        save (['A:\WorkingSet\ThetaBurst_DOSE_ST\SCD-SCS\' SubjectNamestr '_' cond '_SGC-DLPFC_scouts_activation_' date '.mat'],'timedata_avg');
        % 'J_sum','SCD_sum', 'SCS_sum',
        clearvars -except i fileNames pathName Subloc iProtocol showgraphics ProtocolName cond chanlocs62 Channel brainsto_data conditions baseline sigtime AtlasScts
    end
    
end

brainstorm stop
%   clearvars -except i fileNames pathName Subloc iProtocol showgraphics ProtocolName cond chanlocs62 Channel brainsto_data conditions

close all
fclose('all')
%clear all
%% stats - repeated ANOVA
aa=readtable('A:\WorkingSet\ThetaBurst_DOSE_ST\test.xlsx')
aaa=aa(:,2:7);
WithinStructure = table(categorical([1 2 1 2 1 2]'),categorical([1 1 2 2 3 3]'),'VariableNames',{'PrePost','Dose'});
rm = fitrm(aaa, 'pre6,post6,pre12,post12,pre18,post18~1','WithinDesign',WithinStructure);
ranovatable = ranova(rm,'WithinModel','PrePost*Dose');
figure; plotprofile(rm,'Dose'); hold off
post = multcompare(rm,'Dose')
figure; boxplot(aaa{:,:}); hold off
%% Theta dose GRAND averaging
cd 'A:\WorkingSet\ThetaBurst_DOSE_ST\SCD-SCS'
[fileNames, pathName]=Z_getSetsFileNames('mat'); % A:\WorkingSet\ThetaBurst_DOSE_ST\SCD-SCS
datatable=table;
for uuu=1:size(fileNames,1)
    data=load([pathName fileNames{uuu}])
    dat2=struct2cell(data);
    datatable{uuu,1}={fileNames{uuu}(1:2)};
    floc=strfind(fileNames{uuu},'_');
    datatable{uuu,2}={fileNames{uuu}(floc(1)+1:floc(2)-1)};
    datatable{uuu,3}={fileNames{uuu}(floc(2)+1:floc(3)-1)};
    datatable{uuu,4:24}=dat2(1:21)';
    if uuu==1
        fns = fieldnames(data);
        datatable.Properties.VariableNames(1)={'Subject'};
        datatable.Properties.VariableNames(2)={'Stage'};
        datatable.Properties.VariableNames(3)={'Dose'};
        datatable.Properties.VariableNames(4:24)=fns(1:21)
    end
end

datatable=sortrows(datatable,{'Stage' 'Dose' 'Subject'});
save(['A:\WorkingSet\ThetaBurst_DOSE_ST\SCD-SCS\GRANDtable_J_SCD_SCS_SGC_DLPFC_' date '.mat'],"datatable")
%% timecourse figures load files
% load('D:\WORKINGset_D\Trial_compare_SCD_SCS_Yeo_Networks\3trials_stats_struct.mat')
% load('D:\WORKINGset_D\Trial_compare_SCD_SCS_Yeo_Networks\conditions_run_loop_YeoNets.mat'); comp=comp_nets;
% load('D:\WORKINGset_D\Trial_compare_SCD_SCS_Yeo_Networks\brainstorm_YeoNets_time_course_vec.mat')
% load('D:\WORKINGset_D\Trial_compare_SCD_SCS_Yeo_Networks\YeoNet_scout_atlas.mat')
% datatable.filt=true(size(datatable,1),1); datatable=movevars(datatable, 'filt', 'Before', 'J_L_DLPFC');
%  cd 'D:\OneDrive\DATA\Theta-Burst-Dose_tremblay\SCD-SCS';
%  load('timevec.mat'); load('GRANDtable_J_SCD_SCS_SGC_DLPFC_01-Jul-2021.mat')
cd 'A:\WorkingSet\ThetaBurst_DOSE_ST\SCD-SCS'
load('A:\WorkingSet\ThetaBurst_DOSE_ST\SCD-SCS\timevec.mat')
load('A:\WorkingSet\ThetaBurst_DOSE_ST\SCD-SCS\GRANDtable_J_SCD_SCS_SGC_DLPFC_01-Jul-2021.mat')
%% timecourse figures

pretab=contains(datatable{:,2},'pre','IgnoreCase',true);
posttab=contains(datatable{:,2},'post','IgnoreCase',true);
dose600=contains(datatable{:,3},'600','IgnoreCase',true);
dose1200=contains(datatable{:,3},'1200','IgnoreCase',true);
dose1800=contains(datatable{:,3},'1800','IgnoreCase',true);
%pretab posttab dose600 dose1200 dose1800 jtab scdtab scstab tailtt tail1 tail2 LDLPFC RDLPFC SGC ROI
jtab=contains(datatable.Properties.VariableNames,'J_','IgnoreCase',true);
scdtab=contains(datatable.Properties.VariableNames,'SCD','IgnoreCase',true);
scstab=contains(datatable.Properties.VariableNames,'SCS','IgnoreCase',true);
tailtt=contains(datatable.Properties.VariableNames,'tt_','IgnoreCase',true);
tail1=contains(datatable.Properties.VariableNames,'1_','IgnoreCase',true);
tail2=contains(datatable.Properties.VariableNames,'2_','IgnoreCase',true);
%l
tail1=contains(datatable.Properties.VariableNames,'1_','IgnoreCase',true);
tail2=contains(datatable.Properties.VariableNames,'2_','IgnoreCase',true);
tailtt=contains(datatable.Properties.VariableNames,'tt_','IgnoreCase',true);
%
LDLPFC=contains(datatable.Properties.VariableNames,'L_DLPFC','IgnoreCase',true);
RDLPFC=contains(datatable.Properties.VariableNames,'R_DLPFC','IgnoreCase',true);
SGC=contains(datatable.Properties.VariableNames,'SGC','IgnoreCase',true);
%% plotting J SCD SCS
remoutlie=0; smoo_method='movmedian'; smoo_winsize=4; % movmean movmedian % within='o';
sctail=tail1; % tail1 tail2 tailtt
measures=[2];% {'J'} {'SCD'} {'SCS'}
ROIs=[2]; %{'DLPFC_L'} {'DLPFC_R'} {'SGC'}
figplot=1; statstab=1;
stattim=[25 90]; showtim=[0 550]; fillmis=1; fillmiss='movmean';% 'movmean'  movmedian
ROI= {'DLPFC_L' 'DLPFC_R' 'SGC'; LDLPFC RDLPFC SGC};
grouping={}; grouping={'Pre_600' 'Post_600' 'Pre_1200' 'Post_1200' 'Pre_1800' 'Post_1800' ; [pretab & dose600] [posttab & dose600] ...
    [pretab & dose1200] [posttab & dose1200] [pretab & dose1800] [posttab & dose1800]};
mm={'J' 'SCD' 'SCS'; jtab  [scdtab & sctail] [scstab & sctail] ;...
    [strcat('current J (\muA/mm', '^{', '2', '}', ')')]  [strcat('SCD (\muA/mm', '^{', '2', '}', ')')] ['SCS (mm)']}; %
grphtime=timevec.*1000; [~, fix1] = min(abs(grphtime-(showtim(1)))); [~, fix2] = min(abs(grphtime-(showtim(2))));
[~, sta1] = min(abs(grphtime-(stattim(1)))); [~, sta2] = min(abs(grphtime-(stattim(2))));
co=brewermap(4,'RdBu');  c_pre=co(1,:); c_post=co(4,:); colors={c_pre, c_post} ; %RdBu PuOr
clearvars -except datatable fillmis figplot statstab fillmiss stattim ROIs timevec sta1 sta2 fix1 fix2 TimeVec comp2 measures stat_method within pcor remoutlie smoo_method smoo_winsize figg timeind grphtime UU colors co timefunc timevect mm grouping ROI comp grouptab pretab posttab dose600 dose1200 dose1800 jtab scdtab scstab tailtt tail1 tail2 LDLPFC RDLPFC SGC

for r=ROIs %1:size(ROIs,2)
    for j=measures % 1:size(mm(:,measures),2)
        pl=1;
        if figplot==1 %;
            figure('Position', [300 80 700 350]), hold on;
        end
        for i=1:2:size(grouping,2)
            colvec=mm{2,j} & ROI{2,r};
            rowvec1=[grouping{2,i} & datatable.filt]; rowvec2=[grouping{2,i+1} & datatable.filt];
            dat{i}=cell2mat(datatable{rowvec1,colvec}); dat{i+1}=cell2mat(datatable{rowvec2,colvec});
            if remoutlie==1  %%% REMOVE outliers - 2 methods: 1. replace by nan 2. replace by close value
                dat{i}=filloutliers(dat{i},'clip',smoo_method,smoo_winsize,'ThresholdFactor',3); %nearest previous 'movmedian' movmean
                dat{i+1}=filloutliers(dat{i+1},'clip',smoo_method,smoo_winsize,'ThresholdFactor',3); %nearest previous 'movmedian' movmean
            end
            if figplot==1
                subplo{pl}=subplot(2,3,pl); hold on; pl=pl+1;
                f{i}=fill([stattim stattim(2) stattim(1)],[0 0 0 0],'-k','facealpha',.1,'edgecolor','none');
                N{j,r,1}=num2str(size(dat{i},1)); N{j,r,2}=num2str(size(dat{i+1},1));
                avg1=smoothdata(nanmean(dat{i}(:,fix1:fix2),1),'movmean',20); avg2=smoothdata(nanmean(dat{i+1}(:,fix1:fix2),1),'movmean',20);
                stndE1=smoothdata(ste(dat{i}(:,fix1:fix2),1),'movmean',20); stndE2=smoothdata(ste(dat{i+1}(:,fix1:fix2),1),'movmean',20);
                patch([grphtime(fix1:fix2) fliplr(grphtime(fix1:fix2))], [avg1-stndE1 fliplr(avg1+stndE1)],co(4,:),'EdgeColor',co(3,:),'LineWidth',0.05,'FaceVertexAlphaData',0,'FaceAlpha',0.15); %,'EdgeColor',co(1,:)
                patch([grphtime(fix1:fix2) fliplr(grphtime(fix1:fix2))], [avg2-stndE2 fliplr(avg2+stndE2)],co(2,:),'EdgeColor',co(2,:),'LineWidth',0.05,'FaceVertexAlphaData',0,'FaceAlpha',0.15); %,'EdgeColor',co(4,:)
                jplot(j,r,1)=plot(grphtime(fix1:fix2),avg1,'color',[co(4,:) 0.6],'LineWidth',1.8);
                jplot(j,r,2)=plot(grphtime(fix1:fix2),avg2,'color',[co(1,:) 1],'LineWidth',2);
                title(strrep([grouping{1,i}(4:end) '_pulse_' ROI{1,r} ' ' mm{1,j} ],'_',' '))
            end
            clear colvec rowvec avg1 stndE1 avg2 stndE2
        end
        if figplot==1
            linkaxes([subplo{1:3}],'y');
            yy=max(get(gca,'YLim')); subplo{1}.YLim=[0 yy]; subplo{2}.YLim=[0 yy]; subplo{3}.YLim=[0 yy];
            f{1}.YData=[yy yy 0 0]; f{3}.YData=[yy yy 0 0]; f{5}.YData=[yy yy 0 0];
            legend([jplot(j,r,1) jplot(j,r,2)],{['PRE n=' N{j,r,1}] ['POST n=' N{j,r,2}]})
        end
        clear aaa WithinStructure rm ranovatable post
        aaa=cell2mat(cellfun(@(x) mean(x(:,sta1:sta2),2,'omitnan'),dat,'UniformOutput',false)); aaa(aaa==0)=nan;
        if fillmis==1, aaa=fillmissing(aaa,fillmiss,10) ; end %fillmiss='movmedian' 'movmean'
        aaa=array2table(aaa,'VariableNames', grouping(1,:));
        WithinStructure = table(categorical({'Pre';'Post';'Pre';'Post';'Pre';'Post'}),...
            categorical([600;600;1200;1200;1800;1800]),'VariableNames',{'PrePost','Dose'});
        rm = fitrm(aaa, 'Pre_600,Post_600,Pre_1200,Post_1200,Pre_1800,Post_1800~1','WithinDesign',WithinStructure);
        ranovatable = ranova(rm,'WithinModel','PrePost*Dose');
        if figplot==1
            subplo{pl}=subplot(2,9,pl+6:pl+8); hold on; pl=pl+3;
            plotprofile(rm,'PrePost','Group','Dose') ; legend('FontSize',5,'Location','best'); title(strrep([ROI{1,r} ' ' mm{1,j} ],'_',' ')); hold off
            subplo{pl}=subplot(2,9,pl+6:pl+7); hold on; pl=pl+2;
            plotprofile(rm,'Dose'); ylabel('') ; title(strrep([ROI{1,r} ' ' mm{1,j} ],'_',' '));  hold off
            %post = multcompare(rm,'Dose');
            subplo{pl}=subplot(2,9,pl+6:pl+9); hold on; pl=pl+2;
            boxplot(aaa{:,:},'Labels',grouping(1,:)); title(strrep([ROI{1,r} ' ' mm{1,j} ],'_',' ')); hold off
        end
        if statstab==1
            ranovatable(:,[5:end])
            %ufig=uifigure('Position', [21 446 689 211]);
            %tabl=uitable(ufig,'data',ranovatable{:,[5:end]},'Units', 'Normalized','Position',[0, 0, 1, 1],...
            %'ColumnName',ranovatable.Properties.VariableNames(5:end),'RowName',ranovatable.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);%
            %s = uistyle('BackgroundColor',[0.8500 0.3250 0.0980]); addStyle(tabl,s,'column',5)
        end
        %figfig('dose_Sig_current')
        clear subplo
    end
end

%% ISP
remoutlie=0; smoo_method='movmedian'; smoo_winsize=4; % movmean movmedian % within='o';
sctail=tail1; % tail1 tail2 tailtt
measures=[1];% {'J'} {'SCD'} {'SCS'}
ROIs=[2]; %{'DLPFC_L'} {'DLPFC_R'} {'SGC'}
figplot=1; statstab=1;
stattim=[25 80]; showtim=[0 400]; fillmis=1; fillmiss='movmean';% 'movmean'  movmedian
ROI= {'DLPFC_L' 'DLPFC_R' 'SGC'; LDLPFC RDLPFC SGC};
grouping={}; grouping={'Pre_600' 'Post_600' 'Pre_1200' 'Post_1200' 'Pre_1800' 'Post_1800' ; [pretab & dose600] [posttab & dose600] ...
    [pretab & dose1200] [posttab & dose1200] [pretab & dose1800] [posttab & dose1800]};
mm={'J' 'SCD' 'SCS'; jtab  [scdtab & sctail] [scstab & sctail] ;...
    [strcat('current J (\muA/mm', '^{', '2', '}', ')')]  [strcat('SCD (\muA/mm', '^{', '2', '}', ')')] ['SCS (mm)']}; %
grphtime=timevec.*1000; [~, fix1] = min(abs(grphtime-(showtim(1)))); [~, fix2] = min(abs(grphtime-(showtim(2))));
[~, sta1] = min(abs(grphtime-(stattim(1)))); [~, sta2] = min(abs(grphtime-(stattim(2))));
co=brewermap(4,'RdBu');  c_pre=co(1,:); c_post=co(4,:); colors={c_pre, c_post} ; %RdBu PuOr
clearvars -except datatable fillmis figplot statstab fillmiss stattim ROIs timevec sta1 sta2 fix1 fix2 TimeVec comp2 measures stat_method within pcor remoutlie smoo_method smoo_winsize figg timeind grphtime UU colors co timefunc timevect mm grouping ROI comp grouptab pretab posttab dose600 dose1200 dose1800 jtab scdtab scstab tailtt tail1 tail2 LDLPFC RDLPFC SGC

for r=ROIs %{'DLPFC_L'} {'DLPFC_R'} {'SGC'}
    for j=measures % 1:size(mm(:,measures),2)
        pl=1;
        if figplot==1 %;
            figure('Position', [300 80 700 350]), hold on;
        end
        for i=1:2:size(grouping,2)
            colvecL=mm{2,j} & ROI{2,1}; colvecR=mm{2,j} & ROI{2,2};
            rowvec1=[grouping{2,i} & datatable.filt]; rowvec2=[grouping{2,i+1} & datatable.filt];
            dat{i}=cell2mat(datatable{rowvec1,colvecL}); dat{i+1}=cell2mat(datatable{rowvec2,colvec});
            if remoutlie==1  %%% REMOVE outliers - 2 methods: 1. replace by nan 2. replace by close value
                dat{i}=filloutliers(dat{i},'clip',smoo_method,smoo_winsize,'ThresholdFactor',3); %nearest previous 'movmedian' movmean
                dat{i+1}=filloutliers(dat{i+1},'clip',smoo_method,smoo_winsize,'ThresholdFactor',3); %nearest previous 'movmedian' movmean
            end
            if figplot==1
                subplo{pl}=subplot(2,3,pl); hold on; pl=pl+1;
                f{i}=fill([stattim stattim(2) stattim(1)],[0 0 0 0],'-k','facealpha',.1,'edgecolor','none');
                N{j,r,1}=num2str(size(dat{i},1)); N{j,r,2}=num2str(size(dat{i+1},1));
                avg1=smoothdata(nanmean(dat{i}(:,fix1:fix2),1),'movmean',20); avg2=smoothdata(nanmean(dat{i+1}(:,fix1:fix2),1),'movmean',20);
                stndE1=smoothdata(ste(dat{i}(:,fix1:fix2),1),'movmean',20); stndE2=smoothdata(ste(dat{i+1}(:,fix1:fix2),1),'movmean',20);
                patch([grphtime(fix1:fix2) fliplr(grphtime(fix1:fix2))], [avg1-stndE1 fliplr(avg1+stndE1)],co(4,:),'EdgeColor',co(3,:),'LineWidth',0.05,'FaceVertexAlphaData',0,'FaceAlpha',0.15); %,'EdgeColor',co(1,:)
                patch([grphtime(fix1:fix2) fliplr(grphtime(fix1:fix2))], [avg2-stndE2 fliplr(avg2+stndE2)],co(2,:),'EdgeColor',co(2,:),'LineWidth',0.05,'FaceVertexAlphaData',0,'FaceAlpha',0.15); %,'EdgeColor',co(4,:)
                jplot(j,r,1)=plot(grphtime(fix1:fix2),avg1,'color',[co(4,:) 0.6],'LineWidth',1.8);
                jplot(j,r,2)=plot(grphtime(fix1:fix2),avg2,'color',[co(1,:) 1],'LineWidth',2);
                title(strrep([grouping{1,i}(4:end) '_pulse_' ROI{1,r} ' ' mm{1,j} ],'_',' '))
            end
            clear colvec rowvec avg1 stndE1 avg2 stndE2
        end
        if figplot==1
            linkaxes([subplo{1:3}],'y');
            yy=max(get(gca,'YLim')); subplo{1}.YLim=[0 yy]; subplo{2}.YLim=[0 yy]; subplo{3}.YLim=[0 yy];
            f{1}.YData=[yy yy 0 0]; f{3}.YData=[yy yy 0 0]; f{5}.YData=[yy yy 0 0];
            legend([jplot(j,r,1) jplot(j,r,2)],{['PRE n=' N{j,r,1}] ['POST n=' N{j,r,2}]})
        end
        clear aaa WithinStructure rm ranovatable post
        aaa=cell2mat(cellfun(@(x) mean(x(:,sta1:sta2),2,'omitnan'),dat,'UniformOutput',false)); aaa(aaa==0)=nan;
        if fillmis==1, aaa=fillmissing(aaa,fillmiss,10) ; end %fillmiss='movmedian' 'movmean'
        aaa=array2table(aaa,'VariableNames', grouping(1,:));
        WithinStructure = table(categorical({'Pre';'Post';'Pre';'Post';'Pre';'Post'}),...
            categorical([600;600;1200;1200;1800;1800]),'VariableNames',{'PrePost','Dose'});
        rm = fitrm(aaa, 'Pre_600,Post_600,Pre_1200,Post_1200,Pre_1800,Post_1800~1','WithinDesign',WithinStructure);
        ranovatable = ranova(rm,'WithinModel','PrePost*Dose');
        if figplot==1
            subplo{pl}=subplot(2,9,pl+6:pl+8); hold on; pl=pl+3;
            plotprofile(rm,'PrePost','Group','Dose') ; legend('FontSize',5,'Location','best'); title(strrep([ROI{1,r} ' ' mm{1,j} ],'_',' ')); hold off
            subplo{pl}=subplot(2,9,pl+6:pl+7); hold on; pl=pl+2;
            plotprofile(rm,'Dose'); ylabel('') ;title(strrep([ROI{1,r} ' ' mm{1,j} ],'_',' '));  hold off
            %post = multcompare(rm,'Dose');
            subplo{pl}=subplot(2,9,pl+6:pl+7); hold on; pl=pl+2;
            boxplot(aaa{:,:},'Labels',grouping(1,:)); title(strrep([ROI{1,r} ' ' mm{1,j} ],'_',' ')); hold off
        end
        if statstab==1
            ufig=uifigure('Position', [21 446 689 211]);
            tabl=uitable(ufig,'data',ranovatable{:,[5:end]},'Units', 'Normalized','Position',[0, 0, 1, 1],...
            'ColumnName',ranovatable.Properties.VariableNames(5:end),'RowName',ranovatable.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);%
            %s = uistyle('BackgroundColor',[0.8500 0.3250 0.0980]); addStyle(tabl,s,'column',5)
        end
        %figfig('dose_Sig_current')
        clear subplo
    end
end
%%

for c=1 %trials %1:size(comp,2) %size(comp,2)%%[4]%1:3%size(comp,2)
    %    dat2=datatable(datatable.filt, [ROI{2, r} & mm{2,j}]); clear dat; N={};
    for i=1:size(grouping,2)
        %choosing a grouping based on different columns
        %or based on different  rows
        %if find(grouping{2,i})==4 & sum(grouping{2,i})==1
        %   colvec=posttab & mm{2,j}; %& ROI{2,r};
        %   if i==1
        %       rowvec={};
        %       rowvec{i}=~cellfun('isempty',dat2{:,colvec}) & comp{2,c}.filt & comp{2,c}{:,grouping{2,i}};
        %   elseif i==2
        %       rowvec{i}=~cellfun('isempty',dat2{:,colvec}) & comp{2,c}.filt & ~comp{2,c}{:,grouping{2,i}};
        %   end
        %else
        
        if strcmp(within,'on')
            colvec=mm{2,j} & ROI{2,r};
            rowvec=[grouping{2,i} & datatable.filt];
            %rowvec{i}=~cellfun('isempty',dat2{:,(grouping{2,1} & mm{2,j} )})...    %& ROI{2,r}
            %  & ~cellfun('isempty',dat2{:,(grouping{2,2} & mm{2,j} )}) & comp{2,c}.filt; %& ROI{2,r}
            %rowvec{i}=~cellfun('isempty',dat2{rowveci,colvec}); %& ROI{2,r}
            %%% conditions for different trials:
            %if c==3
            %    rowvec{i}=rowvec{i} & comp{2,c}.Real1_Sham0;
            %end
        else
            colvec=mm{2,j} & ROI{2,r};
            rowvec=[grouping{2,i} & datatable.filt];
            %rowvec{i}=~cellfun('isempty',dat2{:,colvec});
            %rowvec{i}=~cellfun('isempty',dat2{rowveci,colvec});
            %%% conditions for different trials:
            %if c==3 && i==2
            %    rowvec{i}=rowvec{i} & comp{2,c}.Real1_Sham0;
            %end
        end
        %end
        %dat{i}=cell2mat(cellfun(@(x) nanmean(x([r.*2-1 r.*2],:),1), dat2{rowvec{i},colvec}, 'UniformOutput',false));
        %dat{i}=nanmean(cell2mat(dat2{:,:}),1);
        dat{i}=cell2mat(datatable{rowvec,colvec});
        %%% REMOVE outliers - 2 methods: 1. replace by nan 2. replace by close value
        if remoutlie==1
            dat{i}=filloutliers(dat{i},'clip',smoo_method,smoo_winsize,'ThresholdFactor',3); %nearest previous 'movmedian' movmean
        end
        if strcmpi(stat_method,'admean')
            % adaptive mean
            [~, maxloc]=max(dat{i}(:,timeind),[],2);
            datsig_t=[];
            for rang=1:size(maxloc,1)
                minrang=find(timeind,1)+maxloc(rang)-sigwin; if minrang<1, minrang=1; end
                maxrang=find(timeind,1)+maxloc(rang)+sigwin; if maxrang>size(dat{i},2), maxrang=size(dat{i},2); end
                datsig_t(rang)=nanmean(dat{i}(rang,minrang:maxrang));
                clear minrang maxrang
            end
            datsig{r,c,i}=datsig_t;
            clear datsig_t
        end
    end
    % stats
    % if r==rroi(1) && j==measures(1)
    %    stats{1,c}=table;
    %    stats{1,c}=dat2(:,[1:3]);
    %end
    
    %             stats{1,c}.stats1=nan(size(rowvec{1},1),1);
    %             stats{1,c}.stats2=nan(size(rowvec{1},1),1);
    %             stats{1,c}.statsDIFF=nan(size(rowvec{1},1),1);
    %             stats{1,c}.statsRAT=nan(size(rowvec{1},1),1);
    stats{1,1}.stats1=nan(size(dat{1,1},1),1);
    stats{1,1}.stats2=nan(size(dat{1,2},1),1);
    stats{1,1}.stats3=nan(size(dat{1,3},1),1);
    stats{1,1}.stats4=nan(size(dat{1,4},1),1);
    stats{1,1}.stats5=nan(size(dat{1,5},1),1);
    stats{1,1}.stats6=nan(size(dat{1,6},1),1);
    %stats{1,1}.statsDIFF=nan(size(rowvec,1),1);
    %stats{1,1}.statsRAT=nan(size(rowvec,1),1);
    clear pval
    if strcmpi(stat_method,'admean')
        %                 stats{1,c}.stats1(find(rowvec{1}))=datsig{r,c,1};
        %                 stats{1,c}.stats2(find(rowvec{2}))=datsig{r,c,2};
        stats{1,c}.stats1(rowvec{1})=datsig{r,c,1};
        stats{1,c}.stats2(rowvec{2})=datsig{r,c,2};
        %[~, pval]=ttest2(datsig{r,c,1},datsig{r,c,2},'tail','both'); p=timeind;
        pval=ranksum(datsig{r,c,1},datsig{r,c,2}); p=timeind;
    elseif strcmpi(stat_method,'mean') || strcmpi(stat_method,'detect') || strcmpi(stat_method,'detectboot')
        %                 stats{1,c}.stats1(find(rowvec{1}))=nanmean(dat{1}(:,timeind),2);
        %                 stats{1,c}.stats2(find(rowvec{2}))=nanmean(dat{2}(:,timeind),2);
        stats{1,1}.stats1(rowvec{1})=nanmean(dat{1}(:,timeind),2);
        stats{1,1}.stats2(rowvec{2})=nanmean(dat{2}(:,timeind),2);
        %[~, pval]=ttest2(stats{1,c}.stats1,stats{1,c}.stats2,'tail','both'); p=timeind;
        pval=ranksum(stats{1,c}.stats1,stats{1,c}.stats2); p=timeind;
    end
    stats{1,c}.statsDIFF=stats{1,c}.stats1-stats{1,c}.stats2;
    stats{1,c}.statsRAT=stats{1,c}.stats2./stats{1,c}.stats1;
    stats{r+1,c}=pval;
    stats{1,c}.Properties.VariableNames(end-3)={strrep([grouping{1,1} '_' ROI{r,2} '_' mm{1,j} stat_method '_tim' num2str([minval maxval])],'  ','_')};
    stats{1,c}.Properties.VariableNames(end-2)={strrep([grouping{1,2} '_' ROI{r,2} '_' mm{1,j} stat_method '_tim' num2str([minval maxval])],'  ','_')};
    stats{1,c}.Properties.VariableNames(end-1)={strrep(['PrePost_DIFF_' ROI{r,2} '_' mm{1,j} stat_method '_tim' num2str([minval maxval])],'  ','_')};
    stats{1,c}.Properties.VariableNames(end)={strrep(['PrePost_RATIO_' ROI{r,2} '_' mm{1,j} stat_method '_tim' num2str([minval maxval])],'  ','_')};
    if strcmpi(stat_method,'detect')
        % mark significant times on timecourse
        p=[]; h=[];
        [h, p]=ttest2(dat{1},dat{2},'dim',1,'tail','both');  %
        if strcmp(pcor,'on')
            [~, p]=fdr(p,0.05);
        else
            h(isnan(p))=0; p=logical(h);
        end
    end
    if strcmpi(stat_method,'detectboot')
        % bootstrap pvalue
        clear dat11 dat1 dat22 dat2 dat111 dat222
        sig_tim=250;
        dat111=dat{1,1}; dat11=permute(dat111,[2 1]); dat1(1,:,:)=dat11(1:sig_tim,:);
        dat222=dat{1,2}; dat22=permute(dat222,[2 1]); dat2(1,:,:)=dat22(1:sig_tim,:);
        [Ty,diff,CI,pval_b,tcrit,df,ori_mask,boot_mask]=limo_yuen_ttest_boot(dat1,dat2);
        p=[boot_mask zeros(1,[length(grphtime)-sig_tim])];
    end
    %Figure
    if figg==1
        % figure('units','normalized','outerposition',[0.01 0.05 0.6 size(trials,2)/3.2]), hold on,
        subplot(size(measures,2),size(trials,2),pl),hold on,
        pl=pl+1;
        % ploting significance patch
        patch([grphtime fliplr(grphtime)], [zeros(1,length(grphtime))-p*10000 fliplr(zeros(1,length(grphtime))+p*10000)],[0.1 0.1 0.1],'EdgeColor','none','FaceVertexAlphaData',0,'FaceAlpha',0.15);
        for i=1:size(grouping,2)
            N{i}=num2str(size(dat{i},1));
            avg=smoothdata(nanmean(dat{i},1),'movmean',10);
            yoyoyo=dat;
            stndE=smoothdata(ste(dat{i},1),'movmean',10);
            patch([grphtime fliplr(grphtime)], [avg-stndE fliplr(avg+stndE)],colors{i},'EdgeColor',co(i+1,:),'FaceVertexAlphaData',0,'FaceAlpha',0.15);
            jplot(j,i)=plot(grphtime,avg,'color',colors{i},'LineWidth',2);
            miny(i)=min(avg-stndE); maxy(i)=max(avg+stndE);
            clear colvec rowvec avg stndE
        end
        %                if exist('pval') %smoo_method smoo_winsize
        if strcmp(stat_method,'mean')
            title(strrep([ROI{r,2} ' ' comp{1,c} ' ' mm{1,j} char newline stat_method '_time:' num2str([minval maxval]) newline ' smoothing ' smoo_method num2str(smoo_winsize)  ' p=' num2str(pval)],'_',' '))
        elseif strcmp(stat_method,'admean')
            title(strrep([ROI{r,2} ' ' comp{1,c} ' ' mm{1,j} char newline num2str(sigwin) stat_method  '_time:' num2str([minval maxval]) newline ' smoothing ' smoo_method num2str(smoo_winsize)  ' p=' num2str(pval)],'_',' '))
        end
        %                     title(strrep([ROI{1,r} ' ' comp{1,c} ' ' mm{1,j} char newline stat_method '_time:' num2str([minval maxval])  newline ' smoothing ' smoo_method smoo_winsize],'_',' '))
        %                end
        ylabel(mm{3,j}),xlabel('time (ms)');
        legend([jplot(j,:)]',  strcat(strrep(grouping(1,:),'_',' '),{' N=' ' N='},N));
        xlim([0 400]);
        %ylim([floor(min(miny)) ceil(max(maxy)+0.1*max(maxy))]);
        ylim([0 [max(maxy)+0.1*max(maxy)]]);
        set(gca,'FontName','Helvetica Neue','FontSize',8) ;
        
    end
end

% export_fig ECT_MST_rTMS_SCD_SCS_time-dynamics_sgACC_hipp.pdf -q101 -painters -append
% append_pdfs.m
%export_fig ECT_MST_rTMS_SCS_insula_NOfdr.pdf -q101 -painters -append

%stats([size(trials,2)+1],trials+1)=comp(1,trials)
%stats(2:6,5)={'G_Ins_lg_and_S_cent_ins L' 'G_Ins_lg_and_S_cent_ins R' 'Lat_Fis-post R' 'S_circular_insula_ant L' 'S_circular_insula_sup R'}
%clear -vars pl r j c i N avg stndE jplot dat2 dat colvec rowvec jplot pval p datsig minrang maxrang maxloc datsig_t
clearvars -except comp_nets comp2 comp1 c_pre c_post stats TimeVec timevec measures yoyoyo trials stat_method within pcor remoutlie smoo_method smoo_winsize figg timeind grphtime timevect sigwin minval maxval UU colors co timefunc timevect mm grouping ROI comp grouptab pretab posttab jtab scdtab scstab SGCtab Hipptab insula1tab insula2tab insula3tab insula4tab insula5tab tail2
%figfig('')
% export_fig ECT_responders_SCD_SCS_sgACC.pdf -q101 -painters -transparent -append
%figfig('')
% export_fig MST_hipp-SCS_timecourse.PNG -q101 -painters -transparent -append
% export_fig ECT_Hipp_SGC_timecourse.pdf -q101 -painters -transparent -append
%% correlations
% load('A:\WORKINGSET-A\Trial_compare_SCD_SCS_SGC_hipp\3trials_stats_struct.mat')

comp1=comp2;
tt='ect' %'mst' 'rtms'
if strcmpi(tt,'mst')
    % MST
    statsc=table;
    comp1.MST.HAMD.Properties.VariableNames(2:end)=strcat('HAMD_',comp1.MST.HAMD.Properties.VariableNames(2:end));
    statsc=outerjoin(stats{1,4},comp1.MST.HAMD(:,[1 4:47])  ,'MergeKeys',true,'Keys', [1]);
    comp1.MST.QIDS.Properties.VariableNames(2:end)=strcat('QIDS_',comp1.MST.QIDS.Properties.VariableNames(2:end));
    statsc=outerjoin(statsc,comp1.MST.QIDS(:,[1 2:9])  ,'MergeKeys',true,'Keys', [1]);
    comp1.MST.SSI.Properties.VariableNames(2:end)=strcat('SSI_',comp1.MST.SSI.Properties.VariableNames(2:end));
    statsc=outerjoin(statsc,comp1.MST.SSI(:,[1:5])  ,'MergeKeys',true,'Keys', [1]);
    comp1.MST.MOCA.Properties.VariableNames(2:end)=strcat('MOCA_',comp1.MST.MOCA.Properties.VariableNames(2:end));
    statsc=outerjoin(statsc,comp1.MST.MOCA(:,[1:38])  ,'MergeKeys',true,'Keys', [1]);
    %comp1.MST.MATRICS_diff.Properties.VariableNames(2:end)=strcat('MATRICS_',comp1.MST.MATRICS_diff.Properties.VariableNames(2:end));
    %statsc=outerjoin(statsc,comp1.MST.MATRICS_diff  ,'MergeKeys',true,'Keys', [1]);
    statsc = removevars(statsc, [sum(isnan(statsc{:,:}))>size(statsc,1)-3]); % remove vars  with too much missing data
elseif strcmpi(tt,'ect')
    %ECT
    statsc=table;
    comp1.ECT.HAMD.Properties.VariableNames(2:end)=strcat('HAMD_',comp1.ECT.HAMD.Properties.VariableNames(2:end));
    statsc=outerjoin(stats{1,1},comp1.ECT.HAMD(:,[1 3:7 9:15])  ,'MergeKeys',true,'Keys', [1]);
    comp1.ECT.BDI.Properties.VariableNames(2:end)=strcat('BDI_',comp1.ECT.BDI.Properties.VariableNames(2:end));
    statsc=outerjoin(statsc,comp1.ECT.BDI(:,[1 2 3])  ,'MergeKeys',true,'Keys', [1]);
    comp1.ECT.MADRS.Properties.VariableNames(2:end)=strcat('MADRS_',comp1.ECT.MADRS.Properties.VariableNames(2:end));
    statsc=outerjoin(statsc,comp1.ECT.MADRS(:,[1 2 3])  ,'MergeKeys',true,'Keys', [1]);
    comp1.ECT.MOCA.Properties.VariableNames(2:end)=strcat('MOCA_',comp1.ECT.MOCA.Properties.VariableNames(2:end));
    statsc=outerjoin(statsc,comp1.ECT.MOCA(:,[1 3:6])  ,'MergeKeys',true,'Keys', [1]);
    statsc = removevars(statsc, [sum(isnan(statsc{:,:}))>size(statsc,1)-3]); % remove vars  with too much missing data
elseif strcmpi(tt,'rtms')
    %rTMS
    statsc=table;
    comp1.rTMS.Clinical.Properties.VariableNames(2:end)=strcat('Clinical_',comp1.rTMS.Clinical.Properties.VariableNames(2:end));
    statsc=outerjoin(stats{1,3},comp1.rTMS.Clinical(:,[1 20:35])  ,'MergeKeys',true,'Keys', [1]);
    statsc = removevars(statsc, [sum(isnan(statsc{:,:}))>size(statsc,1)-3]); % remove vars  with too much missing data
end
% %
vars=[10 22];% [10 22] [14 86] [22 94] [10 103] [18 48] [18 154]
alpha=0.05;
cortype= 'Spearman' ;% 'Spearman' 'Pearson' 'Kendall'
disregard_zeros='o';
disregard_num=-1;
rem_outlie='o';
% ECT: 23 41 Hipp_SCD_ratio/bdi ratio; 10 35 sgc_scd_diff/hamd_diff; 22 35 hipp_scd_diff/hamd_diff;
% MST: 10 22 sgc_SCS/hamd_diff; 14 86 scd diff/moca total18 diff
[correlations]=Corr_simp(statsc, vars,alpha,cortype, disregard_zeros,disregard_num,rem_outlie);
yo=gca; yo.Title.String=[{[upper(tt) ' trial']}; yo.Title.String];
% export_fig ECT_SGC_Hipp_correlations.PDF -q101 -painters -transparent -append
%figfig('')
%% binarized table variables
statsc{:,117}=statsc{:,91}>0;
%% ROC curve
vars=[18];
resp=37;%Final4(:,108);
clsfun='linear'; %logistic   linear
ROCcurve(statsc,vars,resp,clsfun);
% export_fig MST_ROC.pdf -q101 -painters -transparent -append
%figfig('')
%%
% comparing MOCA pre-post values
moca_out=[~isoutlier(statsc.MOCA_DIFF_Moca_TotalT18)];
premoca=statsc.MOCA_Moca_TotalBaseline(moca_out);
postmoca=statsc.MOCA_Moca_TotalT18(moca_out);
Nn= isfinite(premoca) & isfinite(postmoca);
mean(premoca(Nn)); nanstd(premoca(Nn));
mean(postmoca(Nn)); std(postmoca(Nn))
[h,p] = kstest(premoca(Nn));
[h,p] = kstest(postmoca(Nn))
hist([premoca(Nn) postmoca(Nn)]);
%[h,p]=ttest(premoca(Nn),postmoca(Nn))
[p, ~, stats]=ranksum(premoca(Nn),postmoca(Nn))
% [R, P]=corr(statsc{:,vars},'rows','pairwise', 'type', 'Spearman')
%corrplot(statsc(:,[[6:4:26] 36 38]),'type','Spearman','testR','on')
%corrplot(statsc(:,[6 10 21 22 24 26 28 30]),'type','Spearman','testR','on')
%corrplot(statsc(:,[14 18 84:95]),'type','Spearman','testR','on')
%corrplot(statsc(:,[[6:4:26] 29 30]))
zzz=[]; zzz=[statsc,comp1.ECT.HAMD(:,[2:7 13 14 15])  ]
ZZZ=zzz{:,:}; ZZZ(isinf(ZZZ))=nan; zzz{:,:}=ZZZ;
zzz1=zzz(~isnan(zzz{:,6}),:)

zzz=[]; zzz=[statsc comp1.MST.HAMD(:,[5:9])  ];
ZZZ=zzz{:,:}; ZZZ(isinf(ZZZ))=nan; zzz{:,:}=ZZZ;
zzz1=zzz(~isnan(zzz{:,6}),:)
% export_fig MST_SCD_SCS_sgACC.pdf -q101 -painters -transparent -append
%figfig('')
% export_fig ECT_responders_SCD_SCS_sgACC.pdf -q101 -painters -transparent -append
%figfig('')
% vars=[86 87];
% alpha=0.05;
% disregard_zeros='on';
% rem_outlie='on';
% cortype= 'Pearson' ;% 'Spearman' 'Pearson' 'Kendall'
% %[
% [corMatADHDF, corMatHealthyF, corMatF] = Corr_stats(Final_pre_post,vars, 'sham', 'real', alpha, cortype,disregard_zeros, rem_outlie );
%% Table means
nanmean(statsc.Pre_SGC_SCSadmean_tim15_85)
nanstd(statsc.Pre_SGC_SCSadmean_tim15_85)

nanmean(statsc.Post_SGC_SCSadmean_tim15_85)
nanstd(statsc.Post_SGC_SCSadmean_tim15_85)

[p, ~, stats]=ranksum(statsc.Pre_SGC_SCSadmean_tim15_85,statsc.Post_SGC_SCSadmean_tim15_85);


nanmean(statsc.Pre_Hipp_SCDmean_tim15_85)
nanstd(statsc.Pre_Hipp_SCDmean_tim15_85)

nanmean(statsc.Post_Hipp_SCDmean_tim15_85)
nanstd(statsc.Post_Hipp_SCDmean_tim15_85)

[p, ~, stats]=ranksum(statsc.Pre_Hipp_SCDmean_tim15_85,statsc.Post_Hipp_SCDmean_tim15_85);
%% gramm violin plots
%ya=statsc.Pre_SGC_SCSadmean_tim15_85;
%yb=statsc.Post_SGC_SCSadmean_tim15_85;
ya=statsc{:,4};
yb=statsc{:,5};
co=brewermap(4,'PuOr');
colorsg=[co(1,:); co(4,:)] ;
y=[ya;yb];
x=[cellstr(repmat('Pre-Post',size(y,1),1))];
group=[cellstr(repmat('pre',size(ya,1),1)) ; cellstr(repmat('post',size(yb,1),1))];
clear g
g=gramm('x',x,'y',y,'color',group,'group',group);%,'color',bo,
g.stat_violin('fill','transparent' ,'width',1 ,'dodge',0.4 );%'normalization','count', 'npoints',15,
%g.stat_boxplot('width', 0.6 )
g.stat_summary('type','sem' ,'geom','lines','dodge', 0.4 );%'quartile' '95percentile' 'std' 'fitnormalci'
%g.geom_jitter('width', 0.2,'dodge', 0.4 )
%g.geom_vline('xintercept',1)
g.set_color_options('map',colorsg)
g.set_names('x','treatment stage','y','SCD');
g.set_title('pre-post MST');
figure('Position',[100 100 400 500]);
%g.coord_flip();
g.set_order_options('color',-1)
g.draw();

%g.export('file_name','MST_violin_pre-post_Hipp_SCD','file_type','pdf')
%g.export('file_name','MST_violin_pre-post_SGC_SCS','file_type','pdf')
%%%% gramm violin plots
%ya=statsc.Pre_SGC_SCSadmean_tim15_85;
%yb=statsc.Post_SGC_SCSadmean_tim15_85;
ya=statsc{:,4};
yb=statsc{:,5};
co=brewermap(4,'PuOr');
colorsg=[co(1,:); co(4,:)] ;
y=[ya;yb];
x=[cellstr(repmat('Pre-Post',size(y,1),1))];
group=[cellstr(repmat('pre',size(ya,1),1)) ; cellstr(repmat('post',size(yb,1),1))];
clear g
g=gramm('x',x,'y',y,'color',group,'group',group);%,'color',bo,
g.stat_violin('fill','transparent' ,'width',1 ,'dodge',0.4 );%'normalization','count', 'npoints',15,
%g.stat_boxplot('width', 0.6 )
g.stat_summary('type','sem' ,'geom','lines','dodge', 0.4 );%'quartile' '95percentile' 'std' 'fitnormalci'
%g.geom_jitter('width', 0.2,'dodge', 0.4 )
%g.geom_vline('xintercept',1)
g.set_color_options('map',colorsg)
g.set_names('x','treatment stage','y','SCD');
g.set_title('pre-post MST');
figure('Position',[100 100 400 500]);
%g.coord_flip();
g.set_order_options('color',-1)
g.draw();
%% DIFF gramm violin plots
% comp  {'PRE_POST_ECT'}    {'PRE_POST_MST'}    {'PRE_POST_rTMS'}    {'PRE_POST_MST_notAron'}
figure;  bar([nanmean(stats{1,1}.(6)) nanmean(stats{1,3}.(6)) nanmean(stats{1,4}.(6))])
figure;  bar([nanmean(stats{1,1}.(7)) nanmean(stats{1,3}.(7)) nanmean(stats{1,4}.(7))])
figure;  bar([nanmean(stats{1,1}.(10)) nanmean(stats{1,3}.(10)) nanmean(stats{1,4}.(10))])
figure;  bar([nanmean(stats{1,1}.(11)) nanmean(stats{1,3}.(11)) nanmean(stats{1,4}.(11))])

nanmean(stats{1,1}.(10))
%ya=statsc.Pre_SGC_SCSadmean_tim15_85;
%yb=statsc.Post_SGC_SCSadmean_tim15_85;
ya=statsc{:,4};
yb=statsc{:,5};
co=brewermap(4,'PuOr');
colorsg=[co(1,:); co(4,:)] ;
y=[ya;yb];
x=[cellstr(repmat('Pre-Post',size(y,1),1))];
group=[cellstr(repmat('pre',size(ya,1),1)) ; cellstr(repmat('post',size(yb,1),1))];
clear g
g=gramm('x',x,'y',y,'color',group,'group',group);%,'color',bo,
g.stat_violin('fill','transparent' ,'width',1 ,'dodge',0.4 );%'normalization','count', 'npoints',15,
%g.stat_boxplot('width', 0.6 )
g.stat_summary('type','sem' ,'geom','lines','dodge', 0.4 );%'quartile' '95percentile' 'std' 'fitnormalci'
%g.geom_jitter('width', 0.2,'dodge', 0.4 )
%g.geom_vline('xintercept',1)
g.set_color_options('map',colorsg)
g.set_names('x','treatment stage','y','SCD');
g.set_title('pre-post MST');
figure('Position',[100 100 400 500]);
%g.coord_flip();
g.set_order_options('color',-1)
g.draw();

%g.export('file_name','MST_violin_pre-post_Hipp_SCD','file_type','pdf')
%g.export('file_name','MST_violin_pre-post_SGC_SCS','file_type','pdf')
%%%% gramm violin plots
%ya=statsc.Pre_SGC_SCSadmean_tim15_85;
%yb=statsc.Post_SGC_SCSadmean_tim15_85;
ya=statsc{:,4};
yb=statsc{:,5};
co=brewermap(4,'PuOr');
colorsg=[co(1,:); co(4,:)] ;
y=[ya;yb];
x=[cellstr(repmat('Pre-Post',size(y,1),1))];
group=[cellstr(repmat('pre',size(ya,1),1)) ; cellstr(repmat('post',size(yb,1),1))];
clear g
g=gramm('x',x,'y',y,'color',group,'group',group);%,'color',bo,
g.stat_violin('fill','transparent' ,'width',1 ,'dodge',0.4 );%'normalization','count', 'npoints',15,
%g.stat_boxplot('width', 0.6 )
g.stat_summary('type','sem' ,'geom','lines','dodge', 0.4 );%'quartile' '95percentile' 'std' 'fitnormalci'
%g.geom_jitter('width', 0.2,'dodge', 0.4 )
%g.geom_vline('xintercept',1)
g.set_color_options('map',colorsg)
g.set_names('x','treatment stage','y','SCD');
g.set_title('pre-post MST');
figure('Position',[100 100 400 500]);
%g.coord_flip();
g.set_order_options('color',-1)
g.draw();
%% ploting spesific scout/region activation in brainstorm
% - need to work on

verthemis = []; % Left and right hemisphere labelling


for k = 1:size(head.GridLoc,1)
    if find(cortex.Atlas(1, 7).Scouts(1, 1).Vertices == k) % left cortex
        verthemis(k) = 1;
    else
        verthemis(k) = 2;
    end
end

% Atlas variables
% 2 = Destrieux, 148 scouts
% 3 = Desikan-Killiany, 68 scouts

atlas = 2;
vertlabel = [];

for s = 1:size(cortex.Atlas(1,atlas).Scouts,2) %loop through all the scouts
    scoutnames{s} = cortex.Atlas(1, atlas).Scouts(1, s).Label;
end
scoutnames = scoutnames';

for k = 1:size(head.GridLoc,1) %loop through all sources
    for s = 1:size(cortex.Atlas(1,atlas).Scouts,2) %loop through all the scouts
        if find(cortex.Atlas(1, atlas).Scouts(1, s).Vertices == k)
            vertlabel(k) = s;
            vertlabelnames{k} = cortex.Atlas(1, atlas).Scouts(1, s).Label;
        end
    end
end
%% ploting figures individually



clear -vars pl r j c i N avg stndE jplot dat2 dat colvec rowvec jplot pval p datsig minrang maxrang maxloc datsig_t stats
for r=1:size(ROI,2)
    %figure('units','normalized','outerposition',[0 0 1 1]), hold on,
    
    %pl=1;
    for j=3%1:size(mm,2)
        for c=1:size(comp,2)%%[4]%1:3%size(comp,2)
            
            %subplot(3,3,pl),hold on,
            %pl=pl+1;  N={};
            figure,hold on,
            dat2=comp{2,c}; clear dat
            for i=1:size(grouping,2)
                colvec=grouping{2,i} & mm{2,j} & ROI{2,r};
                if c==3 && i==2
                    % this condition if for the rTMS trial only - includes the sham in the pre condition but exludes it in the post condition
                    rowvec{i}=~cellfun('isempty',dat2{:,colvec}) & comp{2,c}.filt & comp{2,c}.Real1_Sham0;
                elseif [c==2 && i==2] | [c==4 && i==2]
                    rowvec{i}=~cellfun('isempty',dat2{:,colvec}) & comp{2,c}.filt & ~comp{2,c}.DROPOUTS;
                else
                    rowvec{i}=~cellfun('isempty',dat2{:,colvec})& comp{2,c}.filt;
                end
                dat{i}=cell2mat(dat2{rowvec{i},colvec});
                if r==1
                    stats=dat2(:,[1:3]);
                end
                % adaptive mean
                [~, maxloc]=max(dat{i}(:,timeind),[],2);
                datsig_t=[];
                for rang=1:size(maxloc,1)
                    minrang=find(timeind,1)+maxloc(rang)-sigwin; if minrang<1, minrang=1; end
                    maxrang=find(timeind,1)+maxloc(rang)+sigwin; if maxrang>size(dat{i},2), maxrang=size(dat{i},2); end
                    datsig_t(rang)=nanmean(dat{i}(rang,minrang:maxrang));
                end
                datsig{i}=datsig_t;
                %REMOVE outliers - 2 methods: 1. replace by nan 2. replace by close value
                %dat{i}(isoutlier(dat{i},'movmedian',4,1))=nan; %'movmean''movmedian'
                %dat{i}=filloutliers(dat{i},'clip','movmean',3,'ThresholdFactor',4); %nearest previous 'movmedian'
            end
            % adaptive mean
            %[~, pval]=ttest2(datsig{1},datsig{2},'tail','both'); p=timeind;
            % time window ttest pvalue
            stats.stats1=nan(size(rowvec{1},1),1);
            stats.stats2=nan(size(rowvec{1},1),1);
            stats.stats1(find(rowvec{1}))=mean(dat{1}(:,timeind),2);
            stats.stats2(find(rowvec{2}))=mean(dat{2}(:,timeind),2);
            stats.statsDIFF=stats.stats1-stats.stats2;
            stats.statsRAT=stats.stats2./stats.stats1;
            [~, pval]=ttest2(stats.stats1,stats.stats2,'tail','both'); p=timeind;
            stats.Properties.VariableNames(end-3)={strrep([grouping{1,1} '_' ROI{1,r} '_' mm{1,j} '_tim' num2str([minval maxval])],'  ','_')};
            stats.Properties.VariableNames(end-2)={strrep([grouping{1,2} '_' ROI{1,r} '_' mm{1,j} '_tim' num2str([minval maxval])],'  ','_')};
            stats.Properties.VariableNames(end-1)={strrep(['PrePost_DIFF_' ROI{1,r} '_' mm{1,j} '_tim' num2str([minval maxval])],'  ','_')};
            stats.Properties.VariableNames(end)={strrep(['PrePost_RATIO_' ROI{1,r} '_' mm{1,j} '_tim' num2str([minval maxval])],'  ','_')};
            
            % Significance shading
            % regular ttest
            %[p, ~]=ttest2(dat{1},dat{2},'dim',1,'tail','both'); p(isnan(p))=0; p=logical(p);
            % bootstrap pvalue
            %                 clear dat11 dat1 dat22 dat2 dat111 dat222
            %                 sig_tim=270;
            %                 dat111=dat{1,1}; dat11=permute(dat111,[2 1]); dat1(1,:,:)=dat11(1:sig_tim,:);
            %                 dat222=dat{1,2}; dat22=permute(dat222,[2 1]); dat2(1,:,:)=dat22(1:sig_tim,:);
            %                 [Ty,diff,CI,pval,tcrit,df,ori_mask,boot_mask]=limo_yuen_ttest_boot(dat1,dat2);
            %                 p=[boot_mask zeros(1,[length(grphtime)-sig_tim])];
            
            %
            
            % ploting significance patch
            patch([grphtime fliplr(grphtime)], [zeros(1,length(grphtime))-p*10000 fliplr(zeros(1,length(grphtime))+p*10000)],[0.1 0.1 0.1],'EdgeColor','none','FaceVertexAlphaData',0,'FaceAlpha',0.15);
            for i=1:size(grouping,2)
                N{i}=num2str(size(dat{i},1));
                avg=smoothdata(nanmean(dat{i},1),'movmean',10);
                stndE=smoothdata(ste(dat{i},1),'movmean',10);
                patch([grphtime fliplr(grphtime)], [avg-stndE fliplr(avg+stndE)],colors{i},'EdgeColor',co(i+1,:),'FaceVertexAlphaData',0,'FaceAlpha',0.15);
                jplot(j,i)=plot(grphtime,avg,'color',colors{i},'LineWidth',2);
                
                miny(i)=min(avg-stndE); maxy(i)=max(avg+stndE);
                clear colvec rowvec avg stndE
            end
            if exist('pval')
                title(strrep([comp{1,c} ' ' mm{1,j} char newline ROI{1,r} ' p=' num2str(pval)],'_',' '))
            else
                title(strrep([comp{1,c} ' ' mm{1,j} char newline ROI{1,r}],'_',' '))
            end
            ylabel(mm{3,j}),xlabel('time (ms)');
            legend([jplot(j,:)]',  strcat(grouping(1,:),{' N=' ' N='},N));
            xlim([0 400]);
            %ylim([floor(min(miny)) ceil(max(maxy)+0.1*max(maxy))]);
            ylim([0 [max(maxy)+0.1*max(maxy)]]);
            
        end
    end
    % export_fig ECT_MST_rTMS_SCD_SCS_time-dynamics_sgACC_hipp.pdf -q101 -painters -transparent -append
end
clear -vars pl r j c i N avg stndE jplot dat2 dat colvec rowvec jplot pval p datsig minrang maxrang maxloc datsig_t


