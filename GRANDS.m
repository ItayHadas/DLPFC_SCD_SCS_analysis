% matlab -nosplash -nodesktop -r "workspace"
 
clear all
if (ispc)
    sep='\';
    not_sep='/';
    rep_space = ' ';
    GITS='D:\GITs\';
    %Path = dir('\\ad.ucsd.edu\ahs\apps\INTERPSYC\DATA\Wellcome_Leap_802232\Neurophysiology_Data\**\*SPD_*.cdt');
    %Path = dir('A:\WorkingSet\TEPs_GRANDS')
    outdir='A:\WorkingSet\TEPs_GRANDS';
    eeglabdir='D:\MATLAB\EEGLAB'; 
elseif (ismac || isunix)
    sep='/';
    not_sep='\';
    rep_space = '\ ';
    GITS='/media/ipp/DATA/GITs';
    %Path = dir('/mnt/INTERPSYC/DATA/Bipolar_800601/Neurophysiology_Data/**/*SPD_*.cdt');
    outdir='/media/ipp/DATA/EEG_DATA/Bipolar_TEP'; 
    mkdir outdir
    eeglabdir='/media/ipp/DATA/Documents/MATLAB/eeglab-2023.1/'; 
    FTdir='/home/ipp/Documents/Documents/MATLAB/fieldtrip/'
    addpath('/media/ipp/DATA/GITs/Localization/','/media/ipp/DATA/GITs/TMS-EEG/' )
end
run([FTdir sep 'ft_defaults'])
addpath([GITS sep 'TESA']); 
addpath(genpath([GITS sep 'TMS-EEG']),genpath([GITS sep 'AARATEPPipeline']));
addpath(eeglabdir); %addpath('/media/ipp/DATA/Documents/MATLAB/eeglab-2023.1/plugins/ICLabel/viewprops/')
eeglab nogui


%% MST
pathName=[outdir sep];
fileNames=[{'pre_MST_MDD_SP_grand_GrandAverage.set'},...
     {'post_MST_MDD_SP_grand_GrandAverage.set'}];
EEG{1}=pop_loadset('filepath' , pathName, 'filename',fileNames{1});
EEG{2}=pop_loadset('filepath' , pathName, 'filename',fileNames{2});
%% welcomeleap
EEG={};
pathName=['/media/ipp/DATA/EEG_DATA/WL_TEP_7RUNb/'];
fileNames=[{'WL_bl_x_GrandAverage.set'}, {'WL_post_x_GrandAverage.set'},...
     {'WL_bl_y_GrandAverage.set'}, ...
     {'WL_post_y_GrandAverage.set'}, {'WL_ALL_GrandAverage.set'}];
EEG{1}=pop_loadset('filepath' , pathName, 'filename',fileNames{1});
% EEG{1} = eeg_checkset( EEG{1} )
% data{1} = eeglab2fieldtrip(EEG{1}, 'raw')
% cfg = struct(options{:});
% cfg.elec = data{1}.elec;
% cfg.method = 'finite'
% [tmpdata] = ft_scalpcurrentdensity(cfg, data);
% ==EEGCD{1}=pop_currentdensity(EEG{1},'method','finite');
EEG{2}=pop_loadset('filepath' , pathName, 'filename',fileNames{2});
EEG{3}=pop_loadset('filepath' , pathName, 'filename',fileNames{3});
EEG{4}=pop_loadset('filepath' , pathName, 'filename',fileNames{4});
EEG{5}=pop_loadset('filepath' , pathName, 'filename',fileNames{5});

%removing patients that don't have pre AND post
[aa bb cc]=intersect(EEG{1}.subjects',EEG{2}.subjects'); 
EEG{1}.subjects=EEG{1}.subjects(bb); EEG{1}.data=EEG{1}.data(:,:,bb); 
EEG{1}.event=EEG{1}.event(1:size(EEG{1}.data,3)); 
EEG{1}.urevent=EEG{1}.urevent(1:size(EEG{1}.data,3));
EEG{1}.epoch=EEG{1}.epoch(1:size(EEG{1}.data,3)); EEG{1}.trials=size(EEG{1}.data,3);
EEG{2}.subjects=EEG{2}.subjects(cc); EEG{2}.data=EEG{2}.data(:,:,cc);
EEG{2}.event=EEG{2}.event(1:size(EEG{2}.data,3)); 
EEG{2}.urevent=EEG{2}.urevent(1:size(EEG{2}.data,3));
EEG{2}.epoch=EEG{2}.epoch(1:size(EEG{2}.data,3)); EEG{2}.trials=size(EEG{2}.data,3);
EEG{1} = eeg_checkset( EEG{1} ) ; EEG{2} = eeg_checkset( EEG{2} )
[aa2 bb2 cc2]=intersect(EEG{3}.subjects',EEG{4}.subjects'); 
EEG{3}.subjects=EEG{3}.subjects(bb2); EEG{3}.data=EEG{3}.data(:,:,bb2);
EEG{3}.event=EEG{3}.event(1:size(EEG{3}.data,3)); 
EEG{3}.urevent=EEG{3}.urevent(1:size(EEG{3}.data,3));
EEG{3}.epoch=EEG{3}.epoch(1:size(EEG{3}.data,3)); EEG{3}.trials=size(EEG{3}.data,3);
EEG{4}.subjects=EEG{4}.subjects(cc2); EEG{4}.data=EEG{4}.data(:,:,cc2);
EEG{4}.event=EEG{4}.event(1:size(EEG{4}.data,3)); 
EEG{4}.urevent=EEG{4}.urevent(1:size(EEG{4}.data,3));
EEG{4}.epoch=EEG{4}.epoch(1:size(EEG{4}.data,3)); EEG{4}.trials=size(EEG{4}.data,3);
%% LMFPs

color=brewermap(2,'set1'); color=[color(2,:); color(1,:)];
%figure; topoplot([],EEG{2}.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG{2}.chaninfo)
    
dat=[3 4]
elec={'fc3' 'fc4' 'f3' 'f5' 'f1' 'fc1' } %'c1' 'c3' 'c5'
timel=[-300 800];
LMFP=[]; LMFPste=[];
for ii=1:length(EEG)
timeInd=EEG{ii}.times>=timel(1) & EEG{ii}.times<=timel(2);
[elecnum, elecname]=elecName(EEG{ii},elec);
LMFP(:,:,ii)=double(squeeze(mean(std(EEG{ii}.data(elecnum,timeInd,:),0,1),3,"omitnan")));
LMFP(:,:,ii)=smoothdata(LMFP(:,:,ii),'gaussian',10);
LMFPste(:,:,ii)=double(squeeze(ste(std(EEG{ii}.data(elecnum,timeInd,:),0,1),3)));
LMFPste(:,:,ii)=smoothdata(LMFPste(:,:,ii),'gaussian',10);
end
% elec={    'fc3' }
% [elecnum, elecname]=elecName(EEG{ii},elec);
% figure; plot(EEG{ii}.times(timeInd),double(squeeze(mean(EEG{ii}.data(elecnum,timeInd,:),3))));


figfig1=figure('position', [500 200 1200 400]); hold on;
pp1=patch([EEG{dat(1)}.times(timeInd) fliplr(EEG{dat(1)}.times(timeInd))],[LMFP(:,:,dat(1))-LMFPste(:,:,dat(1))...
 fliplr(LMFP(:,:,dat(1))+LMFPste(:,:,dat(1)))],color(1,:),'EdgeColor',[color(1,:).*0.7],'LineWidth',0.2,'FaceVertexAlphaData',0.1,'FaceAlpha',0.3) %
p1=plot(EEG{dat(1)}.times(timeInd),LMFP(:,:,dat(1)),'LineWidth',2);
p1.Color=color(1,:);
pp2=patch([EEG{dat(2)}.times(timeInd) fliplr(EEG{dat(2)}.times(timeInd))],[LMFP(:,:,dat(2))-LMFPste(:,:,dat(2))...
 fliplr(LMFP(:,:,dat(2))+LMFPste(:,:,dat(2)))],color(2,:),'EdgeColor',[color(2,:).*0.7],'LineWidth',0.2,'FaceVertexAlphaData',0.1,'FaceAlpha',0.3)%,'FaceVertexAlphaData',0.1,'FaceAlpha',0.3
p2=plot(EEG{dat(2)}.times(timeInd),LMFP(:,:,dat(2)),'LineWidth',2);
p2.Color=color(2,:);
yy=max(get(gca,'YLim')); %yy=yy.*0.1+yy;
%f=fill([statwin statwin(2) statwin(1)],[yy yy 0 0],'-k','facealpha',.15,'edgecolor','none'); %'facealpha',.15,,'-k','LineWidth',2
chH = get(gca,'Children');
set(gca,'Children',[chH(2:end) ;chH(1)])
%fb=fill([mintrim maxtrim maxtrim mintrim],[yy yy 0 0],'w','edgecolor','none');
%axis([timeWin 0 amp]);
ax=gca; ax.XAxisLocation='origin' ; ax.FontName='Helvetica Neue' ; ax.FontSize=11; ax.YAxisLocation = 'origin'; ax.XTick =[min(ax.XLim):50:max(ax.XLim)];% 0 15 30 60 80 100 120 140 160 180 200 250
legend('show'); clear title; 
% title(['Local Mean Field Power ' strjoin([elecname]) '  time:'...
%     strrep(num2str(statwin),'  ','_')], 'interpreter', 'none'); 
title(['Local Mean Field Power ' strjoin([elecname]) ...
    ], 'interpreter', 'none'); 
ylabel ('\muv'); xlabel ('milliseconds');
legend([p1 p2],{strrep([EEG{dat(1)}.condition ' n=' num2str(size(EEG{dat(1)}.subjects,2))],'_',' '),...
strrep([EEG{dat(2)}.condition ' n=' num2str(size(EEG{dat(2)}.subjects,2))],'_',' ')})
ylim([min(get(gca,'YLim')) round(max([pp1.YData ; pp2.YData]))*1.1])
set(ax ,'Layer', 'Top')
% ['ADHD  n=' num2str(size(corrsADHD(NnA,:),1))]
hold off;
clearvars f p1 p2 ax


%% TEP
% pathName='/media/ipp/Data/EEG_DATA/WellcomeLeap_SCS_2ndrun/',  fileName='WL_bl_GrandAverage.set'
% EEG{1}=pop_loadset('filepath' , pathName, 'filename',fileName);
color=brewermap(2,'set1'); color=[color(2,:); color(1,:)];
%figure; topoplot([],EEG{2}.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG{2}.chaninfo)
    
dat=[1]
elec={'f3' 'f5' } %'fc3' 'fc4' 'f3' 'f5' 'f1' 'fc1' 'c1' 'c3' 'c5'
timel=[-300 800];
statim=[55 125];
TEP=[]; TEPste=[];
for ii=1:length(EEG)
timeInd=EEG{ii}.times>=timel(1) & EEG{ii}.times<=timel(2);
[~, fix1] = min(abs(grphtime-(statim(1))));
[~, fix2] = min(abs(grphtime-(statim(2))));
[elecnum, elecname]=elecName(EEG{ii},elec);
TEP(:,:,ii)=double(squeeze(mean(mean(EEG{ii}.data(elecnum,timeInd,:),3,"omitnan"),1,"omitnan")));
TEP(:,:,ii)=smoothdata(TEP(:,:,ii),'gaussian',10);
TEPste(:,:,ii)=double(squeeze(mean(ste(EEG{ii}.data(elecnum,timeInd,:),3),1,"omitnan")));
TEPste(:,:,ii)=smoothdata(TEPste(:,:,ii),'gaussian',10);
end


figfig1=figure('position', [500 200 1200 400]); hold on;
pp1=patch([EEG{dat(1)}.times(timeInd) fliplr(EEG{dat(1)}.times(timeInd))],[TEP(:,:,dat(1))-TEPste(:,:,dat(1))...
 fliplr(TEP(:,:,dat(1))+TEPste(:,:,dat(1)))],color(1,:),'EdgeColor',[color(1,:).*0.7],'LineWidth',0.2,'FaceVertexAlphaData',0.1,'FaceAlpha',0.3) %
p1=plot(EEG{dat(1)}.times(timeInd),TEP(:,:,dat(1)),'LineWidth',2);
p1.Color=color(1,:);
pp2=patch([EEG{dat(2)}.times(timeInd) fliplr(EEG{dat(2)}.times(timeInd))],[TEP(:,:,dat(2))-TEPste(:,:,dat(2))...
 fliplr(TEP(:,:,dat(2))+TEPste(:,:,dat(2)))],color(2,:),'EdgeColor',[color(2,:).*0.7],'LineWidth',0.2,'FaceVertexAlphaData',0.1,'FaceAlpha',0.3)%,'FaceVertexAlphaData',0.1,'FaceAlpha',0.3
p2=plot(EEG{dat(2)}.times(timeInd),TEP(:,:,dat(2)),'LineWidth',2);
p2.Color=color(2,:);
yy=max(get(gca,'YLim')); %yy=yy.*0.1+yy;
%f=fill([statwin statwin(2) statwin(1)],[yy yy 0 0],'-k','facealpha',.15,'edgecolor','none'); %'facealpha',.15,,'-k','LineWidth',2
chH = get(gca,'Children');
set(gca,'Children',[chH(2:end) ;chH(1)])
%fb=fill([mintrim maxtrim maxtrim mintrim],[yy yy 0 0],'w','edgecolor','none');
%axis([timeWin 0 amp]);
ax=gca; ax.XAxisLocation='origin' ; ax.FontName='Helvetica Neue' ; ax.FontSize=11; ax.YAxisLocation = 'origin'; ax.XTick =[min(ax.XLim):50:max(ax.XLim)];% 0 15 30 60 80 100 120 140 160 180 200 250
legend('show'); clear title; 
% title(['Local Mean Field Power ' strjoin([elecname]) '  time:'...
%     strrep(num2str(statwin),'  ','_')], 'interpreter', 'none'); 
title(['TEP ' strjoin([elecname]) ...
    ], 'interpreter', 'none'); 
ylabel ('\muv'); xlabel ('milliseconds');
legend([p1 p2],{strrep([EEG{dat(1)}.condition ' n=' num2str(size(EEG{dat(1)}.subjects,2))],'_',' '),...
strrep([EEG{dat(2)}.condition ' n=' num2str(size(EEG{dat(2)}.subjects,2))],'_',' ')})
% legend([p1],{strrep([EEG{dat(1)}.condition ' n=' num2str(size(EEG{dat(1)}.subjects,2))],'_',' ')})

ylim([min(get(gca,'YLim')) round(max([pp1.YData ; pp2.YData]))*1.1])
set(ax ,'Layer', 'Top')
% ['ADHD  n=' num2str(size(corrsADHD(NnA,:),1))]
hold off;
clearvars f p1 p2 ax



%% Current Density
clear EEGCD
dat=[5];
EEGCD{1}=pop_currentdensity(EEG{1},'method','raw');
EEGCD{2}=pop_currentdensity(EEG{2},'method','raw');
EEGCD{3}=pop_currentdensity(EEG{3},'method','raw');
EEGCD{4}=pop_currentdensity(EEG{4},'method','raw');
EEGCD{5}=pop_currentdensity(EEG{5},'method','finite');

color=brewermap(2,'set1');
    
ii=5; dat=5
elec={'fc3' 'fc4' 'f3' 'f5'  } % 'f1' 'fc1'
timel=[-200 800];
SCD=[]; SCDste=[]; SCDg=[];
for ii=1:length(EEGCD)
timeInd=EEGCD{ii}.times>=timel(1) & EEGCD{ii}.times<=timel(2);
[elecnum, elecname]=elecName(EEGCD{ii},elec);
SCD(:,:,ii)=double(squeeze(mean(mean(EEGCD{ii}.data(elecnum,timeInd,:),1,"omitnan"),3,"omitnan")));
SCD(:,:,ii)=smoothdata(SCD(:,:,ii),'gaussian',10);
SCDste(:,:,ii)=double(squeeze(ste(mean(EEGCD{ii}.data(elecnum,timeInd,:),1,"omitnan"),3)));
SCDste(:,:,ii)=smoothdata(SCDste(:,:,ii),'gaussian',10);
SCDg(:,:,ii)=double(squeeze(mean(EEGCD{ii}.data(elecnum,timeInd,:),1,"omitnan")));
SCDg(:,:,ii)=smoothdata(SCDg(:,:,ii),'gaussian',10);
end
figfig2=figure('position', [500 200 900 400]); hold on;
p5=plot(EEGCD{dat(1)}.times(timeInd),SCDg(:,:,dat(1)),'color', [color(1,:) 0.2],'LineWidth',2);
p6=plot(EEGCD{dat(1)}.times(timeInd),SCD(:,:,dat(1)),'color', [color(2,:) 0.7],'LineWidth',2)

figfig1=figure('position', [500 200 900 400]); hold on;
pp1=patch([EEGCD{dat(1)}.times(timeInd) fliplr(EEGCD{dat(1)}.times(timeInd))],...
    [SCD(:,:,dat(1))-SCDste(:,:,dat(1)) fliplr(SCD(:,:,dat(1))+SCDste(:,:,dat(1)))],color(1,:),'EdgeColor',[color(1,:).*0.7],'LineWidth',0.2,'FaceVertexAlphaData',0.1,'FaceAlpha',0.3) %
p1=plot(EEGCD{dat(1)}.times(timeInd),SCD(:,:,dat(1)),'LineWidth',2);
p1.Color=color(1,:);
pp2=patch([EEGCD{dat(2)}.times(timeInd) fliplr(EEGCD{dat(2)}.times(timeInd))],[SCD(:,:,dat(2))-SCDste(:,:,dat(2)) fliplr(SCD(:,:,dat(2))+SCDste(:,:,dat(2)))],color(2,:),'EdgeColor',[color(2,:).*0.7],'LineWidth',0.2,'FaceVertexAlphaData',0.1,'FaceAlpha',0.3)%,'FaceVertexAlphaData',0.1,'FaceAlpha',0.3
p2=plot(EEGCD{dat(2)}.times(timeInd),SCD(:,:,dat(2)),'LineWidth',2);
p2.Color=color(2,:);
yy=max(get(gca,'YLim')); %yy=yy.*0.1+yy;
chH = get(gca,'Children');
set(gca,'Children',[chH(2:end) ;chH(1)])
ax=gca; ax.XAxisLocation='origin' ; ax.FontName='Helvetica Neue' ; ax.FontSize=11; ax.YAxisLocation = 'origin'; ax.XTick =[min(ax.XLim):50:max(ax.XLim)];% 0 15 30 60 80 100 120 140 160 180 200 250
legend('show'); clear title; 
title(['Scalp Current Density ' strjoin([elecname]) ...
    ], 'interpreter', 'none'); 
ylabel ('Amp'); xlabel ('milliseconds');
legend([p1 p2],{strrep([EEGCD{dat(1)}.setname ' n=' num2str(size(EEGCD{dat(1)}.subjects,2))],'_',' '),...
strrep([EEGCD{dat(2)}.setname ' n=' num2str(size(EEGCD{dat(2)}.subjects,2))],'_',' ')})
ylim([min(get(gca,'YLim')) round(max([pp1.YData ; pp2.YData]))+1])
set(ax ,'Layer', 'Top')
hold off;
clearvars f p1 p2 ax
%% total
color=brewermap(2,'set1'); color=[color(2,:); color(1,:)];
%figure; topoplot([],EEG{2}.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG{2}.chaninfo)
    
dat=[5]
elec={'fc3' 'fc4' 'f3' 'f5' 'f1' 'fc1' } %'c1' 'c3' 'c5'
timel=[-300 800];
clear LMFP LMFPste LMFPg TEP TEPste TEPg
for ii=5 %1:length(EEG)
timeInd=EEG{ii}.times>=timel(1) & EEG{ii}.times<=timel(2);
[elecnum, elecname]=elecName(EEG{ii},elec);
LMFP(:,:,ii)=double(squeeze(mean(std(EEG{ii}.data(elecnum,timeInd,:),0,1),3,"omitnan")));
LMFP(:,:,ii)=smoothdata(LMFP(:,:,ii),'gaussian',10);
LMFPg(:,:,ii)=double(squeeze(std(EEG{ii}.data(elecnum,timeInd,:),0,1)));
LMFPg(:,:,ii)=smoothdata(LMFPg(:,:,ii),'gaussian',10);
LMFPste(:,:,ii)=double(squeeze(ste(std(EEG{ii}.data(elecnum,timeInd,:),0,1),3)));
LMFPste(:,:,ii)=smoothdata(LMFPste(:,:,ii),'gaussian',10);
TEP(:,:,ii)=double(squeeze(mean(mean(EEG{ii}.data(elecnum,timeInd,:),3,"omitnan"),1,"omitnan")));
TEP(:,:,ii)=smoothdata(TEP(:,:,ii),'gaussian',10);
TEPste(:,:,ii)=double(squeeze(mean(ste(EEG{ii}.data(elecnum,timeInd,:),3),1,"omitnan")));
TEPste(:,:,ii)=smoothdata(TEPste(:,:,ii),'gaussian',10);
TEPg(:,:,ii)=double(squeeze(mean(EEG{ii}.data(elecnum,timeInd,:),1,"omitnan")));
TEPg(:,:,ii)=smoothdata(TEPg(:,:,ii),'gaussian',10);
end
%lmfp
figfig1=figure('position', [500 200 1200 400]); hold on;
%pp1=patch([EEG{dat(1)}.times(timeInd) fliplr(EEG{dat(1)}.times(timeInd))],[LMFP(:,:,dat(1))-LMFPste(:,:,dat(1))...
 %fliplr(LMFP(:,:,dat(1))+LMFPste(:,:,dat(1)))],color(1,:),'EdgeColor',[color(1,:).*0.7],'LineWidth',0.2,'FaceVertexAlphaData',0.1,'FaceAlpha',0.3) %
p1=plot(EEG{dat(1)}.times(timeInd),LMFP(:,:,dat(1)),'LineWidth',2);
p1.Color=color(1,:);
p2=plot(EEG{dat(1)}.times(timeInd),LMFPg(:,:,dat(1)),'color',[color(2,:) 0.1], 'LineWidth',2)
yy=max(get(gca,'YLim')); %yy=yy.*0.1+yy;
%f=fill([statwin statwin(2) statwin(1)],[yy yy 0 0],'-k','facealpha',.15,'edgecolor','none'); %'facealpha',.15,,'-k','LineWidth',2
chH = get(gca,'Children');
set(gca,'Children',[chH(2:end) ;chH(1)])
ax=gca; ax.XAxisLocation='origin' ; ax.FontName='Helvetica Neue' ; ax.FontSize=11; ax.YAxisLocation = 'origin'; ax.XTick =[min(ax.XLim):50:max(ax.XLim)];% 0 15 30 60 80 100 120 140 160 180 200 250
legend('show'); clear title; 
title(['Local Mean Field Power ' strjoin([elecname]) ...
    ], 'interpreter', 'none'); 
ylabel ('\muv'); xlabel ('milliseconds');
legend([p1],{strrep([EEG{dat(1)}.condition ' n=' num2str(size(EEG{dat(1)}.subjects,2))],'_',' ')})
%ylim([min(get(gca,'YLim')) round(max([pp1.YData]))*1.1])
set(ax ,'Layer', 'Top')
hold off;
clearvars f p1 p2 ax

%tep
elec={ 'f3' 'f5' } %'c1' 'c3' 'c5' 'fc3' 'fc4' 'f5' 'f1' 'fc1'
[elecnum, elecname]=elecName(EEG{ii},elec);
timel=[-500 400];
figfig1=figure('position', [500 200 1200 400]); hold on;
pp1=patch([EEG{dat(1)}.times(timeInd) fliplr(EEG{dat(1)}.times(timeInd))],[TEP(:,:,dat(1))-TEPste(:,:,dat(1))...
 fliplr(TEP(:,:,dat(1))+TEPste(:,:,dat(1)))],color(1,:),'EdgeColor',[color(1,:).*0.7],'LineWidth',0.2,'FaceVertexAlphaData',0.1,'FaceAlpha',0.3) %
p3=plot(EEG{dat(1)}.times(timeInd),TEP(:,:,dat(1)),'LineWidth',2);
%p1.Color=color(1,:);
%p2=plot(EEG{dat(1)}.times(timeInd),TEPg(:,:,dat(1)),'color',[color(2,:) 0.1], 'LineWidth',2)
%yy=max(get(gca,'YLim')); %yy=yy.*0.1+yy;
%f=fill([statwin statwin(2) statwin(1)],[yy yy 0 0],'-k','facealpha',.15,'edgecolor','none'); %'facealpha',.15,,'-k','LineWidth',2
chH = get(gca,'Children');
set(gca,'Children',[chH(2:end) ;chH(1)])
ax=gca; ax.XAxisLocation='origin' ; ax.FontName='Helvetica Neue' ; ax.FontSize=11; ax.YAxisLocation = 'origin'; ax.XTick =[min(ax.XLim):50:max(ax.XLim)];% 0 15 30 60 80 100 120 140 160 180 200 250
legend('show'); clear title; 
title(['TEP ' strjoin([elecname]) ...
    ], 'interpreter', 'none'); 
ylabel ('\muv'); xlabel ('milliseconds');
legend([p3],{strrep([EEG{dat(1)}.condition ' n=' num2str(size(EEG{dat(1)}.subjects,2))],'_',' ')})
%ylim([min(get(gca,'YLim')) round(max([pp1.YData]))*1.1])
xlim(timel)
set(ax ,'Layer', 'Top')
hold off;
clearvars f p1 p2 ax


figure('position', [500 200 1200 400]); hold on;
plot(EEG{dat(1)}.times(timeInd),LMFPg(:,:,dat(1)),'color',[color(2,:) 0.1], 'LineWidth',2)

figure('position', [500 200 1200 400]); hold on;
plot(EEG{dat(1)}.times(timeInd),TEPg(:,:,dat(1)),'color',[color(2,:) 0.1], 'LineWidth',2)




pp2=patch([EEG{dat(2)}.times(timeInd) fliplr(EEG{dat(2)}.times(timeInd))],[LMFP(:,:,dat(2))-LMFPste(:,:,dat(2))...
 fliplr(LMFP(:,:,dat(2))+LMFPste(:,:,dat(2)))],color(2,:),'EdgeColor',[color(2,:).*0.7],'LineWidth',0.2,'FaceVertexAlphaData',0.1,'FaceAlpha',0.3)%,'FaceVertexAlphaData',0.1,'FaceAlpha',0.3
p2=plot(EEG{dat(2)}.times(timeInd),LMFP(:,:,dat(2)),'LineWidth',2);
p2.Color=color(2,:);
yy=max(get(gca,'YLim')); %yy=yy.*0.1+yy;
%f=fill([statwin statwin(2) statwin(1)],[yy yy 0 0],'-k','facealpha',.15,'edgecolor','none'); %'facealpha',.15,,'-k','LineWidth',2
chH = get(gca,'Children');
set(gca,'Children',[chH(2:end) ;chH(1)])
%fb=fill([mintrim maxtrim maxtrim mintrim],[yy yy 0 0],'w','edgecolor','none');
%axis([timeWin 0 amp]);
ax=gca; ax.XAxisLocation='origin' ; ax.FontName='Helvetica Neue' ; ax.FontSize=11; ax.YAxisLocation = 'origin'; ax.XTick =[min(ax.XLim):50:max(ax.XLim)];% 0 15 30 60 80 100 120 140 160 180 200 250
legend('show'); clear title; 
% title(['Local Mean Field Power ' strjoin([elecname]) '  time:'...
%     strrep(num2str(statwin),'  ','_')], 'interpreter', 'none'); 
title(['Local Mean Field Power ' strjoin([elecname]) ...
    ], 'interpreter', 'none'); 
ylabel ('\muv'); xlabel ('milliseconds');
legend([p1 p2],{strrep([EEG{dat(1)}.condition ' n=' num2str(size(EEG{dat(1)}.subjects,2))],'_',' '),...
strrep([EEG{dat(2)}.condition ' n=' num2str(size(EEG{dat(2)}.subjects,2))],'_',' ')})
ylim([min(get(gca,'YLim')) round(max([pp1.YData ; pp2.YData]))*1.1])
set(ax ,'Layer', 'Top')
% ['ADHD  n=' num2str(size(corrsADHD(NnA,:),1))]
hold off;
clearvars f p1 p2 ax




%[~, t1] = min(abs(EEG.times-(timel(1))));
%[~, t2] = min(abs(EEG.times-(timel(2))));

%% GRANDAVERAGES
%% GRANDAVERAGES
% Channel TEP
%Path = dir('/mnt/BACKUP14/EEG_DATA/MST/EEG_data/**/*DLPFC*_SP2.cnt');
adjustelec=true;
fileNames={Path.name}';
pathi='/media/ipp/Data/EEG_DATA/WellcomeLeap_SCS_2ndrun/';
folder=dir([pathi '**/*.set']);
load([GITS '/TMS-EEG/chanlocs_compumedics_flexnet_64_FULL.mat']);
chan_interp='on'; chanlocs=chanlocs69; STDcalc=0
prepost={}; prepost={'Pre' 'Post'  ; [pretab] [posttab]};
grouping={}; grouping={'x' 'y'  ; [groupx] [groupy]};
fileNames={folder.name}';
bl=contains(fileNames,'bl','IgnoreCase',true) ;
post=contains(fileNames,'post','IgnoreCase',true) ;
xx=contains(blindcode.code,'x','IgnoreCase',true);
x=contains(fileNames,blindcode(xx,:).Subject);
yy=contains(blindcode.code,'y','IgnoreCase',true);
y=contains(fileNames,blindcode(yy,:).Subject);
GRANDAVERAGES_NEW('WL_bl',[1:6],pathi,fileNames([bl],:), chan_interp,chanlocs) 

GRANDAVERAGES_NEW('WL_bl_x',[1:6],pathi,fileNames([x & bl],:), chan_interp,chanlocs) 
GRANDAVERAGES_NEW('WL_bl_y',[1:6],pathi,fileNames([y & bl],:), chan_interp,chanlocs) 
GRANDAVERAGES_NEW('WL_post_x',[1:6],pathi,fileNames([x & post],:), chan_interp,chanlocs) 
GRANDAVERAGES_NEW('WL_post_y',[1:6],pathi,fileNames([y & post],:), chan_interp,chanlocs) 
GRANDAVERAGES_NEW('WL_ALL',[1:6],pathi,fileNames(:,:), chan_interp,chanlocs) 
%% Maintenance

adjustelec=true;
fileNames={Path.name}';
pathi='/media/ipp/DATA/EEG_DATA/Maintenance/TEP/';
folder=dir([pathi '**/*.set']);
load([GITS '/TMS-EEG/chanlocs_compumedics_flexnet_64_FULL.mat']);
chan_interp='on'; chanlocs=chanlocs69; STDcalc=0

fileNames={folder.name}';
bl=contains(fileNames,'s1','IgnoreCase',true) ;
post=contains(fileNames,'s3','IgnoreCase',true) ;
outlier=contains(fileNames,'IND037','IgnoreCase',true) & contains(fileNames,'s1','IgnoreCase',true) ;

GRANDAVERAGES_NEW('Maintenance_bl',[1:6],pathi,fileNames([bl & ~outlier],:), chan_interp,chanlocs) 






% co=brewermap(4,'set1');
load('/media/ipp/DATA/GITs/TMS-EEG/chanlocs_compumedics_flexnet_64_FULL.mat');
pathi='/media/ipp/Data/EEG_DATA/WellcomeLeap_TEP/'
chan_interp='on'
chanlocs=chanlocs69;
STDcalc=0
Z2_grand_average('WL_SP_all',[1:6],pathi, chan_interp,chanlocs,STDcalc) 

      
