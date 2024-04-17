% tcsh BrainStorm.csh
% matlab -nodesktop -r "workspace"

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
    outdir='/media/ipp/Data/EEG_DATA/WellcomeLeap_SCS_2ndrun/';
    %outdir='/media/ipp/Data/EEG_DATA/MST_SCS/';
    eeglabdir='/media/ipp/DATA/Documents/MATLAB/eeglab-2023.1/'; 
    preprocessdir='/media/ipp/DATA/EEG_DATA/WL_TEP/';
    mkdir outdir
    mkdir preprocessdir
    Path = dir([preprocessdir '/**/*.set']); % {Path.name}'
end  
addpath('/media/ipp/DATA/GITs/Localization/','/media/ipp/DATA/GITs/TMS-EEG/' )
addpath(eeglabdir); eeglab nogui
%condi='post';
%fileNames={Path(contains({Path.name},condi,'IgnoreCase',true)).name}';
% fileNames={Path.name};
% cd(outdir)
% brainsto_data=['/media/ipp/DATA/Brainstorm_database'];
% ProtocolName='WellcomeLeap_inter_SCS';
% load([GITS '/TMS-EEG/chanlocs_compumedics_flexnet_64_FULL.mat']);
% chanlocs62=chanlocs69;
% Subloc=[1:6];
% baseline=[-0.95, -0.35];
% sigtime=[0.015, 0.615];

%% analysis
% clear all
addpath([GITS '/Localization/'],[GITS '/TMS-EEG/'],[GITS '/TMS-EEG/BrewerMap/'] )
%outdir='/media/ipp/Data/EEG_DATA/MST_SCS/';
%outdir='/media/ipp/Data/EEG_DATA/rTMS_SCS/'; 
%outdir='/media/ipp/DATA/EEG_DATA/Maintenance/SCS/'
%outdir='/media/ipp/Data/EEG_DATA/ECT_SCS/'
outdir='/media/ipp/Data/EEG_DATA/WellcomeLeap_SCS_2ndrun/'
cd(outdir)
load([outdir 'grand_average_SCD_SCS.mat']);
%outdir='/media/ipp/Data/EEG_DATA/rTMS_SCS/'; 
% outdir='/media/ipp/Data/EEG_DATA/ECT_SCS/'


timevec=datatable.timevec(1,:);
datatable.filt=true(size(datatable,1),1);
pretab=contains(datatable{:,2},'pre','IgnoreCase',true)...
| contains(datatable{:,2},'bl','IgnoreCase',true);
posttab=contains(datatable{:,2},'post','IgnoreCase',true);
jtab=contains(datatable.Properties.VariableNames,'J_','IgnoreCase',true);
scdtab=contains(datatable.Properties.VariableNames,'SCD','IgnoreCase',true);
scstab=contains(datatable.Properties.VariableNames,'SCS','IgnoreCase',true);
LDLPFC=contains(datatable.Properties.VariableNames,'L_DLPFC','IgnoreCase',true);
RDLPFC=contains(datatable.Properties.VariableNames,'R_DLPFC','IgnoreCase',true);
SGC=contains(datatable.Properties.VariableNames,'SGC','IgnoreCase',true);

% timecourse figures
% plotting J SCD SCS
measures=[1];% {'J'} {'SCD'} {'SCS'}
ROIs=[2]; %{'DLPFC_L'} {'DLPFC_R'} {'SGC'}
stattim=[60 120]; showtim=[-50 500]; fillmis=1; fillmiss='movmean';% 'movmean'  movmedian
ROI= {'DLPFC_L' 'DLPFC_R' 'SGC'; LDLPFC RDLPFC SGC};
remoutlie=0; smoo_method='movmedian'; smoo_winsize=10; % movmean movmedian % within='o';
figplot=1; statstab=1;
grouping={}; grouping={'Pre' 'Post'  ; [pretab] [posttab]}; 
trial=char(datatable{1,1}); trial=trial(1:3); %'Maintenance' %datatable{1,1}(1:3);
mm={'J' 'SCD' 'SCS'; jtab  [scdtab] [scstab] ;...
    [strcat('current J (\muA/mm', '^{', '2', '}', ')')]...
     [strcat('SCD (\muA/mm', '^{', '2', '}', ')')] ['SCS (mm)']}; %
grphtime=timevec.*1000; [~, fix1] = min(abs(grphtime-(showtim(1)))); [~, fix2] = min(abs(grphtime-(showtim(2))));
[~, sta1] = min(abs(grphtime-(stattim(1)))); [~, sta2] = min(abs(grphtime-(stattim(2))));
co=brewermap(4,'RdBu');  c_pre=co(1,:); c_post=co(4,:); colors={c_pre, c_post} ; %RdBu PuOr
clearvars -except outdir datatable trial fillmis figplot statstab fillmiss stattim ROIs timevec sta1 sta2 fix1 fix2 TimeVec comp2 measures stat_method within pcor remoutlie smoo_method smoo_winsize figg timeind grphtime UU colors co timefunc timevect mm grouping ROI comp grouptab pretab posttab dose600 dose1200 dose1800 jtab scdtab scstab tailtt tail1 tail2 LDLPFC RDLPFC SGC

for r=ROIs %1:size(ROIs,2)
    for j= measures %1:size(mm(:,measures),2)
        for i=1 %:(grouping,2)
            colvec=mm{2,j} & ROI{2,r};
            rowvec1=[grouping{2,i} & datatable.filt]; rowvec2=[grouping{2,i+1} & datatable.filt];
            dat{i}=datatable{rowvec1,colvec}; dat{i+1}=datatable{rowvec2,colvec};
            if remoutlie==1  %%% REMOVE outliers - 2 methods: 1. replace by nan 2. replace by close value
                dat{i}=filloutliers(dat{i},'clip',smoo_method,smoo_winsize,'ThresholdFactor',3); %nearest previous 'movmedian' movmean
                dat{i+1}=filloutliers(dat{i+1},'clip',smoo_method,smoo_winsize,'ThresholdFactor',3); %nearest previous 'movmedian' movmean
            end
            if remoutlie  
                outlierpre(i)=sum(isoutlier(nanmean(dat{i,1}(:,sta1:sta2),2)));
                outlierpost(i)=sum(isoutlier(nanmean(dat{i,2}(:,sta1:sta2),2)));
                dat{i,1}=dat{i,1}(~isoutlier(nanmean(dat{i,1}(:,sta1:sta2),2)),:);
                dat{i,2}=dat{i,2}(~isoutlier(nanmean(dat{i,2}(:,sta1:sta2),2)),:);
                % figure; plot(dat{:,1}'), figure; plot(dat{:,2}')
            end
            if remzero 
                zerospre(i)=sum(nanmean(dat{i,1}(:,sta1:sta2),2)==0);
                zerospost(i)=sum(nanmean(dat{i,2}(:,sta1:sta2),2)==0);
                dat{i,1}=dat{i,1}(~(nanmean(dat{i,1}(:,sta1:sta2),2))==0,:);
                dat{i,2}=dat{i,2}(~(nanmean(dat{i,2}(:,sta1:sta2),2))==0,:);
                % figure; plot(dat{:,1}'), figure; plot(dat{:,2}')
            end 
            if figplot==1
                figure('Position', [300 80 700 350]), hold on;
                f{i}=fill([stattim stattim(2) stattim(1)],[0 0 0 0],'-k','facealpha',.1,'edgecolor','none');
                N{j,r,1}=num2str(size(dat{i},1)); N{j,r,2}=num2str(size(dat{i+1},1));
                avg1=smoothdata(nanmean(dat{i}(:,fix1:fix2),1),'movmean',smoo_winsize); avg2=smoothdata(nanmean(dat{i+1}(:,fix1:fix2),1),'movmean',smoo_winsize);
                stndE1=smoothdata(ste(dat{i}(:,fix1:fix2),1),'movmean',smoo_winsize); stndE2=smoothdata(ste(dat{i+1}(:,fix1:fix2),1),'movmean',smoo_winsize);
                patch([grphtime(fix1:fix2) fliplr(grphtime(fix1:fix2))], [avg1-stndE1 fliplr(avg1+stndE1)],co(4,:),'EdgeColor',co(3,:),'LineWidth',0.05,'FaceVertexAlphaData',0,'FaceAlpha',0.15); %,'EdgeColor',co(1,:)
                patch([grphtime(fix1:fix2) fliplr(grphtime(fix1:fix2))], [avg2-stndE2 fliplr(avg2+stndE2)],co(2,:),'EdgeColor',co(2,:),'LineWidth',0.05,'FaceVertexAlphaData',0,'FaceAlpha',0.15); %,'EdgeColor',co(4,:)
                jplot(j,r,1)=plot(grphtime(fix1:fix2),avg1,'color',[co(4,:) 0.6],'LineWidth',1.8);
                jplot(j,r,2)=plot(grphtime(fix1:fix2),avg2,'color',[co(1,:) 1],'LineWidth',2);
                title(strrep([trial ' ' ROI{1,r} ' ' mm{1,j} ],'_',' '))
                legend([jplot(j,r,1) jplot(j,r,2)],{['PRE n=' N{j,r,1}] ['POST n=' N{j,r,2}]})
            end
            clear colvec rowvec avg1 stndE1 avg2 stndE2
        end
        figfig(['/media/ipp/Data/EEG_DATA/' trial '_' ROI{1,r} '_' mm{1,j} '.png'])
        clear subplo
    end
end

%% scatter plot
statim=[70 120];
measures=[2];j=measures;% {'J'} {'SCD'} {'SCS'}
ROIs=[3];r=ROIs; %{'DLPFC_L'} {'DLPFC_R'} {'SGC'}
remoutlie=1; scat=1;
grphtime=timevec.*1000; [~, fix1] = min(abs(grphtime-(statim(1)))); [~, fix2] = min(abs(grphtime-(statim(2))));
ROI= {'DLPFC_L' 'DLPFC_R' 'SGC'; LDLPFC RDLPFC SGC};
grouping={}; grouping={'Pre' 'Post'  ; [pretab] [posttab]};
mm={'J' 'SCD' 'SCS'; jtab  [scdtab] [scstab] ;...
    [strcat('current J (\muA/mm', '^{', '2', '}', ')')]  [strcat('SCD (\muA/mm', '^{', '2', '}', ')')] ['SCS (mm)']}; %
colvec=mm{2,j} & ROI{2,r};
rowvec1=[grouping{2,i} & datatable.filt]; rowvec2=[grouping{2,i+1} & datatable.filt];
bl_table=[table(string(datatable{rowvec1,[1]}),'VariableNames',{'subject'}) ...
    table(nanmean(datatable{rowvec1,colvec}(:,fix1:fix2),2),...
    'VariableNames',{[mm{1,j} '_' ROI{1,r} '_' strrep(num2str(statim),'  ','_') '_pre']})];
post_table=[table(string(datatable{rowvec2,[1]}),'VariableNames',{'subject'}) ...
    table(nanmean(datatable{rowvec2,colvec}(:,fix1:fix2),2),...
    'VariableNames',{[mm{1,j} '_' ROI{1,r} '_' strrep(num2str(statim),'  ','_') '_post']})];
pre_post=outerjoin(bl_table,post_table,'MergeKeys',true,'Keys',{'subject'});
if remoutlie==1 
pre_post(:,2)=filloutliers(pre_post(:,2),NaN,"mean");
pre_post(:,3)=filloutliers(pre_post(:,3),NaN,"mean");
end
pre_post.preminuspost=pre_post{:,2}-pre_post{:,3};
depscor=readtable([outdir 'WL_prepost.xlsx']); depscor.Properties.VariableNames(1)=pre_post.Properties.VariableNames(1);
full_table=outerjoin(pre_post,depscor,'MergeKeys',true,'Keys',{'subject'});
full_table.MADRSdiff1wk=full_table.MADRS_Pre-full_table.x1_WeekMADRS_Post;
full_table.MADRSdiff4wk=full_table.MADRS_Pre-full_table.x4_WeekMADRS_Post;
full_table.HAMDdiff1wk=full_table.HAMD_Pre-full_table.x1_WeekHAMD_Post;
full_table.HAMDdiff4wk=full_table.HAMD_Pre-full_table.x4_WeekHAMD_Post;

vars=[2 3]
if scat==1
figure('Renderer', 'painters', 'Position', [100 100 450 500]); hold on, 
scatter(ones(size(full_table,1),1),full_table{:,vars(1)},80,...
    'MarkerFaceColor',co(4,:),'MarkerEdgeColor',co(4,:),'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',.95); 
scatter(ones(size(full_table,1),1)+1,full_table{:,vars(2)},80,...
    'MarkerFaceColor',co(1,:),'MarkerEdgeColor',co(1,:),'MarkerFaceAlpha',.6,'MarkerEdgeAlpha',0.95)
plot([ones(size(full_table,1),1)';[ones(size(full_table,1),1)+1]'], [full_table{:,vars(1)}';full_table{:,vars(2)}'],...
    '-', 'Color',[0 0 0 0.1], 'LineWidth',1.5)
xlim([0.5 2.5]); set(gca,'XTick',[1 2],...
    'xticklabel',[{'pre'} {'post'}],'FontSize',15); ylabel(strrep([ROI{1,r} '  ' mm{3,j}],'_',' '))
[~, P , ~, desc]=ttest2(double(full_table{:,vars(1)}),double(full_table{:,vars(2)}))
cohensD = meanEffectSize(double(full_table{:,vars(1)}),double(full_table{:,vars(2)}),Effect="cohen")
title({['p=' num2str(P)]},{['Cohen''s D=' num2str(cohensD.Effect)]})
%figfig('scatter_SCD.png')
end
%% historgams
co2=brewermap(5,'RdYlBu'); %PuOr RdBu%PuOr RdBu
vars=[5 6 7]; %16 17 18 26 37
a=full_table{:,vars(1)}; a=a(isfinite(a)); 
b=full_table{:,vars(2)}; b=b(isfinite(b)); 
c=full_table{:,vars(3)}; c=c(isfinite(c));
% a = bootstrp(1000,@(x)[mean(x) std(x)],a)
%a = bootstrp(100,@nanmean,a); b=bootstrp(100,@nanmean,b);
%pd = makedist('Normal','mu',nanmean(a),'sigma',nanstd(a))
[fia,xia] = ksdensity(a);
[fib,xib] = ksdensity(b);
[fic,xic] = ksdensity(c);
figure; hold on;
%plot(xia,fia); plot(xib,fib);
%hist(a) gca.Healthycolor; hist(b,ADHDcolor)
fill(xia,fia,co2(1,:),'EdgeColor',[co2(1,:).*0.6],'FaceVertexAlphaData',0.1,'FaceAlpha',0.3); 
fill(xib,fib,co2(3,:),'EdgeColor',[co2(3,:).*0.6],'FaceVertexAlphaData',0.1,'FaceAlpha',0.3);
fill(xic,fic,co2(5,:),'EdgeColor',[co2(5,:).*0.6],'FaceVertexAlphaData',0.1,'FaceAlpha',0.3);
legend(strrep(full_table.Properties.VariableNames(vars),'_',' ')); legend boxoff; %title(strrep(full_table.Properties.VariableNames(vars),'_',' '));
hold off
figfig('clinical_scores_histograms.png')
%a = datasample(a,1000); b = datasample(b,1000)
%binsa=ceil(sqrt(sum(isfinite(a)))); binsb=ceil(sqrt(sum(isfinite(b))));
%figure; hold on;
%generatePDF(a,Healthycolor,'area',binsa); generatePDF(b,ADHDcolor,'area',binsb)
%% clasififer
tab=table
var=[2 3]
tab.var1=[full_table{:,var(1)} ; full_table{:,var(2)}]
tab.Properties.VariableNames(1)= full_table.Properties.VariableNames(var(1)); %[full_table.Properties.VariableNames(var(1)) full_table.Properties.VariableNames(var(2))];
tab.resp=[repmat({'pre'},size(full_table{:,var(1)},1),1) ; repmat({'post'},size(full_table{:,var(2)},1),1)]

clsfun='linear'; %logistic   linear
[CM] = ROCcurve(tab,1,2,clsfun)
figfig('linear_discriminant.png')

%% Clinical Scores
co1=brewermap(4,'PuOr'); %PuOr RdBu

vars=[8 9 10]
figure('Renderer', 'painters', 'Position', [100 100 650 500]); hold on, 
scatter(ones(size(full_table,1),1),full_table{:,vars(1)},80,...
    'MarkerFaceColor',co1(3,:),'MarkerEdgeColor',co1(4,:),'MarkerFaceAlpha',.75,'MarkerEdgeAlpha',.95); 
scatter(ones(size(full_table,1),1)+1,full_table{:,vars(2)},80,...
    'MarkerFaceColor',co1(3,:),'MarkerEdgeColor',co1(4,:),'MarkerFaceAlpha',.75,'MarkerEdgeAlpha',0.95)
scatter(ones(size(full_table,1),1)+2,full_table{:,vars(3)},80,...
    'MarkerFaceColor',co1(3,:),'MarkerEdgeColor',co1(4,:),'MarkerFaceAlpha',.75,'MarkerEdgeAlpha',0.95)
plot([ones(size(full_table,1),1)';[ones(size(full_table,1),1)+1]';[ones(size(full_table,1),1)+2]'], [full_table{:,vars(1)}';full_table{:,vars(2)}';full_table{:,vars(3)}'],...
    '-', 'Color',[0 0 0 0.1], 'LineWidth',1.5)
xlim([0.5 3.5]); set(gca,'XTick',[1 2 3],...
    'xticklabel',[strrep(full_table.Properties.VariableNames(vars(1)),'_',' ')...
    strrep(full_table.Properties.VariableNames(vars(2)),'_',' ') ...
    strrep(full_table.Properties.VariableNames(vars(3)),'_',' ')],'FontSize',15); 
figfig('WL_clinical_scores.png')
%ylabel(strrep([ROI{1,r} '  ' mm{3,j}],'_',' '))
%
%% Correlation matrix for all table
[cormat.R cormat.P]=corrcoef(table2array(full_table(:,2:end)),'rows','pairwise');
corMat=nan(size(cormat.R));
corMat(cormat.P<0.05)=cormat.R(cormat.P<0.05);
corMat=array2table(corMat,'VariableNames',full_table.Properties.VariableNames(2:end),'RowNames',full_table.Properties.VariableNames(2:end)')

%% Correlation plot-  one group 
vars=[4 14];
alpha=0.05;
cortype= 'Pearson' ;% 'Spearman' 'Pearson' 'Kendall'
disregard_zeros='off';
disregard_num=0; %should be 0 to disable
rem_outlie='on';

[correlations]=Corr_simp(full_table, vars,alpha,cortype,disregard_zeros,disregard_num, rem_outlie);
figfig('correlation.png')


%%


%%%%%%% JUNK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2-groups correlations
vars=[1 48];
alpha=0.05;
disregard_zeros='off';
disregard_num=0; %should be 0 to disable
rem_outlie='on';
cortype= 'Spearman' ;% 'Spearman' 'Pearson' 'Kendall'
[corMatADHDF, corMatHealthyF, corMatF] = Corr_stats(full_table,vars, 'Healthy', 'MDD', alpha, cortype,disregard_zeros, disregard_num, rem_outlie );


scatter(pre_post,"subject",[string([mm{1,j} '_' ROI{1,r} '_pre']), string([mm{1,j} '_' ROI{1,r} '_post'])])

figure; plot([pre_post{:,[mm{1,j} '_' ROI{1,r} '_pre']} pre_post{:,[mm{1,j} '_' ROI{1,r} '_post']}])
figure; plot([pre_post{:,[mm{1,j} '_' ROI{1,r} '_pre']} pre_post{:,[mm{1,j} '_' ROI{1,r} '_post']}])
figure; plot([1 2 1 2 1 2],[2 5 3 4 5 6])
gplotmatrix(X,[],group)

swarmchart(pre_post,[ones(size(pre_post,1),1) ones(size(pre_post,1),1)+1],[string([mm{1,j} '_' ROI{1,r} '_pre']), string([mm{1,j} '_' ROI{1,r} '_post'])])
figure; hold on; scatter(ones(size(pre_post,1),1),pre_post{:,2}); scatter(ones(size(pre_post,1),1)+2,pre_post{:,3})
figure; swarmchart(ones(size(pre_post,1),1),pre_post{:,3})


