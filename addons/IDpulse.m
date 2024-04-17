function [EEG] = IDpulse (EEG,method,DoublePulseINT,electrodes,newtrigs,evname)
%    IDpulse(EEG,method,DoublePulseINT,impotant,newtrigs,evname)
% if DoublePulseINT=0 then the file is treated as single pulse
% evname= '128'  % evname= '128' 
% electrodes=impotant
% newtrigs=1; %newtrigs=0
%output=table;
%method=2


% evname= 'Tms';
%     newtrigs=1; %newtrigs=0
%     method=1;
%     DoublePulseINT=0; % paired-pulse interval, Zero (0) in case of single pulse
%     DoublePulseINT=2;
%     [impotant ~]=elecName(EEG,{'f5' 'f3' 'f1' 'fc5' 'fc3' 'fc1' 'fcz' 'fc2' 'c1' 'cz' 'c2' 'c4' 'c6' 'cp5' 'cp3'});
%     electrodes=impotant
%  

CriterionORG=40; % starting detection criteria
minISI=1000*EEG.srate/1000; %minimum inter-stimulation interval in samples
stim_no=120; % number of TMS pulses recorded in the data
%electrodes=[7 8 9 16 17 18 ]; %region of stimulation 25 26 27
%electrodes=[7 8 9 16 17 18 3 25 39 45]; %  ST theta burst project electrodes for region of stimulation 3 25 39 45
%%%%%%%%
%electrodes=[3 6 7 8]; %  for brain products demo
%%%%%%%%%

   if method==1 
    for ff=1:size(electrodes,2)
        e=electrodes(ff);
        x=EEG.data(e,:);
        
        %  the "derivative" of the signal
        dx = x(2:end)-x(1:end-1);
        Criterion=CriterionORG;
        ii{ff}=find(abs(dx)>mean(abs(dx)) +Criterion*std(abs(dx)) );
        
        while size(ii{ff},2)<20
            % 'reducing detection Criterion'
            Criterion=Criterion-1;
            ii{ff}=find(abs(dx)>mean(abs(dx)) +Criterion*std(abs(dx)) );
        end
        % cluster the jumps by looking at where they occur
        di = ii{ff}(2:end)-ii{ff}(1:end-1);
        % i.e., look for jumps which comes after long intervals
        kk{ff}=find(di>minISI & di>10*median(di))+1;
        if size(ii{ff},2)>di(1)
            kk{ff}(1:end+1)=[di(1) kk{ff}];% preserve the first pulse
        end
        clear -vars di dx x ev Criterion
        
    end
    
    for yy=1:size(kk,2)
        df(yy)=std(diff(kk{yy}));
        %dff(yy)=ste(diff(kk{yy}));
        dff(yy)=std(diff(kk{yy}))./sqrt(length(kk{yy})-1);
    end
    
    if sum(df<10 & dff<5)>=1 %sum(df<5 & dff<1.5)>=1
        cellsz = cell2mat(cellfun(@size,kk,'uni',false));  cellsz = cellsz(:,2:2:end);
        %atleast 2 epochs at most 150 epochs and at most 5 std in lags
        %conds = find(cellsz>2 & cellsz<150 & df<5 & dff<1); %atleast 2 epochs at most 150 epochs and at most 5 std in lags
        if DoublePulseINT>0
            el = find(cellsz>2 & cellsz<stim_no*2+30 & df<10 & dff<5);
            elseif DoublePulseINT==0
            el = find(cellsz>2 & cellsz<stim_no+30 & df<10 & dff<5);
            end
        el= el(dff(el)==min(dff(el)));%find(df(conds))  ./cellsz(conds)==min(df(conds)./cellsz(conds)));
        
        
    else
        error('The data Triggers are highlly iregular - please look at the data')
        
        %out.original_events_diff=99999;
        %return
    end
    
    
    % Handeling The EEG.events
    
    %el=6
    ii=ii{el};
    kk=kk{el};
    % pop_eegplot( EEG, 1, 1, 1);
    
    % inserting new identified events.
    if newtrigs==1 || size(EEG.event,2)<5 || size(EEG.event,2)+5<sum(kk~=0)
        event=[];
        for ev=1 : sum(kk~=0)
            event(:, ev)= ii(kk(:,ev));
        end
        event=event(event>1000);
        clear -vars ev kk ii
        % pop_eegplot( EEG, 1, 1, 1);
        % Backuping OLD events
        if size(EEG.event,2)<=1
            oldevents=0;
        else
            [EEG.event(1:end).type]=deal('old');
            oldevents=size(EEG.event,2);
            EEG.events_stats.OLD_total=size(EEG.event,2);
            EEG.events_stats.OLD_mean_latency= mean(diff([EEG.event.latency]));
            
        end
        
        % creating new events
        newevents=1;
        for ev= oldevents+1 : oldevents+size(event,2)
            EEG.event(1, ev).type=evname;
            if DoublePulseINT>0
            EEG.event(1, ev).latency= event(newevents)+DoublePulseINT.*EEG.srate/1000;
            elseif DoublePulseINT==0
            EEG.event(1, ev).latency= event(newevents);
            end
            EEG.event(1, ev).urevent= newevents;
            newevents=newevents+1;
                     
        end
        EEG.events_stats.NEW_total= sum(strcmp(cellstr(char(EEG.event.type)),evname));
        EEG.events_stats.NEW_mean_latency= mean(diff([EEG.event(strcmp(cellstr(char(EEG.event.type)),evname)).latency]));
       
        % pop_eegplot( EEG, 1, 1, 1);
        
        
        clear -var di dx x yy ff stdDx Channel Criterion event ev df el sz newevents oldevents e cellsz conds df dff
    end
    
    %clear -var di dx x yy ii kk ff stdDx Channel Criterion event ev el sz newevents oldevents e cellsz conds df dff
   elseif method==2
       
      x=EEG.data(electrodes,:);
       
        [pjj pl1]=   findpeaks(x,'MinPeakProminence',700,'MinPeakDistance',490);
      
        [pjj pl2]=   findpeaks(abs(x),'MinPeakProminence',700,'MinPeakDistance',490);
       
        if  sum(pl1(1:3))>sum(pl2(1:3))
            pl=pl2;
        elseif sum(pl2(1:3))>sum(pl1(1:3))
            pl=pl1;
        end
       
        if size(EEG.event,2)<=1
            oldevents=0;
        else
            [EEG.event(1:end).type]=deal('old');
            oldevents=size(EEG.event,2);
        end
        
        % creating new events
        newevents=1;
        for ev= oldevents+1 : oldevents+size(pl,2)
            EEG.event(1, ev).type=evname;
            if DoublePulseINT>0
            EEG.event(1, ev).latency= pl(newevents)+DoublePulseINT;
            elseif DoublePulseINT==0
            EEG.event(1, ev).latency= pl(newevents);
            end
            EEG.event(1, ev).urevent= newevents;
            newevents=newevents+1;
            
        end
              disp(['created ' num2str(size(pl,2)) ' events.' ])
    end
    