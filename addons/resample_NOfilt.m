function [EEG] = resample_NOfilt(EEG,freq)
%freq=1000;
EEG.data=double(EEG.data);
oldpnts  = EEG.pnts; preSrate=EEG.srate;
[p,rate] = rat(freq./preSrate, 1e-12);
rate1=floor(rate./p);
[~, ind_t0] = min(abs(EEG.times-(0)));
T=EEG.times; %D=EEG.data;
Dd=[smoothdata(EEG.data(:,[ind_t0-1]:-1:1,:),2,'movmean',rate1) smoothdata(EEG.data(:,ind_t0:end,:),2,'movmean',rate1)]; %movmean gaussian
EEG.times=[]; EEG.times= [flip(T(ind_t0:-rate1:1)) T(ind_t0+rate1:rate1:end)];
EEG.data=[]; EEG.data=[flip(Dd(:,ind_t0:-rate1:1,:),2) Dd(:,ind_t0+rate1:rate1:end,:)];
EEG.srate=preSrate./rate1;
EEG.pnts=size(EEG.data,2);
EEG.xmax=EEG.xmin + (EEG.pnts-1)/EEG.srate;
EEG.times=linspace(EEG.xmin*1000, EEG.xmax*1000, EEG.pnts);
% recompute all event latencies
% -----------------------------
if isfield(EEG.event, 'latency')
    if EEG.trials > 1
        for iEvt = 1:length( EEG.event )
            EEG.event( iEvt ).latency = ( EEG.event( iEvt ).latency - ( EEG.event(iEvt).epoch - 1 ) * oldpnts - 1 ) * p / rate + ( EEG.event(iEvt).epoch - 1 ) * EEG.pnts + 1;
            if isfield(EEG.event, 'duration') && ~isempty(EEG.event(iEvt).duration)
                EEG.event( iEvt ).duration = EEG.event( iEvt ).duration * p/rate;
            end
        end

    end
end
EEG.urevent=[];
if exist('evname')
    EEG.urevent=EEG.event(strcmp({EEG.event.type},evname));
else
    EEG.urevent=EEG.event(strcmp({EEG.event.type},'128'));
end
if isfield(EEG.urevent, 'urevent'), EEG.urevent = rmfield(EEG.urevent, 'urevent'); end
if isfield(EEG.urevent, 'epoch'), EEG.urevent = rmfield(EEG.urevent, 'epoch'); end

EEG = eeg_checkset(EEG, 'eventconsistency');
EEG = Z_append(EEG,'_Dec');
EEG.data=double(EEG.data);
EEG = eeg_checkset( EEG );
end
