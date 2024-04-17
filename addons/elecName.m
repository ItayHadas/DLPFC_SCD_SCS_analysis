function [ eNum, eName ] = elecName( EEG, elecname )
%converts electrodes namae to its number in spesific EEG dataset
%elecname=[{'16'}]
% elecname=ROI; elecNum
%elecname={'fp1' 'fpz' 'fp2' 'oo' 'fc4' 'fp2'} % elecname={'cpz'}
e = struct2cell(EEG.chanlocs);
elec = e(1,:,:);%squeeze(e(1,:,:));

eNum = [];
eName = {};

for a = 1:size(elecname,2)
    if isempty(find(strcmpi(elecname{1,a},elec))) & ~contains(class(elecname{1,a}),'double')
        
        warning(['electrode ' elecname{1,a} ...
            ' is not present in the EEG.chanlocs set. Electrode not included ']);
        eName(1,size(eName,2)+1) = {'NONEabab'};
        eNum (1,size(eNum,2)+1) =[nan]  ;
    elseif ~isempty(find(strcmpi(elecname{1,a},elecname(1:a-1))))
        ff=find(strcmpi(elecname{1,a},elec));
        warning(['electrode ' EEG.chanlocs(ff).labels ' on channel ' num2str(EEG.chanlocs(ff).urchan)...
            ' is a duplicate. Electrode not included ']);
        eName(1,size(eName,2)+1) = {'NONEabab'};
        eNum (1,size(eNum,2)+1) =[nan]  ;
    else
        if ~contains(class(elecname{1,a}),'double')
        eNum (1,size(eNum,2)+1) = find(strcmpi(elecname{1,a},elec));
        eName (1,size(eName,2)+1) = elec(strcmpi(elecname{1,a},elec));
        else contains(class(elecname{1,a}),'double')
            eNum (1,size(eNum,2)+1) = elecname{:}(a);
            eName (1,size(eName,2)+1) = elec(elecname{:}(a));
    end
    
end
eName=eName(~strcmp(eName,'NONEabab'));
eNum=eNum(~isnan(eNum));

end

