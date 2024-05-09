%% SNR run over grandaverage - for reviewer comments on clinical neurophysiology

[~, bc1] = min(abs(EEG.times-(-180)));
[~, bc2] = min(abs(EEG.times-(-100)));
[~, bc3] = min(abs(EEG.times-(20)));
[~, bc4] = min(abs(EEG.times-(100)));
figure; hold on; plot(EEG.times,squeeze(EEG.data(elecName(EEG,{'f3'}),:,:))); ylim([-40 40]); xlim([-300 300])
fill([-180 -100 -100 -180],[-40 -40 40 40],'-k','facealpha',.15,'edgecolor','none'); 
fill([20 100 100 20],[-40 -40 40 40],'-k','facealpha',.15,'edgecolor','none'); 
fill([-2 20 20 -2],[-40 -40 40 40],'w','edgecolor','none'); 
snr(squeeze(EEG.data(elecName(EEG,{'f3'}),[bc3:bc4],:)),squeeze(EEG.data(elecName(EEG,{'f3'}),[bc1:bc2],:)))
figfig('signal2noise')