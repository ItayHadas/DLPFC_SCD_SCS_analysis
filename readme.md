## Acute iTBS dose response neurophysiological metrics

###### - Sara Tremblay Project

Complete Matlab code for Sara Tremablay project: measuring TMS-EEG pre/post acute iTBS dose response (Pulse / time).

The repository includes:

1. codes to fix original files channel locations in the EEG.chalocs field - are located in the OLDER_stats folder
2. AARATEPPipeline shell for automatic TMS-EEG preprocessing
   Cline, Christopher C., Molly V. Lucas, Yinming Sun, Matthew Menezes, and Amit Etkin. “Advanced Artifact Removal for Automated TMS-EEG Data Processing.” In  *2021 10th International IEEE/EMBS Conference on Neural Engineering (NER)* , 1039–42, 2021. [https://doi.org/10.1109/NER49283.2021.9441147](https://doi.org/10.1109/NER49283.2021.9441147).
3. significant current density and current scattering (SCD an SCS) preprocessing script, brainstorm based.
4. Analysis, ploting and repeated measure ANOVA scripts for SCD and SCS analysis
5. added a folder with brainstorm anatomy atlases - deep stractures, and yeo network atlas (I think there are newer network atlases out there)

refs:

primary outcome paper, including TEP, LICI and SICI including spectral decomposition analysis:

Desforges, Manon, Itay Hadas, Brian Mihov, Yan Morin, Mathilde Rochette Braün, Pantelis Lioumis, Reza Zomorrodi, et al. “Dose-Response of Intermittent Theta Burst Stimulation of the Prefrontal Cortex: A TMS-EEG Study.” *Clinical Neurophysiology* 136 (January 20, 2022): 158–72. [https://doi.org/10.1016/j.clinph.2021.12.018](https://doi.org/10.1016/j.clinph.2021.12.018).
