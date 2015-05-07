memo_rep
========

We hope to make this repo a source for electrophysiological recording analysis scripts. Currently, only iEEG data analysis scripts are present, and they rely heavily on Fieldtrip functions and other external scripts. It is desirable, however, that the scripts would eventually become self-contained and widely applicable.

Due to the nature of the memory replay project, however, there will be specialized scripts (e.g. alignment script for clin1 and clin2 data of NYU continuous ECoG recordings) uploaded here, in order to facilitate exchange. It is also likely that, given the potential copyright issues, this repository will be made private (i.e. restricted access to contributors only) in the near future.

## Notes for specific functions

### outlier\_chopperC\_nonincl: 
Two of the inputs needed, "data\_out" and "data\_raw", can be obtained from the following files (for NY394 specifically):
- data\_out: /space/mdeh4/1/halgdev/projects/xjiang/data\_rep/NY\_EEG/NY394/NY394\_D4\_morn\_clin1\_HGPsmoothed.mat
- data\_raw: /space/mdeh4/1/halgdev/projects/xjiang/data\_rep/NY\_EEG/NY394/NY394\_D4\_morn\_clin1\_raw\_matrix.mat
 
### outlier\_unifier: 
Similar to the function above, two of the inputs needed can be found in the above paths. "unifier" also outputs a 61x1 cell containing HGP peak outlier indices for peaks that occur within clean time bins. Example (partial) output can be found at the following path, with variable "good\_peaks" representing the HGP peak outlier indices:
- /space/mdeh4/1/halgdev/projects/xjiang/data\_rep/NY\_EEG/NY394/NY394\_D4\_morn\_clin1\_cleaned\_bins.mat

NY394 good channel labels are also in the above folder. Testing, testing
