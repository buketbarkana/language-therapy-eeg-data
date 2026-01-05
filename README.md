# language-therapy-eeg-data
EEG features and analysis code for language therapy style study
# language-therapy-eeg-data

This repository contains EEG features and MATLAB analysis code
used in the manuscript:

“Language therapy style shapes neural engagement across recognition,
encoding, and retrieval stages in early childhood”
(submitted to Nature Communications).

## Repository structure

data/
- Excel files containing participant-level EEG features
  organized by task (FWP, UWE, UWR) and region
  (Frnt, Occ, Par, Tmp)

code/matlab/
- MATLAB scripts used to generate figures and tables in the manuscript

## Data description

Each Excel file contains:
- Column 1: Participant ID
- Column 2: Word index (not analyzed)
- Columns 3–13: EEG features and band power
- Columns 14–78: MFCC features (delta–broadband)

## Reproducibility

To reproduce the main results:
1. Download this repository
2. Open MATLAB
3. Navigate to code/matlab
4. Run the relevant script:
   - make_Figure1_NC.m
   - make_Figure2and3.m
   - make_Table2_NC.m
   - make_Table3_NC.m
   - make_Table4_NC.m

## Ethics and privacy

All data are de-identified and shared in accordance
with IRB approval.
