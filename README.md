# debra
View the pulseq sequence at https://bughht.github.io/seqeyes_plugin/?url=https://github.com/xingwangyong/debra/blob/main/debra.seq . The sequence starts with a short GRE acquisition for sensitivity map calculation. Then the real EPI imaging readout starts. The EPI readout has three shots. The 1st is dummy, the 2nd is b-value=0, and the 3rd is b-value=1000s/mm<sup>2</sup>.
![](overview.jpg)
## Install
Step 2,3,4 is only needed for reconstruction
1. Clone this repo. Please note that if you choose to download the repo as a .zip file, the k-space data would **NOT** get downloaded. 
2. After clone, run `git lfs ls-files` to make sure k-space data is downloaded. It should say something like `data/ScanArchive_781487MR4_20260612_210510504.h5`
3. Download the following 3rd party tools
   - ISMRM water fat toolbox: http://cds.ismrm.org/protected/FatWater12_data/fwtoolbox_v1_code.zip
   - GE orchestra-sdk matlab: https://github.com/GEHC-External/MR-Orchestra-SDK-Matlab
4. In `setuppath.m`, set the path to 3rd party toolboxes

## Usage
- `main_gen_sequence_brain_phantom.m`, generate the single shot pulseq sequence used for brain and phantom scan
- `main_recon.m`, recon the phantom data

The `main_recon.m` would output three images
![](result.jpg)