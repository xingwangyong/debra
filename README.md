# debra
![](overview.jpg)
## Install
Step 2&3 is only needed for reconstruction
1. Download this repo. 
2. Download the following 3rd party tools
   - ISMRM water fat toolbox: http://cds.ismrm.org/protected/FatWater12_data/fwtoolbox_v1_code.zip
   - GE orchestra-sdk matlab
3. In `setuppath.m`, set the path to 3rd party toolboxes

## Usage
- `main_gen_sequence_brain_phantom.m`, generate the single shot pulseq sequence used for brain and phantom scan
- `main_recon.m`, recon the phantom data
