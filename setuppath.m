water_fat_toolbox_path = 'C:\Users\76494/Syncbb/others_toolboxes/fwtoolbox_v1/fwtoolbox_v1_code/';
setup_water_fat_toolbox(water_fat_toolbox_path);
orchestra_path = 'C:\Users\76494\Syncbb/others_toolboxes/orchestra-sdk-2.1-1/orchestra-sdk-2.1-1.matlab';
addpath(orchestra_path);


addpath ./3rd_party/pulseq_v151/matlab
addpath ./3rd_party/BUDA_SLORAKS-main/BUDA_SLORAKS/codes
addpath ./utils

%%
function setup_water_fat_toolbox(tbpath)
BASEPATH = fullfile(tbpath,'/hernando/');
addpath([BASEPATH 'common/']);
addpath([BASEPATH 'graphcut/']);
addpath([BASEPATH 'descent/']);
addpath([BASEPATH 'mixed_fitting/']);
addpath([BASEPATH 'create_synthetic/']);
addpath([BASEPATH 'matlab_bgl/']);



BASEPATH = fullfile(tbpath,'/lu/');
addpath(BASEPATH);
addpath([BASEPATH '/multiResSep']);


BASEPATH = fullfile(tbpath,'/tsao_jiang/');
addpath(BASEPATH);
end