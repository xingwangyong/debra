clear
close all

setuppath

ghost_cor_method = 'phzShift'; 
isEnableGhostCorrection = 1;
num_virtual_chan = 12;
isUseSensFromB0 = 0; % deprecated, always set to 0
coil_compression_method = 'svd'; 

B0map_fn = '';

dat_fn = './data/ScanArchive_781487MR4_20260612_210510504.h5';
seq_fn = './debra.seq';


%--------------------------------------------------------------------------
%% get some definitions in .seq file
%--------------------------------------------------------------------------
seq = mr.Sequence();
seq.read(seq_fn,'detectRFuse');

TE = seq.getDefinition('EchoTime');
Nx = seq.getDefinition('Nx');
Ry = seq.getDefinition('Ry');
Ry_low = seq.getDefinition('Ry_low');
Ny_meas = seq.getDefinition('Ny_meas');
Ny = Nx;
ro_os = seq.getDefinition('ro_os');
Nx = ro_os*Nx;warning('tmp!!!!!!!!!!!!!!')
FOV = seq.getDefinition('FOV');
deltaTE  = seq.getDefinition('deltaTE');
nImagingFrame = seq.getDefinition('nImagingFrame');
partFourierFactor = seq.getDefinition('partFourierFactor');
pislquant = seq.getDefinition('pislquant');
Nslices = seq.getDefinition('Nslices');
numSamples = seq.getDefinition('numSamples');
isGEscanner = seq.getDefinition('isGEscanner');
if isGEscanner
    B0 = 3;
else
    B0=2.89;
end
isRefscan = seq.getDefinition('isRefscan');
isConventionalSegmentSpacing = seq.getDefinition('isConventionalSegmentSpacing');
RSegment = seq.getDefinition('RSegment');
isFatSat = seq.getDefinition('isFatSat');
dkyScale_debra = seq.getDefinition('dkyScale_debra');
distance_rf180_readout = seq.getDefinition('distance_rf180_readout');
if isConventionalSegmentSpacing
    dist_ky = dkyScale_debra*Ry*RSegment; % distance between two ky lines, in the unit of deltak
else
    dist_ky = dkyScale_debra*Ry;
end
ky_shift_echo1_aligned = 0;
ky_shift_echo2_aligned = -1*dist_ky;
ky_shift_echo3_aligned = -2*dist_ky;
Necho = seq.getDefinition('debra_nechoes');
dTE_esp_ratio = seq.getDefinition('dTE_esp_ratio');
all_kyShiftFactors_frame = seq.getDefinition('all_kyShiftFactors_frame');

if dTE_esp_ratio==1 && Necho==3
    ky_shift_echo2_aligned = 1*dist_ky;
    ky_shift_echo3_aligned = 2*dist_ky;
end




bval = seq.getDefinition('bval');
bvalue0_threshold = seq.getDefinition('bvalue0_threshold');
Ny_low_meas = seq.getDefinition('Ny_low_meas');
numSamplesLow = seq.getDefinition('numSamplesLow');
Nx_low_org = seq.getDefinition('Nx_low_org');
Nx_low_org = ro_os*Nx_low_org;warning('tmp!!!!!!!!!!!!!!')
base_shift_between_seg = seq.getDefinition('base_shift_between_seg');
isKyShiftBetweenSegments = seq.getDefinition('isKyShiftBetweenSegments');
if isempty(isKyShiftBetweenSegments)
    warning('tmp!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    isKyShiftBetweenSegments = 1;
end
base_shift_between_seg = base_shift_between_seg.*logical(isKyShiftBetweenSegments);
is2ndEcho = seq.getDefinition('is2ndEcho');
if ~is2ndEcho
    Ny_low_meas = 0;
end
isUseGreACS = seq.getDefinition('isUseGreACS');
EchoSpacing = seq.getDefinition('EchoSpacing');
samples_in_gre_acs = seq.getDefinition('samples_in_gre_acs');
length_y_acs = samples_in_gre_acs/numSamples;%isUseGreACS*prod([num_acs,Nslices]);
num_adcs_total = seq.getDefinition('num_adcs_total');
num_adcs_1stEcho = seq.getDefinition('num_adcs_1stEcho');
num_adcs_2ndEcho = seq.getDefinition('num_adcs_2ndEcho');
num_adcs_gre = seq.getDefinition('num_adcs_gre');
delay_time_pi_2_phase= seq.getDefinition('delay_time_pi_2_phase');
Ny_lines_acs = seq.getDefinition('Ny_lines_acs');
isSpiralInNav = seq.getDefinition('isSpiralInNav');
isEpiNav = seq.getDefinition('isEpiNav');
isNavEnabled = (isSpiralInNav || isEpiNav);
assert(~(isSpiralInNav && isEpiNav),'Both isSpiralInNav and isEpiNav are enabled. Only one navigator type is supported.');
adcSamplesNav = seq.getDefinition('adcSamplesNav');
nNavAdc = seq.getDefinition('nNavAdc');
numSamplesNavEach = seq.getDefinition('numSamplesNavEach');
Ny_epiNav_meas = seq.getDefinition('Ny_epiNav_meas');
adcSamplesNavDropped = seq.getDefinition('adcSamplesNavDropped');
isNavSegmentedADC = seq.getDefinition('isNavSegmentedADC');

requiredDefNames = {'numSamples','numSamplesLow','Ny_meas','Ny_low_meas','Nslices','RSegment','nImagingFrame', ...
    'num_adcs_1stEcho','num_adcs_2ndEcho','num_adcs_gre','adcSamplesNav','nNavAdc','numSamplesNavEach','isNavSegmentedADC'};
requiredDefVals = {numSamples,numSamplesLow,Ny_meas,Ny_low_meas,Nslices,RSegment,nImagingFrame, ...
    num_adcs_1stEcho,num_adcs_2ndEcho,num_adcs_gre,adcSamplesNav,nNavAdc,numSamplesNavEach,isNavSegmentedADC};
for iDef = 1:numel(requiredDefNames)
    if isempty(requiredDefVals{iDef})
        error('Missing required seq definition: %s', requiredDefNames{iDef});
    end
end

if isNavEnabled
    if isempty(Ny_epiNav_meas)
        Ny_epiNav_meas = 13;
        warning('Ny_epiNav_meas is missing in .seq. Temporarily hard-coded to 13. TODO: setDefinition(''Ny_epiNav_meas'') in sequence.');
    end
else
    if isempty(Ny_epiNav_meas)
        Ny_epiNav_meas = 0;
    end
end

assert(isNavSegmentedADC==1,'Sequence must use new segmented-nav ADC interface.');
assert(adcSamplesNav==nNavAdc*numSamplesNavEach,'adcSamplesNav mismatch: expected nNavAdc*numSamplesNavEach.');
if isNavEnabled
    assert(adcSamplesNav==Ny_epiNav_meas*numSamplesNavEach, ...
        'adcSamplesNav mismatch: expected Ny_epiNav_meas*numSamplesNavEach (Ny_epiNav_meas=%d, numSamplesNavEach=%d, adcSamplesNav=%d).', ...
        Ny_epiNav_meas, numSamplesNavEach, adcSamplesNav);
end

num_imaging_blocks = numel(seq.blockEvents);
imaging_blocks = cell2mat(seq.blockEvents);
imaging_blocks = reshape(imaging_blocks,[],num_imaging_blocks);
num_adcs = sum(  imaging_blocks(6,:)~=0  ); % these are LARGE ADCs concatenated from small adcs (atom adc)
fprintf('#adc in .seq file is: %d\n',num_adcs);


ind_undo_interleave_epi = tdr_sliceorder(Nslices,2);
ind_undo_interleave_epi = flip(ind_undo_interleave_epi); % 1st to last slice, is from feet to head


%--------------------------------------------------------------------------
%% load data: gre acs and EPI
% first 32x32 phase encodes belong to GRE ACS
%--------------------------------------------------------------------------
isSimulatedData = endsWith(dat_fn,'.mat');
if isGEscanner
    archive = GERecon('Archive.Load', dat_fn);

    assert(num_adcs==archive.FrameCount, 'ADCs in .seq file is different from that in rawdata!')

    nEchoGRE = numel(  seq.getDefinition('TE_gre')  );
    gre_ky_lines = num_adcs_gre-pislquant*nEchoGRE;


    % skip past receive gain calibration TRs (pislquant)
    for n = 1:pislquant*nEchoGRE
        currentControl = GERecon('Archive.Next', archive);
    end

    % read the gre first
    currentControl = GERecon('Archive.Next', archive);
    [nx1 nc] = size(currentControl.Data);
    assert(numSamples==nx1,'ADC samples mismatch!')
    ny1 = nx1;
    d1 = zeros(nx1, nc, gre_ky_lines);
    d1(:,:,1) = currentControl.Data;

    for iy = 2:gre_ky_lines
        currentControl = GERecon('Archive.Next', archive);
        d1(:,:,iy) = currentControl.Data;
    end
    tmp_dat_gre = permute(d1,[1 3 2]);

    tmp_dat_1stEcho = zeros(nx1,nc,num_adcs_1stEcho);
    tmp_dat_2ndEcho = zeros(numSamplesLow, nc, num_adcs_2ndEcho);

    num_repeat_1 = num_adcs_1stEcho/Ny_meas;
    num_repeat_2 = num_adcs_2ndEcho/Ny_low_meas;
    if Ny_low_meas~=0
        assert(num_repeat_1==num_repeat_2,'Mismatch!')
    end
    num_repeat = num_repeat_1;
    counter_1 = 0;
    counter_2 = 0;
    for iter_repeat = 1:num_repeat
        for iter_1stEcho = 1:Ny_meas
            counter_1 = counter_1 + 1;
            currentControl = GERecon('Archive.Next', archive);
            tmp_dat_1stEcho(:,:,counter_1) = currentControl.Data;
        end
        for iter_2ndEcho = 1:Ny_low_meas
            counter_2 = counter_2 + 1;
            currentControl = GERecon('Archive.Next', archive);
            tmp_dat_2ndEcho(:,:,counter_2) = currentControl.Data;
        end
    end
    % disp('finished reading GE data')

    tmp_dat_1stEcho = permute(tmp_dat_1stEcho,[1 3 2]);
    tmp_dat_2ndEcho = permute(tmp_dat_2ndEcho,[1 3 2]);


    %--------------------------------------------------------------------------
    %% parse data: gre acs
    % first 32x32 phase encodes belong to GRE ACS
    %--------------------------------------------------------------------------

    % num_acs = [32];
    Nreadoutlen = numSamples;
    num_chan = nc;


    data_acs = tmp_dat_gre;
    if isUseGreACS
        TE_gre = seq.getDefinition('TE_gre');

        kspace_acs = reshape(data_acs, [Nreadoutlen,numel(TE_gre), Ny_lines_acs, Nslices, num_chan]);
        kspace_acs = permute(kspace_acs,[1 3 2 4 5]);

        % only keep the 1st echo for now
        kspace_acs = squeeze( kspace_acs(:,:,1,:,:) );
        sz = size(kspace_acs);
        kspace_acs = reshape(kspace_acs,sz(1),sz(2),Nslices,sz(end));

        img_acs = ifft2call(kspace_acs);

        % imagesc3d2( rsos(img_acs,4), s(img_acs)/2, 1, [90,90,-90], [0,1e-2]), setGcf(.5)
        ims(rsos(img_acs,4));title 'acs'
    end


    %--------------------------------------------------------------------------
    %% parse data: epi
    %--------------------------------------------------------------------------
    % data after acs scan
    % data_epi_only = data_raw(:,1+length_y_acs:end,:);

    % data without gy/gz blips for phase correction
    data_epi_refscan = tmp_dat_1stEcho(:,1:Ny_meas*Nslices,:);
    data_epi_refscan = reshape(data_epi_refscan, [Nreadoutlen, Ny_meas,Nslices,1,num_chan]);

    img_epi_refscan = ifft2call(data_epi_refscan);


    % imaging data
    if isRefscan
        data_epi_imaging = tmp_dat_1stEcho(:,1+Ny_meas*Nslices:end,:);
    else
        data_epi_imaging = tmp_dat_1stEcho(:,1:end,:);
    end
    data_epi_imaging = reshape(data_epi_imaging, [Nreadoutlen, (Ny_meas), Nslices,RSegment,nImagingFrame,num_chan]);


    if is2ndEcho
        if isRefscan
            data_epi_refscan_low = tmp_dat_2ndEcho(:,1:Ny_low_meas*Nslices,:);
            data_epi_refscan_low = reshape(data_epi_refscan_low, [numSamplesLow, (Ny_low_meas),Nslices,1,num_chan]);

            data_epi_imaging_low = tmp_dat_2ndEcho(:,1+Ny_low_meas *Nslices:end,:);
            data_epi_imaging_low = reshape(data_epi_imaging_low, [numSamplesLow, (Ny_low_meas), Nslices,RSegment,nImagingFrame,num_chan]);
        else
            errro('TODO')
        end
    end





elseif isSimulatedData
    raw_data_struct = load(dat_fn);
    if isfield(raw_data_struct,'signature_value')
        dat_sig = raw_data_struct.signature_value;
    else
        dat_sig = raw_data_struct.DEF.signature;
    end
    seq_sig = seq.signatureValue;
    assert(isequal(dat_sig,seq_sig),'Signature mismatch')

    raw = raw_data_struct.raw;
    Nc = size(raw,2);

    samples_gre = num_adcs_gre*numSamples;
    samples_1stEcho = num_adcs_1stEcho*numSamples;
    samples_2ndEcho = num_adcs_2ndEcho*numSamplesLow;
    samples_total_expected = samples_gre + samples_1stEcho + samples_2ndEcho;

    % assert(size(raw,1)==samples_total_expected,'Mismatch!')

    % first split data into gre and EPI
    tmp_dat_gre = raw(1:samples_gre,:);
    tmp_dat_nav_epi = raw(samples_gre+1:end,:);


    num_repeat_1 = num_adcs_1stEcho/Ny_meas;
    assert(abs(num_repeat_1-round(num_repeat_1))<1e-9,'num_adcs_1stEcho/Ny_meas must be integer.');
    num_repeat = round(num_repeat_1);
    if is2ndEcho
        num_repeat_2 = num_adcs_2ndEcho/Ny_low_meas;
        assert(abs(num_repeat_2-round(num_repeat_2))<1e-9,'num_adcs_2ndEcho/Ny_low_meas must be integer.');
        assert(num_repeat==round(num_repeat_2),'Mismatch between repeats in 1st and 2nd echo.')
    end

    samples_nav_per_repeat = isNavEnabled*adcSamplesNav;
    samples_epi1_per_repeat = Ny_meas*numSamples;
    samples_epi2_per_repeat = Ny_low_meas*numSamplesLow;
    samples_per_repeat = samples_nav_per_repeat + samples_epi1_per_repeat + samples_epi2_per_repeat;
    assert(samples_per_repeat>0,'Invalid samples_per_repeat');
    assert(size(tmp_dat_nav_epi,1)==num_repeat*samples_per_repeat,'Raw EPI samples do not match repeat budget.');

    tmp_dat_nav_epi = reshape(tmp_dat_nav_epi,samples_per_repeat,num_repeat,Nc);
    if isNavEnabled
        tmp_dat_spiNav = tmp_dat_nav_epi(1:samples_nav_per_repeat,:,:);
    else
        tmp_dat_spiNav = zeros(0,num_repeat,Nc,'like',tmp_dat_nav_epi);
    end
    tmp_dat_epi1 = tmp_dat_nav_epi(samples_nav_per_repeat+1:samples_nav_per_repeat+samples_epi1_per_repeat,:,:);
    tmp_dat_epi2 = tmp_dat_nav_epi(samples_nav_per_repeat+samples_epi1_per_repeat+1:end,:,:);

    tmp_dat_1stEcho = reshape(tmp_dat_epi1,numSamples,[],Nc);
    tmp_dat_gre = reshape(tmp_dat_gre,numSamples,[],Nc);
    if is2ndEcho
        tmp_dat_2ndEcho = reshape(tmp_dat_epi2,numSamplesLow,[],Nc);
    else
        tmp_dat_2ndEcho = zeros(numSamplesLow,0,Nc,'like',tmp_dat_1stEcho);
    end
    if isNavEnabled
        tmp_dat_spiNav = reshape(tmp_dat_spiNav,adcSamplesNav,[],Nc);
    end

    assert(size(tmp_dat_1stEcho,2)==num_adcs_1stEcho,'Parsed 1st echo samples do not match num_adcs_1stEcho.');
    if is2ndEcho
        assert(size(tmp_dat_2ndEcho,2)==num_adcs_2ndEcho,'Parsed 2nd echo samples do not match num_adcs_2ndEcho.');
    end


    %--------------------------------------------------------------------------
    %% parse data: gre acs
    % first 32x32 phase encodes belong to GRE ACS
    %--------------------------------------------------------------------------

    % num_acs = [32];
    Nreadoutlen = numSamples;
    num_chan = Nc;


    data_acs = tmp_dat_gre;
    if isUseGreACS
        TE_gre = seq.getDefinition('TE_gre');

        kspace_acs = reshape(data_acs, [Nreadoutlen,numel(TE_gre), Ny_lines_acs, Nslices, num_chan]);
        kspace_acs = permute(kspace_acs,[1 3 2 4 5]);

        % only keep the 1st echo for now
        kspace_acs = squeeze( kspace_acs(:,:,1,:,:) );
        sz = size(kspace_acs);
        kspace_acs = reshape(kspace_acs,sz(1),sz(2),Nslices,sz(end));

        img_acs = ifft2call(kspace_acs);

        % imagesc3d2( rsos(img_acs,4), s(img_acs)/2, 1, [90,90,-90], [0,1e-2]), setGcf(.5)
        ims(rsos(img_acs,4));title 'acs'
    end


    %--------------------------------------------------------------------------
    %% parse data: epi
    %--------------------------------------------------------------------------
    % data after acs scan
    % data_epi_only = data_raw(:,1+length_y_acs:end,:);

    % data without gy/gz blips for phase correction
    data_epi_refscan = tmp_dat_1stEcho(:,1:Ny_meas*Nslices,:);
    data_epi_refscan = reshape(data_epi_refscan, [Nreadoutlen, Ny_meas,Nslices,1,num_chan]);

    img_epi_refscan = ifft2call(data_epi_refscan);


    % imaging data
    if isRefscan
        data_epi_imaging = tmp_dat_1stEcho(:,1+Ny_meas*Nslices:end,:);
    else
        data_epi_imaging = tmp_dat_1stEcho(:,1:end,:);
    end
    data_epi_imaging = reshape(data_epi_imaging, [Nreadoutlen, (Ny_meas), Nslices,RSegment,nImagingFrame,num_chan]);


    if is2ndEcho
        if isRefscan
            data_epi_refscan_low = tmp_dat_2ndEcho(:,1:Ny_low_meas*Nslices,:);
            data_epi_refscan_low = reshape(data_epi_refscan_low, [numSamplesLow, (Ny_low_meas),Nslices,1,num_chan]);

            data_epi_imaging_low = tmp_dat_2ndEcho(:,1+Ny_low_meas *Nslices:end,:);
            data_epi_imaging_low = reshape(data_epi_imaging_low, [numSamplesLow, (Ny_low_meas), Nslices,RSegment,nImagingFrame,num_chan]);
        else
            errro('TODO')
        end
    end


else % is siemens    

    [twix_cells, twix_meta] = mapVBVD_mixed_readout(dat_fn);

    % sanity check
    twix_sig = twix_cells{1}{2}.hdr.Dicom.tSequenceVariant;%dt_main_last.hdr.Dicom.tSequenceVariant;
    seq_sig = seq.signatureValue;
    if ~isequal(twix_sig,seq_sig) && ~contains(twix_sig,seq_sig) % ideally these two should be equal, but somehow sometimes the former has some extra trailing chars
        % warning('tmp, disabled signature check!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
        error('Signature mismatch')
    end


    nx_list = twix_meta.nx_list;
    assert(~isempty(nx_list),'mapVBVD_mixed_readout returned no image groups.');

    idx_main = find(nx_list==numSamples,1,'first');
    idx_low = find(nx_list==numSamplesLow,1,'first');
    idx_nav = find(nx_list==numSamplesNavEach,1,'first');

    if isempty(idx_main) || (is2ndEcho && isempty(idx_low)) || (isNavEnabled && isempty(idx_nav))
        error(['Missing required Nx groups from mapVBVD_mixed_readout. ' ...
            'Expected main=%d, low=%d, nav=%d (if nav enabled), found nx_list=%s'], ...
            numSamples, numSamplesLow, numSamplesNavEach, mat2str(nx_list));
    end

    dt_main = twix_cells{idx_main};
    if iscell(dt_main)
        dt_main_last = dt_main{end};
    else
        dt_main_last = dt_main;
    end



    data_raw = dt_main_last.image.unsorted();
    data_raw = permute(data_raw, [1,3,2]);
    assert(size(data_raw,1)==numSamples, ...
        'Main-readout Nx mismatch: expected numSamples=%d, got %d.', numSamples, size(data_raw,1));
    fprintf('Number of coils in raw data: %d\n',size(data_raw,3));


    noise_data = [];
    try
        if iscell(dt_main)
            if isfield(dt_main{1},'noise')
                noise_data = double(dt_main{1}.noise());
            end
        else
            if isfield(dt_main,'noise')
                noise_data = double(dt_main.noise());
            end
        end
    catch
        noise_data = [];
    end
    if ~isempty(noise_data)
        %noise_data is has a size of 512-cha-128-1-1-2
        noise_data = mean(noise_data,6);
        noise_data = permute(noise_data,[1,3,2]);%512*128*16
        noise_data = reshape(noise_data, [size(noise_data,1)*size(noise_data,2), size(noise_data,3)]);%65536*16
    else
        warning('No noise data found in Twix object. Continue without noise pre-whitening.');
    end



    if is2ndEcho
        dt_low = twix_cells{idx_low};
        if iscell(dt_low)
            dt_low_last = dt_low{end};
        else
            dt_low_last = dt_low;
        end
        data_raw_low = dt_low_last.image.unsorted();
        data_raw_low = permute(data_raw_low, [1,3,2]);
        assert(size(data_raw_low,1)==numSamplesLow, ...
            'Low-readout Nx mismatch: expected numSamplesLow=%d, got %d.', numSamplesLow, size(data_raw_low,1));
    else
        data_raw_low = zeros(numSamplesLow,0,size(data_raw,3),'like',data_raw);
    end

    if isNavEnabled
        dt_nav = twix_cells{idx_nav};
        if iscell(dt_nav)
            dt_nav_last = dt_nav{end};
        else
            dt_nav_last = dt_nav;
        end
        data_raw_nav = dt_nav_last.image.unsorted();
        data_raw_nav = permute(data_raw_nav, [1,3,2]);
        assert(size(data_raw_nav,1)==numSamplesNavEach, ...
            'Navigator Nx mismatch: expected numSamplesNavEach=%d, got %d.', numSamplesNavEach, size(data_raw_nav,1));
    else
        data_raw_nav = [];
    end

    %--------------------------------------------------------------------------
    %% parse data: gre acs
    % first 32x32 phase encodes belong to GRE ACS
    %--------------------------------------------------------------------------

    % num_acs = [32];
    Nreadoutlen = numSamples;
    num_chan = size(data_raw,3);


    data_acs = data_raw(:,1:length_y_acs,:);
    if isUseGreACS
        TE_gre = seq.getDefinition('TE_gre');

        kspace_acs = reshape(data_acs, [Nreadoutlen,numel(TE_gre), Ny_lines_acs, Nslices, num_chan]);
        kspace_acs = permute(kspace_acs,[1 3 2 4 5]);

        % only keep the 1st echo for now
        kspace_acs = squeeze( kspace_acs(:,:,1,:,:) );
        sz = size(kspace_acs);
        kspace_acs = reshape(kspace_acs,sz(1),sz(2),Nslices,sz(end));

        img_acs = ifft2call(kspace_acs);

        % imagesc3d2( rsos(img_acs,4), s(img_acs)/2, 1, [90,90,-90], [0,1e-2]), setGcf(.5)
        ims(rsos(img_acs,4));title 'acs'
    end

    %--------------------------------------------------------------------------
    %% parse data: epi
    %--------------------------------------------------------------------------
    % data after acs scan
    data_epi_only = data_raw(:,1+length_y_acs:end,:);

    num_repeat_1 = num_adcs_1stEcho/Ny_meas;
    assert(abs(num_repeat_1-round(num_repeat_1))<1e-9,'num_adcs_1stEcho/Ny_meas must be integer.');
    num_repeat = round(num_repeat_1);
    if is2ndEcho
        num_repeat_2 = num_adcs_2ndEcho/Ny_low_meas;
        assert(abs(num_repeat_2-round(num_repeat_2))<1e-9,'num_adcs_2ndEcho/Ny_low_meas must be integer.');
        assert(num_repeat==round(num_repeat_2),'Mismatch between repeats in 1st and 2nd echo.')
    end

    assert(size(data_epi_only,2)==Ny_meas*num_repeat, ...
        'Siemens 1st echo line budget mismatch: got %d total lines, expected Ny_meas*num_repeat=%d (Ny_meas=%d, num_repeat=%d)', ...
        size(data_epi_only,2), Ny_meas*num_repeat, Ny_meas, num_repeat);

    tmp_dat_1stEcho = reshape(data_epi_only, Nreadoutlen, Ny_meas*num_repeat, num_chan);
    assert(size(tmp_dat_1stEcho,2)==num_adcs_1stEcho,'Parsed Siemens 1st echo lines do not match num_adcs_1stEcho.');

    if isNavEnabled
        assert(size(data_raw_nav,2)==Ny_epiNav_meas*num_repeat, ...
            'Siemens nav line budget mismatch: got %d nav lines, expected Ny_epiNav_meas*num_repeat=%d (Ny_epiNav_meas=%d, num_repeat=%d)', ...
            size(data_raw_nav,2), Ny_epiNav_meas*num_repeat, Ny_epiNav_meas, num_repeat);
        assert(size(data_raw_nav,3)==num_chan, ...
            'Navigator coil count mismatch: nav has %d, main readout has %d.', size(data_raw_nav,3), num_chan);
        tmp_dat_spiNav = reshape(data_raw_nav, numSamplesNavEach, Ny_epiNav_meas, num_repeat, num_chan);
        tmp_dat_spiNav = reshape(tmp_dat_spiNav, adcSamplesNav, num_repeat, num_chan);
        assert(size(tmp_dat_spiNav,1)==adcSamplesNav,'Parsed Siemens nav samples do not match adcSamplesNav.');
    else
        tmp_dat_spiNav = zeros(0,num_repeat,num_chan,'like',data_epi_only);
    end

    if is2ndEcho
        data_epi_only_low = data_raw_low;
        assert(size(data_epi_only_low,2)==Ny_low_meas*num_repeat, ...
            'Siemens 2nd echo line budget mismatch: got %d total lines, expected Ny_low_meas*num_repeat=%d (Ny_low_meas=%d, num_repeat=%d)', ...
            size(data_epi_only_low,2), Ny_low_meas*num_repeat, Ny_low_meas, num_repeat);
        tmp_dat_2ndEcho = reshape(data_epi_only_low, numSamplesLow, Ny_low_meas*num_repeat, num_chan);
        assert(size(tmp_dat_2ndEcho,2)==num_adcs_2ndEcho,'Parsed Siemens 2nd echo lines do not match num_adcs_2ndEcho.');
    else
        tmp_dat_2ndEcho = zeros(numSamplesLow,0,num_chan,'like',data_epi_only);
    end

    % data without gy/gz blips for phase correction
    data_epi_refscan = tmp_dat_1stEcho(:,1:Ny_meas*Nslices,:);
    data_epi_refscan = reshape(data_epi_refscan, [Nreadoutlen, Ny_meas,Nslices,1,num_chan]);

    img_epi_refscan = ifft2call(data_epi_refscan);


    % imaging data
    if isRefscan
        data_epi_imaging = tmp_dat_1stEcho(:,1+Ny_meas*Nslices:end,:);
    else
        data_epi_imaging = tmp_dat_1stEcho(:,1:end,:);
    end
    expected_numel_imaging = Nreadoutlen*Ny_meas*Nslices*RSegment*nImagingFrame*num_chan;
    assert(numel(data_epi_imaging)==expected_numel_imaging, ...
        'EPI imaging reshape mismatch: got %d elems, expected %d (Nreadoutlen=%d, Ny_meas=%d, Nslices=%d, RSegment=%d, nImagingFrame=%d, num_chan=%d)', ...
        numel(data_epi_imaging), expected_numel_imaging, Nreadoutlen, Ny_meas, Nslices, RSegment, nImagingFrame, num_chan);
    data_epi_imaging = reshape(data_epi_imaging, [Nreadoutlen, Ny_meas, Nslices,RSegment,nImagingFrame,num_chan]);


    if is2ndEcho
        if isRefscan
            data_epi_refscan_low = tmp_dat_2ndEcho(:,1:Ny_low_meas*Nslices,:);
            data_epi_refscan_low = reshape(data_epi_refscan_low, [numSamplesLow, Ny_low_meas,Nslices,1,num_chan]);

            data_epi_imaging_low = tmp_dat_2ndEcho(:,1+Ny_low_meas*Nslices:end,:);
            expected_numel_imaging_low = numSamplesLow*Ny_low_meas*Nslices*RSegment*nImagingFrame*num_chan;
            assert(numel(data_epi_imaging_low)==expected_numel_imaging_low, ...
                'EPI low imaging reshape mismatch: got %d elems, expected %d (numSamplesLow=%d, Ny_low_meas=%d, Nslices=%d, RSegment=%d, nImagingFrame=%d, num_chan=%d)', ...
                numel(data_epi_imaging_low), expected_numel_imaging_low, numSamplesLow, Ny_low_meas, Nslices, RSegment, nImagingFrame, num_chan);
            data_epi_imaging_low = reshape(data_epi_imaging_low, [numSamplesLow, Ny_low_meas, Nslices,RSegment,nImagingFrame,num_chan]);
        else
            errro('TODO')
        end
    end
end






%--------------------------------------------------------------------------
%% load trajectory for regridding
%--------------------------------------------------------------------------
traj_recon_delay = 0e-6;
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP('trajectory_delay', traj_recon_delay);

start_ind = 1 + Nreadoutlen*(length_y_acs);

ktraj_adc_spiralInNav = ktraj_adc(:,start_ind:end);
nsegments_in_ref = 1;
tmp_repeat = Nslices*nsegments_in_ref + Nslices*RSegment*nImagingFrame;
ktraj_adc_spiralInNav = reshape(ktraj_adc_spiralInNav,3,[],tmp_repeat);

if isNavEnabled
    start_ind2 = 1 + adcSamplesNav;
else
    start_ind2 = 1;
end
ktraj_adc_epi = ktraj_adc_spiralInNav(:,start_ind2:end,:); % the 2nd dim is all samples in one slice
ktraj_adc_epi = reshape(ktraj_adc_epi, size(ktraj_adc_epi,1),[]);
ktraj_adc_spiNav = ktraj_adc_spiralInNav(:,1:start_ind2-1,:);



nsegments_in_ref = 1;
number_of_elements_in_ref = Nslices*nsegments_in_ref*Ny_meas*numSamples;
if is2ndEcho
    number_of_elements_in_ref = number_of_elements_in_ref + Nslices*nsegments_in_ref*Ny_low_meas*numSamplesLow;
end
ktraj_adc_epi_refonly = ktraj_adc_epi(:,1:number_of_elements_in_ref);
ktraj_adc_epi_refonly_high_low = reshape(ktraj_adc_epi_refonly,3,[],Nslices);
ktraj_adc_epi_refonly_high = ktraj_adc_epi_refonly_high_low(:,1:Ny_meas*numSamples,:,:,:);
ktraj_adc_epi_refonly_high = reshape(ktraj_adc_epi_refonly_high,3,numSamples,Ny_meas,Nslices);
ktraj_adc_epi_refonly_low  = ktraj_adc_epi_refonly_high_low(:,1+Ny_meas*numSamples:end,:,:,:);
ktraj_adc_epi_refonly_low = reshape(ktraj_adc_epi_refonly_low,3,numSamplesLow,Ny_low_meas,Nslices);




if isRefscan
    ktraj_adc_epi_high_low = ktraj_adc_epi(:,1+number_of_elements_in_ref:end);
else
    ktraj_adc_epi_high_low = ktraj_adc_epi(:,1:end);
end
ktraj_adc_epi_high_low = reshape(ktraj_adc_epi_high_low,3,[],Nslices,RSegment,nImagingFrame);
ktraj_adc_epi_high = ktraj_adc_epi_high_low(:,1:Ny_meas*numSamples,:,:,:);
ktraj_adc_epi_high = reshape(ktraj_adc_epi_high,3,numSamples,Ny_meas,Nslices,RSegment,nImagingFrame);


ktraj_adc_epi_low = ktraj_adc_epi_high_low(:,1+Ny_meas*numSamples:end,:,:,:);
ktraj_adc_epi_low = reshape(ktraj_adc_epi_low,3,numSamplesLow,Ny_low_meas,Nslices,RSegment,nImagingFrame);
% figure;plot(vec(ktraj_adc_epi_low(1,:,:,1,:,:)/(1/FOV(1))),vec(ktraj_adc_epi_low(2,:,:,1,:,:))/(1/1/FOV(1)),'.')



slice_ind = Nslices;
seg_ind = 1:RSegment;
% frame_ind = isRefscan+nImagingFrame;


% figure;
% plot(vec(ktraj_adc_epi_tmp(1,:,:,slice_ind,seg_ind,frame_ind)),vec(ktraj_adc_epi_tmp(2,:,:,slice_ind,seg_ind,frame_ind)),'-o')
%
% figure;
% plot(vec(ktraj_adc_epi_tmp(1,:,:,slice_ind,seg_ind,frame_ind)),'-o')
%
% figure;
% plot(vec(ktraj_adc_epi_tmp(2,:,:,slice_ind,seg_ind,frame_ind)),'-o')
%
% figure;
% plot(unique(vec(ktraj_adc_epi_tmp(2,:,:,slice_ind,seg_ind,frame_ind))),'-o')


ky_traj = squeeze(  ktraj_adc_epi_high(2,1,:,slice_ind,seg_ind,end-nImagingFrame+1:end)  );
% assert(numel(unique((ky_traj)))==RSegment*Ny_meas);
kymax = max(ky_traj(:));
kymin = min(ky_traj(:));
Nymax = round( kymax/(1/FOV(2)) );
Nymin = round( kymin/(1/FOV(2)) );
Ny_actual = 2*max(abs(Nymax),abs(Nymin));


ky_traj_low = squeeze(  ktraj_adc_epi_low(2,1,:,slice_ind,seg_ind,end-nImagingFrame+1:end)  );
kylowmax = max(ky_traj_low(:));
kylowmin = min(ky_traj_low(:));
Nylowmax = round( kylowmax/(1/FOV(2)) );
Nylowmin = round( kylowmin/(1/FOV(2)) );
Ny_actual_low = 2*max(abs(Nylowmax),abs(Nylowmin));



figure;
plot(ktraj_adc_epi(1,:),ktraj_adc_epi(2,:),'.')
title 'traj, epi-only'
figure;
plot(ktraj_adc(1,:),ktraj_adc(2,:),'.')


kx_ind = 1;
ky = squeeze(   ktraj_adc_epi_high(2,kx_ind,:,slice_ind,:,end-nImagingFrame+1:end)   );
ky = round(  ky/(1/FOV(2))  );
ky_seg1 = ky(1:Ny_meas);
ky_seg1_echo1 = ky_seg1(1:Necho:end);%assert(Necho==3,'currently only support necho=3')
ky_shift_to_center = abs(ky_seg1_echo1(end)) - Ny_actual/2; % later, when shifting the k-space, we assume that 1st echo of 1st seg is at the edge
% ky_shift_to_center = 1 + ky_shift_to_center; % in matlab, array of length 2N, DC is at N+1

if is2ndEcho
ky_low = squeeze(   ktraj_adc_epi_low(2,kx_ind,:,slice_ind,:,end-nImagingFrame+1:end)   );
ky_low = round(  ky_low/(1/FOV(2))  );
ky_seg1 = ky_low(1:Ny_low_meas);
ky_shift_to_center_low = abs(ky_seg1(end)) - Ny_actual_low/2;
end

% figure;
% plot(ky(:),'-o');
% hold on;
% yline(0,'m--')
% xlabel 'index'
% ylabel 'ky'




%--------------------------------------------------------------------------
%% SVD coil compression and espirit coil sensitivity estimation
%--------------------------------------------------------------------------
if isUseSensFromB0
    

elseif isUseGreACS
    if size(kspace_acs,4)/num_virtual_chan >2
        msgbox('Coil compression might be aggressive')
    end

    if num_virtual_chan==size(kspace_acs,4)% no comression
        kspace_acs_svd = kspace_acs;
        cmp_mtx = eye(num_virtual_chan);
    else
        if strcmp(coil_compression_method,'svd')
            [kspace_acs_svd, cmp_mtx] = svd_compress3d(kspace_acs, num_virtual_chan, 1);           
        else
            error(' ');
        end
    end

    if exist('noise_data','var') && ~strcmp(coil_compression_method,'gcc') && ~isempty(noise_data)
        noise_data_compressed = noise_data*cmp_mtx;

        % Demean the compressed noise data
        noise_mean_compressed = mean(noise_data_compressed, 1);
        temp_compressed = noise_data_compressed - noise_mean_compressed;

        % Calculate the standard, unbiased noise covariance matrix.
        % This is the most direct and theoretically sound way.
        noise_matrix_compressed = cov(temp_compressed);

        % Now proceed to calculate sigma_noise from this standard covariance matrix.
        channel_variances = diag(noise_matrix_compressed);
        mean_variance = mean(channel_variances);
        sigma_noise = sqrt(mean_variance);
    end

    rmse_acs = rmse(rsos(kspace_acs_svd,4), rsos(kspace_acs,4))

    sz = size(data_epi_imaging);
    nc = sz(end);
    if strcmp(coil_compression_method,'svd')
        tmp = reshape(data_epi_imaging,sz(1),sz(2),[],nc);
        data_epi_imaging_svd = svd_apply3d(tmp, cmp_mtx);
    else
        error();
    end
    data_epi_imaging_svd = reshape(data_epi_imaging_svd,[sz(1:end-1),num_virtual_chan]);


    sz = size(data_epi_refscan);
    tmp = reshape(data_epi_refscan,sz(1),sz(2),[],nc);
    data_epi_refscan_svd = svd_apply3d(tmp, cmp_mtx);
    data_epi_refscan_svd = reshape(data_epi_refscan_svd,[sz(1:end-1),num_virtual_chan]);

    if is2ndEcho
        sz = size(data_epi_imaging_low);
        tmp = reshape(data_epi_imaging_low,sz(1),sz(2),[],nc);
        data_epi_imaging_low_svd = svd_apply3d(tmp, cmp_mtx);
        data_epi_imaging_low_svd = reshape(data_epi_imaging_low_svd,[sz(1:end-1),num_virtual_chan]);
    end

    % zero pad acs to epi's matrix size
    kspace_gre_svd_zpad = zpad(kspace_acs_svd,size(kspace_acs_svd,1),Ny_actual,Nslices,num_virtual_chan);
    img_acs_pad = ifft2call(kspace_gre_svd_zpad);

    num_acs = 32;
    kernel_size = [6,6];
    eigen_thresh = 0.8;

    sens = zeross(size(img_acs_pad));

    % c = parcluster('local');    % build the 'local' cluster object

    % total_cores = c.NumWorkers;
    % parpool(ceil(total_cores/2))

    tic
    parfor slc_select = 1:size(img_acs_pad,3)
        % disp(num2str(slc_select))

        [maps, weights] = ecalib_soft( fft2c( sq(img_acs_pad(:,:,slc_select,:)) ), num_acs, kernel_size, eigen_thresh );

        sens(:,:,slc_select,:) = permute(dot_mult(maps, weights >= eigen_thresh ), [1,2,4,3]);
    end
    toc

    % delete(gcp)

    % save(strcat(data_path, 'sens_', file_name_epi, '.mat'), 'sens', '-V7.3');
    sens_orig = sens;
    sens = crop(sens, [Nx,Ny_actual,Nslices,num_virtual_chan]);
    sens = sens(:,:,ind_undo_interleave_epi,:);

    sens = flip(sens,2); warning('flip sens')




    if is2ndEcho                                    
    kspace_gre_svd_zpad = zpad(kspace_acs_svd,size(kspace_acs_svd,1),Ny_actual_low,Nslices,num_virtual_chan);
    % kspace_gre_svd_zpad = crop(kspace_acs_svd,Nx_low_org,Ny_actual_low,Nslices,num_virtual_chan);
    img_acs_pad = ifft2call(kspace_gre_svd_zpad);
    sens_low = zeross(size(img_acs_pad));
    parfor slc_select = 1:size(img_acs_pad,3)
        % disp(num2str(slc_select))

        [maps, weights] = ecalib_soft( fft2c( sq(img_acs_pad(:,:,slc_select,:)) ), num_acs, kernel_size, eigen_thresh );

        sens_low(:,:,slc_select,:) = permute(dot_mult(maps, weights >= eigen_thresh ), [1,2,4,3]);
    end
    sens_low = crop(sens_low, [Nx_low_org,Ny_actual_low,Nslices,num_virtual_chan]);
    sens_low = sens_low(:,:,ind_undo_interleave_epi,:);

    sens_low = flip(sens_low,2); warning('flip sens')
    end


    % get a 32*32 sens map
    Nx_low_spi = seq.getDefinition('Nx_low_spi');
    kspace_gre_svd_tmp = crop(kspace_acs_svd,Nx_low_spi,Nx_low_spi,Nslices,num_virtual_chan);
    sz = size(kspace_gre_svd_tmp);
    wind1 = hamming(sz(1));
    wind2 = wind1(:)*wind1(:).';
    kspace_gre_svd_tmp = kspace_gre_svd_tmp.*wind2;
    img_acs_pad = ifft2call(kspace_gre_svd_tmp);
    sens_low_for_spi = zeross(size(img_acs_pad));
    for slc_select = 1:size(img_acs_pad,3)
        [maps, weights] = ecalib_soft( fft2c( sq(img_acs_pad(:,:,slc_select,:)) ), num_acs, kernel_size, eigen_thresh );

        sens_low_for_spi(:,:,slc_select,:) = permute(dot_mult(maps, weights >= eigen_thresh ), [1,2,4,3]);
    end
    sens_low_for_spi = sens_low_for_spi(:,:,ind_undo_interleave_epi,:);
    sens_low_for_spi = flip(sens_low_for_spi,2); warning('flip sens')



    if numel(TE_gre)>1
        delta_TE_gre = diff(TE_gre);
        kspace_acs = reshape(data_acs, [Nreadoutlen,numel(TE_gre), Ny_lines_acs, Nslices, num_chan]);
        kspace_acs = permute(kspace_acs,[1 3 2 4 5]);


        sz = size(kspace_acs);
        tmp = reshape(kspace_acs,sz(1),sz(2),[],nc);
        kspace_acs_svd = svd_apply3d(tmp, cmp_mtx);
        kspace_acs = reshape(kspace_acs_svd,[sz(1:end-1),num_virtual_chan]);



        sz = size(kspace_acs);
        sz(2) = Ny_actual;
        kspace_acs_fieldmapping = zpad(kspace_acs,sz);
        kspace_acs_fieldmapping = kspace_acs_fieldmapping(:,:,:,ind_undo_interleave_epi,:);


        gre_img = ifft2call(kspace_acs_fieldmapping);
        sz = size(gre_img);
        sz(1) = Nx;
        gre_img = crop(gre_img,sz);
        gre_img = flip(gre_img,2);
        gre_img = permute(gre_img,[1 2 4 5 3]);
        gre_img = sum(gre_img.*conj(sens),4);

        % phase_diff = angle(  gre_img(:,:,:,:,2).*conj(gre_img(:,:,:,:,1)) );
        % phase_diff = squeeze(phase_diff);
        % b0map_gre_2echo_pulseq = phase_diff/2/pi/delta_TE_gre;


        %%
        sz = size(sens);
        b0map_gre_2echo_pulseq = zeros( sz(1:end-1) );
        for iter_slice = 1:size(gre_img,3)
            mask = logical(  sos(  sens(:,:,iter_slice,:)  ));
            ph = angle( squeeze( gre_img(:,:,iter_slice,1,:))  );
            [ b0map_gre_2echo_pulseq(:,:,iter_slice),fit_error,phase ] = dB_fitting_JumpCorrect( ph ,TE_gre,mask,1);
        end

        %%
        if numel(TE_gre)>2
            if isGEscanner
                B0 = 3;
            else
                B0=2.89;
            end
            method = 'wenmiao'; % wenmiao, hernando_graphcut

            sz = size(gre_img);
            water_image_from_gre = zeros(sz(1:3));
            fat_image_from_gre = water_image_from_gre;
            field_map_from_gre = water_image_from_gre;
            if isGEscanner
                PrecessionIsClockwise__ = -1;
            else
                PrecessionIsClockwise__ = 1;
            end
            parfor iter_slice = 1:size(gre_img,3)
                fprintf('dixon for gre slice %d/%d\n', iter_slice, size(gre_img,3));
                data = [];
                data.TE = TE_gre;
                data.FieldStrength = B0;
                data.PrecessionIsClockwise = 1;
                data.images = gre_img(:,:,iter_slice,:,:);   % the actual data, array of dimensions nx * ny * nz * ncoils * nte
                data.DO_MIXED_FIT = 0;
                switch method
                    case 'hernando_graphcut'
                        outParams = func_fw_i2cm1i_3pluspoint_hernando_graphcut(data);
                    case 'MultiSeedRegionGrowing'
                        outParams = func_MultiSeedRegionGrowing(data); % very bad
                    case 'wenmiao'
                        outParams = func_Synthetic_wenmiao_111219(data);
                    case 'RGideal'
                        sz = size(data.images);
                        data.images = zpad(data.images,[N N sz(3:end)]);
                        outParams = func_fw_i2cm0c_3pluspoint_RGmulticoil(data);
                    otherwise
                        error()
                end
                water_image_from_gre(:,:,iter_slice) = outParams.species(1).amps;
                fat_image_from_gre(:,:,iter_slice) = outParams.species(2).amps;
                field_map_from_gre(:,:,iter_slice) = outParams.fieldmap;
            end
            disp('finished dixon for gre')
        end


    end
    rmse_epi = rmse(rsos(data_epi_imaging_svd,6), rsos(data_epi_imaging,6))
    rmse_ref = rmse(rsos(data_epi_refscan_svd,5), rsos(data_epi_refscan,5))
else
    data_epi_imaging_svd = data_epi_imaging;
    data_epi_refscan_svd = data_epi_refscan;
    num_virtual_chan = num_chan;
end


%--------------------------------------------------------------------------
%% apply ramp sampling correction
%--------------------------------------------------------------------------
slc_ind = 1;
kx_adc_ref_seq = ktraj_adc_epi_refonly_high(1,:,:,slc_ind);
kx_adc_ref_seq = reshape(kx_adc_ref_seq, Nreadoutlen,Ny_meas);

kxmin=min(kx_adc_ref_seq(:));
kxmax=max(kx_adc_ref_seq(:));
kxmax1=kxmax/(Nx/2-1)*(Nx/2); % this compensates for the non-symmetric center definition in FFT
kmaxabs=max(kxmax1, -kxmin);

kxx= ((-Nx/2):(Nx/2-1))/(Nx/2)*kmaxabs; % kx-sample positions

if is2ndEcho
    kx_adc_ref_seq_low = ktraj_adc_epi_refonly_low(1,:,:,slc_ind);
    kx_adc_ref_seq_low = reshape(kx_adc_ref_seq_low,numSamplesLow,Ny_low_meas);
    kxmin=min(kx_adc_ref_seq_low(:));
    kxmax=max(kx_adc_ref_seq_low(:));
    kxmax1=kxmax/(Nx_low_org/2-1)*(Nx_low_org/2); % this compensates for the non-symmetric center definition in FFT
    kmaxabs=max(kxmax1, -kxmin);
    kxx_low= ((-Nx_low_org/2):(Nx_low_org/2-1))/(Nx_low_org/2)*kmaxabs; % kx-sample positions
end


data_epi_regrid = zeross([Nx,Ny_meas*Nslices*RSegment*nImagingFrame,num_virtual_chan]);
data_ref_regrid = zeross([Nx,Ny_meas*Nslices*1*1            ,num_virtual_chan]);
data_epi_regrid_low = zeross([Nx_low_org,Ny_low_meas*Nslices*RSegment*nImagingFrame,num_virtual_chan]);

tmp_data_epi_imaging_svd = reshape(data_epi_imaging_svd,Nreadoutlen,Ny_meas*Nslices*RSegment*nImagingFrame,num_virtual_chan);
tmp_data_epi_refscan_svd = reshape(data_epi_refscan_svd,Nreadoutlen,Ny_meas*Nslices*1*1            ,num_virtual_chan);
ktraj_use1 = repmat(kx_adc_ref_seq,[1,Nslices*RSegment*nImagingFrame]);
for cc = 1:num_virtual_chan
    for aa = 1:Ny_meas*Nslices*RSegment*nImagingFrame
        data_epi_regrid(:,aa,cc) = interp1(ktraj_use1(:,aa), tmp_data_epi_imaging_svd(:,aa,cc), kxx, 'spline', 0);
    end
end
ktraj_use1 = repmat(kx_adc_ref_seq,[1,Nslices]);
for cc = 1:num_virtual_chan
    for aa = 1:Ny_meas*Nslices
        data_ref_regrid(:,aa,cc) = interp1(ktraj_use1(:,aa), tmp_data_epi_refscan_svd(:,aa,cc), kxx, 'spline', 0);
    end
end
data_epi_regrid = reshape(data_epi_regrid,Nx,Ny_meas,Nslices,RSegment,nImagingFrame,num_virtual_chan);
data_ref_regrid = reshape(data_ref_regrid,Nx,Ny_meas,Nslices,1,              num_virtual_chan);

data_epi_regrid = data_epi_regrid(:,:,ind_undo_interleave_epi,:,:,:);
data_ref_regrid = data_ref_regrid(:,:,ind_undo_interleave_epi,:,:);


img_epi_regrid = ifft2call(data_epi_regrid);
img_ref_regrid = ifft2call(data_ref_regrid);






if is2ndEcho
    tmp_data_epi_imaging_low_svd = reshape(data_epi_imaging_low_svd,numSamplesLow,Ny_low_meas*Nslices*RSegment*nImagingFrame,num_virtual_chan);
    ktraj_use1_low = repmat(kx_adc_ref_seq_low,[1,Nslices*RSegment*nImagingFrame]);
    for cc = 1:num_virtual_chan
        for aa = 1:Ny_low_meas*Nslices*RSegment*nImagingFrame
            data_epi_regrid_low(:,aa,cc) = interp1(ktraj_use1_low(:,aa), tmp_data_epi_imaging_low_svd(:,aa,cc), kxx_low, 'spline', 0);
        end
    end
    data_epi_regrid_low = reshape(data_epi_regrid_low,Nx_low_org,Ny_low_meas,Nslices,RSegment,nImagingFrame,num_virtual_chan);
    data_epi_regrid_low = data_epi_regrid_low(:,:,ind_undo_interleave_epi,:,:,:);
    img_epi_regrid_low = ifft2call(data_epi_regrid_low);
end



% imagesc3d2( rsos(img_acs_pad,4), s(img_acs_pad)/2, 11, [90,90,-90], [0,1.5e-3], [], 'gre acs pad'), setGcf(.5)
% imagesc3d2( rsos(img_epi_regrid,4), s(img_epi_imaging)/2, 12, [90,90,-90], [0,1.5e-3], [], 'regridded epi'), setGcf(.5)
% imagesc3d2( rsos(img_epi_imaging,4), s(img_epi_imaging)/2, 13, [90,90,-90], [0,1.5e-3], [], 'before regridding'), setGcf(.5)


%--------------------------------------------------------------------------
%% k-space based ghost correction
%--------------------------------------------------------------------------

if isEnableGhostCorrection&&isRefscan
    % if isConventionalSegmentSpacing
    %     AccY = Ry*RSegment;
    %     AccY_low = Ry_low*RSegment;
    % else
    %     AccY = Ry;
    %     AccY_low = Ry_low;
    % end

    AccY = 1;
    AccY_low = 1;

    navCor_kspace = 1*data_epi_regrid;
    navCor_kspace_low = 1 * data_epi_regrid_low;
    navCor_ref = data_ref_regrid;

    switch ghost_cor_method
        case 'phzShift'
            parfor iter_diff_frame = 1:nImagingFrame
                for iter_slice = 1:Nslices
                    for iter_seg = 1:RSegment
                        s=squeeze(data_ref_regrid(:,2:4,iter_slice,1,:));
                        k_1 = squeeze(data_epi_regrid(:,:,iter_slice,iter_seg,iter_diff_frame,:));

                        S0 = fif(mean(s(:,[1 3],:),2));
                        S1 = fif(mean(s(:,[ 2 ],:),2));
                        y = zeros(2, size(S0,3));
                        for cnt=1:size(S0,3)
                            [~,y(:,cnt)]=phzshift( S1(:,:,cnt).', S0(:,:,cnt).', {'nofft'});
                        end
                        tmp = phzapply( permute(k_1(:,1:end,:),[2 1 3]), y);
                        navCor_kspace(:,:,iter_slice,iter_seg,iter_diff_frame,:) = permute(tmp,[2 1 3]);
                    end
                end
            end

            if is2ndEcho

                parfor iter_diff_frame = 1:nImagingFrame
                    for iter_slice = 1:Nslices
                        for iter_seg = 1:RSegment
                            s=squeeze(data_ref_regrid(:,2:4,iter_slice,1,:));
                            k_1 = squeeze(data_epi_regrid_low(:,:,iter_slice,iter_seg,iter_diff_frame,:));

                            S0 = fif(mean(s(:,[1 3],:),2));
                            S1 = fif(mean(s(:,[ 2 ],:),2));
                            y = zeros(2, size(S0,3));
                            for cnt=1:size(S0,3)
                                [~,y(:,cnt)]=phzshift( S1(:,:,cnt).', S0(:,:,cnt).', {'nofft'});
                            end
                            tmp = phzapply( permute(k_1(:,1:end,:),[2 1 3]), y);
                            navCor_kspace_low(:,:,iter_slice,iter_seg,iter_diff_frame,:) = permute(tmp,[2 1 3]);
                        end
                    end
                end

            end
        otherwise 
            error()
    end
    navCor_img = sq((ifft2call(navCor_kspace)));
    navCor_img_low = sq((ifft2call(navCor_kspace_low)));

    % imagesc3d2( sq(rsos(navCor_img,3)), [Nx, mtx_size(2:3)]./[1,Ry,1]/2, 4, [90,90,-90], [0,1.5e-3], [], 'regrid & ghost corr'), setGcf(.5)
    % imagesc3d2( rsos(img_epi_regrid,4), s(img_epi_imaging)/2, 5, [90,90,-90], [0,1.5e-3], [], 'regridded epi'), setGcf(.5)

else
    navCor_kspace = data_epi_regrid;
end



%% re-arrange the k-space to align different echoes, segments, diffusion frames
ksp_full_shifted = zeros(Nx,Ny_actual,Nslices,Necho*RSegment,nImagingFrame,num_virtual_chan);
% all_kyShiftFactors_frame;

for iter_slice = 1:Nslices
    for iter_diff_frame = 1:nImagingFrame
        kyShiftFactor_frame = all_kyShiftFactors_frame(iter_diff_frame);


        fig_num = 10;

        gen_ind = @(x0,step,num_elements) x0 + (0:num_elements-1) * step;

        if isConventionalSegmentSpacing
            ky_dist_sameEcho_sameSeg = Necho*Ry*RSegment;
        else
            ky_dist_sameEcho_sameSeg = Necho*Ry;
        end

        switch Necho
            case 1
                tmp = navCor_kspace(:,1:Necho:end,iter_slice,:,iter_diff_frame,:);
                ind = gen_ind(1,ky_dist_sameEcho_sameSeg,size(tmp,2));
                tmpksp = zeros(Nx,max(ind),RSegment,num_virtual_chan);
                tmpksp(:,ind,:,:) = squeeze(tmp);
                pre_pad = Ny_actual-size(tmpksp,2);
                ksp_echo1 = padarray(tmpksp,[0 pre_pad],'pre');

                ksp_shifted_echo1 = circshift(ksp_echo1,[0,ky_shift_echo1_aligned]);

                ksp_shifted_echo1 = circshift(ksp_shifted_echo1,[0,ky_shift_to_center]);

                % undo shift between different segments
                for iter_seg = 1:RSegment
                    shift =  (iter_seg-1)*base_shift_between_seg;
                    ksp_shifted_echo1(:,:,iter_seg,:) = circshift( ksp_shifted_echo1(:,:,iter_seg,:), [0 shift] );
                end

                ksp_allSeg_allEcho = ksp_shifted_echo1;
            case 2
                tmp = navCor_kspace(:,1:Necho:end,iter_slice,:,iter_diff_frame,:);
                ind = gen_ind(1,ky_dist_sameEcho_sameSeg,size(tmp,2));
                tmpksp = zeros(Nx,max(ind),RSegment,num_virtual_chan);
                tmpksp(:,ind,:,:) = squeeze(tmp);
                pre_pad = Ny_actual-size(tmpksp,2);
                ksp_echo1 = padarray(tmpksp,[0 pre_pad],'pre');

                tmp = navCor_kspace(:,2:Necho:end,iter_slice,:,iter_diff_frame,:);
                ind = gen_ind(1,ky_dist_sameEcho_sameSeg,size(tmp,2));
                tmpksp = zeros(Nx,max(ind),RSegment,num_virtual_chan);
                tmpksp(:,ind,:,:) = squeeze(tmp);
                ksp_echo2 = padarray(tmpksp,[0 pre_pad],'pre');
                post_pad = Ny_actual - size(ksp_echo2,2);
                ksp_echo2 = padarray(ksp_echo2,[0 post_pad],'post');

                ksp_shifted_echo1 = circshift(ksp_echo1,[0,ky_shift_echo1_aligned]);
                ksp_shifted_echo2 = circshift(ksp_echo2,[0,ky_shift_echo2_aligned]);

                ksp_shifted_echo1 = circshift(ksp_shifted_echo1,[0,ky_shift_to_center]);
                ksp_shifted_echo2 = circshift(ksp_shifted_echo2,[0,ky_shift_to_center]);

                % undo shift between different segments
                for iter_seg = 1:RSegment
                    shift =  (iter_seg-1)*base_shift_between_seg;
                    ksp_shifted_echo1(:,:,iter_seg,:) = circshift( ksp_shifted_echo1(:,:,iter_seg,:), [0 shift] );
                    ksp_shifted_echo2(:,:,iter_seg,:) = circshift( ksp_shifted_echo2(:,:,iter_seg,:), [0 shift] );
                end

                ksp_allSeg_allEcho = cat(3,ksp_shifted_echo1,ksp_shifted_echo2);
            case 3
                tmp = navCor_kspace(:,1:Necho:end,iter_slice,:,iter_diff_frame,:);
                ind = gen_ind(1,ky_dist_sameEcho_sameSeg,size(tmp,2));
                tmpksp = zeros(Nx,max(ind),RSegment,num_virtual_chan);
                tmpksp(:,ind,:,:) = squeeze(tmp);
                pre_pad = Ny_actual-size(tmpksp,2);
                ksp_echo1 = padarray(tmpksp,[0 pre_pad],'pre');

                tmp = navCor_kspace(:,2:Necho:end,iter_slice,:,iter_diff_frame,:);
                ind = gen_ind(1,ky_dist_sameEcho_sameSeg,size(tmp,2));
                tmpksp = zeros(Nx,max(ind),RSegment,num_virtual_chan);
                tmpksp(:,ind,:,:) = squeeze(tmp);
                ksp_echo2 = padarray(tmpksp,[0 pre_pad],'pre');
                post_pad = Ny_actual - size(ksp_echo2,2);
                ksp_echo2 = padarray(ksp_echo2,[0 post_pad],'post');

                tmp = navCor_kspace(:,3:Necho:end,iter_slice,:,iter_diff_frame,:);
                ind = gen_ind(1,ky_dist_sameEcho_sameSeg,size(tmp,2));
                tmpksp = zeros(Nx,max(ind),RSegment,num_virtual_chan);
                tmpksp(:,ind,:,:) = squeeze(tmp);
                ksp_echo3 = padarray(tmpksp,[0 pre_pad],'pre');
                post_pad = Ny_actual - size(ksp_echo3,2);
                ksp_echo3 = padarray(ksp_echo3,[0 post_pad],'post');

                ksp_shifted_echo1 = circshift(ksp_echo1,[0,ky_shift_echo1_aligned]);
                ksp_shifted_echo2 = circshift(ksp_echo2,[0,ky_shift_echo2_aligned]);
                ksp_shifted_echo3 = circshift(ksp_echo3,[0,ky_shift_echo3_aligned]);

                ksp_shifted_echo1 = circshift(ksp_shifted_echo1,[0,ky_shift_to_center]);
                ksp_shifted_echo2 = circshift(ksp_shifted_echo2,[0,ky_shift_to_center]);
                ksp_shifted_echo3 = circshift(ksp_shifted_echo3,[0,ky_shift_to_center]);

                % undo shift between different segments
                for iter_seg = 1:RSegment
                    shift =  (iter_seg-1)*base_shift_between_seg;
                    ksp_shifted_echo1(:,:,iter_seg,:) = circshift( ksp_shifted_echo1(:,:,iter_seg,:), [0 shift] );
                    ksp_shifted_echo2(:,:,iter_seg,:) = circshift( ksp_shifted_echo2(:,:,iter_seg,:), [0 shift] );
                    ksp_shifted_echo3(:,:,iter_seg,:) = circshift( ksp_shifted_echo3(:,:,iter_seg,:), [0 shift] );
                end

                ksp_allSeg_allEcho = cat(3,ksp_shifted_echo1,ksp_shifted_echo2,ksp_shifted_echo3);

            otherwise
                error('undefined');
        end

        % undo shift between diff volumes
        ksp_full_shifted(:,:,iter_slice,:,iter_diff_frame,:) = circshift(ksp_allSeg_allEcho,[0 kyShiftFactor_frame]);
    end
end


ksp_full_shifted_echoSegSeparated = reshape(ksp_full_shifted,Nx,Ny_actual,Nslices,RSegment,Necho,nImagingFrame,num_virtual_chan);
clear ksp_full_shifted


%% 2nd echo recon
if is2ndEcho
if isConventionalSegmentSpacing
    dist = Rynmapsnmaps_low*RSegment;
else
    dist = Ry_low;
end
ind = gen_ind(1,dist,Ny_low_meas);
tmpksp = zeros(Nx_low_org, max(ind),Nslices,RSegment,nImagingFrame, num_virtual_chan);
tmpksp(:,ind,:,:,:,:) = navCor_kspace_low;
post_pad = Ny_actual_low-size(tmpksp,2);
tmpksp = padarray(tmpksp,[0 post_pad],'post'); % note that in the acq part, partial fourier is in the other side of the kspace

% undo shift between different segments
ksp_shifted_low = 0*tmpksp;
for iter_seg = 1:RSegment
    shift =  1 * (iter_seg-1)*Ry_low;
    ksp_shifted_low(:,:,:,iter_seg,:,:) = circshift( tmpksp(:,:,:,iter_seg,:,:), [0 shift] );
end

% use bart to recon low res
sz = size(ksp_shifted_low);
img_low_water = zeros(sz(1:end-1));
img_low_fat = img_low_water;
img_low_water_fat_combined = img_low_water;
isSoftSense = 0;
if isSoftSense
    nmaps = 2; % see bart ecalib -h,   -m maps         % Number of maps to compute.
else
    nmaps = 1;
end
phase_diff_calc_method = 'direct';% esnails, direct
ecalib_tool = 'matlab'; % matlab, bart

assert( ~isSoftSense || ~strcmp(phase_diff_calc_method, 'direct'), ...
    'SoftSense requires phase_diff_calc_method = ''esnails''' )
sens_nav = zeros(Nx,Ny_actual,Nslices,num_virtual_chan,nmaps,RSegment,nImagingFrame);
% parfor iter_slice = 1:Nslices
for iter_slice = 21%1:Nslices
    for iter_diff_frame = 1:nImagingFrame
        imgs_nav_temp = zeros(size(sens_low,1), size(sens_low,2), RSegment);
        for iter_seg = 1:RSegment
            ksp_shifted_low_tmp = squeeze(ksp_shifted_low(:,:,iter_slice,iter_seg,iter_diff_frame,:));
            ksp_shifted_low_tmp = reshape(ksp_shifted_low_tmp, [size(ksp_shifted_low_tmp,1), size(ksp_shifted_low_tmp,2), 1, size(ksp_shifted_low_tmp,3)]);
            if 1 % simple sense
                imnav = bart('pics -S -R Q:1e-2', ksp_shifted_low_tmp, sens_low(:,:,iter_slice,:)  );
                imgs_nav_temp(:,:,iter_seg) = imnav;
                img_low_water(:,:,iter_slice,iter_seg,iter_diff_frame) = imnav;
            else % sense, also with fat
                params = [];
                params.ESP = EchoSpacing;
                gyr = 42.576e6;
                fatChemShift_ppm = 3.5e-6; % ppm
                water_fat_freq_diff_Hz = fatChemShift_ppm * B0 * gyr; % Fat frequency difference in Hz
                delta_phi_per_line = 2 * pi * water_fat_freq_diff_Hz * EchoSpacing/Ry;
                % phase_ramp_vector_pe = exp(-1i * (1:Ny_actual) * delta_phi_per_line);
                phase_ramp_vector_pe = exp(-1i * (1:size(sens_low,2)) * delta_phi_per_line);
                params.phase_mod_kspace_fat = repmat(phase_ramp_vector_pe,size(sens_low,1),1);
                params.masks = ksp_shifted_low_tmp~=0;

                params.N1 = size(sens_low ,1); % Assuming Nx is PE direction
                params.N2 = size(sens_low ,2); % Assuming Ny is RO direction
                params.sens =  sens_low(:,:,iter_slice,:) ;

                lambda_tik = 0.001;

                Afor = @(in) apply_forward_adjoint_kspace_model2(in, 'notransp', params );
                Aadj = @(in) apply_forward_adjoint_kspace_model2(in, 'transp', params );
                AhA=@(in)Aadj(Afor(in));                
                Ahd = Aadj( ksp_shifted_low_tmp(:) );

                M_operator = @(in_x_vec) Aadj(Afor(in_x_vec)) + lambda_tik * in_x_vec;

                res = pcg(M_operator, Ahd);
                tmp1= reshape(res,size(sens_low,1),size(sens_low,2),2);


                if 0
                    imgs_nav_temp(:,:,iter_seg) = tmp1(:,:,1);
                else
                    imgs_nav_temp(:,:,iter_seg) = sum(tmp1,3);
                end
                img_low_water(:,:,iter_slice,iter_seg,iter_diff_frame)              = tmp1(:,:,1);
                img_low_fat(:,:,iter_slice,iter_seg,iter_diff_frame)                = tmp1(:,:,2);
                img_low_water_fat_combined(:,:,iter_slice,iter_seg,iter_diff_frame) = sum(tmp1,3);              
            end
        end



        % % esnails
        if strcmp(phase_diff_calc_method,'esnails')
            tmpim = squeeze( imgs_nav_temp );
            tmpim_sens = tmpim .* sens_low(:,:,iter_slice,:);
            sz = size(tmpim_sens);
            tmpim_sens = reshape(tmpim_sens,sz(1),sz(2),1,sz(3)*sz(4));
            tmpksp = xyrecon.fft2c(tmpim_sens);

            wind1 = hamming(sz(1));
            wind2 = wind1(:)*wind1(:).';
            tmpksp = tmpksp.*wind2;

            tmpksp = bart(sprintf('resize -c 0 %d 1 %d',Nx,Ny_actual),tmpksp);
            if strcmp(ecalib_tool,'bart')
                cmd = sprintf('ecalib -d5  -m %d -c 0.01 -r 24', nmaps);
                [ecal_spiral,eigvals] = bart(cmd,tmpksp);
                %
                % ecal_spiral = bart('slice 4 0', ecal_spiral);
            elseif strcmp(ecalib_tool,'matlab')
                num_acs = 24;
                kernel_size = [6,6];
                % kernel_size = [8,8];
                eigen_thresh = 0.8;			% for mask size
                [maps, weights] = ecalib_soft_hack( squeeze(tmpksp), num_acs, kernel_size, eigen_thresh );
                [ecal_spiral] = dot_mult(maps, weights >= eigen_thresh);
                eigvals = weights;
            else
                error('error')
            end
            ecal_spiral = reshape(ecal_spiral,Nx,Ny_actual,1,RSegment,num_virtual_chan,nmaps);
            ecal_spiral = permute(ecal_spiral,[1 2 3 5 6 4]);
            sens_nav(:,:,iter_slice,:,:,:,iter_diff_frame) = ecal_spiral .* logical( sos(sens(:,:,iter_slice,:),4) );


        elseif strcmp(phase_diff_calc_method,'direct')
            tmpksp = xyrecon.fft2c(imgs_nav_temp);

            sz = size(tmpksp);
            wind1 = hamming(sz(1));
            wind2 = hamming(sz(2));
            wind  = wind1 * wind2.';
            tmpksp = tmpksp.*wind;

            % tmpksp = xyrecon.zpad(tmpksp,Nx,Ny_actual,RSegment);
            tmpksp = bart(sprintf('resize -c 0 %d 1 %d',Nx,Ny_actual),tmpksp);
            im = xyrecon.ifft2c(tmpksp);
            if 0
                %                             phs_diff = im .* conj(  im(:,:,1)  );
                phs_diff = conj(im) .* im(:,:,1);
                phs_diff = exp(1i*angle(phs_diff));
            else
                phs_diff = 0 * im;
                for iter_diff_calc = 1:size(im,3)
                    phs_diff(:,:,iter_diff_calc) = exp(1i * medfilt2(angle( conj(im(:,:,iter_diff_calc)) .* im(:,:,1) ),[15,15],'symmetric'));
                end
            end


            phs_diff = reshape(phs_diff,Nx,Ny_actual,1,1,1,RSegment);


            sens_nav(:,:,iter_slice,:,:,:,iter_diff_frame) = phs_diff.*sens(:,:,iter_slice,:);

        else
            error()
        end




    end
end
disp('nav calc finished')
end



%% get TE for the echoes
if ismatrix(ky)
    sz = size(ky);
    ky = reshape(ky,sz(1),RSegment,sz(2));
end
sz = size(ky);
assert(sz(2)==RSegment)
assert(sz(3)==nImagingFrame)
sz(1) = Necho;
dist_start_of_readout_to_center = zeros(sz);
for iter2 = 1:sz(2)
    for iter3 = 1:sz(3)
        dist_start_of_readout_to_center(:,iter2,iter3)=get_dur_to_center(ky(:,iter2,iter3),Necho,EchoSpacing);
    end
end

TE = seq.getDefinition('EchoTime');
assert(~isempty(TE));

tmp_distance_rf180_readout = distance_rf180_readout(:).';
tmp_distance_rf180_readout = repmat(tmp_distance_rf180_readout,Necho,1,sz(3));warning('tmp');


TEs_for_dixon =  tmp_distance_rf180_readout + dist_start_of_readout_to_center - TE/2;




%%
ky_allSeg_allDiffFrames = reshape(ky,Ny_meas,RSegment,nImagingFrame);
isSingleChannelLoraks = 1;


recon_strategy = 'perSegPerDir'; % perSegPerDir, jointOverSegments, jointOverSegmentsAndDirections



if isGEscanner
    B0 = 3.0;
else
    B0 = 2.89;
end


if ~isempty(B0map_fn)
    t = load(B0map_fn);
    B0map = t.B0map;
    B0map = imresize3(B0map,[Nx,Ny_actual,Nslices]);
    B0map_use = B0map;
end

b0map_gre_2echo_pulseq_filtered = 0*b0map_gre_2echo_pulseq;
for iter_slice = 1:Nslices
    b0map_gre_2echo_pulseq_filtered(:,:,iter_slice) = imgaussfilt(b0map_gre_2echo_pulseq(:,:,iter_slice),3);
end

B0map_use = b0map_gre_2echo_pulseq_filtered;


imgs_vc_mussels = zeros(Nx,Ny_actual,Nslices,Necho,RSegment,nImagingFrame);
imgs_S_loraks = imgs_vc_mussels;
imgs_S_loraks_single_channel = imgs_vc_mussels;
imgs_sep_recon = imgs_S_loraks;
imgs_MUSE = imgs_S_loraks_single_channel;


bval_for_imaging_frame = bval(end-nImagingFrame+1:end);
diff_vol_ind_b0 = bval_for_imaging_frame<=bvalue0_threshold;
diff_vol_ind_bdiff = ~diff_vol_ind_b0;
img_ind_b0 = 1:sum(diff_vol_ind_b0);
img_ind_bdiff = sum(diff_vol_ind_b0)+1:nImagingFrame;

temp_b0_mussels = cell(Nslices,RSegment);
temp_bdiff_mussels = temp_b0_mussels;
temp_b0_sloraks = temp_b0_mussels;
temp_bdiff_sloraks = temp_b0_mussels;



water_oneByOne = zeros(Nx,Ny_actual,Nslices,numel(unique(bval)));
fat_oneByOne = water_oneByOne;
fieldmap_oneByOne = water_oneByOne;

water_avgOverSeg = water_oneByOne;
fat_avgOverSeg = water_oneByOne;
fieldmap_avgOverSeg = water_oneByOne;

water_avgOverSegOverDirection = water_oneByOne;
fat_avgOverSegOverDirection = water_oneByOne;
fieldmap_avgOverSegOverDirection = water_oneByOne;


water_avgOverDirection = water_oneByOne;
fat_avgOverDirection = water_oneByOne;
fieldmap_avgOverDirection = water_oneByOne;

water_individualTE_avgOverSameTE = water_oneByOne;
fat_individualTE_avgOverSameTE = water_oneByOne;
fieldmap_individualTE_avgOverSameTE = water_oneByOne;

water_individualTE_jointSeg_jointBvol = water_oneByOne;
fat_individualTE_jointSeg_jointBvol =  water_oneByOne;
fieldmap_individualTE_jointSeg_jointBvol =  water_oneByOne;

water_individualTE_iterSeg_iterBvol = zeros(Nx,Ny_actual,Nslices,RSegment,nImagingFrame);
fat_individualTE_iterSeg_iterBvol = zeros(Nx,Ny_actual,Nslices,RSegment,nImagingFrame);
fieldmap_individualTE_iterSeg_iterBvol = zeros(Nx,Ny_actual,Nslices,RSegment,nImagingFrame);

water_individualTE_jointSeg_iterBvol = zeros(Nx,Ny_actual,Nslices,nImagingFrame);
fat_individualTE_jointSeg_iterBvol = zeros(Nx,Ny_actual,Nslices,nImagingFrame);
fieldmap_individualTE_jointSeg_iterBvol = zeros(Nx,Ny_actual,Nslices,nImagingFrame);

water_individualTE_b0fmPrior_iterSeg_iterBvol = zeros(Nx,Ny_actual,Nslices,RSegment,nImagingFrame);
fat_individualTE_b0fmPrior_iterSeg_iterBvol = zeros(Nx,Ny_actual,Nslices,RSegment,nImagingFrame);
fieldmap_individualTE_b0fmPrior_iterSeg_iterBvol = zeros(Nx,Ny_actual,Nslices,RSegment,nImagingFrame);
resnorm_individualTE_b0fmPrior_iterSeg_iterBvol = nan(Nx,Ny_actual,Nslices,RSegment,nImagingFrame);
residual_rms_individualTE_b0fmPrior_iterSeg_iterBvol = nan(Nx,Ny_actual,Nslices,RSegment,nImagingFrame);

water_individualTE_b0fmPrior_jointSeg_iterBvol = zeros(Nx,Ny_actual,Nslices,nImagingFrame);
fat_individualTE_b0fmPrior_jointSeg_iterBvol = zeros(Nx,Ny_actual,Nslices,nImagingFrame);
fieldmap_individualTE_b0fmPrior_jointSeg_iterBvol = zeros(Nx,Ny_actual,Nslices,nImagingFrame);
resnorm_individualTE_b0fmPrior_jointSeg_iterBvol = nan(Nx,Ny_actual,Nslices,nImagingFrame);
residual_rms_individualTE_b0fmPrior_jointSeg_iterBvol = nan(Nx,Ny_actual,Nslices,nImagingFrame);


if strcmp(recon_strategy,'perSegPerDir')
    seg_ind_ranges = num2cell(1:RSegment);
    diff_ind_ranges = num2cell(1:nImagingFrame);
elseif strcmp(recon_strategy,'jointOverSegments')
    seg_ind_ranges = {1:RSegment};
    diff_ind_ranges = num2cell(1:nImagingFrame);
elseif strcmp(recon_strategy,'jointOverSegmentsAndDirections')
    seg_ind_ranges = {1:RSegment};
    unqbval = unique(bval);
    diff_ind_ranges = cell(1,numel(unqbval));
    for iter_tmp = 1:numel(unqbval)
        diff_ind_ranges{iter_tmp} = find( unqbval(iter_tmp)==bval );
    end
else
    error('undefined');
end


num_iter_joint_grappa = 4;
% parfor iter_slice = 1:Nslices
isCompositeSens = 0;
for iter_slice = 13;%1:Nslices
    imgs_vc_mussels_slc = zeros(Nx,Ny_actual,Necho,RSegment,nImagingFrame);
    imgs_S_loraks_single_channel_slc = zeros(Nx,Ny_actual,Necho,RSegment,nImagingFrame);
    imgs_MUSE_slc = zeros(Nx,Ny_actual,Necho,RSegment,nImagingFrame);
    water_individualTE_iterSeg_iterBvol_slc = zeros(Nx,Ny_actual,RSegment,nImagingFrame);
    fat_individualTE_iterSeg_iterBvol_slc = zeros(Nx,Ny_actual,RSegment,nImagingFrame);
    fieldmap_individualTE_iterSeg_iterBvol_slc = zeros(Nx,Ny_actual,RSegment,nImagingFrame);
    water_individualTE_jointSeg_iterBvol_slc = zeros(Nx,Ny_actual,nImagingFrame);
    fat_individualTE_jointSeg_iterBvol_slc = zeros(Nx,Ny_actual,nImagingFrame);
    fieldmap_individualTE_jointSeg_iterBvol_slc = zeros(Nx,Ny_actual,nImagingFrame);
    water_individualTE_jointSeg_jointBvol_slc = zeros(Nx,Ny_actual,numel(unique(bval)));
    fat_individualTE_jointSeg_jointBvol_slc = zeros(Nx,Ny_actual,numel(unique(bval)));
    fieldmap_individualTE_jointSeg_jointBvol_slc = zeros(Nx,Ny_actual,numel(unique(bval)));
    water_individualTE_b0fmPrior_iterSeg_iterBvol_slc = zeros(Nx,Ny_actual,RSegment,nImagingFrame);
    fat_individualTE_b0fmPrior_iterSeg_iterBvol_slc = zeros(Nx,Ny_actual,RSegment,nImagingFrame);
    fieldmap_individualTE_b0fmPrior_iterSeg_iterBvol_slc = zeros(Nx,Ny_actual,RSegment,nImagingFrame);
    resnorm_individualTE_b0fmPrior_iterSeg_iterBvol_slc = nan(Nx,Ny_actual,RSegment,nImagingFrame);
    residual_rms_individualTE_b0fmPrior_iterSeg_iterBvol_slc = nan(Nx,Ny_actual,RSegment,nImagingFrame);
    water_individualTE_b0fmPrior_jointSeg_iterBvol_slc = zeros(Nx,Ny_actual,nImagingFrame);
    fat_individualTE_b0fmPrior_jointSeg_iterBvol_slc = zeros(Nx,Ny_actual,nImagingFrame);
    fieldmap_individualTE_b0fmPrior_jointSeg_iterBvol_slc = zeros(Nx,Ny_actual,nImagingFrame);
    resnorm_individualTE_b0fmPrior_jointSeg_iterBvol_slc = nan(Nx,Ny_actual,nImagingFrame);
    residual_rms_individualTE_b0fmPrior_jointSeg_iterBvol_slc = nan(Nx,Ny_actual,nImagingFrame);
    water_avgOverDirection_slc = zeros(Nx,Ny_actual,numel(unique(bval)));
    fat_avgOverDirection_slc = zeros(Nx,Ny_actual,numel(unique(bval)));
    fieldmap_avgOverDirection_slc = zeros(Nx,Ny_actual,numel(unique(bval)));
    water_avgOverSegOverDirection_slc = zeros(Nx,Ny_actual,numel(unique(bval)));
    fat_avgOverSegOverDirection_slc = zeros(Nx,Ny_actual,numel(unique(bval)));
    fieldmap_avgOverSegOverDirection_slc = zeros(Nx,Ny_actual,numel(unique(bval)));

    for iter_diff_frame_ind = 1:numel(diff_ind_ranges)

        iter_diff_frame = diff_ind_ranges{iter_diff_frame_ind};

        if bval(iter_diff_frame)<=bvalue0_threshold
            isb0 = true;
        else
            isb0 = false;
        end

        % get acs from b-value=0 volumes
        ind = bval<=bvalue0_threshold;
        ind_bvol = ind(1); % get from the 1st b-value=0 volume
        ksp_bval0_allSeg_allEcho =  squeeze(   ksp_full_shifted_echoSegSeparated(:,:,iter_slice,:,:,ind_bvol,:)  );
        ksp_acs_allEcho = squeeze(sum(ksp_bval0_allSeg_allEcho,3)); % sum over all segments (shots)

        % separate grappa for each echo
        num_acs = [24,24];
        kernel_size = [3,3];            % odd kernel size
        Rz = 1;
        lambda_tik = 1e-10;
        lambda_Tik = 1e-10;
        subs = 0;
        Img_Grappa_J_allseg = zeros(Nx,Ny_actual,num_virtual_chan,Necho,RSegment    );


        for iter_seg_ind = 1:numel(seg_ind_ranges)
            iter_seg = seg_ind_ranges{iter_seg_ind};

            ksp_allSeg_allEcho =  squeeze(   ksp_full_shifted_echoSegSeparated(:,:,iter_slice,:,:,iter_diff_frame,:)  );

            Delz_allEchoes = zeros(num_virtual_chan,Necho);
            Dely_allEchoes = Delz_allEchoes;
            indcies_1st_non_zero = zeros(1,Necho);



            ksp_for_mussels_loraks = squeeze(ksp_full_shifted_echoSegSeparated(:,:,iter_slice,iter_seg,:,iter_diff_frame,:)) ;
            sz=size(ksp_for_mussels_loraks );
            ksp_for_mussels_loraks = reshape(ksp_for_mussels_loraks,sz(1),sz(2),[],sz(end));

          
            %--------------------------------------------------------------------------
            %% S-loraks, single channel
            %--------------------------------------------------------------------------
            if isSingleChannelLoraks
                % cg_sense2 settings
                fig_num = 30;



                nx = Nx;  % nx of cEPI image
                % ny = Ny;  % ny of cEPI image
                ny = size(ksp_for_mussels_loraks,2);

                if 1 % old simple sens
                    params = [];
                    params.N1 = nx;
                    params.N2 = ny;
                    params.masks = ksp_for_mussels_loraks~=0;
                    params.Necho = size(ksp_for_mussels_loraks,3);                    
                    if isCompositeSens
                        error('temp, not supported due to parfor');
                    else
                        params.sens = sens(:,:,iter_slice,:); % no squeeze(), use bart format here!
                        A_for = @(in)apply_sense_nechoes_noreshape(in,'notransp',params);
                        A_adj = @(in)apply_sense_nechoes_noreshape(in,'transp',params);
                    end

                    AHA = @(in) A_adj(A_for(in));
                    AtA = AHA;
                    Ahd = A_adj(ksp_for_mussels_loraks);


                    if 0 % muse style recon
                        lambda_l2 = 0.01;
                        AHA_with_L2 = @(in) A_adj(A_for(in)) + lambda_l2 * in;
                        AHA_vec = @(x_vec) reshape( AHA_with_L2( reshape(x_vec, [nx, ny, params.Necho]) ), [], 1 );
                        [x_opt_vec, flag, relres, pcg_iter] = pcg(AHA_vec, Ahd(:));
                        img_recon = reshape(x_opt_vec, [nx, ny, params.Necho]);

                        sz = size(img_recon);
                        tmp = reshape(img_recon,sz(1),sz(2),numel(iter_seg),Necho,numel(iter_diff_frame));
                        tmp1 = permute( squeeze( tmp ), [1 2 4 3 5] );
                        imgs_MUSE_slc(:,:,:,iter_seg,iter_diff_frame) = tmp1;
                    end



                else
                    params_all_shots = cell(1, size(ksp_for_mussels_loraks,3));
                    Ahd_allS = zeros(Nx * Ny_actual, size(ksp_for_mussels_loraks,3));
                    for iter_s = 1:size(ksp_for_mussels_loraks,3)
                        % data is seg1echo1, seg1echo2, seg1echo3, seg2echo1,..
                        current_seg_ind = ceil(iter_s/Necho);
                        current_echo_ind  = mod(iter_s-1,Necho) + 1;
                        ky_current_seg = ky_allSeg_allDiffFrames(:,current_seg_ind,iter_diff_frame);


                        params = [];
                        params.Nx = Nx;
                        params.sens = sens(:,:,iter_slice,:);
                        params.EchoSpacing = Necho*EchoSpacing;
                        params.numCha = size(sens,4);
                        params.Ny_actual = Ny_actual;
                        params.ky_acquired = ky_current_seg(current_echo_ind:Necho:end);
                        params.FieldStrength_T = B0;
                        params.B0map_Hz = B0map_use(:,:,iter_slice);%imgaussfilt(b0map_gre_2echo_pulseq(:,:,iter_slice),3);
                        params.time =  distance_rf180_readout + ((1:numel(params.ky_acquired))-0.5  )*EchoSpacing - TE/2; % at spin-echo time TE, time is 0

                        params_all_shots{iter_s} = params;
                        A_adj = @(in)apply_sense_B0_informed(in, 'transp', params);
                        kdat = squeeze(navCor_kspace(:,current_echo_ind:Necho:end,iter_slice,current_seg_ind,iter_diff_frame ,:));
                        Ahd_allS(:, iter_s) = A_adj(kdat);
                    end

                    AtA_multishot = @(x_vec) apply_ata_B0informed_multishot(x_vec, params_all_shots, size(ksp_for_mussels_loraks,3));


                    AtA = AtA_multishot;
                    Ahd = Ahd_allS;
                end



                % im_tmp = CG_Recon(AHA,Ahd,Niter_cgsense2,tol_cgsense2,Display_cgsense2);

                % LORAKS
                disp('LORAKS Reconstruction');
                ns = size(ksp_for_mussels_loraks,3);
                LORAKS_type = 1; % S-matrix (Support + Phase + parallell imaging)

                % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$$$$$$$$$$$$$$$$$$$$$$$####################
                % S-LORAKS parameters

                % z = im_tmp;
                z = Ahd;zeros(Nx,Ny_actual,ns);%im_vc_mussels;
                % z = 0*z;warning('tmp')
                Radius = 3;
                rank = 60;
                lambda = 0.001;%0.003;%0.018;
                tol = 1e-4;
                max_iter = 15;  % change this based on necessary # iterations of S-loraks

                [in1,in2] = meshgrid(-Radius:Radius,-Radius:Radius);
                idx = find(in1.^2+in2.^2<=Radius^2);
                patchSize = numel(idx);

                coil_sens = params.sens;
                % coil_sens 1,= permute(coil_sens,[1,2,4,3]);

                %
                B = @(x) ft2(x); %
                Bh = @(x)  ift2(x); %

                %         P_M = @(x) LORAKS_operators(x,Nx,Ny,ns*num_virtual_chan,Radius,LORAKS_type,[]);
                %         Ph_M = @(x) LORAKS_operators(x,Nx,Ny,ns*num_virtual_chan,Radius,-LORAKS_type,[]);
                P_M = @(x) LORAKS_operators(x,nx,ny,ns,Radius,LORAKS_type,[]);
                Ph_M = @(x) LORAKS_operators(x,nx,ny,ns,Radius,-LORAKS_type,[]);

                N1 = nx;
                N2 = ny;
                Nc = ns;

                ZD = @(x) padarray(reshape(x,[N1 N2 Nc]),[2*Radius, 2*Radius], 'post');
                ZD_H = @(x) x(1:N1,1:N2,:,:);

                %
                tic
                for iter = 1:max_iter

                    z_cur = z;
                    pz = B(z);
                    MM = P_M(pz);
                    l = sum(MM,1);
                    l = find(abs(l)>0);
                    MM = MM(:,l);

                    [Um,s] = my_svd_left(MM);
                    nmm = Um(:,rank+1:end)'; % null space
                    Bhr = 0;


                    if 0
                        % [~,s] = svd(MM);
                        f = figure(105);
                        clf(f);
                        plot(s./max(s));
                        hold on;
                        xline(rank,'k')
                        drawnow
                    end

                    if LORAKS_type == 1 % S

                        nf = size(nmm,1);
                        nmm = reshape(nmm,[nf, patchSize, 2*Nc]);
                        nss_h = reshape(nmm(:,:,1:2:end)+1j*nmm(:,:,2:2:end),[nf, patchSize*Nc]);
                        Nis = filtfilt(nss_h,'C',N1,N2,Nc,Radius);
                        Nis2 = filtfilt(nss_h,'S',N1,N2,Nc,Radius);

                        L1 = @(x) ZD_H(ifft2(squeeze(sum(Nis.*repmat(fft2(ZD(B(x))),[1 1 1 Nc]),3))));
                        L2 = @(x) ZD_H(ifft2(squeeze(sum(Nis2.*repmat(conj(fft2(ZD(B(x)))),[1 1 1 Nc]),3))));
                        LhL = @(x) 2*Bh(reshape(L1(x)-L2(x),[nx ny ns 1]));

                    end

                    % data fitting
                    % M = @(x) AtA(x) + lambda*LhL(x);
                    reshape_input = @(x) reshape(x,[nx ny ns]);
                    M = @(x) vec(  vec(AtA(reshape_input(x))) + vec(  lambda*LhL(reshape_input(x)) )  ); % for pcg
                    [z,~] = pcg(M, Ahd(:) + lambda*Bhr);
                    z = reshape_input(z);

                    t = (norm(z_cur(:)-z(:))/norm(z(:)));

                    % display the status
                    if ~rem(iter,1)
                        disp(['iter ' num2str(iter) ', relative change in solution: ' num2str(t)]);
                        img_all_shots = reshape(z, [Nx, Ny_actual,  ns]);
                        img_all_shots_disp = img_all_shots;
                        % mosaic(img_all_shots_disp, 1, ns , 104, strcat('iter ',num2str(iter)), quantile(abs(Ahd(:)),[0.01,0.99]), 90);drawnow;
                        % figure(104); imagesc(imtile(imrotate(abs(img_all_shots_disp), 90)), quantile(abs(Ahd(:)),[0.01, 0.99])); axis image off; colormap gray; title(['iter ', num2str(iter)], 'Color', 'r', 'FontSize', 48); drawnow;
                    end

                    % check for convergence
                    if t < tol
                        disp('Convergence tolerance met: change in solution is small');
                        break;
                    end




                    im_S_LORAKS = z;


                end
                % imgs_S_loraks_single_channel(:,:,iter_slice,:,iter_diff_frame) = im_S_LORAKS;


                sz = size(im_S_LORAKS);
                tmp = reshape(im_S_LORAKS,sz(1),sz(2),numel(iter_seg),Necho,numel(iter_diff_frame));
                tmp1 = permute( squeeze( tmp ), [1 2 4 3 5] );
                imgs_S_loraks_single_channel_slc(:,:,:,iter_seg,iter_diff_frame) = tmp1;


                toc
            end % of single cha loraks
        end % of iter seg
    end % of iter_diff_





    %% start dixon
    isDixon = 1;
    isDxion_jointSeg_jointBvol = strcmp('jointOverSegmentsAndDirections',recon_strategy);
    isDxion_jointSeg_iterBvol = strcmp('jointOverSegmentsAndDirections',recon_strategy) || strcmp('jointOverSegments',recon_strategy);

    if isGEscanner
        PrecessionIsClockwise__ = -1; % 1 or -1
    else
        PrecessionIsClockwise__ = 1;
    end

    isOneFatPeak = 0;
    isManualFatPeakRatio = 0;



    manualFatPeakPpm = [-3.80, -3.40, -2.60, -1.94, -0.39, 0.60];
    manualFatRelAmps = [0.087 0.693 0.128 0.004 0.039 0.048];


    % % this is from leg scan
    %    manualFatPeakPpm = [-0.39, 0.60];
    % manualFatRelAmps = [0.214745762165254          0.785254237834746]; 

    
    


    manualFatRelAmps = manualFatRelAmps / sum(manualFatRelAmps);
    if isOneFatPeak
        hernando_dixon_func = @(x)func_fw_i2cm1i_3pluspoint_hernando_graphcut_oneFatPeak(x);
        wenmiao_dixon_func = @(x)func_Synthetic_wenmiao_111219_oneFatPeak(x);
        tsaojiang_dixon_func = @(x)func_fw_i2cm0c_3pluspoint_tsaojiang_wrapper(x, [-3.40], 1);
        fat_model_spec_for_ls = struct('name', 'custom_onepeak', ...
            'frequency_ppm', 0.60, 'relAmps', 1);
    else
        if isManualFatPeakRatio
            hernando_dixon_func = @(x)func_fw_i2cm1i_3pluspoint_hernando_graphcut_manual_ratio(x, manualFatPeakPpm, manualFatRelAmps);
            wenmiao_dixon_func = @(x)func_Synthetic_wenmiao_111219_manual_ratio(x, manualFatPeakPpm, manualFatRelAmps);
            tsaojiang_dixon_func = @(x)func_fw_i2cm0c_3pluspoint_tsaojiang_wrapper(x, manualFatPeakPpm, manualFatRelAmps);
            fat_model_spec_for_ls = struct('name', 'manual_ratio', ...
                'frequency_ppm', manualFatPeakPpm, 'relAmps', manualFatRelAmps);
        elseif isFatSat
            hernando_dixon_func = @(x)func_fw_i2cm1i_3pluspoint_hernando_graphcut_newPeakProportion(x);
            wenmiao_dixon_func = @(x)func_Synthetic_wenmiao_111219_newPeakProportion(x);
            tsaojiang_dixon_func = @(x)func_fw_i2cm0c_3pluspoint_tsaojiang_wrapper(x, [-3.40, -2.60, -0.39, 0.60], []);
            fat_model_spec_for_ls = struct('name', 'fatsat_4peak', ...
                'frequency_ppm', [-3.40, -2.60, -0.39, 0.60], ...
                'relAmps', ones(1,4) / 4);
        else
            hernando_dixon_func = @(x)func_fw_i2cm1i_3pluspoint_hernando_graphcut(x);
            wenmiao_dixon_func = @(x)func_Synthetic_wenmiao_111219(x);
            tsaojiang_dixon_func = @(x)func_fw_i2cm0c_3pluspoint_tsaojiang_wrapper(x, [-3.40, -2.60, -0.39, 0.60], []);
            fat_model_spec_for_ls = struct('name', 'hernando_6peak', ...
                'frequency_ppm', [-3.80, -3.40, -2.60, -1.94, -0.39, 0.60], ...
                'relAmps', [0.087 0.693 0.128 0.004 0.039 0.048]);
        end
    end




    if isDixon



        ECHO_DIM = 4;
        SEG_DIM = 5;
        DIFF_DIM = 6;
        reshape_func = @(x) reshape(x, [size(x,1), size(x,2), 1, 1, size(x,3)]);
        reshape_func2 = @(x) reshape(x, size(x,1), size(x,2), 1, 1, prod(size(x,3:ndims(x))));
        bvalunq = unique(bval);
        method = 'hernando_graphcut';
        method = 'wenmiao';
        % method = 'tsaojiang_Hierarchical_IDEAL';
        method_b0fmPrior = 'linear_ls_known_fieldmap';
        linear_ls_tikhonov_lambda = 1e-6;
        DO_MIXED_FIT = 0;


        %% calibration
        if 0

            isAverageSameTEData_precalib = 1;
            


            algoParams = [];
            algoParams.species(1).name = 'water';
            algoParams.species(1).frequency = 0;
            algoParams.species(1).relAmps = 1;
            fat_peak_ppm = [ -0.39, 0.60];
            for iter_peak = 1:numel(fat_peak_ppm)
                species_ind = iter_peak + 1;
                algoParams.species(species_ind).name = sprintf('fat_%+.2fppm', fat_peak_ppm(iter_peak));
                algoParams.species(species_ind).frequency = fat_peak_ppm(iter_peak);
                algoParams.species(species_ind).relAmps = 1;
            end


            for iter_bvol_dixon = 1:1%nImagingFrame % only use b-value=0
                img_slc = imgs_S_loraks_single_channel_slc(:,:,:, :,iter_bvol_dixon);
                img4dixon = squeeze(  img_slc  );

                data = [];
                tmpTE = vec(TEs_for_dixon(:,:,iter_bvol_dixon)).';
                tmpImages = reshape_func2(img4dixon);
                [tmpImages, tmpTE] = prepare_same_te_dixon_input(tmpImages, tmpTE, isAverageSameTEData_precalib);
                data.TE = tmpTE;
                data.FieldStrength = B0;
                data.PrecessionIsClockwise = PrecessionIsClockwise__;
                data.images = tmpImages;
                data.DO_MIXED_FIT = 0;

                outParams = fw_i2cm0c_3pluspoint_tsaojiang(data, algoParams);

                tsao_water = outParams.species(1).amps;
                tsao_fat_species = cat(3, outParams.species(2:end).amps);
                tsao_fat_species_abs = abs(tsao_fat_species);
                tsao_fat = sum(tsao_fat_species, 3);
                tsao_fat_sum = sum(tsao_fat_species_abs, 3);
                tsao_peak_ratio_maps = tsao_fat_species_abs ./ max(tsao_fat_sum, eps);
                tsao_peak_names = {algoParams.species(2:end).name};
                tsao_phasemap = outParams.phasemap;
                tsao_r2starmap = outParams.r2starmap;
                tsao_fiterror = outParams.fiterror;

                tsao_dTE = median(diff(sort(unique(data.TE))));
                tsao_fieldmap_hz = angle(tsao_phasemap * exp(-1i * 2 * pi / 3)) / (2 * pi * tsao_dTE);

                % Example if you want cropped outputs:
                % tmp_water_individualTE_oneByOne2(:,:,iter_slice,iter_seg_dixon,iter_bvol_dixon) = crop(tsao_water, Nx, Ny_actual);
                % tmp_fat_individualTE_oneByOne2(:,:,iter_slice,iter_seg_dixon,iter_bvol_dixon) = crop(tsao_fat, Nx, Ny_actual);
                % tmp_fieldmap_individualTE_oneByOne2(:,:,iter_slice,iter_seg_dixon,iter_bvol_dixon) = crop(tsao_fieldmap_hz, Nx, Ny_actual);
            end

            figure; imagesc(abs(tsao_fat_sum)); axis image; colormap gray; colorbar;
            title('Draw fat-rich ROI on tsao\_fat\_sum, then double click');
            tsao_fat_mask = roipoly;

            tsao_peak_amp_global = zeros(1, size(tsao_fat_species_abs, 3));
            for iter_peak = 1:size(tsao_fat_species_abs, 3)
                tmp_peak = tsao_fat_species_abs(:,:,iter_peak);
                tsao_peak_amp_global(iter_peak) = mean(tmp_peak(tsao_fat_mask));
            end
            tsao_peak_ratio_global = tsao_peak_amp_global / max(sum(tsao_peak_amp_global), eps);

            disp('tsao_peak_names =');
            disp(tsao_peak_names);
            disp('tsao_peak_amp_global =');
            disp(tsao_peak_amp_global);
            disp('tsao_peak_ratio_global =');
            disp(tsao_peak_ratio_global);
        end        



        



        %% toolbox dixon, shot by shot, volume by volume
        if 1
            for iter_bvol_dixon = 1:nImagingFrame
                % parfor iter_seg_dixon = 1:RSegment
                for iter_seg_dixon = 1:RSegment
                    img_slc = imgs_S_loraks_single_channel_slc(:,:,:,iter_seg_dixon,iter_bvol_dixon);
                    img4dixon = squeeze(  img_slc  );

                    data = [];
                    data.TE = vec(TEs_for_dixon(:,iter_seg_dixon,iter_bvol_dixon)).';
                    data.FieldStrength = B0;
                    data.PrecessionIsClockwise = PrecessionIsClockwise__;
                    data.images = reshape_func(   img4dixon   );
                    data.DO_MIXED_FIT = DO_MIXED_FIT;
                    outParams = [];
                    switch method
                        case 'hernando_graphcut'
                            if data.DO_MIXED_FIT
                                [~,outParams] = hernando_dixon_func(data);
                            else
                                outParams = hernando_dixon_func(data);
                            end
                        case 'MultiSeedRegionGrowing'
                            outParams = func_MultiSeedRegionGrowing(data); % very bad
                        case 'wenmiao'
                            outParams = wenmiao_dixon_func(data);
                        case 'tsaojiang_Hierarchical_IDEAL'
                            outParams = tsaojiang_dixon_func(data);
                        case 'RGideal'
                            sz = size(data.images);
                            data.images = zpad(data.images,[npad1 npad2 sz(3:end)]);
                            outParams = func_fw_i2cm0c_3pluspoint_RGmulticoil(data);
                        otherwise
                            error()
                    end

                    water_individualTE_iterSeg_iterBvol_slc(:,:,iter_seg_dixon,iter_bvol_dixon) = crop(outParams.species(1).amps,Nx,Ny_actual);
                    fat_individualTE_iterSeg_iterBvol_slc(:,:,iter_seg_dixon,iter_bvol_dixon) = crop(outParams.species(2).amps, Nx,Ny_actual);
                    fieldmap_individualTE_iterSeg_iterBvol_slc(:,:,iter_seg_dixon,iter_bvol_dixon) = crop(outParams.fieldmap, Nx,Ny_actual);
                end
            end
        end


        %% toolbox dixon, volume by volume, joint over segments
        if isDxion_jointSeg_iterBvol
            for iter_bvol_dixon = 1:nImagingFrame
            % for iter_bvol_dixon = 1%:nImagingFrame
                img_slc = imgs_S_loraks_single_channel_slc(:,:,:, :,iter_bvol_dixon);
                img4dixon = squeeze(  img_slc  );
                sz = size(img4dixon);
                img4dixon = reshape(img4dixon,sz(1),sz(2),[]);

                if 0 
                    warning('tmp!!!!!!!!!!!!!')                    
                    valid_TE_ind = vec(TEs_for_dixon(:,:,iter_bvol_dixon)).' > 10e-6;
                    tmpTE = vec(TEs_for_dixon(:,:,iter_bvol_dixon)).';
                    tmpDixonTE = tmpTE(valid_TE_ind);
                    img4dixon = img4dixon(:,:,valid_TE_ind);
                end

                data = [];
                data.TE = vec(TEs_for_dixon(:,:,iter_bvol_dixon)).';
                data.FieldStrength = B0;
                data.PrecessionIsClockwise = PrecessionIsClockwise__;
                data.images = reshape_func(   img4dixon   );
                data.DO_MIXED_FIT = DO_MIXED_FIT;
                outParams = [];
                switch method
                    case 'hernando_graphcut'
                        if data.DO_MIXED_FIT
                            [~,outParams] = hernando_dixon_func(data);
                        else
                            outParams = hernando_dixon_func(data);
                        end
                    case 'MultiSeedRegionGrowing'
                        outParams = func_MultiSeedRegionGrowing(data); % very bad
                    case 'wenmiao'
                        outParams = wenmiao_dixon_func(data);
                    case 'tsaojiang_Hierarchical_IDEAL'
                        outParams = tsaojiang_dixon_func(data);
                    case 'RGideal'
                        sz = size(data.images);
                        data.images = zpad(data.images,[npad1 npad2 sz(3:end)]);
                        outParams = func_fw_i2cm0c_3pluspoint_RGmulticoil(data);
                    otherwise
                        error()
                end

                water_individualTE_jointSeg_iterBvol_slc(:,:,iter_bvol_dixon) = crop(outParams.species(1).amps,Nx,Ny_actual);
                fat_individualTE_jointSeg_iterBvol_slc(:,:,iter_bvol_dixon) = crop(outParams.species(2).amps, Nx,Ny_actual);
                fieldmap_individualTE_jointSeg_iterBvol_slc(:,:,iter_bvol_dixon) = crop(outParams.fieldmap, Nx,Ny_actual);
            end
        end


        %% dixon on b-value=0 first, get field map, then use this field map to do Dixon on b-value>0
        % when b-value=0, dixon uses data from all shots, just like water_individualTE_jointSeg_iterBvol
        % when b-value>0, dixon can either use data from a single shot, just like water_individualTE_iterSeg_iterBvol, or use data from all shots, just like water_individualTE_jointSeg_iterBvol        
        if 1
            flag_use_all_shots_when_bval_larger_than_0 = 0;

            
            b0_frames = find(bval<=bvalue0_threshold);
            assert(~isempty(b0_frames), 'Need at least one b=0 frame for the b0-fieldmap-prior Dixon branch.');
            ref_b0_frame = b0_frames(1);
            fieldmap_prior_b0 = fieldmap_individualTE_jointSeg_iterBvol_slc(:,:,ref_b0_frame);

            for iter_bvol_dixon = 1:nImagingFrame
                if bval(iter_bvol_dixon) <= bvalue0_threshold
                    water_individualTE_b0fmPrior_jointSeg_iterBvol_slc(:,:,iter_bvol_dixon) = ...
                        water_individualTE_jointSeg_iterBvol_slc(:,:,iter_bvol_dixon);
                    fat_individualTE_b0fmPrior_jointSeg_iterBvol_slc(:,:,iter_bvol_dixon) = ...
                        fat_individualTE_jointSeg_iterBvol_slc(:,:,iter_bvol_dixon);
                    fieldmap_individualTE_b0fmPrior_jointSeg_iterBvol_slc(:,:,iter_bvol_dixon) = ...
                        fieldmap_individualTE_jointSeg_iterBvol_slc(:,:,iter_bvol_dixon);

                    for iter_seg_dixon = 1:RSegment
                        water_individualTE_b0fmPrior_iterSeg_iterBvol_slc(:,:,iter_seg_dixon,iter_bvol_dixon) = ...
                            water_individualTE_jointSeg_iterBvol_slc(:,:,iter_bvol_dixon);
                        fat_individualTE_b0fmPrior_iterSeg_iterBvol_slc(:,:,iter_seg_dixon,iter_bvol_dixon) = ...
                            fat_individualTE_jointSeg_iterBvol_slc(:,:,iter_bvol_dixon);
                        fieldmap_individualTE_b0fmPrior_iterSeg_iterBvol_slc(:,:,iter_seg_dixon,iter_bvol_dixon) = ...
                            fieldmap_individualTE_jointSeg_iterBvol_slc(:,:,iter_bvol_dixon);
                    end
                    continue
                end

                if flag_use_all_shots_when_bval_larger_than_0
                    img_slc = imgs_S_loraks_single_channel_slc(:,:,:, :,iter_bvol_dixon);
                    img4dixon = squeeze(img_slc);
                    sz = size(img4dixon);
                    img4dixon = reshape(img4dixon,sz(1),sz(2),[]);

                    data = [];
                    data.TE = vec(TEs_for_dixon(:,:,iter_bvol_dixon)).';
                    data.FieldStrength = B0;
                    data.PrecessionIsClockwise = PrecessionIsClockwise__;
                    data.images = reshape_func(img4dixon);
                    data.DO_MIXED_FIT = DO_MIXED_FIT;
                    data.fieldmap_init = fieldmap_prior_b0;
                    data.r2starmap_init = zeros(size(fieldmap_prior_b0));
                    data.linear_ls_tikhonov_lambda = linear_ls_tikhonov_lambda;

                    outParams = run_dixon_with_optional_fieldmap_prior( ...
                        data, method_b0fmPrior, hernando_dixon_func, wenmiao_dixon_func, tsaojiang_dixon_func, fat_model_spec_for_ls);

                    water_individualTE_b0fmPrior_jointSeg_iterBvol_slc(:,:,iter_bvol_dixon) = crop(outParams.species(1).amps,Nx,Ny_actual);
                    fat_individualTE_b0fmPrior_jointSeg_iterBvol_slc(:,:,iter_bvol_dixon) = crop(outParams.species(2).amps,Nx,Ny_actual);
                    fieldmap_individualTE_b0fmPrior_jointSeg_iterBvol_slc(:,:,iter_bvol_dixon) = crop(outParams.fieldmap,Nx,Ny_actual);
                    resnorm_individualTE_b0fmPrior_jointSeg_iterBvol_slc(:,:,iter_bvol_dixon) = crop(outParams.resnorm, Nx, Ny_actual);
                    residual_rms_individualTE_b0fmPrior_jointSeg_iterBvol_slc(:,:,iter_bvol_dixon) = crop(outParams.residual_rms, Nx, Ny_actual);
                else
                    for iter_seg_dixon = 1:RSegment
                        img_slc = imgs_S_loraks_single_channel_slc(:,:,:,iter_seg_dixon,iter_bvol_dixon);
                        img4dixon = squeeze(img_slc);

                        data = [];
                        data.TE = vec(TEs_for_dixon(:,iter_seg_dixon,iter_bvol_dixon)).';
                        data.FieldStrength = B0;
                        data.PrecessionIsClockwise = PrecessionIsClockwise__;
                        data.images = reshape_func(img4dixon);
                        data.DO_MIXED_FIT = DO_MIXED_FIT;
                        data.fieldmap_init = fieldmap_prior_b0;
                        data.r2starmap_init = zeros(size(fieldmap_prior_b0));
                        data.linear_ls_tikhonov_lambda = linear_ls_tikhonov_lambda;

                        outParams = run_dixon_with_optional_fieldmap_prior( ...
                            data, method_b0fmPrior, hernando_dixon_func, wenmiao_dixon_func, tsaojiang_dixon_func, fat_model_spec_for_ls);

                        water_individualTE_b0fmPrior_iterSeg_iterBvol_slc(:,:,iter_seg_dixon,iter_bvol_dixon) = crop(outParams.species(1).amps,Nx,Ny_actual);
                        fat_individualTE_b0fmPrior_iterSeg_iterBvol_slc(:,:,iter_seg_dixon,iter_bvol_dixon) = crop(outParams.species(2).amps,Nx,Ny_actual);
                        fieldmap_individualTE_b0fmPrior_iterSeg_iterBvol_slc(:,:,iter_seg_dixon,iter_bvol_dixon) = crop(outParams.fieldmap,Nx,Ny_actual);
                        resnorm_individualTE_b0fmPrior_iterSeg_iterBvol_slc(:,:,iter_seg_dixon,iter_bvol_dixon) = crop(outParams.resnorm, Nx, Ny_actual);
                        residual_rms_individualTE_b0fmPrior_iterSeg_iterBvol_slc(:,:,iter_seg_dixon,iter_bvol_dixon) = crop(outParams.residual_rms, Nx, Ny_actual);
                    end
                end
            end
        end
  
    end % of isDixon

    imgs_vc_mussels(:,:,iter_slice,:,:,:) = imgs_vc_mussels_slc;
    imgs_S_loraks_single_channel(:,:,iter_slice,:,:,:) = imgs_S_loraks_single_channel_slc;
    imgs_MUSE(:,:,iter_slice,:,:,:) = imgs_MUSE_slc;
    water_individualTE_iterSeg_iterBvol(:,:,iter_slice,:,:) = water_individualTE_iterSeg_iterBvol_slc;
    fat_individualTE_iterSeg_iterBvol(:,:,iter_slice,:,:) = fat_individualTE_iterSeg_iterBvol_slc;
    fieldmap_individualTE_iterSeg_iterBvol(:,:,iter_slice,:,:) = fieldmap_individualTE_iterSeg_iterBvol_slc;
    water_individualTE_jointSeg_iterBvol(:,:,iter_slice,:) = water_individualTE_jointSeg_iterBvol_slc;
    fat_individualTE_jointSeg_iterBvol(:,:,iter_slice,:) = fat_individualTE_jointSeg_iterBvol_slc;
    fieldmap_individualTE_jointSeg_iterBvol(:,:,iter_slice,:) = fieldmap_individualTE_jointSeg_iterBvol_slc;
    water_individualTE_jointSeg_jointBvol(:,:,iter_slice,:) = water_individualTE_jointSeg_jointBvol_slc;
    fat_individualTE_jointSeg_jointBvol(:,:,iter_slice,:) = fat_individualTE_jointSeg_jointBvol_slc;
    fieldmap_individualTE_jointSeg_jointBvol(:,:,iter_slice,:) = fieldmap_individualTE_jointSeg_jointBvol_slc;
    water_individualTE_b0fmPrior_iterSeg_iterBvol(:,:,iter_slice,:,:) = water_individualTE_b0fmPrior_iterSeg_iterBvol_slc;
    fat_individualTE_b0fmPrior_iterSeg_iterBvol(:,:,iter_slice,:,:) = fat_individualTE_b0fmPrior_iterSeg_iterBvol_slc;
    fieldmap_individualTE_b0fmPrior_iterSeg_iterBvol(:,:,iter_slice,:,:) = fieldmap_individualTE_b0fmPrior_iterSeg_iterBvol_slc;
    resnorm_individualTE_b0fmPrior_iterSeg_iterBvol(:,:,iter_slice,:,:) = resnorm_individualTE_b0fmPrior_iterSeg_iterBvol_slc;
    residual_rms_individualTE_b0fmPrior_iterSeg_iterBvol(:,:,iter_slice,:,:) = residual_rms_individualTE_b0fmPrior_iterSeg_iterBvol_slc;
    water_individualTE_b0fmPrior_jointSeg_iterBvol(:,:,iter_slice,:) = water_individualTE_b0fmPrior_jointSeg_iterBvol_slc;
    fat_individualTE_b0fmPrior_jointSeg_iterBvol(:,:,iter_slice,:) = fat_individualTE_b0fmPrior_jointSeg_iterBvol_slc;
    fieldmap_individualTE_b0fmPrior_jointSeg_iterBvol(:,:,iter_slice,:) = fieldmap_individualTE_b0fmPrior_jointSeg_iterBvol_slc;
    resnorm_individualTE_b0fmPrior_jointSeg_iterBvol(:,:,iter_slice,:) = resnorm_individualTE_b0fmPrior_jointSeg_iterBvol_slc;
    residual_rms_individualTE_b0fmPrior_jointSeg_iterBvol(:,:,iter_slice,:) = residual_rms_individualTE_b0fmPrior_jointSeg_iterBvol_slc;
    water_avgOverDirection(:,:,iter_slice,:) = water_avgOverDirection_slc;
    fat_avgOverDirection(:,:,iter_slice,:) = fat_avgOverDirection_slc;
    fieldmap_avgOverDirection(:,:,iter_slice,:) = fieldmap_avgOverDirection_slc;
    water_avgOverSegOverDirection(:,:,iter_slice,:) = water_avgOverSegOverDirection_slc;
    fat_avgOverSegOverDirection(:,:,iter_slice,:) = fat_avgOverSegOverDirection_slc;
    fieldmap_avgOverSegOverDirection(:,:,iter_slice,:) = fieldmap_avgOverSegOverDirection_slc;

end % of iter slice

disp done



slc_ind_show = 13;
vol_ind_show = 2; % b1k
src_img_to_show = imgs_S_loraks_single_channel(:,:,slc_ind_show,1,vol_ind_show);
water_to_show = squeeze(water_individualTE_iterSeg_iterBvol(:,:,slc_ind_show,1,vol_ind_show));
fat_to_show = squeeze(fat_individualTE_iterSeg_iterBvol(:,:,slc_ind_show,1,vol_ind_show));


ims(src_img_to_show);title 'source image';clim = get(gca, 'CLim');
ims(water_to_show,clim);title 'water';colorbar off
ims(fat_to_show,clim);title 'fat';colorbar off
return

%%
[p,n] = fileparts(dat_fn);
res_img_path = fullfile(p,'res_img',[recon_strategy,'_','isCompositeSens',num2str(isCompositeSens),'_',n,'_',datestr(now,30)]);
[~, res_img_prefix] = fileparts(seq_fn);
res_name = fullfile(res_img_path,res_img_prefix);


if ~isfolder(res_img_path)
    mkdir(res_img_path)
end

save(res_name,...
    'imgs_S_loraks_single_channel',...
    'fieldmap_individualTE_jointSeg_iterBvol',...
    'water_individualTE_iterSeg_iterBvol',...
    'fat_individualTE_iterSeg_iterBvol',...
    'fieldmap_individualTE_iterSeg_iterBvol',...
    'water_individualTE_b0fmPrior_iterSeg_iterBvol',...
    'fat_individualTE_b0fmPrior_iterSeg_iterBvol',...
    'fieldmap_individualTE_b0fmPrior_iterSeg_iterBvol',...
        'gre_img',...
    'fat_image_from_gre',...
    'water_image_from_gre',...
    '-v7.3');


filename = mfilename('fullpath');
copyfile([filename,'.m'], fileparts(res_name));



return
% end



function res = apply_sense_nechoes_noreshape(in, tflag, params )
% apply sense for each shot, i.e. each shot gives an image


N1 = params.N1;
N2 = params.N2;
sens = params.sens;
masks = params.masks;
Necho = params.Necho;

sz = size(sens);
num_chan = sz(end);


if strcmp(tflag,'transp')

    % Transposed SENSE operator:
    % IFFT coil k-space, multiply by conjugate of coil sensitivities, then
    % sum across channels


    kspace = reshape(in, [N1,N2, Necho, num_chan]);
    ksp = kspace .* masks;
    % img_coils = bart('fft -i -u 3',ksp);
    img_coils = xyrecon.ifft2c(ksp);
    res = sum(img_coils .* conj(sens), 4);

    % res = res(:);

elseif strcmp(tflag,'notransp')

    % Forward SENSE operator:
    % multiply by coil sensitivities, take undersampled FFT

    img = reshape(in, N1,N2,Necho) .* sens;
    % kspace = bart('fft -u 3',img);
    kspace = xyrecon.fft2c(img);
    res = kspace .* masks;
    % res = res(:);
else
    error('Syntax error, un-recognized 2nd param: %s', tflag);
end
end



function [U, s] = my_svd_left(A, r)
% Returns the left singular vectors U and singular values s (in vector form) simultaneously
if nargin < 2
    % --- Compute all singular values ---
    [U, E] = eig(A * A');

    % Get eigenvalues and sort in descending order
    [sorted_eigenvalues, idx] = sort(abs(diag(E)), 'descend');

    % Reorder eigenvectors to match
    U = U(:,idx);

    % Compute singular values s (take square root)
    s = sqrt(sorted_eigenvalues);
else
    % --- Compute only the r largest singular values ---
    [U, E] = eigs(double(A * A'), r);

    % Compute singular values s (eigs returns the largest by default, no need to sort)
    s = sqrt(diag(E));

    % Convert type
    U = cast(U, 'like', A);
end
end



function [water_image, fat_image, resnorm_map, residual_rms_map] = dixon_least_squares(images, TEs, B0map_Hz, FieldStrength_T, fat_model_choice, phase_flag, lambda_tikhonov)
%DIXON_LEAST_SQUARES Separates water and fat from multi-echo images using a linear least-squares fit.
%
% [water_image, fat_image] = dixon_least_squares(images, TEs, B0map_Hz, FieldStrength_T, fat_model_choice)
%
% This function solves the simplified Dixon model on a pixel-by-pixel basis,
% assuming the B0 field map is known. The model is:
%
%   s(t_n) = exp(i*2*pi*f_B*t_n) * (rho_w + rho_f * C_f(t_n))
%
% Rearranging this gives a linear system for each pixel:
%   s(t_n) * exp(-i*2*pi*f_B*t_n) = rho_w + rho_f * C_f(t_n)
%
% This can be written as a linear system A*x = b, where:
%   x = [rho_w; rho_f]
%   b_n = s(t_n) * exp(-i*2*pi*f_B*t_n)
%   A_n = [1, C_f(t_n)]
%
% INPUTS:
%   images           - Complex-valued multi-echo images [Nx, Ny, N_TE].
%   TEs              - Vector of echo times in seconds [1, N_TE].
%   B0map_Hz         - B0 field map in Hz [Nx, Ny].
%   FieldStrength_T  - Main field strength in Tesla (e.g., 3.0).
%   fat_model_choice - String specifying the fat model. Options:
%                      'hernando_6peak', 'berglund_9peak', 'tsao_jiang_6peak',
%                      'tsao_jiang_1peak', 'lu_6peak'.
%
% OUTPUTS:
%   water_image      - Separated complex water image [Nx, Ny].
%   fat_image        - Separated complex fat image [Nx, Ny].
%

%% --- 1. Input Validation and Parameter Setup ---
if nargin < 5
    error('At least 5 input arguments are required: images, TEs, B0map_Hz, FieldStrength_T, fat_model_choice.');
end

if nargin < 6 || isempty(phase_flag)
    phase_flag = 1;
end

if nargin < 7 || isempty(lambda_tikhonov)
    lambda_tikhonov = 0;
end

[Nx, Ny, N_TE] = size(images);

if length(TEs) ~= N_TE
    error('Number of TEs must match the 3rd dimension of the images array.');
end

if any(size(B0map_Hz) ~= [Nx, Ny])
    error('B0 map dimensions must match the image dimensions.');
end

% Ensure TEs is a column vector for matrix operations
TEs = TEs(:); % [N_TE, 1]

%% --- 2. Define Multi-Peak Fat Model ---
% This section is adapted from your provided code.
persistent fat_model;
if ~isstruct(fat_model_choice) && (isempty(fat_model) || ~isfield(fat_model, fat_model_choice))

    % --- Hernando, 6-peak model ---
    spec.hernando_6peak.name = 'Fat (Hernando, 6 peaks)';
    spec.hernando_6peak.frequency_ppm = [-3.80, -3.40, -2.60, -1.94, -0.39, 0.60];
    spec.hernando_6peak.relAmps = [0.087 0.693 0.128 0.004 0.039 0.048];

    % --- Berglund, 9-peak model ---
    spec.berglund_9peak.name = 'Fat (Berglund, 9 peaks)';
    spec.berglund_9peak.frequency_ppm = [0.90, 1.30, 1.60, 2.02, 2.24, 2.75, 4.20, 5.19, 5.29] - 4.70;
    spec.berglund_9peak.relAmps = [88 642 58 62 58 6 39 10 37];
    spec.berglund_9peak.relAmps = spec.berglund_9peak.relAmps / sum(spec.berglund_9peak.relAmps); % Normalize amps

    % --- Tsao/Jiang, 6-peak model ---
    spec.tsao_jiang_6peak.name = 'Fat (Tsao/Jiang, 6 peaks)';
    spec.tsao_jiang_6peak.frequency_ppm = [0.9000, 1.3000, 2.1000, 2.7600, 4.3100, 5.3000] - 4.7;
    spec.tsao_jiang_6peak.relAmps =   [0.0871, 0.6937, 0.1281, 0.0040, 0.0390, 0.0480];

    % --- Tsao/Jiang, single-peak model ---
    spec.tsao_jiang_1peak.name = 'Fat (Tsao/Jiang, 1 peak)';
    spec.tsao_jiang_1peak.frequency_ppm = [1.3000] - 4.7;
    spec.tsao_jiang_1peak.relAmps =   [1.0000];

    % --- Lu, 6-peak model (same as Hernando) ---
    spec.lu_6peak = spec.hernando_6peak;
    spec.lu_6peak.name = 'Fat (Lu, 6 peaks)';

    fat_model = spec;
end

% Select the active model
if isstruct(fat_model_choice)
    active_model = fat_model_choice;
else
    active_model = fat_model.(fat_model_choice);
end
fat_freq_ppm = active_model.frequency_ppm;
fat_rel_amps = active_model.relAmps;
fat_rel_amps = fat_rel_amps / sum(fat_rel_amps);

% Convert fat frequencies from ppm to Hz
gyro = 42.576e6; % Gyromagnetic ratio for protons in Hz/T
fat_freq_Hz = FieldStrength_T * gyro * fat_freq_ppm / 1e6;

% Ensure dimensions are correct for matrix operations
fat_freq_Hz = fat_freq_Hz(:).';  % Ensure [1, N_peaks]
fat_rel_amps = fat_rel_amps(:).'; % Ensure [1, N_peaks]


%% --- 3. Demodulate the B0 field effect from the signal ---
% This prepares the 'b' vector in A*x = b.

% Reshape for vectorized operations
images_vec = reshape(images, [Nx*Ny, N_TE]);
B0map_vec = reshape(B0map_Hz, [Nx*Ny, 1]);

% B0 phase term: exp(-i*2*pi*f_B*t_n)
b0_phase_demod = exp(-phase_flag*1i * 2*pi * B0map_vec .* TEs'); % [N_pixels, N_TE]

% Demodulated signal (our 'b' vector for each pixel)
b = images_vec .* b0_phase_demod; % [N_pixels, N_TE]


%% --- 4. Construct the System Matrix A ---
% The matrix A is the same for every pixel.
% It has N_TE rows and 2 columns.
% Column 1 is for water (all ones).
% Column 2 is for fat (the complex modulator C_f(t_n)).

% Calculate the phase for each TE and each fat peak using an outer product.
% TEs is [N_TE x 1]
% fat_freq_Hz is [1 x N_peaks]
% Result is an [N_TE x N_peaks] matrix.
phase_matrix = exp(phase_flag*1i * 2*pi * (TEs * fat_freq_Hz));

% Calculate the fat complex modulator C_f(t_n) for all TEs.
% This is a weighted sum over the peaks for each TE.
% 'fat_rel_amps' is [1 x N_peaks].
% The result of the matrix multiplication will be a column vector [N_TE x 1].
fat_modulator = phase_matrix * fat_rel_amps.';

% Construct the system matrix A
A = zeros(N_TE, 2);
A(:, 1) = 1;             % Column for water
A(:, 2) = fat_modulator; % Column for fat


%% --- 5. Solve the Linear System for All Pixels ---
% Tikhonov-regularized least squares:
% x = (A^H * A + lambda * I)^-1 * A^H * b
A_reg = (A' * A + lambda_tikhonov * eye(2)) \ A';

% Solve for all pixels at once using matrix multiplication.
% (The 'b' matrix needs to be transposed for the multiplication)
% x_all_pixels will have size [2, N_pixels]
x_all_pixels = A_reg * b.';

% Residual diagnostics
residual_all_pixels = b.' - A * x_all_pixels; % [N_TE, N_pixels]
resnorm_vec = sum(abs(residual_all_pixels).^2, 1);
residual_rms_vec = sqrt(mean(abs(residual_all_pixels).^2, 1));


%% --- 6. Reshape the output back to images ---
water_image_vec = x_all_pixels(1, :);
fat_image_vec = x_all_pixels(2, :);

water_image = reshape(water_image_vec, [Nx, Ny]);
fat_image = reshape(fat_image_vec, [Nx, Ny]);
resnorm_map = reshape(resnorm_vec, [Nx, Ny]);
residual_rms_map = reshape(residual_rms_vec, [Nx, Ny]);

end
function [ dB,fit_error,phase ] = dB_fitting_JumpCorrect( phase,TEs,mask,unwrap_flag)
% Fitting of the B0 inhomogeneity field (dB) using multi-echo phase data
% Fuyixue Wang, 2018

% updated, clean up
% Zijing Dong, 2019

if(unwrap_flag==1)
    phase=unwrap(phase,[],3);
end
phase=phase./(2*pi);
[nx,ny,nt]=size(phase);
A = [ones(nt,1),TEs(1:nt)];
dB=zeros(nx,ny);
fit_error=zeros(nx,ny);
for ii=1:nx
    for jj=1:ny
        if mask(ii,jj) == 1
            % fit voxel
            signal=double(squeeze(phase(ii,jj,:)));
            %            signal_diff=signal(2:end)-signal(1:end-1);
            %            Jump_index=find(abs(signal_diff)>0.3);
            %            if (~isempty(Jump_index))&&(length(Jump_index)<3)
            %            end
            param = A\signal;
            phase_fitting=A*param;
            error = abs(phase_fitting-signal);
            %            mask_temp=(error > thershold*mean(error(:)));
            %            jump_num=sum(mask_temp(:));
            %            if (jump_num>0)&&(jump_num<2)
            %                mask_temp=logical(1-mask_temp);
            %                param = A(mask_temp,:)\signal(mask_temp);
            %                phase_fitting=A(mask_temp,:)*param;
            %                error = abs(phase_fitting-signal(mask_temp));
            %            end
            dB(ii,jj) = param(2);
            fit_error(ii,jj)=sqrt(sum(error.^2)/size(error,1));
        end
    end % end jj
end % end ii
dB(isnan(dB)) = 0;
% dB=dB-mean(dB(mask(:)));
dB = dB.*mask;
phase=phase.*(2*pi);
end



% function outParamsMixed = func_fw_i2cm1i_3pluspoint_hernando_graphcut(data)
function [outParams,outParamsMixed] = func_fw_i2cm1i_3pluspoint_hernando_graphcut(data)
%% test_hernando_111110
%%
%% Test fat-water algorithms on (Peter Kellman's) data
%%
%% Author: Diego Hernando
%% Date created: August 18, 2011
%% Date last modified: February 29, 2011

% Add to matlab path
% if ispc
%     BASEPATH = '';
% else
%     BASEPATH = '';
% end
% 
% addpath([BASEPATH 'common/']);
% addpath([BASEPATH 'graphcut/']);
% addpath([BASEPATH 'descent/']);
% addpath([BASEPATH 'mixed_fitting/']);
% addpath([BASEPATH 'create_synthetic/']);
% addpath([BASEPATH 'matlab_bgl/']);



%% Load some data
% foldername = [BASEPATH '../../fwtoolbox_v1_data/kellman_data/'];
% fn = dir([foldername '*.mat']);
% file_index = ceil(rand*length(fn));
% disp([foldername fn(file_index).name]);
% load([foldername fn(file_index).name]);
imDataParams = data;
imDataParams.images = double(data.images);
% $$$ imDataParams.FieldStrength = 1.5;
% $$$ imDataParams.PrecessionIsClockwise = -1;

%% Set recon parameters
% General parameters
algoParams.species(1).name = 'water';
algoParams.species(1).frequency = 0;
algoParams.species(1).relAmps = 1;
algoParams.species(2).name = 'fat';
algoParams.species(2).frequency = [-3.80, -3.40, -2.60, -1.94, -0.39, 0.60];
algoParams.species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048];

% Algorithm-specific parameters
algoParams.size_clique = 1; % Size of MRF neighborhood (1 uses an 8-neighborhood, common in 2D)
algoParams.range_r2star = [0 100]; % Range of R2* values
algoParams.NUM_R2STARS = 11; % Numbre of R2* values for quantization
algoParams.range_fm = [-400 400]; % Range of field map values
algoParams.NUM_FMS = 801; % Number of field map values to discretize
algoParams.NUM_ITERS = 40; % Number of graph cut iterations
algoParams.SUBSAMPLE = 1; % Spatial subsampling for field map estimation (for speed)
algoParams.DO_OT = 1; % 0,1 flag to enable optimization transfer descent (final stage of field map estimation)
algoParams.LMAP_POWER = 2; % Spatially-varying regularization (2 gives ~ uniformn resolution)
algoParams.lambda = 0.05; % Regularization parameter
algoParams.LMAP_EXTRA = 0.05; % More smoothing for low-signal regions
algoParams.TRY_PERIODIC_RESIDUAL = 0;
THRESHOLD = 0.01;

if isfield(data,'fieldmap_init') && ~isempty(data.fieldmap_init)
    algoParams.fieldmap = sanitize_init_map(data.fieldmap_init, 0, [-400 400]);
    if isfield(data,'r2starmap_init') && ~isempty(data.r2starmap_init)
        algoParams.r2starmap = sanitize_init_map(data.r2starmap_init, 0, [0 200]);
    else
        algoParams.r2starmap = zeros(size(algoParams.fieldmap));
    end
    algoParams.NUM_MAGN = 1;
    algoParams.THRESHOLD = 0.04;
    algoParams.range_r2star = [0 200];
    try
        outParams = fw_i2xm1c_3pluspoint_hernando_mixedfit( imDataParams, algoParams );
        outParamsMixed = outParams;
        return
    catch ME
        warning('Mixed-fit with fieldmap prior failed in hernando_graphcut: %s. Falling back to graph-cut initialization.', char(ME.message));
    end
end

%% Recon -- graph cut
%% (Hernando D, Kellman P, Haldar JP, Liang ZP. Robust water/fat separation in the presence of large
%% field inhomogeneities using a graph cut algorithm. Magn Reson Med. 2010 Jan;63(1):79-90.)
tic
outParams = fw_i2cm1i_3pluspoint_hernando_graphcut( imDataParams, algoParams );
toc

DO_MIXED_FIT = data.DO_MIXED_FIT;0;
if DO_MIXED_FIT > 0
    try
        %% Recon -- mixed fit for phase error correction
        % Initialize mixed fitting to graph cut solution
        algoParams.fieldmap = outParams.fieldmap;
        algoParams.r2starmap = outParams.r2starmap;
        algoParams.NUM_MAGN = 1;
        algoParams.THRESHOLD = 0.04;
        algoParams.range_r2star = [0 200];

        % Do mixed fitting
        %% (Hernando D, Hines CDG, Yu H, Reeder SB. Addressing phase errors in fat-water imaging
        %% using a mixed magnitude/complex fitting method. Magn Reson Med; 2011.)
        outParamsMixed = fw_i2xm1c_3pluspoint_hernando_mixedfit( imDataParams, algoParams );
    end
end


end


function [outParams,outParamsMixed] = func_fw_i2cm1i_3pluspoint_hernando_graphcut_newPeakProportion(data)
%% test_hernando_111110
%%
%% Test fat-water algorithms on (Peter Kellman's) data
%%
%% Author: Diego Hernando
%% Date created: August 18, 2011
%% Date last modified: February 29, 2011

% Add to matlab path
% if ispc
%     BASEPATH = '';
% else
%     BASEPATH = '';
% end
% 
% addpath([BASEPATH 'common/']);
% addpath([BASEPATH 'graphcut/']);
% addpath([BASEPATH 'descent/']);
% addpath([BASEPATH 'mixed_fitting/']);
% addpath([BASEPATH 'create_synthetic/']);
% addpath([BASEPATH 'matlab_bgl/']);



%% Load some data
% foldername = [BASEPATH '../../fwtoolbox_v1_data/kellman_data/'];
% fn = dir([foldername '*.mat']);
% file_index = ceil(rand*length(fn));
% disp([foldername fn(file_index).name]);
% load([foldername fn(file_index).name]);
imDataParams = data;
imDataParams.images = double(data.images);
% $$$ imDataParams.FieldStrength = 1.5;
% $$$ imDataParams.PrecessionIsClockwise = -1;

%% Set recon parameters
% General parameters
algoParams.species(1).name = 'water';
algoParams.species(1).frequency = 0;
algoParams.species(1).relAmps = 1;
algoParams.species(2).name = 'fat';
algoParams.species(2).frequency = [-3.80, -3.40, -2.60, -1.94, -0.39, 0.60];
algoParams.species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048];
algoParams.species(2).relAmps = [0.0252    0.0935    0.0514    0.0030    0.0386    0.0479];
algoParams.species(2).relAmps = [0.0252    0.0935    0.0514    0.0030    0.0386    0.0479] / sum([0.0252    0.0935    0.0514    0.0030    0.0386    0.0479]);

% Algorithm-specific parameters
algoParams.size_clique = 1; % Size of MRF neighborhood (1 uses an 8-neighborhood, common in 2D)
algoParams.range_r2star = [0 100]; % Range of R2* values
algoParams.NUM_R2STARS = 11; % Numbre of R2* values for quantization
algoParams.range_fm = [-400 400]; % Range of field map values
algoParams.NUM_FMS = 301; % Number of field map values to discretize
algoParams.NUM_ITERS = 40; % Number of graph cut iterations
algoParams.SUBSAMPLE = 2; % Spatial subsampling for field map estimation (for speed)
algoParams.DO_OT = 1; % 0,1 flag to enable optimization transfer descent (final stage of field map estimation)
algoParams.LMAP_POWER = 2; % Spatially-varying regularization (2 gives ~ uniformn resolution)
algoParams.lambda = 0.05; % Regularization parameter
algoParams.LMAP_EXTRA = 0.05; % More smoothing for low-signal regions
algoParams.TRY_PERIODIC_RESIDUAL = 0;
THRESHOLD = 0.01;

%% Recon -- graph cut
%% (Hernando D, Kellman P, Haldar JP, Liang ZP. Robust water/fat separation in the presence of large
%% field inhomogeneities using a graph cut algorithm. Magn Reson Med. 2010 Jan;63(1):79-90.)
tic
outParams = fw_i2cm1i_3pluspoint_hernando_graphcut( imDataParams, algoParams );
toc

DO_MIXED_FIT = data.DO_MIXED_FIT;0;
if DO_MIXED_FIT > 0
    % try
        %% Recon -- mixed fit for phase error correction
        % Initialize mixed fitting to graph cut solution
        algoParams.fieldmap = outParams.fieldmap;
        algoParams.r2starmap = outParams.r2starmap;
        algoParams.NUM_MAGN = 1;
        algoParams.THRESHOLD = 0.04;
        algoParams.range_r2star = [0 200];

        % Do mixed fitting
        %% (Hernando D, Hines CDG, Yu H, Reeder SB. Addressing phase errors in fat-water imaging
        %% using a mixed magnitude/complex fitting method. Magn Reson Med; 2011.)
        outParamsMixed = fw_i2xm1c_3pluspoint_hernando_mixedfit( imDataParams, algoParams );
    % end
end


end



function [echo_te] = get_dur_to_center(ky, N_echo, esp)
% Calculate the time from the start of readout to the center (ky==0) of each echo
%
% Inputs:
%   ky      - k-space sampling sequence (in acquisition order)
%   N_echo  - number of echoes
%   esp     - echo spacing (time interval between consecutive k-space samples)
%
% Output:
%   echo_te - TE corresponding to ky=0 for each echo

t = (1:length(ky)) * esp;
t = t - 0.5*esp;

echo_te = zeros(1, N_echo);

for echo_idx = 1:N_echo
    % Extract ky and corresponding time for current echo (sampling with step N_echo)
    ky_echo = ky(echo_idx:N_echo:end);
    t_echo = t(echo_idx:N_echo:end);

    % Check if current echo contains ky=0
    if any(ky_echo == 0)
        % Directly get the time when ky=0
        echo_te(echo_idx) = t_echo(ky_echo == 0);
    else
        % Linear interpolation to find time corresponding to ky=0
        echo_te(echo_idx) = interp1(ky_echo, t_echo, 0, 'linear');
    end
end

end



function outParams2 = func_Synthetic_wenmiao_111219(data)
%% testSynthetic_hernando_111110
%%
%% Test fat-water algorithms on synthetic data
%%
%% Author: Wenmiao Lu
%% Date created: August 19, 2011
%% Date last modified: November 10, 2011

% Add to matlab path
% if ispc
%     BASEPATH = '';
% else
%     BASEPATH = '';
% end
% % BASEPATH = './';
% addpath(BASEPATH);
% addpath([BASEPATH '/multiResSep']);
% addpath ../hernando/create_synthetic/

%% Create some synthetic "true" data, and set acquisition parameters
sx = 128;sy=sx;
imtest = phantom(sx);
threshold = 0.8;
N = 3; TEinit = 1.2e-3; dTE = 2e-3;
TE = TEinit + [0:N-1]*dTE;
trueParams.species(1).amps = imtest.*(imtest<threshold); % Water
trueParams.species(2).amps = imtest.*(imtest>=threshold); % Fat

[X,Y] = meshgrid(linspace(-1,1,sx),linspace(-1,1,sy));
trueParams.fieldmap = 20*randn*ones(sx,sy) + 40*randn*X + 40*randn*Y + 100*randn*X.^2 + 100*randn*Y.^2  + 100*randn*X.*Y.^2 + 100*randn*X.^3  + 100*randn*Y.^3;
trueParams.r2starmap = 0 - 0*imtest;

imDataParams0.TE = TE;
imDataParams0.FieldStrength = 1.5;
imDataParams0.PrecessionIsClockwise = -1;
imDataParams0 = data;

%% Set recon parameters
% General parameters
algoParams.species(1).name = 'water';
algoParams.species(1).frequency = 0;
algoParams.species(1).relAmps = 1;
algoParams.species(2).name = 'fat';
algoParams.species(2).frequency = -[-3.80, -3.40, -2.60, -1.94, -0.39, 0.60];
algoParams.species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048];


% Simulate data
% imDataParams = createSynthetic_imageSpace( imDataParams0, algoParams, trueParams );
imDataParams = imDataParams0;

% Algorithm-specific parameters
algoParams.size_clique = 1; % Size of MRF neighborhood (1 uses an 8-neighborhood, common in 2D)
algoParams.range_r2star = [0 120]; % Range of R2* values
algoParams.NUM_R2STARS = 11; % Numbre of R2* values for quantization
algoParams.range_fm = [-400 400]; % Range of field map values
algoParams.NUM_FMS = 301; % Number of field map values to discretize
algoParams.NUM_ITERS = 40; % Number of graph cut iterations
algoParams.SUBSAMPLE = 2; % Spatial subsampling for field map estimation (for speed)
algoParams.DO_OT = 1; % 0,1 flag to enable optimization transfer descent (final stage of field map estimation)
algoParams.LMAP_POWER = 2; % Spatially-varying regularization (2 gives ~ uniformn resolution)
algoParams.lambda = 0.05; % Regularization parameter
algoParams.LMAP_EXTRA = 0.05; % More smoothing for low-signal regions
algoParams.TRY_PERIODIC_RESIDUAL = 0;
THRESHOLD = 0.01;

%% Recon -- graph cut
%% (Hernando D, Kellman P, Haldar JP, Liang ZP. Robust water/fat separation in the presence of large
%% field inhomogeneities using a graph cut algorithm. Magn Reson Med. 2010 Jan;63(1):79-90.)
outParams2 = fw_3point_wm_goldSect( imDataParams, algoParams );
%outParams = fw_i2cm1i_3pluspoint_hernando_graphcut( imDataParams, algoParams );

% $$$ %% Recon -- mixed fit for phase error correction
% $$$ % Initialize mixed fitting to graph cut solution
% $$$ algoParams.fieldmap = outParams.fieldmap;
% $$$ algoParams.r2starmap = outParams.r2starmap;
% $$$ algoParams.NUM_MAGN = 1;
% $$$ algoParams.THRESHOLD = 0.04;
% $$$ algoParams.range_r2star = [0 200];
% $$$
% $$$ % Do mixed fitting
% $$$ outParamsMixed = fw_i2xm1c_3pluspoint_hernando_mixedfit( imDataParams, algoParams )




end


function outParams2 = func_Synthetic_wenmiao_111219_newPeakProportion(data)
%% testSynthetic_hernando_111110
%%
%% Test fat-water algorithms on synthetic data
%%
%% Author: Wenmiao Lu
%% Date created: August 19, 2011
%% Date last modified: November 10, 2011

% Add to matlab path
% if ispc
%     BASEPATH = '';
% else
%     BASEPATH = '';
% end
% % BASEPATH = './';
% addpath(BASEPATH);
% addpath([BASEPATH '/multiResSep']);
% addpath ../hernando/create_synthetic/

%% Create some synthetic "true" data, and set acquisition parameters
sx = 128;sy=sx;
imtest = phantom(sx);
threshold = 0.8;
N = 3; TEinit = 1.2e-3; dTE = 2e-3;
TE = TEinit + [0:N-1]*dTE;
trueParams.species(1).amps = imtest.*(imtest<threshold); % Water
trueParams.species(2).amps = imtest.*(imtest>=threshold); % Fat

[X,Y] = meshgrid(linspace(-1,1,sx),linspace(-1,1,sy));
trueParams.fieldmap = 20*randn*ones(sx,sy) + 40*randn*X + 40*randn*Y + 100*randn*X.^2 + 100*randn*Y.^2  + 100*randn*X.*Y.^2 + 100*randn*X.^3  + 100*randn*Y.^3;
trueParams.r2starmap = 0 - 0*imtest;

imDataParams0.TE = TE;
imDataParams0.FieldStrength = 1.5;
imDataParams0.PrecessionIsClockwise = -1;
imDataParams0 = data;

%% Set recon parameters
% General parameters
algoParams.species(1).name = 'water';
algoParams.species(1).frequency = 0;
algoParams.species(1).relAmps = 1;
algoParams.species(2).name = 'fat';
algoParams.species(2).frequency = -[-3.80, -3.40, -2.60, -1.94, -0.39, 0.60];
algoParams.species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048];
algoParams.species(2).relAmps = [ 0.0252    0.0935    0.0514    0.0030    0.0386    0.0479];
algoParams.species(2).relAmps = [0.0252    0.0935    0.0514    0.0030    0.0386    0.0479] / sum([0.0252    0.0935    0.0514    0.0030    0.0386    0.0479]);


% Simulate data
% imDataParams = createSynthetic_imageSpace( imDataParams0, algoParams, trueParams );
imDataParams = imDataParams0;

% Algorithm-specific parameters
algoParams.size_clique = 1; % Size of MRF neighborhood (1 uses an 8-neighborhood, common in 2D)
algoParams.range_r2star = [0 120]; % Range of R2* values
algoParams.NUM_R2STARS = 11; % Numbre of R2* values for quantization
algoParams.range_fm = [-400 400]; % Range of field map values
algoParams.NUM_FMS = 301; % Number of field map values to discretize
algoParams.NUM_ITERS = 40; % Number of graph cut iterations
algoParams.SUBSAMPLE = 2; % Spatial subsampling for field map estimation (for speed)
algoParams.DO_OT = 1; % 0,1 flag to enable optimization transfer descent (final stage of field map estimation)
algoParams.LMAP_POWER = 2; % Spatially-varying regularization (2 gives ~ uniformn resolution)
algoParams.lambda = 0.05; % Regularization parameter
algoParams.LMAP_EXTRA = 0.05; % More smoothing for low-signal regions
algoParams.TRY_PERIODIC_RESIDUAL = 0;
THRESHOLD = 0.01;

%% Recon -- graph cut
%% (Hernando D, Kellman P, Haldar JP, Liang ZP. Robust water/fat separation in the presence of large
%% field inhomogeneities using a graph cut algorithm. Magn Reson Med. 2010 Jan;63(1):79-90.)
outParams2 = fw_3point_wm_goldSect( imDataParams, algoParams );
%outParams = fw_i2cm1i_3pluspoint_hernando_graphcut( imDataParams, algoParams );

% $$$ %% Recon -- mixed fit for phase error correction
% $$$ % Initialize mixed fitting to graph cut solution
% $$$ algoParams.fieldmap = outParams.fieldmap;
% $$$ algoParams.r2starmap = outParams.r2starmap;
% $$$ algoParams.NUM_MAGN = 1;
% $$$ algoParams.THRESHOLD = 0.04;
% $$$ algoParams.range_r2star = [0 200];
% $$$
% $$$ % Do mixed fitting
% $$$ outParamsMixed = fw_i2xm1c_3pluspoint_hernando_mixedfit( imDataParams, algoParams )




end



function x = vec(x)
x = x(:);
end
function outParams = func_fw_i2cm0c_3pluspoint_tsaojiang_wrapper(data, fatPeakPpm, fatRelAmps)
if nargin < 2 || isempty(fatPeakPpm)
    fatPeakPpm = [-3.40, -2.60, -0.39, 0.60];
end

algoParams = [];
algoParams.species(1).name = 'water';
algoParams.species(1).frequency = 0;
algoParams.species(1).relAmps = 1;

if nargin >= 3 && ~isempty(fatRelAmps)
    fatRelAmps = fatRelAmps(:).';
    fatRelAmps = fatRelAmps / sum(fatRelAmps);
    algoParams.species(2).name = 'fat';
    algoParams.species(2).frequency = fatPeakPpm;
    algoParams.species(2).relAmps = fatRelAmps;
else
    for iter_peak = 1:numel(fatPeakPpm)
        species_ind = iter_peak + 1;
        algoParams.species(species_ind).name = sprintf('fat_%+.2fppm', fatPeakPpm(iter_peak));
        algoParams.species(species_ind).frequency = fatPeakPpm(iter_peak);
        algoParams.species(species_ind).relAmps = 1;
    end
end

outParamsTsao = fw_i2cm0c_3pluspoint_tsaojiang(data, algoParams);

outParams = outParamsTsao;
outParams.species(1).amps = outParamsTsao.species(1).amps;
if numel(outParamsTsao.species) >= 2
    outParams.species(2).amps = sum(cat(3, outParamsTsao.species(2:end).amps), 3);
else
    outParams.species(2).amps = zeros(size(outParamsTsao.species(1).amps));
end
dTE_tsao = median(diff(sort(unique(data.TE))));
outParams.fieldmap = angle(outParamsTsao.phasemap * exp(-1i * 2 * pi / 3)) / (2 * pi * dTE_tsao);
end

function outParams = run_dixon_with_optional_fieldmap_prior(data, method, hernando_dixon_func, wenmiao_dixon_func, tsaojiang_dixon_func, fat_model_spec_for_ls)

switch method
    case 'hernando_graphcut'
        if data.DO_MIXED_FIT
            [~,outParams] = hernando_dixon_func(data);
        else
            outParams = hernando_dixon_func(data);
        end
    case 'MultiSeedRegionGrowing'
        outParams = func_MultiSeedRegionGrowing(data);
    case 'wenmiao'
        if isfield(data,'fieldmap_init') && ~isempty(data.fieldmap_init)
            warning('Fieldmap prior is not wired into wenmiao Dixon here; falling back to the standard wenmiao solver.');
        end
        outParams = wenmiao_dixon_func(data);
    case 'tsaojiang_Hierarchical_IDEAL'
        if isfield(data,'fieldmap_init') && ~isempty(data.fieldmap_init)
            warning('Fieldmap prior is not wired into tsaojiang_Hierarchical_IDEAL here; falling back to the standard Tsao-Jiang solver.');
        end
        outParams = tsaojiang_dixon_func(data);
    case 'linear_ls_known_fieldmap'
        outParams = func_dixon_least_squares_known_fieldmap(data, fat_model_spec_for_ls);
    otherwise
        error('Unsupported Dixon method: %s', method);
end
end

function outParams = func_dixon_least_squares_known_fieldmap(data, fat_model_spec)
if ~isfield(data,'fieldmap_init') || isempty(data.fieldmap_init)
    error('linear_ls_known_fieldmap requires data.fieldmap_init.');
end

images = squeeze(data.images);
if ndims(images) ~= 3
    error('Expected squeezed data.images to have size [Nx, Ny, N_TE], got ndims=%d.', ndims(images));
end

fieldmap_use = sanitize_init_map(data.fieldmap_init, 0, []);
phase_flag = 1;
if isfield(data,'PrecessionIsClockwise') && ~isempty(data.PrecessionIsClockwise)
    phase_flag = sign(data.PrecessionIsClockwise);
    if phase_flag == 0
        phase_flag = 1;
    end
end

if isfield(data,'linear_ls_tikhonov_lambda') && ~isempty(data.linear_ls_tikhonov_lambda)
    lambda_tikhonov = data.linear_ls_tikhonov_lambda;
else
    lambda_tikhonov = 0;
end

[water_image, fat_image, resnorm_map, residual_rms_map] = dixon_least_squares( ...
    images, data.TE, fieldmap_use, data.FieldStrength, fat_model_spec, phase_flag, lambda_tikhonov);

outParams = [];
outParams.species(1).amps = water_image;
outParams.species(2).amps = fat_image;
outParams.fieldmap = fieldmap_use;
outParams.resnorm = resnorm_map;
outParams.residual_rms = residual_rms_map;
end

function out = sanitize_init_map(in, fill_value, clamp_range)
out = double(in);
out(~isfinite(out)) = fill_value;
if nargin >= 3 && ~isempty(clamp_range)
    out = min(max(out, clamp_range(1)), clamp_range(2));
end
end
