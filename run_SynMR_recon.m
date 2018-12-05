


% Requirements:
% 1) DICOM MR images reconstructed by the MR console are required to obtain MR
% image coordinates, use "spm_dicom2nifti.m" to convert DICOM images into Nifti
% 
% 2) NIFTI and SPM(12) packages


addpath(genpath('C:\SynMR-Recon\'));
addpath('c:\nifti\')
addpath('c:\spm12\')

io.rawdata_dir = 'c:\SynMR-Recon\example_dataset_synMr\';
io.results = [io.rawdata_dir 'results\']; mkdir(io.results);
io.dataset(1).nii = [io.rawdata_dir 'MPRAGE_VD_x3.nii'];
io.dataset(1).dat = [io.rawdata_dir 'meas_MID00453_FID05990_mprage_VD_x3.dat'];
io.dataset(2).nii = [io.rawdata_dir 'MPRAGE_VD_x14.nii'];
io.dataset(2).dat = [io.rawdata_dir 'meas_MID00458_FID05995_mprage_VD_x14.dat'];
io.dataset(3).nii = [io.rawdata_dir 'T2PREP_VD_x3.nii'];
io.dataset(3).dat = [io.rawdata_dir 'meas_MID00451_FID05988_t2prep_VD_x3.dat'];
io.dataset(4).nii = [io.rawdata_dir 'T2PREP_VD_x14.nii'];
io.dataset(4).dat = [io.rawdata_dir 'meas_MID00457_FID05994_t2prep_VD_x14.dat'];
io.dataset(1).MrInfo = getNiftiDataInfo(io.dataset(1).nii);
io.dataset(2).MrInfo = getNiftiDataInfo(io.dataset(2).nii);
io.dataset(3).MrInfo = getNiftiDataInfo(io.dataset(3).nii);
io.dataset(4).MrInfo = getNiftiDataInfo(io.dataset(4).nii);

%% Construct MR objects ---------------------------------
for i = 1:length(io.dataset)
    twixt = mapVBVD(io.dataset(i).dat);
    twixt = twixt{2};
    arg.RemoveOS = 1;
    [data,hdr]=getwixdata(twixt,'image',arg);
    % scale data to ease selection of the regularization parameter
    data.ftKspaceData = 1e5*data.ftKspaceData./max(abs(data.ftKspaceData(:)));
    data.centerOfkSpaceMask = data.ftKspaceData(:,:,:,1)~=0;
    
    sensOpt.null = 0;
    mrObj = mrReconClass(sensOpt,data,hdr);
    % is undersampled
    mrObj.kUnderSampling.is = 1;
    % set reference affine transfromation obtained from nifti files
    mrObj.setRefAffineTransfMatrix(io.dataset(i).MrInfo);
    % calculate coil sensitivity map
    arg.useCenterOfKspace = 15;
    mrObj.getCoilSensitivityMap(arg)
    % save mr objects
    [~,io.dataset(i).name] = fileparts(io.dataset(i).nii);
    io.dataset(i).mat = [io.results io.dataset(i).name '.mat'];
    save(io.dataset(i).mat,'mrObj','-v7.3');
end

save([io.results 'io.mat'],'io');
%% Standard reconstructions ---------------------------------

for i = 1:length(io.dataset)
    load(io.dataset(i).mat)
    % 1- recon Zero_fill sum of square image
    [~, imgZeroFill] = mrObj.reconCoilsImages();
    imgZeroFill_mrs = mrObj.mapNativeSpaceToNiftiRefSpace(abs(imgZeroFill));
    % 2- basic SENSE reconstruction
    sensOpt.niter = 10;
    imgSENSE = mrObj.SENSE_CG(sensOpt);
    imgSENSE_mrs = mrObj.mapNativeSpaceToNiftiRefSpace(abs(imgSENSE));
     % 3-TV regualrization
    optTV.imCropFactor = [7,0,0];
    mrObj.BuildNativeResolutionPrior(optTV)
    optTV.MrADMMPenaltyParameter = 0.1;
    optTV.MrRegularizationParameter = 0.05;
    optTV.ADMM_niter = 5;
    optTV.SENSE_niter = 2;
    imgTV = mrObj.SENSE_TV_ADMM(optTV);
    imgTV_mrs = mrObj.mapNativeSpaceToNiftiRefSpace(abs(imgTV));
    % save
    [MrInfo,MrNifti] = getNiftiDataInfo(io.dataset(i).nii);
    save_nifti(['imgZF_', io.dataset(i).name],1e3*imgZeroFill_mrs, io.results ,MrNifti,MrInfo,0)
    save_nifti(['imgSENSE_', io.dataset(i).name],1e3*imgSENSE_mrs, io.results ,MrNifti,MrInfo,0)
    save_nifti(['imgTV_', io.dataset(i).name],1e3*imgTV_mrs, io.results ,MrNifti,MrInfo,0)
end
%% Synergistic reconstruction of "MPRAGE_VD_x3" "MPRAGE_VD_x14","T2PREP_VD_x14" ---------------------------------
id = [1,2,4];
% Initialise Prior objects
mrObjs = cell(length(id),1);
opt.imCropFactor = [7,0,0];
opt.sWindowSize = 3; % Neighborhood size
for i = 1:length(id)
    load(io.dataset(id(i)).mat)
    mrObj.BuildNativeResolutionPrior(opt);
    mrObjs{i} = mrObj;
    clear mrObj
end

opt.global_niter = 10;
opt.Display = 45;
opt.MrPriorType = 'Quadratic';
opt.MrSigma = [0.03, 0.03,0.03];
opt.SENSE_niter = [4,4,4];
opt.MrRegularizationParameter = [15,15,15];
opt.message = ['"MPRAGE-VD-x3" "MPRAGE-VD-x14","T2PREP-VD-x14"'];

vNew = synRecon(mrObjs,opt);
% save 
for i = 1:length(id)
    [MrInfo,MrNifti] = getNiftiDataInfo(io.dataset(id(i)).nii);
    save_nifti(['Syn_', io.dataset(id(i)).name],1e3*vNew{i}, io.results ,MrNifti,MrInfo,0)
end

%% Compare reconstruction, map images to the 1st dataset's image space ---------------------------------
id = [1,2,3,4];
[imgZF,imgSENSE,imgTV,Syn] = deal(cell(length(id),1));
map = @(img,i)flip(permute(mapSpaceAToSpaceBspm(img,io.dataset(i).MrInfo,io.dataset(1).MrInfo),[2,1,3]),1);
for i = 1:length(id)
    tmp = load_untouch_nii([io.results,'imgZF_', io.dataset(id(i)).name,'.nii']);
    imgZF{i} = map(tmp.img,i);
    tmp = load_nii([io.results,'imgSENSE_', io.dataset(id(i)).name,'.nii']);
    imgSENSE{i} = map(tmp.img,i);
    tmp = load_nii([io.results,'imgTV_', io.dataset(id(i)).name,'.nii']);
    imgTV{i} = map(tmp.img,i);
    if i==3, continue, end %
    tmp = load_nii([io.results,'Syn_', io.dataset(id(i)).name,'.nii']);
    Syn{i} = map(tmp.img,i);
end
%%
figure
p = 30;
L1 = 50; H1 = 1000;
L2 = 50; H2 = 1000;
%T1
ax(1) = subplot(241);imshow(imgSENSE{1}(:,:,p),[L1,H1]), title('SENSE (3x)')
ax(2) = subplot(242);imshow(imgSENSE{2}(:,:,p),[L1,H1]), title('SENSE (14x)')
ax(3) = subplot(243);imshow(imgTV{2}(:,:,p),[L1,H1]), title('TV-SENSE (14x)')
ax(4) = subplot(244);imshow(Syn{2}(:,:,p),[L1,H1]),title('Synergistic (14x)')
%T2
ax(5) = subplot(245);imshow(imgSENSE{3}(:,:,p),[L2,H2])
ax(6) = subplot(246);imshow(imgSENSE{4}(:,:,p),[L2,H2])
ax(7) = subplot(247);imshow(imgTV{4}(:,:,p),[L2,H2])
ax(8) = subplot(248);imshow(Syn{4}(:,:,p),[L2,H2])
linkaxes(ax,'xy')
zoom(1.2)
