

function [Data,hdr]=getwixdata(twixh,dataType,arg)

% dataType:
%     .image:         image scan
%     .noise:         for noise scan
%     .phasecor:      phase correction scan
%     .phasestab:     phase stabilization scan
%     .phasestabRef0: phase stab. ref. (MDH_REFPHASESTABSCAN && !MDH_PHASESTABSCAN)
%     .phasestabRef1: phase stab. ref. (MDH_REFPHASESTABSCAN &&  MDH_PHASESTABSCAN)
%     .refscan:       parallel imaging reference scan
%     .refscanPC:     phase correction scan for reference data
%     .refscanPS:     phase stabilization scan for reference data
%     .refscanPSRef0: phase stab. ref scan for reference data
%     .refscanPSRef1: phase stab. ref scan for reference data
%     .RTfeedback:    realtime feedback data
%     .vop:           vop rf data
if nargin==1
    dataType = 'image';
end
arg.null =0;

% opt
opt.RemoveOS = 0;
opt.DoAverage = 1;
opt.AverageReps = 0;
opt.AverageSets = 1;
opt.IgnoreSeg = 1;
opt.RampSampRegrid = 1;
opt.SkipToFirstLine = 0;
opt.DoRawDataCorrect = 1;

opt = getFiledsFromUsersOpt(opt,arg);

if ~isfield(twixh,dataType)
    error('%s is unavailable\n',dataType);
end

%%
twixh.(dataType).flagRemoveOS = opt.RemoveOS;
twixh.(dataType).flagDoAverage = opt.DoAverage;
twixh.(dataType).flagAverageReps = opt.AverageReps;
twixh.(dataType).flagAverageSets = opt.AverageSets;
twixh.(dataType).flagIgnoreSeg = opt.IgnoreSeg;
if ~isempty(twixh.(dataType).rampSampTrj)
    twixh.(dataType).flagRampSampRegrid = opt.RampSampRegrid;
end
twixh.(dataType).flagSkipToFirstLine = 0;
if ~isempty(twixh.(dataType).RawDataCorrectionFactors)
    twixh.(dataType).flagDoRawDataCorrect = opt.DoRawDataCorrect;
end

data = squeeze(twixh.(dataType)());
% check if multi-channel, then permute dimension to Col, Lin, Par,...
if twixh.(dataType).NCha>1
    data = permute(data,[1,3,4,2,5:ndims(data)]);
end

hdr.img.matrixSize = [twixh.hdr.Meas.NImageCols,twixh.hdr.Meas.NImageLins,twixh.hdr.Meas.NImagePar];
hdr.img.FOV = [twixh.hdr.MeasYaps.sSliceArray.asSlice{1}.dReadoutFOV,twixh.hdr.MeasYaps.sSliceArray.asSlice{1}.dPhaseFOV,twixh.hdr.MeasYaps.sSliceArray.asSlice{1}.dThickness];
hdr.img.voxelSize = [hdr.img.FOV./hdr.img.matrixSize];
sPosition = twixh.hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition;
position = [nan, nan, nan];
if isfield(sPosition,'dSag'), position(1) = sPosition.dSag; end
if isfield(sPosition,'dCor'), position(2) = sPosition.dCor; end
if isfield(sPosition,'dTra'), position(3) = sPosition.dTra; end
hdr.img.position = position;

hdr.data.ftSize = [twixh.hdr.Config.NColMeas, twixh.hdr.Config.NPeftLen,twixh.hdr.Config.NPaftLen];
hdr.data.sqzSize = twixh.(dataType).sqzSize([1,3,4,2,5:ndims(data)]);
hdr.data.sqzDim = twixh.(dataType).sqzDims([1,3,4,2,5:ndims(data)]);
if isfield(twixh.hdr.MeasYaps.sKSpace,'dSliceOversamplingForDialog')
    hdr.data.sliceOverSampling = twixh.hdr.MeasYaps.sKSpace.dSliceOversamplingForDialog;
end
hdr.data.centerOfkSpace = [twixh.(dataType).centerCol(1),twixh.(dataType).centerLin(1),twixh.(dataType).centerPar(1)];
hdr.data.NCol = twixh.(dataType).NCol;
hdr.data.NLin = twixh.(dataType).NLin;
hdr.data.NCha = twixh.(dataType).NCha;
hdr.data.NPar = twixh.(dataType).NPar;
hdr.data.NRep = twixh.(dataType).NRep;
hdr.data.NSeg = twixh.(dataType).NSeg;
hdr.data.NSli = twixh.(dataType).NSli;
hdr.data.NAve = twixh.(dataType).NAve;
hdr.data.NAcq = twixh.(dataType).NAcq;
hdr.data.is3d = twixh.(dataType).NPar>1;
hdr.data.RemoveOS = opt.RemoveOS;

[data, centerOfkSpaceMask] = centrelizeKspace(data,hdr.data.sqzSize,hdr.data.ftSize,hdr.data.centerOfkSpace);
Data.ftKspaceData = data;
Data.centerOfkSpaceMask =centerOfkSpaceMask;


function [kSpaceDataCentre, centerOfkSpaceMask] = centrelizeKspace(kSpaceData,sqzSzie,FtMatrixSize,centerOfkSpace)

% sqzSzie: size of k-space data returned by twix
% FtMatrixSize: size of k-space for FT trasnfrom
% centerOfkSpace: center of k-space data retured by twix

centerOfkSpaceMask = true(FtMatrixSize(1:3));
delta = floor(FtMatrixSize(1:3)/2+1- centerOfkSpace) ;
if all(delta ==0)
    kSpaceDataCentre = kSpaceData;
    return;
end

kSpaceDataCentre = zeros([FtMatrixSize,sqzSzie(4:end)],'single');

if delta(1)>0
    centerOfkSpaceMask(1:delta(1),:,:) = 0;
elseif delta(1)<0
    centerOfkSpaceMask(FtMatrixSize(1)+delta(1)+1:end,:,:) = 0;
end
if delta(2)>0
    centerOfkSpaceMask(:,1:delta(2),:) = 0;
elseif delta(2)<0
    centerOfkSpaceMask(:,FtMatrixSize(2)+delta(2)+1:end,:) = 0;
end
if delta(3)>0
    centerOfkSpaceMask(:,:,1:delta(3)) = 0;
elseif delta(3)<0
    centerOfkSpaceMask(:,:,FtMatrixSize(3)+delta(3)+1:end) = 0;
end
% %

kSpaceDataCentre(1:sqzSzie(1),1:sqzSzie(2),1:sqzSzie(3),:,:,:,:) = kSpaceData;

for i = 1:length(delta)
    if delta(i)~=0
        kSpaceDataCentre = circshift(kSpaceDataCentre,delta(i),i);
    end
end


