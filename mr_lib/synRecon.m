
function vNew = synRecon(mrObjs,arg)
%{
Synergistic reconstruction of multi-contrast undersampled MRI data 

Reference:
Mehranian et al, "Multi?modal synergistic PET and MR reconstruction 
using mutually weighted quadratic priors", Magn Reson Med. 2018;1–15.

Abolfazl Mehranian @ King’s College London, 2018
%}
numberOfDatasets = length(mrObjs);
opt.global_niter = 20;
opt.Display = 45;
opt.gaussianKernel = 1.4; %for reducing noise in joint coefficients
opt.MrPriorType = 'Quadratic';
opt.MrPreCompWeights = 1;
opt.message = [];
opt.ReconUnderSampledkSpace = ones(numberOfDatasets,1);
opt.SENSE_niter = 2*ones(numberOfDatasets,1);
opt.MrRegularizationParameter = 0*ones(numberOfDatasets,1);
opt.MrSigma = 0.03*ones(numberOfDatasets,1);

opt = getFiledsFromUsersOpt (opt,arg);

senseOpts = struct('SENSE_niter',[],'MrRegularizationParameter',[],'MrPreCompWeights',[]);

[RHS,wMr,vOld,vNew] = deal(cell(numberOfDatasets,1));
for i = 1:numberOfDatasets
    senseOpts(i).SENSE_niter = opt.SENSE_niter(i);
    senseOpts(i).MrRegularizationParameter = opt.MrRegularizationParameter(i);    
    RHS{i} = mrObjs{i}.FH(mrObjs{i}.ftKspaceData);
    wMr{i} = 1;
    vOld{i} = 0*RHS{i};
end

GaussFilt = @(x) gauss3DFilter(abs(x),[1,1,1],opt.gaussianKernel,1);
if opt.Display, figure, end   

for it = 1: opt.global_niter
    % reconstruct
    for i = 1:numberOfDatasets
        senseOpts(i).MrPreCompWeights = wMr{i};
        vNew{i} = mrObjs{i}.SENSE_CG(senseOpts(i),vOld{i},RHS{i});
    end
    % map native spaces
    vNewMapped = cell(numberOfDatasets,numberOfDatasets-1);
    ids = 1:numberOfDatasets;
    for i = 1:numberOfDatasets
        idx = ids(ids~=i);
        for ii = 1:length(idx)
             vNewMapped{i,ii} = mrObjs{idx(ii)}.mapNativeSpaceAToNativeSpaceB(mrObjs{i},vNew{idx(ii)});
        end
    end
    % calculate joint weighting coefficients
    for i = 1:numberOfDatasets
        tmp = GaussFilt(vNew{i});
        wMri = mrObjs{i}.Prior.W_JointEntropy(tmp./max(tmp(:)),opt.MrSigma(i));
        idx = ids(ids~=i);
        for ii = 1:length(idx)
            tmp = GaussFilt(vNewMapped{i,ii});
            wMri = wMri.*mrObjs{i}.Prior.W_JointEntropy(tmp./max(tmp(:)),opt.MrSigma(idx(ii)));
        end
        wMr{i} = wMri./repmat(sum(wMri,2),[1,mrObjs{i}.Prior.nS]);
    end  
    % display
    if opt.Display
        wMri = mrObjs{i}.Prior.UndoImCrop(reshape(sum(wMri,2),mrObjs{i}.Prior.CropedImageSize));
        for i = 1:numberOfDatasets
            subplot(1,numberOfDatasets+1,i),imshow(abs(vNew{i}(:,:,opt.Display)),[]), title(['#' num2str(i)]), drawnow
        end
        subplot(1,numberOfDatasets+1,i+1),imshow(wMri(:,:,opt.Display),[]), title('Joint weights'), drawnow
        suptitle([opt.message ', Iter: ' num2str(it)]), drawnow
    end
    clear vNewMapped wMri
    vOld = vNew;
    fprintf(2,'Iteration: %d\n',it);
end

% map into Nifti reference space
for i = 1:numberOfDatasets
    vNew{i} = mrObjs{i}.mapNativeSpaceToNiftiRefSpace(abs(vNew{i}));
end
