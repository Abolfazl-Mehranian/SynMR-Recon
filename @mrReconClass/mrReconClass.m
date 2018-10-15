
classdef mrReconClass < matlab.mixin.Copyable
    
    
    properties (SetAccess = public)
        
        isSimulation
        is3D
        hdr
        kScaleFactor
        kUnderSampling
        CentralizingMask
        PhaseFoV
        ReadFoV
        dThickness
        Coil
        Prior
        nReps
        SR
        ftKspaceData
        ftKspaceDataSize
        ftImgSize
        reconImgSize
        reconVoxelSize
        refAffineTransfMatrix
    end
    
    methods
        function ObjMRI = mrReconClass(varargin)
            % call the constructor
            ObjMRI.isSimulation     = 0;
            ObjMRI.is3D             = [];
            ObjMRI.kScaleFactor     = 1;
            ObjMRI.CentralizingMask             = 1;
            ObjMRI.kUnderSampling.is            = 0;
            ObjMRI.kUnderSampling.trajectory    = 1;
            ObjMRI.kUnderSampling.factor        =0;
            ObjMRI.kUnderSampling.method        = 'radial';
            ObjMRI.Coil.N                       = 1;
            ObjMRI.Coil.sensitivityMap           = [];
            ObjMRI.Coil.estimationMethod         = 'rsos';%'Walsh'
            ObjMRI.Coil.supportMask              = [];          
            ObjMRI.SR.is = 0;
            ObjMRI.SR.imgRefH =[];
            ObjMRI.SR.imgRefL =[];
            ObjMRI.SR.method = 'matlab-imref3d'; %matlab
            ObjMRI.SR.interp = 'linear';
            ObjMRI.SR.PSF = [0 0 0]; %mm anisotric Gaussian
           
            if ~isempty(varargin) && isstruct(varargin{1})
                % update object's properties from user's options
                ObjMRI = getFiledsFromUsersOpt(ObjMRI,varargin{1});
            else
                error('input is required')
            end
            if ObjMRI.isSimulation
                
            else
                data  = varargin{2};
                hdr  = varargin{3};
                getwixMRdata(ObjMRI,data,hdr);
            end
        end
    end
    
    
    
    methods (Access = public)
        % ///////////////////////////////////////////////////////////////////////////////
        %                           DATA PROCESSING FUNCTIONS
        % ///////////////////////////////////////////////////////////////////////////////
        function ObjMRI = getwixMRdata(ObjMRI,data,hdr)            
            ObjMRI.hdr = hdr;
            ObjMRI.Coil.N = hdr.data.NCha;
            ObjMRI.is3D    = hdr.data.is3d;
            ObjMRI.nReps = hdr.data.NRep;
            ObjMRI.ftKspaceData  = data.ftKspaceData;
            ObjMRI.ftKspaceDataSize = hdr.data.ftSize;
            ObjMRI.ftImgSize = ObjMRI.ftKspaceDataSize(1:3);
            ObjMRI.reconImgSize = hdr.img.matrixSize;
            ObjMRI.reconVoxelSize  = hdr.img.voxelSize;
            ObjMRI.CentralizingMask = 1;
            ObjMRI.kUnderSampling.trajectory = data.centerOfkSpaceMask;
        end
        
        function scaledData = kScale(ObjMRI,data)
            ObjMRI.kScaleFactor  = ObjMRI.kScaleFactor ./(max(abs(data(:))));
            scaledData = ObjMRI.kScaleFactor * data;
        end
        
        function doRestropetiveKspaceUndersampling(ObjMRI,arg)
            
            arg.null = 0;
            opt.undersamplingMethod = 'radial'; % 'cartesian', 'spiral', 'vd-random'
            opt.radial_nspokes = 40;
            opt.cartesian_R_pe = 2; % phase encoding
            opt.cartesian_R_se = 1; % slice encoding
            opt = getFiledsFromUsersOpt(opt,arg);
            
            if strcmpi(opt.undersamplingMethod, 'cartesian')%
                ObjMRI.kUnderSampling.cartesian_R_pe = opt.cartesian_R_pe;
                ObjMRI.kUnderSampling.cartesian_R_se = opt.cartesian_R_se;
                Xi = false([ObjMRI.ftKspaceDataSize]);
                
                Xi(:,1:opt.cartesian_R_pe:end,1:opt.cartesian_R_se:end) = 1;
                
            elseif strcmpi(opt.undersamplingMethod, 'radial')% stack-of-radials
                ObjMRI.kUnderSampling.radial_nspokes = opt.radial_nspokes;
                % stack of 2D radial trajectories
                Col = ObjMRI.ftKspaceDataSize(1);
                Lin = ObjMRI.ftKspaceDataSize(2);
                n = max(Col,Lin);
                Theta = linspace(0,pi,opt.radial_nspokes+1);Theta(end) = [];
                Xi = zeros(n,n);
                for theta = Theta
                    t = linspace(-1,1,3*n)*n;
                    X = round(t.*cos(theta)) + n/2+1;
                    Y = round(t.*sin(theta)) + n/2+1;
                    I = find(X>0 & X<=n & Y>0 & Y<=n);
                    X = X(I);
                    Y = Y(I);
                    Xi(X+(Y-1)*n) = 1;
                end
                dx = Col - Lin;
                if dx>0
                    subset = (abs(dx)/2+1):(Col-abs(dx)/2);
                    Xi = Xi(:,subset);
                elseif dx<0
                    subset = (abs(dx)/2+1):(Lin-abs(dx)/2);
                    Xi = Xi(subset,:);
                end
                Xi = repmat(Xi,[1,1,ObjMRI.ftKspaceDataSize(3)]);
            elseif strcmpi(opt.undersamplingMethod, 'spiral')%
                error('not impemented yet')
            elseif strcmpi(opt.undersamplingMethod, 'vd-random')%
                error('not impemented yet')
            else
                error('unknown')
            end
            ObjMRI.kUnderSampling.is = 1;
            ObjMRI.kUnderSampling.method = opt.undersamplingMethod;
            ObjMRI.kUnderSampling.trajectory = Xi;
            
            for rep = 1:ObjMRI.nReps
                for coil = 1:ObjMRI.Coil.N
                    ObjMRI.ftKspaceData(:,:,:,coil,rep) = ObjMRI.ftKspaceData(:,:,:,coil,rep).*ObjMRI.kUnderSampling.trajectory;
                end
            end
            ObjMRI.kUnderSampling.factor = prod(ObjMRI.ftKspaceDataSize)/numel(find(ObjMRI.kUnderSampling.trajectory));
        end
        
        % ///////////////////////////////////////////////////////////////////////////////
        %                           FOURIER ENCODING FUNCTIONS
        % ///////////////////////////////////////////////////////////////////////////////
        
        function m = F0(ObjMRI,x,mask)
            if ObjMRI.isSimulation
                m = fftn(x).*mask;
            else
                m = ifftshift(fftn(fftshift(x))).*mask;%
            end
        end
        
        function x = FH0(ObjMRI,m,mask)
            if ObjMRI.isSimulation
                x = ifftn(m.*mask);
            else
                x = fftshift(ifftn(ifftshift(m.*mask)));
            end
        end
        
        function m = F(ObjMRI,x)
            % perfroms MRI forward operator
            
            if ObjMRI.kUnderSampling.is
                mask = ObjMRI.kUnderSampling.trajectory;
            else
                mask = 1;%ObjMRI.CentralizingMask;
            end
            if ObjMRI.SR.is
                x = ObjMRI.DownSampling(x);
            end
            m = zeros([ObjMRI.ftKspaceDataSize, ObjMRI.Coil.N]);
            
            for i = 1:ObjMRI.Coil.N
                if isempty(ObjMRI.Coil.sensitivityMap)
                    Sen = 1;
                else
                    if ObjMRI.is3D
                        Sen = ObjMRI.Coil.sensitivityMap(:,:,:,i);
                    else
                        Sen = ObjMRI.Coil.sensitivityMap(:,:,i);
                    end
                end
                K = ObjMRI.F0(x.*Sen,mask);
                if ObjMRI.is3D
                    m(:,:,:,i) = K;
                else
                    m(:,:,i) = K;
                end
                
            end
        end
        
        function x = FH(ObjMRI,m)
            % perfroms MRI adjoint operator
            
            if ObjMRI.kUnderSampling.is
                mask = ObjMRI.kUnderSampling.trajectory;
            else
                mask = 1;%ObjMRI.CentralizingMask;
            end
            x = zeros(ObjMRI.ftKspaceDataSize);
            
            for i = 1: ObjMRI.Coil.N
                if isempty(ObjMRI.Coil.sensitivityMap)
                    Sen = 1;
                else
                    if ObjMRI.is3D
                        Sen = ObjMRI.Coil.sensitivityMap(:,:,:,i);
                        K = m(:,:,:,i);
                    else
                        Sen = ObjMRI.Coil.sensitivityMap(:,:,i);
                        K = m(:,:,i);
                    end
                end
                x = x + conj(Sen).*ObjMRI.FH0(K,mask);
            end
            if ObjMRI.SR.is
                x = ObjMRI.UpSampling(x);
            end
        end
        
        function x = FHF(ObjMRI,m)
            % perfroms MRI Gram Matrix opertor
            x = ObjMRI.FH(ObjMRI.F(m));
        end
        
        % ///////////////////////////////////////////////////////////////////////////////
        %                           RECONSTRUCTION FUNCTIONS
        % ///////////////////////////////////////////////////////////////////////////////
        
        function [ImgCoils, ImgCoilsSOS] = reconCoilsImages(ObjMRI,RepKspaceData,mask)
            if nargin==1 || (nargin==3 && isempty(RepKspaceData))
                %for 3D data, if there are 5th dimension,
                %RepKspaceData=1,2,3,.... else all
                RepKspaceData = 1:ObjMRI.nReps;
            end
            if nargin<=2
                mask = 1;
            end
            % ImgCoils from fully-sampled k-space, so use CentralizingMask
            ImgCoils = zeros([ObjMRI.ftKspaceDataSize, ObjMRI.Coil.N,length(RepKspaceData)]);
            for ii = 1:length(RepKspaceData)
                for i = 1:ObjMRI.Coil.N
                    if ObjMRI.is3D
                        ImgCoils(:,:,:,i,ii) = ObjMRI.FH0(ObjMRI.ftKspaceData(:,:,:,i,RepKspaceData(ii)),ObjMRI.CentralizingMask.*mask);
                    else
                        ImgCoils(:,:,i) = ObjMRI.FH0(ObjMRI.ftKspaceData(:,:,i),ObjMRI.CentralizingMask.*mask);
                    end
                end
            end
            if ObjMRI.is3D
                ImgCoilsSOS = zeros([ObjMRI.ftKspaceDataSize, length(RepKspaceData)]);
                for ii = 1:length(RepKspaceData)
                    ImgCoilsSOS(:,:,:,ii) = ObjMRI.RSS(ImgCoils(:,:,:,:,RepKspaceData(ii)));
                end
                if isempty(ObjMRI.Coil.supportMask)
                    ObjMRI.Coil.supportMask = ObjMRI.MakeCoilSupportMask(ImgCoilsSOS(:,:,:,1));
                end
            else
                ImgCoilsSOS = ObjMRI.RSS(ImgCoils);
                if isempty(ObjMRI.Coil.supportMask)
                    ObjMRI.Coil.supportMask = ObjMRI.MakeCoilSupportMask(ImgCoilsSOS);
                end
            end
            
        end
        
        function setPointSpreadFunction(ObjMRI, psf)
            if nargin==1
                psf = 0;
            end
            ObjMRI.SR.PSF = psf;
        end
        
        % ///////////////////////////////////////////////////////////////////////////////
        %                           HELPER FUNCTIONS
        % ///////////////////////////////////////////////////////////////////////////////
        
        function Y = RSS(ObjMRI,X)
            % Calculates root of sum of squares
            if ndims(X)==2 %#ok<ISMAT> % Gradient vectors
                dims = 2;
            else % images
                if ObjMRI.is3D
                    dims = 4;
                else
                    dims = 3;
                end
            end
            Y  = sqrt(sum(abs(X).^2,dims));
        end
        
        function x = Magnitude(~,x)
            x = sqrt(sum(abs(x(:)).^2));
        end
        
        function Y = softThreshold(ObjMRI,Norm,Z,rho,lambda,w)
            
            Norm = repmat(Norm + 1e-5, [1,ObjMRI.Prior.nS] );
            Y =  max(0, Norm - w./(rho/lambda)).*(Z./Norm);
        end
        
        function display(ObjMRI)
            disp(ObjMRI)
            methods(ObjMRI)
        end
        
        function ObjMRI = Revise(ObjMRI,opt)
            % to revise the properties of a given object without
            % re-instantiation
            vfields = fieldnames(opt);
            prop = properties(ObjMRI);
            for i = 1:length(vfields)
                field = vfields{i};
                if sum(strcmpi(prop, field )) > 0
                    ObjMRI.(field) = opt.(field);
                end
            end
        end
        % ///////////////////////////////////////////////////////////////////////////////
        %                           COIL FUNCTIONS
        % ///////////////////////////////////////////////////////////////////////////////
        function setCoilSensitivityMap(ObjMRI,CSM)
            % in the case of a pre-calculated coil map
            ObjMRI.Coil.sensitivityMap = CSM;
        end
        
        function getCoilSensitivityMap(ObjMRI,arg)
            arg.null = 0;
            opt.coilEstimationMethod = 'rsos';
            opt.useCenterOfKspace = 0; %to be done
            opt = getFiledsFromUsersOpt(opt,arg);
            
            dataSize = size(ObjMRI.ftKspaceData);
            if opt.useCenterOfKspace
                
                [x0,y0] = meshgrid(-dataSize(3)/2+1:dataSize(3)/2,-dataSize(2)/2+1:dataSize(2)/2);
                centralMask = (x0.^2+y0.^2)<opt.useCenterOfKspace^2;% 15 is relatively optimal
                centralMask = repmat(centralMask,[1,1,dataSize(1)]);
                centralMask = permute(centralMask,[3,1,2]);
            else
                centralMask = 1;
            end
            [ImgCoils, ImgCoilsSOS] = reconCoilsImages(ObjMRI,[],centralMask);
            
            if strcmpi(opt.coilEstimationMethod,'rsos')
                if ndims(ImgCoils)==5
                    CoilEst = 0*ImgCoils;
                    for i = 1:size(ImgCoils,5)
                        CoilEst(:,:,:,:,i) =  ImgCoils(:,:,:,:,i)./repmat(ImgCoilsSOS(:,:,:,i),[1,1,1,ObjMRI.Coil.N]);
                    end
                    CoilEst = mean(CoilEst,5);
                else
                    CoilEst = ImgCoils./repmat(ImgCoilsSOS,[1,1,1,ObjMRI.Coil.N]);
                end
                CoilEstFilt = 0*CoilEst;
                for i = 1:ObjMRI.Coil.N
                    CoilEstFilt(:,:,:,i) = gauss3DFilter(CoilEst(:,:,:,i),[1,1,1],2);
                end
                CoilEst = CoilEstFilt;
                % mask the CSMs
                for i = 1:ObjMRI.Coil.N
                    CoilEst(:,:,:,i) = CoilEst(:,:,:,i).*ObjMRI.Coil.supportMask;
                end
                ObjMRI.Coil.sensitivityMap = (CoilEst);
            else
                
            end
        end
        
        function mask = MakeCoilSupportMask(ObjMRI,imgCoilSOS)
            
            if ObjMRI.is3D
                mask = imgaussfilt3(imgCoilSOS./max(imgCoilSOS(:)),8);
            else
                mask = imgaussfilt(imgCoilSOS./max(imgCoilSOS(:)),8);
            end
            
            level = graythresh(mask);
            if level==0
                fprintf('MakeCoilSupportMask:: zero threshold level !!\n')
                mask = 1;
                return
            end
            mask = mask >(level*0.25);
            [x,y,z] = meshgrid(-1:1:1);
            r = (x.^2+y.^2+z.^2)<2;
            
            for i = 1:5
                mask = imdilate(single(mask),r);
            end
            if ObjMRI.is3D
                for i=1:size(mask,3)
                    mask(:,:,i) = imfill(mask(:,:,i),'holes');
                end
            else
                mask= imfill(mask,'holes');
            end
            mask = logical(mask);
        end
        
        % ///////////////////////////////////////////////////////////////////////////////
        %                           PRIOR OBJECT FUNCTIONS
        % ///////////////////////////////////////////////////////////////////////////////
        
        function setSuperResolutionImgRef(ObjMRI,imgRefH,imgRefL,opt)
            
            % 1) imgRefH and imgRefH: nifti (.nii) file address of the high
            % resolution image (e.i T1-weighted MR) and the low-resolution
            % image (i.e. one of ASL time serie image), use SPM's nifti
            % function to get info, sampling method is set to 'matlab-spm'
            
            % 2) imgRefH and imgRefH: spatial structre generated by
            % getNifiDataInfo (which calls SPM's nifti function, sampling
            % method is set to 'matlab-spm'
            
            % 3) imgRefH/imgRefH: matlab's imref3d (Reference 3-D image to
            % world coordinate object), sampling method is set to 'matlab-imref3d'
            
            opt.null = 0;
            if ischar(imgRefH) % --> imgRefH and imgRefL are address files
                imgRefH = ObjMRI.getNiftiDataInfo(imgRefH);
                imgRefL = ObjMRI.getNiftiDataInfo(imgRefL);
                ObjMRI.SR.method = 'matlab-spm';
                if isfield(opt,'interp') %0: nearest neighbours, 1: Trilinear
                    ObjMRI.SR.interp = opt.interp;
                else
                    ObjMRI.SR.interp = 0;
                end
            elseif isstruct(imgRefH)
                ObjMRI.SR.method = 'matlab-spm';
                if isfield(opt,'interp') %0: nearest neighbours, 1: Trilinear
                    ObjMRI.SR.interp = opt.interp;
                else
                    ObjMRI.SR.interp = 0;
                end
            elseif strcmpi(class(imgRefH),'imref3d')
                ObjMRI.SR.method = 'matlab-imref3d';
                if isfield(opt,'interp')
                    ObjMRI.SR.interp = opt.interp;
                else
                    ObjMRI.SR.interp = 'linear';
                end
            end
            ObjMRI.SR.imgRefH = imgRefH;
            ObjMRI.SR.imgRefL = imgRefL;
        end
        
        function BuildSuperResolutionPrior(ObjMRI,arg)
            arg.null = 0;
            opt.imCropFactor = [8,8,8];
            opt.sWindowSize = 3;
            opt.lWindowSize = 1;
            opt.imCropHandle = [];
            
            opt = getFiledsFromUsersOpt(opt,arg);
            
            ObjMRI.SR.is = 1;
            if isempty(ObjMRI.SR.imgRefH)
                error('First set setSuperResolutionImgRef');
            end
            if strcmpi(ObjMRI.SR.method,'matlab-imref3d')
                opt.ImageSize = ObjMRI.SR.imgRefH.ImageSize;
                ObjMRI.reconVoxelSize = [ObjMRI.SR.imgRefH.PixelExtentInWorldX, ObjMRI.SR.imgRefH.PixelExtentInWorldY, ObjMRI.SR.imgRefH.PixelExtentInWorldZ];
            elseif strcmpi(ObjMRI.SR.method,'matlab-spm')
                opt.ImageSize = ObjMRI.SR.imgRefH.dim;
                ObjMRI.reconVoxelSize = ObjMRI.SR.imgRefH.pixdim;
            end
            if ~isempty(opt.imCropHandle)
                opt.imCropFactor = 0; % for Siemens mMR
            end
            ObjMRI.Prior = PriorsClass(opt);
        end
        
        function BuildNativeResolutionPrior(ObjMRI,arg)
            ObjMRI.SR.is = 0;
            arg.null = 0;
            opt.imCropFactor = [8,8,0];
            opt.sWindowSize = 3;
            opt.lWindowSize = 1;
            opt.imCropHandle = [];
            opt.isforTV = 0;
            opt = getFiledsFromUsersOpt(opt,arg);
            
            opt.ImageSize = ObjMRI.ftImgSize;
            fprintf('Prior was built for native image size .i.e. %d x %d x %d\n',opt.ImageSize(1),opt.ImageSize(2),opt.ImageSize(3));
            if ~isempty(opt.imCropHandle)
                opt.imCropFactor = 0; % for Siemens mMR
            end
            ObjMRI.Prior = PriorsClass(opt);
        end
        
        function W = getMrGaussianWeights(ObjMRI,MrImage,MrSigma)
            % normalize MrImage to [0,1]
            MrImage = abs(MrImage)./max(abs(MrImage(:)));
            
            W = ObjMRI.Prior.W_GaussianKernel(MrImage,MrSigma);
            W = W./repmat(sum(W,2),[1,ObjMRI.Prior.nS]);
        end
        % ///////////////////////////////////////////////////////////////////////////////
        %                           SPACE MAPPING FUNCTIONS
        % ///////////////////////////////////////////////////////////////////////////////
        
        function imagOut = mapNativeSpaceToNiftiRefSpace(ObjMRI,nativeSpaceImgs)
            % this function maps the reconstcuted images from the kspace into the
            % correposnding nitifi space of the the same images reconstructed by the
            % scanner, need validation for other dataset, currently is hard coded,
            % better to figure out the coordinate system of the images and do
            % resampling
            
            if ObjMRI.hdr.data.RemoveOS
                imagOut = nativeSpaceImgs;
            else
                idx = ObjMRI.hdr.data.NCol/4;
                imagOut = nativeSpaceImgs(idx+1:end-idx,:,:,:);
            end
            nSlice = size(imagOut,3) -ObjMRI.hdr.img.matrixSize(3);
            nPhaseEn = size(imagOut,2) -ObjMRI.hdr.img.matrixSize(2);
            
            imagOut = imagOut(:,nPhaseEn/2+1+(1:ObjMRI.hdr.img.matrixSize(2)),nSlice/2+1+(1:ObjMRI.hdr.img.matrixSize(3)),:);
            imagOut = flip(imagOut,3);
            imagOut = circshift(imagOut,-1,2);
        end
        
        function imagOut = mapNiftiRefSpaceToNativeSpace(ObjMRI,refSpaceImgs)
            
            % inverts the mapNativeSpaceToNiftiRefSpace
            temp = refSpaceImgs;
            temp = circshift(temp,+1,2);
            temp = flip(temp,3);
            
            imagOut = zeros(ObjMRI.hdr.data.ftSize);
            nSlice = size(imagOut,3) -ObjMRI.hdr.img.matrixSize(3);
            nPhaseEn = size(imagOut,2) -ObjMRI.hdr.img.matrixSize(2);
            imagOut(:,nPhaseEn/2+1+(1:ObjMRI.hdr.img.matrixSize(2)),nSlice/2+1+(1:ObjMRI.hdr.img.matrixSize(3)),:) = temp;
            
            if ~ObjMRI.hdr.data.RemoveOS
                idx = ObjMRI.hdr.data.NCol/4;
                temp = imagOut;
                imagOut = zeros(ObjMRI.hdr.data.ftSize);
                imagOut (idx+1:end-idx,:,:,:) = temp;
            end
        end
        
        function imgOut = mapNativeSpaceToMrSpace(ObjMRI, nativeSpaceImgs)
            
            temp = ObjMRI.mapNativeSpaceToNiftiRefSpace(nativeSpaceImgs);
            
            if strcmpi(ObjMRI.SR.method,'matlab-spm')
                imgOut = mapSpaceAToSpaceBspm(temp,ObjMRI.SR.imgRefL,ObjMRI.SR.imgRefH,ObjMRI.SR.interp);
            elseif strcmpi(ObjMRI.SR.method,'matlab-imref3d')
                error('todo');
            end
        end
        
        function imgOut = mapMrSpaceToNativeSpace(ObjMRI, mrImg)
            % inverts mapNativeSpaceToMrSpace
            
            if strcmpi(ObjMRI.SR.method,'matlab-spm')
                temp = mapSpaceAToSpaceBspm(mrImg,ObjMRI.SR.imgRefH,ObjMRI.SR.imgRefL,ObjMRI.SR.interp);
            elseif strcmpi(ObjMRI.SR.method,'matlab-imref3d')
                error('todo');
            end
            imgOut = ObjMRI.mapNiftiRefSpaceToNativeSpace(temp);
            
        end
        
        function imgOut = mapNativeSpaceAToNativeSpaceB(ObjMRI,ObjMRIB,img)
            
            % map native space A into reference nifti A
            imgRefA = ObjMRI.mapNativeSpaceToNiftiRefSpace(abs(img));
            % map nifti space A to nifit space B
            imgARef_RefB = mapSpaceAToSpaceBspm(imgRefA,ObjMRI.refAffineTransfMatrix,ObjMRIB.refAffineTransfMatrix,0);
            % map nifti space B to native space B
            imgOut = ObjMRIB.mapNiftiRefSpaceToNativeSpace(imgARef_RefB);
        end
        
        function setRefAffineTransfMatrix(ObjMRI,t1MrInfo)
            t1MrInfo = rmfield(t1MrInfo,'img');
            ObjMRI.refAffineTransfMatrix = t1MrInfo;
        end
        
        function ImgNew = DownSampling(ObjMRI,Img,interp)
            
            
            if nargin==2
                interp = ObjMRI.SR.interp;
            end
            if strcmpi(ObjMRI.SR.method,'matlab-spm')
                
                if ObjMRI.SR.is == 0 % in case of low-res deconv
                    ImgNew = Img;
                else
                    ImgNew = ObjMRI.mapMrSpaceToNativeSpace(Img);
                end
                
                if any(ObjMRI.SR.PSF)
                    voxelsize = ObjMRI.SR.imgRefL.pixdim;
                    ImgNew = gauss3DFilter(ImgNew,voxelsize,ObjMRI.SR.PSF);
                end
                
            elseif strcmpi(ObjMRI.SR.method,'matlab-imref3d')
                if any(ObjMRI.SR.PSF)
                    voxelsize = [ObjMRI.SR.imgRefH.PixelExtentInWorldX, ObjMRI.SR.imgRefH.PixelExtentInWorldY,ObjMRI.SR.imgRefH.PixelExtentInWorldZ];
                    Img = gauss3DFilter(Img,voxelsize,ObjMRI.SR.PSF);
                end
                ImgNew = ImageResample(Img, ObjMRI.SR.imgRefH, ObjMRI.SR.imgRefL,interp);
                
            end
            
        end
        
        function ImgNew = UpSampling(ObjMRI,Img,interp)
            
            if nargin==2
                interp = ObjMRI.SR.interp;
            end
            
            if ObjMRI.SR.is
                % trim, --< during Usampling and Downsampling there are some
                % very high-intesne planes, need to be removed.
                a = 0*Img;
                pitch = 2;
                if ObjMRI.is3D
                    a(pitch:end-pitch,pitch:end-pitch,pitch:end-pitch) = 1;
                else
                    a(pitch:end-pitch,pitch:end-pitch) = 1;
                end
                Img = a.*Img;
            end
            %
            if strcmpi(ObjMRI.SR.method,'matlab-spm')
                
                if any(ObjMRI.SR.PSF)
                    voxelsize = ObjMRI.SR.imgRefL.pixdim;
                    Img = gauss3DFilter(Img,voxelsize,ObjMRI.SR.PSF);
                end
                
                if ObjMRI.SR.is == 0 % in case of low-res deconv
                    ImgNew = Img;
                else
                    ImgNew = ObjMRI.mapNativeSpaceToMrSpace(Img);
                end
            elseif strcmpi(ObjMRI.SR.method,'matlab-imref3d')
                ImgNew = ImageResample(Img, ObjMRI.SR.imgRefL, ObjMRI.SR.imgRefH,interp);
                if any(ObjMRI.SR.PSF)
                    voxelsize = [ObjMRI.SR.imgRefH.PixelExtentInWorldX, ObjMRI.SR.imgRefH.PixelExtentInWorldY,ObjMRI.SR.imgRefH.PixelExtentInWorldZ];
                    ImgNew = gauss3DFilter(ImgNew,voxelsize,ObjMRI.SR.PSF);
                end
            end
            
            if ObjMRI.SR.is
                a = 0*ImgNew;
                pitch = 3;
                if ObjMRI.is3D
                    a(pitch:end-pitch,pitch:end-pitch,pitch:end-pitch) = 1;
                else
                    a(pitch:end-pitch,pitch:end-pitch) = 1;
                end
                ImgNew = a.*ImgNew;
            end
            
        end
        
        function [out, tmp]= getNiftiDataInfo(ObjMRI,flname)
            % This function requires both NIFIT and SPM libraries
            nifiObj =nifti(flname);
            out.mat = nifiObj.mat;
            
            tmp = load_untouch_nii(flname);
            %             out.img = double(tmp.img);
            out.dim = tmp.hdr.dime.dim(2:tmp.hdr.dime.dim(1)+1);
            out.pixdim = tmp.hdr.dime.pixdim(2:4);
            out.flname = flname;
            % pulling out 4th dimesion,in case of 4D nifti datafile
            out.dim = out.dim(1:3);
            % for compatibility
            out.ImageSize = out.dim;
            out.dataType = tmp.hdr.dime.datatype;
            
        end
        
        % ///////////////////////////////////////////////////////////////////////////////
        %                           OPTIMIZATION ALGORITHMS
        % ///////////////////////////////////////////////////////////////////////////////
        function X = SENSE_CG(ObjMRI,arg,initialEstimate,RHS)
            
            opt.RepKspaceData = 1; %or 'inf' for reconstruction of all Reps
            opt.niter = 10;
            % opt.ReconUnderSampledkSpace = 0;
            opt.MrRegularizationParameter = 0;
            opt.MrPriorType = 'Quadratic'; % 'joint'
            opt.MrPreCompWeights = 1;
            opt.display = 0;
            opt.save = 0;
            
            arg.null = 0;
            
            opt = getFiledsFromUsersOpt(opt,arg);
            
            
            if isinf(opt.RepKspaceData) % if inf, then reconstruct all Reps
                nRe = 1:ObjMRI.nReps;
            else
                nRe = opt.RepKspaceData; % else, reconstruct the specified Rep, default, Rep = 1
            end
            
            if nargin>=3
                img = initialEstimate;
            else
                if ObjMRI.SR.is
                    img = zeros([ObjMRI.SR.imgRefH.ImageSize, length(nRe)],'single');
                else
                    img = zeros([ObjMRI.ftKspaceDataSize length(nRe)],'single');
                end
            end
            
            X = img;
            
            % some reports
            if ObjMRI.SR.is
                fprintf('High-resolution reconstruction\n');
            else
                fprintf('Native-resolution reconstruction\n');
            end
            if ObjMRI.kUnderSampling.is
                fprintf('Undersampled data\n');
            else
                fprintf('Fully sampled\n');
            end
            
            if opt.MrRegularizationParameter, fprintf('%s Regualrization...\n',opt.MrPriorType); end
            
            for i=1:length(nRe)
                fprintf('Rep: #%d\n',nRe(i))
                % K-space data
                Data = ObjMRI.ftKspaceData(:,:,:,:,nRe(i)); % [nCol,nLin,nPar,nCoil,nReps]
                
                % solve Ax = b, where A = FHF + beta*DTwD, b = FH
                if nargin ==4
                    b = RHS(:,:,:,nRe(i));
                else
                    b = ObjMRI.FH(Data);
                end
                
                if opt.MrRegularizationParameter
                    if strcmpi(opt.MrPriorType,'Quadratic')
                        if size(opt.MrPreCompWeights,1) ==1
                            W = 1./ObjMRI.Prior.nS;
                        else
                            W = opt.MrPreCompWeights;
                        end
                        A = @(x, dummy)ObjMRI.FHF(x) + ...
                            opt.MrRegularizationParameter*ObjMRI.Prior.TransGraphGradUndoCrop(W.*ObjMRI.Prior.GraphGradCrop(x));
                    else % for joint weight calculation
                        A = @(x,W)ObjMRI.FHF(x) + ...
                            opt.MrRegularizationParameter*ObjMRI.Prior.TransGraphGradUndoCrop(W.*ObjMRI.Prior.GraphGradCrop(x));
                    end
                else
                    A = @(x,dummy)ObjMRI.FHF(x);
                end
                X(:,:,:,i) = ObjMRI.PCG(b,A,img(:,:,:,i),opt.niter,1,opt);
            end
        end
        
        function [u,report] = SENSE_TV_ADMM(ObjMRI,arg)
            
            if nargin==1, arg.save = 0; end
            
            % default values
            opt.MrADMMPenaltyParameter = 0.1; %\rho
            opt.MrRegularizationParameter = 0.0005; % \lambda
            opt.SigmaParameter = 0;% 0: TV, >0: NCX
            opt.ADMM_niter = 50;
            opt.SENSE_niter = 2;
            opt.display = 0;
            opt.message = [];
            opt.report = 0;
            opt.mask = true(ObjMRI.ftImgSize);
            opt.groundTruth = [];
            
            
            opt = getFiledsFromUsersOpt(opt,arg);
            report = [];
            if opt.report
                Norm = @(x) sqrt(sum(abs(x(:)).^2));
                report.relativeError = zeros(1,opt.ADMM_niter);
                if ~isempty(opt.groundTruth)
                    report.absoluteError = zeros(1,opt.ADMM_niter);
                end
            end
            
            % only for undersampled data
            if ~ObjMRI.kUnderSampling.is
                fprintf('Only undersampled recon is supported\n')
                u = [];
                return
            end
            
            [gamma_u,zu] = deal(zeros(prod(ObjMRI.Prior.CropedImageSize),ObjMRI.Prior.nS,'single'));
            
            FHy = ObjMRI.FH(ObjMRI.ftKspaceData);
            u = zeros(ObjMRI.ftImgSize,'single'); % in native resolution
            
            if opt.display, figure ,end
            display_admm = opt.display;
            opt.display = 0;
            for i = 1:opt.ADMM_niter
                uOld = u;
                RHS = FHy + ObjMRI.Prior.TransGraphGradUndoCrop(opt.MrADMMPenaltyParameter*zu - gamma_u);
                u = ObjMRI.PCG(RHS, @(x, dump) ObjMRI.FH(ObjMRI.F(x)) + ...
                    opt.MrADMMPenaltyParameter*ObjMRI.Prior.TransGraphGradUndoCrop(ObjMRI.Prior.GraphGradCrop(x)), uOld, opt.SENSE_niter, 1,opt);
                
                if opt.SigmaParameter
                    Lu = ObjMRI.RSS(zu);
                    w = exp(-opt.SigmaParameter*Lu./ObjMRI.Magnitude(Lu+eps));
                    w = repmat(w,[1,ObjMRI.Prior.nS]);
                else
                    w = 1;
                end
                
                Du = ObjMRI.Prior.GraphGradCrop(u);
                z_tilde_u = Du + gamma_u/opt.MrADMMPenaltyParameter;
                
                zu = ObjMRI.softThreshold(ObjMRI.RSS(z_tilde_u), z_tilde_u, opt.MrADMMPenaltyParameter,opt.MrRegularizationParameter, w);
                
                gamma_u = gamma_u + opt.MrADMMPenaltyParameter*(Du - zu);
                
                
                if display_admm
                    fprintf('Iteration: %d\n',i);
                    
                    if ObjMRI.is3D
                        drawnow, imshow(abs(u(:,:,display_admm)),[]);
                    else
                        drawnow, imshow(abs(u),[]);
                    end
                    title([opt.message ' Iteration: #' num2str(i)]),drawnow
                end
                
                if opt.report
                    report.relativeError(i) = Norm(u(opt.mask) - uOld(opt.mask))./Norm(uOld(opt.mask));
                    if ~isempty(opt.groundTruth)
                        report.absoluteError(i) = Norm(u(opt.mask) - opt.groundTruth(opt.mask))./Norm(opt.groundTruth(opt.mask));
                    end
                end
            end
        end
        
        function x = PCG(ObjMRI,a,A,x0,Nit,P,arg)
            
            % a: RHS of the normal equation, i.e. FH(mri data)
            % A: FHF() handel object
            % x0: initial image estimate
            % Nit: number of sense iterations
            % P: Precondictioner
            % opt: options
            
            opt.message =[];
            opt.MrPriorType = [];
            opt.save = 0;
            opt.display = 0;
            opt = getFiledsFromUsersOpt(opt,arg);
            
            if opt.display, figure; end
           
            W = 1;
            a = P.*a;
            u = x0./P;
            r = a - P.*A(P.*u,W);
            p = r;
            for i=1:Nit
                uOld = u;
                q = P.*A(P.*p,W);
                alpha = r(:)'*r(:)/(p(:)'*q(:));
                u = uOld + alpha*p;
                rnew = r - alpha*q;
                p = rnew + rnew(:)'*rnew(:)/(r(:)'*r(:))*p;
                r = rnew;
                x = P.*u;
                if opt.display
                    if ObjMRI.is3D
                        drawnow, imshow(abs(u(:,:,opt.display)),[])
                        title([opt.message ' Iteration: #' num2str(i)]),pause(0.1)
                    end
                end
                
            end
        end
        
        function [X, Report] = gradientDescent(ObjMRI,arg,initialEstimate)
            
            arg.null = 0;
            opt.stepSize = 0.2;
            opt.niter = 20;
            opt.RepKspaceData = 1; %or 'inf' for reconstruction of all Reps
            opt.MrRegularizationParameter = 0;
            opt.MrPriorType = 'Quadratic'; % 'joint'
            opt.MrPreCompWeights = 1;
            opt.fSigma = 0.3;
            opt.display = 0;
            opt.save = 0;
            opt.message = [];
            opt.report = 0;
            opt.stepSizeOptimization = 0;
            opt = getFiledsFromUsersOpt(opt,arg);
            
            
            if opt.display, figure; end
            if isinf(opt.RepKspaceData) % if inf, then reconstruct all Reps
                nRe = 1:ObjMRI.nReps;
            else
                nRe = opt.RepKspaceData; % else, reconstruct the specified Rep, default, Rep = 1
            end
            
            if nargin==3
                img = initialEstimate;
            else
                if ObjMRI.SR.is
                    img = zeros([ObjMRI.SR.imgRefH.ImageSize, length(nRe)],'single');
                else
                    img = zeros([ObjMRI.ftKspaceDataSize length(nRe)],'single');
                end
            end
            
            X = img;
            % some reports
            if opt.report
                Report.relativeError = zeros(opt.niter,length(nRe),'single');
                Report.stepSize = zeros(opt.niter,length(nRe),'single');
                Norm = @(x,y) 100*(norm(x(:)-y(:)))/norm(y(:));
                Xp = X;
            else
                Report = [];
            end
            % some report ----------------------------
            if ObjMRI.SR.is
                fprintf('- High-resolution reconstruction.\n');
            else
                fprintf('- Native-resolution reconstruction.\n');
            end
            if ObjMRI.kUnderSampling.is
                fprintf('- Undersampled data.\n');
            else
                fprintf('- Fully sampled.\n');
            end
            if opt.MrRegularizationParameter
                fprintf('- %s Regularization.\n',opt.MrPriorType)
                D = @(x,W) opt.MrRegularizationParameter*ObjMRI.Prior.TransGraphGradUndoCrop(W.*ObjMRI.Prior.GraphGradCrop(x));
            else
                D =@(x,W) 0;
            end
            
            if opt.stepSizeOptimization
                fprintf('- Step size optimization\n');
                col = @(x) x(:);
                H = @(x,W) col(ObjMRI.FHF(x) + D(x,W));
            end
            
            if ~ObjMRI.SR.is, warning('SR is off'); end
            for i=1:length(nRe)
                fprintf('Rep: #%d\n',nRe(i))
                % K-space data
                Data = ObjMRI.ftKspaceData(:,:,:,:,nRe(i)); % [nCol,nLin,nPar,nCoil,nReps]
                
                for j =1:opt.niter
                    % get weights coeffients
                    if strcmpi(opt.MrPriorType , 'Quadratic')
                        W = opt.MrPreCompWeights;
                    elseif strcmpi(opt.MrPriorType , 'Joint')
                        x = X(:,:,:,i);
                        W = opt.MrPreCompWeights .* ObjMRI.Prior.W_GaussianKernel(abs(x)./max(abs(x(:))),opt.fSigma);
                        W = W./repmat(sum(W,2),[1,ObjMRI.Prior.nS]);
                    else
                        W = 1;
                    end
                    % do GD update
                    dPhix =(ObjMRI.FH(ObjMRI.F(X(:,:,:,i))- Data) + D(X(:,:,:,i),W));
                    if opt.stepSizeOptimization
                        stepSize = dot(dPhix(:), dPhix(:))./ dot(dPhix(:), H(dPhix,W));
                        Report.stepSize(j,i) = stepSize;
                    else
                        stepSize = opt.stepSize;
                    end
                    X(:,:,:,i) = X(:,:,:,i) - stepSize*dPhix;
                    
                    % display and report relative L2-norm between iterates
                    if opt.display
                        if ObjMRI.is3D
                            x = X(:,:,:,i);
                            drawnow, imshow(abs(x(:,:,opt.display)),[])
                            title([opt.message ' Iteration: #' num2str(j)]),pause(0.1)
                        end
                    end
                    if opt.report
                        Report.relativeError(j,i) = Norm(X(:,:,:,i),Xp(:,:,:,i));
                        Xp(:,:,:,i) = X(:,:,:,i);
                    end
                end
            end
        end
        
        function [imageOut,Report] = deconvolution(ObjMRI,Imgs,arg)
            arg.null = 0;
            opt.niter = 20;
            opt.RepImageData = 1; %or 'inf' for reconstruction of all Reps
            opt.MrRegularizationParameter = 0;
            opt.MrPriorType = 'Quadratic'; % 'joint'
            opt.optimizationMethod = 'gradientDescent';% 'LucyRichardson-OSL'; conjugateGradient
            opt.MrPreCompWeights = 1;
            opt.fSigma = 0.3;
            opt.display = 0;
            opt.save = 0;
            opt.message = [];
            opt.report = 0;
            opt.imgMask = 1;
            opt.stepSize = 0.2;
            opt.stepSizeOptimization = 0;
            opt = getFiledsFromUsersOpt(opt,arg);
            
            if opt.display, figure; end
            if isinf(opt.RepImageData) % if inf, then reconstruct all Reps
                nRe = 1:size(Imgs,4);
            else
                nRe = opt.RepImageData; % else, reconstruct the specified Rep, default, Rep = 1
            end
            
            % some reports
            if opt.report
                Report.relativeError = zeros(opt.niter,length(nRe),'single');
                Report.stepSize = zeros(opt.niter,length(nRe),'single');
                Norm = @(x,y) 100*(norm(x(:)-y(:)))/norm(y(:));
            else
                Report = [];
            end
            
            if ObjMRI.SR.is
                fprintf('High-resolution deconvolution\n');
                imageOut = zeros([ObjMRI.SR.imgRefH.ImageSize,length(nRe)],'single');
            else
                fprintf('Native-resolution deconvolution\n');
                imageOut = Imgs;
            end
            %
            if opt.MrRegularizationParameter
                fprintf('%s Regularization...\n',opt.MrPriorType)
                D = @(x,W) opt.MrRegularizationParameter*ObjMRI.Prior.TransGraphGradUndoCrop(W.*ObjMRI.Prior.GraphGradCrop(x));
            else
                D =@(x,W) 0;
            end
            
            Ft = @(x) ObjMRI.UpSampling(x);
            F = @(x) ObjMRI.DownSampling(x);
            
            if opt.stepSizeOptimization
                col = @(x) x(:);
                H = @(x,W)Ft(F(x)) + D(x,W);
            end
            
            
            for i = 1: length(nRe)
                fprintf('Rep: #%d\n',nRe(i))
                Ft1 = Ft(ones(size(Imgs(:,:,:,nRe(i)))));
                Ftx = Ft(Imgs(:,:,:,nRe(i)));
                ImageOut = Ftx;
                if opt.report, ImageOutP = ImageOut; end
                for j = 1:opt.niter
                    
                    % get weights coeffients
                    if strcmpi(opt.MrPriorType , 'Quadratic')
                        W = opt.MrPreCompWeights;
                    elseif strcmpi(opt.MrPriorType , 'Joint')
                        W = opt.MrPreCompWeights .* ObjMRI.Prior.W_GaussianKernel(abs(ImageOut)./max(abs(ImageOut(:))),opt.fSigma);
                        W = W./repmat(sum(W,2),[1,ObjMRI.Prior.nS]);
                    else
                        W = 1;
                    end
                    % optimization method
                    if strcmpi(opt.optimizationMethod, 'LucyRichardson')
                        ImageOut =  (ImageOut.*opt.imgMask./(Ft1+D(ImageOut,W)+1e-8)).*Ft(Imgs(:,:,:,nRe(i))./(F(ImageOut)+1e-8));
                    elseif strcmpi(opt.optimizationMethod, 'gradientDescent')
                        dPhix = Ft(F(ImageOut) - Imgs(:,:,:,nRe(i))) + D(ImageOut,W);
                        if opt.stepSizeOptimization
                            stepSize = dot(dPhix(:), dPhix(:))./ dot(dPhix(:), col(H(dPhix,W)));
                            Report.stepSize(j,i) = stepSize;
                        else
                            stepSize = opt.stepSize;
                        end
                        ImageOut = ImageOut - stepSize*dPhix;
                    elseif strcmpi(opt.optimizationMethod, 'conjugateGradient')
                        [ImageOut] = ObjMRI.PCG(Ftx, H,ImageOut,opt.niter,1,opt);
                        % currently no report for PCG
                        Report =[];
                        break
                    else
                        error('unknown OptimizationMethod')
                    end
                    
                    % display and report relative L2-norm between iterates
                    if opt.display
                        if ObjMRI.is3D
                            drawnow, imshow(abs(ImageOut(:,:,opt.display)),[])
                            title([opt.message ' Iteration: #' num2str(j)]),pause(0.1)
                        end
                    end
                    if opt.report
                        Report.relativeError(j,i) = Norm(ImageOut,ImageOutP);
                        ImageOutP = ImageOut;
                    end
                    
                end
                imageOut(:,:,:,i) = ImageOut;
            end
        end
        
        % ///////////////////////////////////////////////////////////////////////////////
        %                           PVC AND DENOISING METHODS
        % ///////////////////////////////////////////////////////////////////////////////
        
        function ImageOut = nonLocalMeans(ObjMRI,Image,PriorImg,GaussianWeightSigma)
            Image = ObjMRI.Prior.imCrop(Image);
            if ~isempty(PriorImg) % if empty [], self-similarities
                % Calculate non-local similarity coeff from a prior image
                PriorImg = ObjMRI.Prior.imCrop(PriorImg);
            else
                % Calculate non-local similarity coeff from the image
                % itself
                PriorImg = Image;
            end
            w = ObjMRI.Prior.W_GaussianKernel(PriorImg,GaussianWeightSigma);
            w = w./repmat(sum(w,2),[1,ObjMRI.Prior.nS]);
            ImageOut = reshape(sum(Image(ObjMRI.Prior.SearchWindow(:,:)).*w,2),ObjMRI.Prior.CropedImageSize);
            ImageOut = ObjMRI.Prior.UndoImCrop(ImageOut);
        end
        
    end
    
end

