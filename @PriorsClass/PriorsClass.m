classdef PriorsClass < handle
    properties (SetAccess = public)
        
        ImageSize
        CropedImageSize
        imCropFactor
        imCropHandle
        sWindowSize % search window size
        lWindowSize % local window size (local neighborhood)
        SearchWindow
        LocalWindow
        chunkSize
        Wd
        nS
        nL
        is3D
        isforTV
    end
    
    methods
        
        function ObjPrior = PriorsClass(varargin)
            ObjPrior.ImageSize = [344,344,1];
            ObjPrior.CropedImageSize = [];
            ObjPrior.sWindowSize = 7;
            ObjPrior.lWindowSize = 3;
            ObjPrior.SearchWindow = [];
            ObjPrior.LocalWindow =[];
            ObjPrior.Wd = [];
            ObjPrior.nS = [];
            ObjPrior.nL = [];
            ObjPrior.chunkSize = 1;
            ObjPrior.is3D = 0;
            ObjPrior.imCropFactor = 4;
            ObjPrior.imCropHandle =[];
            ObjPrior.chunkSize = 5e6;
            ObjPrior.isforTV = 0;
            
            if isempty(varargin{1}.ImageSize)
                error('Image size should be specified');
            end
            
            % get fields from user's input
            ObjPrior = getFiledsFromUsersOpt(ObjPrior,varargin{1});
            
            if length(varargin{1}.ImageSize)==3 && varargin{1}.ImageSize(3)>1
                ObjPrior.is3D = 1;
            end
            InitializePriors(ObjPrior);
        end
    end
    
    methods (Access = private)
        
        function [N,D] = Neighborhood(ObjPrior,w)
            
            % patch-size
            n = ObjPrior.CropedImageSize(1);
            m = ObjPrior.CropedImageSize(2);
            if ObjPrior.is3D
                h = ObjPrior.CropedImageSize(3);
            else
                h = 1;
            end
            
            
            wlen = 2*floor(w/2); % length of neighborhood window
            widx = -wlen/2:wlen/2;
            xidx = widx; yidx = widx;
            
            if h==1
                zidx = 0;
                nN = w*w;
            else
                zidx = widx;
                nN = w*w*w;
            end
            
            % image grid
            [X, Y, Z] = ndgrid(1:n,1:m,1:h);
            Y = single(Y);
            X = single(X);
            Z = single(Z);
            % index and distance
            
            N = zeros(n*m*h, nN,'single');
            D = N;
            l = 1;
            for x = xidx
                Xnew = ObjPrior.setBoundary1(X + x, n);
                for y = yidx
                    Ynew = ObjPrior.setBoundary1(Y + y, m);
                    for z = zidx
                        Znew = ObjPrior.setBoundary1(Z + z, h);
                        N(:,l) = Xnew + (Ynew-1).*n + (Znew-1)*n*m;
                        D(:,l) = sqrt(x^2+y^2+z^2);
                        l = l + 1;
                    end
                end
            end
            
            D = 1./D;
            D(isinf(D))= 0;
            D = D./repmat(sum(D,2),[1,nN]);
            
            % {
            % for local first-order neighborhood 3x3x3
            if ObjPrior.isforTV && ObjPrior.sWindowSize ==3 && ObjPrior.lWindowSize ==1 && size(N,2)>1
                if h>1 %3D
                    nearsetVoxels = [5,11,13,14,15,17,23];
                    ObjPrior.nS = 7;
                else
                    nearsetVoxels = [2,4,5,6,8];
                    ObjPrior.nS = 5;
                end
                N = N(:,nearsetVoxels);
                D = D(:,nearsetVoxels);
            end
            %}
        end
        
        function X = setBoundary1(~,X, n)
            
            idx = X<1;
            X(idx) = 2-X(idx);
            idx = X>n;
            X(idx) = 2*n-X(idx);
            X=X(:);
        end
        
    end
    methods (Access = public)
        
        function InitializePriors(ObjPrior)
            [~,ObjPrior.CropedImageSize] = imCrop(ObjPrior);
            
            if ~rem(ObjPrior.sWindowSize,2)
                error('The size of search window should be odd');
            end
            if ~rem(ObjPrior.lWindowSize,2)
                error('The size of local window should be odd');
            end
            
            if ObjPrior.is3D
                ObjPrior.nS = ObjPrior.sWindowSize^3;
                ObjPrior.nL = ObjPrior.lWindowSize^3;
            else
                ObjPrior.nS = ObjPrior.sWindowSize^2;
                ObjPrior.nL = ObjPrior.lWindowSize^2;
            end
            
            [ObjPrior.SearchWindow, ObjPrior.Wd] = Neighborhood(ObjPrior,ObjPrior.sWindowSize);
            ObjPrior.LocalWindow = Neighborhood(ObjPrior,ObjPrior.lWindowSize);
        end
        
        function [Img,newSize] = imCrop(ObjPrior,Img)
            % 0, [0 0 0]
            % 2,3,...
            % [2, 2, 0]
            if ~isempty(ObjPrior.imCropHandle) %intialized by imcrop3d
                newSize = ObjPrior.imCropHandle.imCropSize;
                if nargin==1
                    Img = [];
                else
                    Img = ObjPrior.imCropHandle.crop(Img);
                end
            else
                if all(ObjPrior.imCropFactor==0)
                    newSize = ObjPrior.ImageSize;
                    if nargin==1
                        Img = [];
                    end
                else
                    if length(ObjPrior.imCropFactor)== 1
                        if ObjPrior.is3D
                            ObjPrior.imCropFactor = ObjPrior.imCropFactor*[1 1 0];
                        else
                            ObjPrior.imCropFactor = ObjPrior.imCropFactor*[1 1];
                        end
                    end
                    
                    J = 0;
                    if ObjPrior.imCropFactor(1)
                        ObjPrior.imCropFactor(1) = max(2.5, ObjPrior.imCropFactor(1));
                        J = floor(ObjPrior.ImageSize(1)/ObjPrior.imCropFactor(1));
                    end
                    
                    I = 0;
                    if ObjPrior.imCropFactor(2)
                        ObjPrior.imCropFactor(2) = max(2.5, ObjPrior.imCropFactor(2));
                        I = floor(ObjPrior.ImageSize(2)/ObjPrior.imCropFactor(2));
                    end
                    
                    if ObjPrior.is3D
                        K = 0;
                        if ObjPrior.imCropFactor(3)
                            ObjPrior.imCropFactor(3) = max(2.5, ObjPrior.imCropFactor(3));
                            K = floor(ObjPrior.ImageSize(3)/ObjPrior.imCropFactor(3));
                        end
                        newSize = [length((J:(ObjPrior.ImageSize(1)-J-1))+1),length((I:(ObjPrior.ImageSize(2)-I-1))+1),length((K:(ObjPrior.ImageSize(3)-K-1))+1)];
                    else
                        newSize = [length((J:(ObjPrior.ImageSize(1)-J-1))+1),length((I:(ObjPrior.ImageSize(2)-I-1))+1)];
                        if length(ObjPrior.ImageSize)==3
                            newSize = [ newSize ,1];
                        end
                    end
                    if nargin==1
                        Img = [];
                    else
                        if ObjPrior.is3D
                            Img = Img((J:(ObjPrior.ImageSize(1)-J-1))+1,(I:(ObjPrior.ImageSize(2)-I-1))+1,(K:(ObjPrior.ImageSize(3)-K-1))+1);
                        else
                            Img = Img((J:(ObjPrior.ImageSize(1)-J-1))+1,(I:(ObjPrior.ImageSize(2)-I-1))+1);
                        end
                    end
                end
            end
        end
        
        function ImgNew = UndoImCrop(ObjPrior,Img)
            
            if ~isempty(ObjPrior.imCropHandle) %intialized by imcrop3d
                ImgNew = ObjPrior.imCropHandle.undoCrop(Img);
            else
                if all(ObjPrior.imCropFactor==0)
                    ImgNew = Img;
                    return
                end
                
                
                ImgNew = zeros(ObjPrior.ImageSize,'single');
                
                S = (ObjPrior.ImageSize - ObjPrior.CropedImageSize)/2;
                J = S(1); I = S(2);
                if ObjPrior.is3D
                    K = S(3);
                    ImgNew((J:(ObjPrior.ImageSize(1)-S(1)-1))+1,(I:(ObjPrior.ImageSize(2)-I-1))+1,(K:(ObjPrior.ImageSize(3)-K-1))+1) = Img;
                else
                    ImgNew((J:(ObjPrior.ImageSize(1)-S(1)-1))+1,(I:(ObjPrior.ImageSize(2)-I-1))+1) = Img;
                end
            end
        end
        
        function ObjPrior = RevisePrior(ObjPrior,opt)
            % to revise the properties of the object
            
            ObjPrior = getFiledsFromUsersOpt(ObjPrior,opt);
            if isfield(opt,'sWindowSize') || isfield(opt,'lWindowSize') ...
                    || isfield(opt,'imCropFactor') || isfield(opt,'ImageSize')
                InitializePriors(ObjPrior);
            end
        end
        
        function imgGrad = GraphGrad(ObjPrior,Img)
            imgGrad = (Img(ObjPrior.SearchWindow)-repmat(Img(:),[1,ObjPrior.nS]));
        end
        
        function imgGrad = GraphGradCrop(ObjPrior,Img)
            Img = imCrop(ObjPrior,(Img));
            imgGrad = (Img(ObjPrior.SearchWindow)-repmat(Img(:),[1,ObjPrior.nS]));
            imgGrad(isnan(imgGrad)) = 0;
        end
        
        function imgDiv = GraphDivCrop(ObjPrior,Img)
            Img = imCrop(ObjPrior,single(Img));
            imgDiv = (Img(ObjPrior.SearchWindow)+repmat(Img(:),[1,ObjPrior.nS]));
        end
        function dP = TransGraphGradUndoCrop(ObjPrior,imgGrad)
            dP = -2* sum(ObjPrior.Wd.*imgGrad,2);
            dP = reshape(dP,ObjPrior.CropedImageSize);
            dP = UndoImCrop(ObjPrior,dP);
            dP(isnan(dP)) = 0;
        end
        
        function imgGradW = TV_weights(ObjPrior,imgGrad,beta)
            Norm = repmat(sqrt(sum(abs(imgGrad).^2,2)+ beta.^2),[1,ObjPrior.nS]);
            imgGradW = imgGrad./Norm/2;
        end
        
        % neighborhood weights
        function Wg = W_GaussianKernel(ObjPrior,Img,KernelSigma)
            
            imgSize = size(Img);
            if ~all(ObjPrior.CropedImageSize(1:2) == imgSize(1:2))
                Img = ObjPrior.imCrop(Img);
            end
            
            nVoxels = prod(ObjPrior.CropedImageSize);
            Wg = zeros(nVoxels,ObjPrior.nS,'single');
            
            normalize = 1;
            for i = 1:ObjPrior.chunkSize:nVoxels
                voxels = i: min(i+ObjPrior.chunkSize-1,nVoxels);
                
                imgPatch = Img(ObjPrior.LocalWindow(ObjPrior.SearchWindow(voxels,:),:));
                imgLocalWindow = Img(ObjPrior.LocalWindow(voxels,:));
                
                
                imgLocalWindow = repmat(imgLocalWindow,[ObjPrior.nS,1]);
                D = reshape(sum( abs(imgPatch - imgLocalWindow).^2, 2 ),length(voxels),ObjPrior.nS);
                clear imgPatch imgLocalWindow
                
                if normalize
                    D_norm = repmat(max(sqrt(sum(D.^2,2)),eps),[1,ObjPrior.nS]);
                else
                    D_norm = 1;
                end
                Wg(voxels,:) = exp( -D/KernelSigma^2 ./D_norm);%
                clear D D_norm
            end
            % check for errors
            sumWg = sum(Wg,2);
            i = sumWg==0 | isinf(sumWg) | isnan(sumWg);
            Wg(i,:) = 1./ObjPrior.nS;
            
        end
        
        function Wb = W_Bowsher(ObjPrior,Img,B)
            % B: user defined number of the neighboring voxels that have the
            % highest similarity on the anatomical image based on thier
            % absoulte intensity differences
            imgSize = size(Img);
            if ~all(ObjPrior.CropedImageSize(1:2) == imgSize(1:2))
                Img = ObjPrior.imCrop(Img);
            end
            if B > ObjPrior.nS
                error ('B can be in maximum %d\n',ObjPrior.nS)
            end
            abs_imgGrad = abs(ObjPrior.GraphGrad(Img));
            Wb = 0*abs_imgGrad;
            
            for i = 1:size(abs_imgGrad,1)
                [~,idx] = sort(abs_imgGrad(i,:));
                
                Wb(i,idx(1:B)) = 1;
            end
            % Wb = Wb./repmat(sum(Wb,2),[1,ObjPrior.nS]);
        end
        function Wje = W_JointEntropy(ObjPrior,Img,sigma)
            imgSize = size(Img);
            if ~all(ObjPrior.CropedImageSize(1:2) == imgSize(1:2))
                Img = ObjPrior.imCrop(Img);
            end
            Wje = ObjPrior.normPDF(ObjPrior.GraphGrad(Img),0,sigma);
        end
        
        function x = normPDF(~,x,y,sigma)
            x = 1./sqrt(2*pi*sigma^2).*exp(-0.5*(x-y).^2./sigma.^2);
        end
        
        function display(ObjPrior)
            disp(ObjPrior)
            methods(ObjPrior)
        end
    end
end
