function v = mapSpaceAToSpaceBspm(rImg,Pr,Ps,interp)

% rImg: image to be resliced
% Pr: reslice image P structure
% Ps: source image P structure
if nargin ==3
    interp = 0; %0: nearest neighbours, 1: Trilinear
end

% assuming rImg is 3D or 4D
nReps = size(rImg,4);
Ps.dim = Ps.dim(1:3);
v = zeros([Ps.dim(1:3) nReps]);
if isreal(rImg)
    for i = 1:nReps
        v(:,:,:,i) = resliceSpm(rImg(:,:,:,i),Pr,Ps,interp);
    end
else
    for i = 1:nReps
        v(:,:,:,i) = resliceSpm(real(rImg(:,:,:,i)),Pr,Ps,interp) +...
            1j*resliceSpm(imag(rImg(:,:,:,i)),Pr,Ps,interp);
%         v(:,:,:,i) = resliceSpm(abs(rImg(:,:,:,i)),Pr,Ps,interp) .*...
%             exp(1j*resliceSpm(angle(rImg(:,:,:,i)),Pr,Ps,interp));        
    end
end

function v = resliceSpm(rImg,Pr,Ps,interp)
% Ps.dim
% Ps.mat




[x1,x2] = ndgrid(1:Ps.dim(1),1:Ps.dim(2));
d = [interp*[1 1 1]' [0 0 0]'];

% C = spm_bsplinc(Pr, d); % similar to imagetoreslice but rotated

v = zeros(Ps.dim);
for x3 = 1:Ps.dim(3)
    [~,y1,y2,y3] = getmask(inv(Ps.mat\Pr.mat),x1,x2,x3,Pr.dim(1:3),[0 0 0]');
    v(:,:,x3)      = spm_bsplins(double(rImg), y1,y2,y3, d);
    % v(~tmp)      = 0;
    
end
v(isnan(v)) = 0;
v(isinf(v)) = 0;

%==========================================================================
%-function [Mask,y1,y2,y3] = getmask(M,x1,x2,x3,dim,wrp)
%==========================================================================
function [Mask,y1,y2,y3] = getmask(M,x1,x2,x3,dim,wrp)
tiny = 5e-2; % From spm_vol_utils.c
y1   = M(1,1)*x1+M(1,2)*x2+(M(1,3)*x3+M(1,4));
y2   = M(2,1)*x1+M(2,2)*x2+(M(2,3)*x3+M(2,4));
y3   = M(3,1)*x1+M(3,2)*x2+(M(3,3)*x3+M(3,4));
Mask = true(size(y1));
if ~wrp(1), Mask = Mask & (y1 >= (1-tiny) & y1 <= (dim(1)+tiny)); end
if ~wrp(2), Mask = Mask & (y2 >= (1-tiny) & y2 <= (dim(2)+tiny)); end
if ~wrp(3), Mask = Mask & (y3 >= (1-tiny) & y3 <= (dim(3)+tiny)); end

