
function gf3d = gauss3DFilter (data, vox3dsz,fwhm,gpu)
%  % 3D aniso/isotropic Gaussian filtering considering voxels size and fwhm

if all(fwhm==0)
    gf3d = data;
    return;
end
if nargin<=3
    gpu = 0;
end


if length(fwhm)==1
    fwhm = fwhm*ones(1,3);
end
sigma=fwhm/sqrt(2^3*log(2));
matsz=ceil(2*fwhm./vox3dsz);
for i=1:size(matsz,2)
    if isequal(mod(matsz(i),2),0), matsz(i)=matsz(i)+1; end
end
padSize = (matsz-1)/2;
bound=padSize.*vox3dsz;
[x,y,z] = meshgrid(-bound(2):vox3dsz(2):bound(2), -bound(1):vox3dsz(1):bound(1), -bound(3):vox3dsz(3):bound(3));
% h = exp(-(x.*x + y.*y + z.*z)/(2*gsigmm*gsigmm));
h = (2*pi)^(-3/2)/sigma(1)/sigma(2)/sigma(3) * exp(-(x.^2/sigma(1)^2/2 + y.^2/sigma(2)^2/2 + z.^2/sigma(3)^2/2));
h = h/sum(h(:));
numDims = length(padSize);
idx = cell(numDims,1);
for k = 1:numDims
    M = size(data,k);
    onesVector = ones(1,padSize(k));
    idx{k} = [onesVector 1:M M*onesVector];
end
b = data(idx{:});
if gpu
    b = gpuArray(b);
    h = gpuArray(h);
    gf3d = gather(convn(b,h, 'valid'));
else
    gf3d = convn(b,h, 'valid');
end

