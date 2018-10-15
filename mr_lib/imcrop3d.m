function H = imcrop3d(img,manual)
if nargin==1
    manual = 0;
end
[nc,nr,ns] = size(img);
if manual
    figure,
    pos=[0 0 nc nr];
    imagesc(max(img,[],3)), axis 'square'
    h = imrect(gca,pos);
    wait(h);
    pos1 = round(getPosition(h));
    
    figure
    pos=[0 0 nc ns];
    imagesc(squeeze(max(img,[],1))), axis 'square'
    h = imrect(gca,pos);
    wait(h);
    pos2 = round(getPosition(h));
    H.indx = [pos1(2),pos1(2)+pos1(4),pos1(1),pos1(1)+pos1(3), pos2(1),pos2(1)+pos2(3)];
else
    a = sum(max(img,[],3));
    id = find(a~=0);
    
    H.indx  = [4,nc-4,[max(1,id(1)-4),min(id(end)+4,nr)], 4, ns - 4];
    
end

H.crop = @(x) x(H.indx(1):H.indx(2),H.indx(3):H.indx(4),H.indx(5):H.indx(6));
H.imgSize = size(img);
H.imCropSize = size(H.crop(img));
H.undoCrop = @(x) undoCrop(x,H);

function y = undoCrop(x,H)
if ~all(size(x)==H.imCropSize)
    error('image size doesnot match')
end
y = zeros(H.imgSize);
y(H.indx(1):H.indx(2),H.indx(3):H.indx(4),H.indx(5):H.indx(6)) = x;


