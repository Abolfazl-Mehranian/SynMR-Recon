
function d = read_shape(siz,format,filename)
if nargin<3
    if nargin==1,
        format = 'float';
        [filename, pathname]= uigetfile('*.*');
        filename =  [ pathname filename];
        
    else
        [filename, pathname]= uigetfile('*.*');
        filename =  [ pathname filename];
    end
end


fid=fopen(filename,'r');

d=fread(fid,prod(siz),format);

d= reshape(d,siz);
fclose(fid);
