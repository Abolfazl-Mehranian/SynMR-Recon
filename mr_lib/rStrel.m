function r = rStrel()
[x,y,z] = meshgrid(-1:1:1);
r = (x.^2+y.^2+z.^2)<2;