function [out, tmp]= getNiftiDataInfo(flname)

nifiObj =nifti(flname);
out.mat = nifiObj.mat;


tmp = load_untouch_nii(flname);
out.img = double(tmp.img);
out.dim = tmp.hdr.dime.dim(2:tmp.hdr.dime.dim(1)+1);
out.pixdim = tmp.hdr.dime.pixdim(2:4);
out.dataType = tmp.hdr.dime.datatype;