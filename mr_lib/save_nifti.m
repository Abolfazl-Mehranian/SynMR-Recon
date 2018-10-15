
function save_nifti(flname,xVar,savePath,nifti_tmp,niftiInfo,zipit)

nifti_tmp.img = convert2NiftiDataType(xVar, niftiInfo.dataType);

if~exist(savePath,'dir'), mkdir(savePath); end
fln = [savePath flname '.nii'];
save_untouch_nii(nifti_tmp,fln );
if zipit
    gzip(fln)
    delete(fln);
end


function x = convert2NiftiDataType(x, dataType)

   switch double(dataType)
   case   1
      precision = 'ubit1';
   case   2
      precision = 'uint8';
   case   4
      precision = 'int16';
   case   8
      precision = 'int32';
   case  16
      precision = 'single';%'float32';
   case  32
      precision = 'single';%'float32';
   case  64
      precision = 'double';%'float64';
   case 128
      precision = 'uint8';
   case 256 
      precision = 'int8';
   case 512 
      precision = 'uint16';
   case 768 
      precision = 'uint32';
   case 1024
      precision = 'int64';
   case 1280
      precision = 'uint64';
   case 1792
      precision = 'double';%'float64';
   otherwise
      error('This datatype is not supported');
   end
   x = cast(x, precision);
   
   