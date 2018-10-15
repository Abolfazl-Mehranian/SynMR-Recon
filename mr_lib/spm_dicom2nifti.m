function  out_dir = spm_dicom2nifti(in_dir)

files = spm_select('FPList',in_dir,'.*');
hdr = spm_dicom_headers(files,1);
out_dir = [in_dir '\nii\'];
mkdir(out_dir);
out = spm_dicom_convert(hdr,'all','flat','nii',out_dir);

out_dir = fileparts(out.files{1});
end