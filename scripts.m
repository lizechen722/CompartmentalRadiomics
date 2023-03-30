%Clear workspace
clear all;
close all;

%Please set 'Current Folder' to the directory with these imaging data
%Import data from NIFTI

data_path = '182864.nii.gz';
mask_path = 'seg182864_mc.nii.gz';


data_key = load_untouch_nii(data_path);
data = data_key.img;
mask_key = load_untouch_nii(mask_path);
mask = mask_key.img;

%pixSize = header.PixelDimensions(1)*header.PixelDimensions(2);

SIZE = size(data);

%data = data *header.raw.scl_slope + header.raw.scl_inter;

%segmask: 3D array with data with just ROI
segmask = zeros (SIZE)-1000;

%1D array for histogram plot
segdata = 0;
n = 1;

%Assigning values to segmask and segdata
for i = 1:SIZE(1)
    for j = 1:SIZE(2)
        for k = 1:SIZE(3)
            if mask (i,j,k) > 0
                if data (i,j,k) > -150 
                    segmask (i,j,k) = data(i,j,k);
                    segdata (n) = data(i,j,k);
                    n = n + 1;
                end
            end
        end
    end
 end

%Compartmentalisation
fat = [-150 -10];
necrosis = [-9 55];
st = [56 110];
iv = [111 250];

fat_cpm = mask;
necrosis_cpm = mask;
st_cpm = mask;
iv_cpm = mask;

dim3 = 0;

for i = 1:SIZE(1)
    for j = 1:SIZE(2)
        for k = 1:SIZE(3)
            fat_cpm(i,j,k) = 0;
            necrosis_cpm(i,j,k) = 0;
            st_cpm(i,j,k) = 0;
            iv_cpm(i,j,k) = 0;
            if (segmask(i,j,k) >= fat(1) &&  segmask(i,j,k) <= fat(2))
                fat_cpm(i,j,k) = 1;
            elseif (segmask(i,j,k) >= necrosis(1) &&  segmask(i,j,k) <= necrosis(2))
                necrosis_cpm(i,j,k) = 1;
            elseif(segmask(i,j,k) >= st(1) && segmask(i,j,k) <= st(2))
                st_cpm(i,j,k) = 1;
            elseif(segmask(i,j,k) >= iv(1) && segmask(i,j,k) <= iv(2))
                iv_cpm(i,j,k) = 1;
            end
        end
    end
end

%Plot histogram
%edges = linspace(-150, 300, 100); % Create 20 bins.
%h = histogram (segdata, 'BinEdges', edges);

%export
mask_key.img = int16(fat_cpm);
save_untouch_nii(mask_key, '182864fat.nii.gz');

mask_key.img = int16(necrosis_cpm);
save_untouch_nii(mask_key, '182864necrosis.nii.gz');

mask_key.img = int16(st_cpm);
save_untouch_nii(mask_key, '182864st.nii.gz');

mask_key.img = int16(iv_cpm);
save_untouch_nii(mask_key, '182864iv.nii.gz');





