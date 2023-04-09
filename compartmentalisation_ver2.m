%Clear workspace
clear all;
close all;

%Please set 'Current Folder' to the directory with these imaging data
%Import data from NIFTI

data_path = '2151.nii.gz';
mask_path = 'seg2151_mc.nii.gz';

data = niftiread(data_path);
mask = niftiread(mask_path);
header = niftiinfo(data_path);
header_mask = niftiinfo(mask_path);

pixSize = header.PixelDimensions(1)*header.PixelDimensions(2);

SIZE = size(data);

data = data *header.raw.scl_slope + header.raw.scl_inter;

%segmask: 3D array with data with just ROI
segmask = zeros (SIZE);

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
edges = linspace(-150, 300, 100); % Create 20 bins.
h = histogram (segdata, 'BinEdges', edges);

%export

info = niftiinfo('seg2151_mc.nii.gz');


niftiwrite(fat_cpm,'test2151fat',info,'Compressed',false, 'Version', NifTI2);
niftiwrite(necrosis_cpm,'new2151necrosis',info,'Compressed',false);
niftiwrite(st_cpm,'new2151st',info,'Compressed',false);
niftiwrite(iv_cpm,'new2151iv',info,'Compressed',false);

