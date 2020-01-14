% ***********************************************
% Copyright (c) 2020 UChicago Argonne, LLC
% See LICENSE file.
% ***********************************************
%
% this code demonstrate how to sum images in the script batch mode

fpath_flist = 'N:\2014-2\wchen_giwaxs_201406p1\flist';
fpath_summed =  'N:\2014-2\wchen_giwaxs_201406p1\summed';
fpath_image = 'N:\2014-2\wchen_giwaxs_201406p1\';

flist = dir(fullfile(fpath_flist,'*.txt'));
flist = {flist().name};

for ii=1:length(flist)
    f = importdata(fullfile(fpath_flist,flist{ii}));
    im = 0;
    for jj=1:size(f,1)
        im = im+imread(fullfile(fpath_image,f{jj}));
    end
    ind = strfind(f{1},'sec');
    t = size(f,1)*str2num(f{1}(ind(end)+3:ind(end)+5));
    fexport = ['summed_',f{ii}]
    info = imfinfo(f{1});
    info.ImageDescription = sprintf('%s%s\n\n',info.ImageDescription,['# Total_exposure: ',num2str(t)] );
   % info.RawImageFiles = f;
    imwrite2tif(im,info,fullfile(fpath_summed,fexport),'int32');
end

