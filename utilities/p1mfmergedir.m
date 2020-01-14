function [ numOK ] = p1mfmergedir( samplename, upscan, downscan, basepath )
% ***********************************************
% Copyright (c) 2020 UChicago Argonne, LLC
% See LICENSE file.
% ***********************************************
%
%p1mfmergedir fills gaps for corresponding images from 2 ccdscans
%   numgapfilledfiles = p1mfmergedir(samplename, upscan, downscan,
%   basepath)

% Note: ffpath, the location of flatfield data, must be edited for your system!

% p1mfmerge_ff_extrapixel_batch_fn    Joe Strzalka 2012.02.08
% recasts p1mfmerge_ff_extrapixel_batch as a function
% use this script to fill gaps in batch mode for data collected in up and
% down positions using the ccdscan command.   Corresponding pairs from 2
% scan directories are combined.
% 
% based on p1mfmerge_ff_extrapixel:
% Compensates for gaps in Pilatus 1M-F by combining data from two files
% taken with a vertical offset of n pixels.  Flatfield correction is 
% applied to each image before merging.  To compensate for detector non-
% uniformity at the borders of modules, replacement of gap pixels is 
% extended by extrapixels number of rows above and below the gap.  
% Suggested best practice: collect images with offset of nshift =
% 17(gapsize) + extrapixel+1 = 23 pixel shift or 3.9560 mm shift.
%
% requires: imwrite2tif.m by Zhang Jiang, part of the gixsgui distribution

% 0) define useful quantities
% number of pixels shifted between images
nshift =23;
% number of extra pixels to add above and below each gap (helps compensate
% for non-uniformity of detector at edge of module
extrapixel=5;
% show plots of data? 1 = yes; 2 = no;
graphicoutput = 0;
% color scale limits for plots
clims = [1 10000];
% EDIT ffpath AND ffname FOR YOUR SYSTEM
% flatfield correction file and path
% next line for the beamline linux system
ffpath = '/home/8-id/MATLAB/GIXSGUI/pilatus_1mf/flatfield';
% ffname = 'FFcorr_p1m0108_E11p9_T5p0_vrf_m0p15.tif';
% ffname = 'FFcorr_p1m108_E11p9_T5p9_vrf_m0p2.tif';
% ffname = 'FFcorr_p1m0108_E11p9_T9p9_vrf_m0p3.tif';
ffname = 'FF_p1m0108_E11p9_T9p9_vrf_m0p3_20111130.tif';

% 1) specify the directories to be combined
%  CHANGE THE INPUTS HERE TO READ IN YOUR DATA
%   define the base path as directory containing all the data
% basepath = 'C:\Documents and Settings\strzalka\My Documents\MATLAB\develop\p1mfmerge\data\';
% define the samplename as the subdirectory for a particular sample
%samplename = 'PBIT1'
%samplename = 'BDTC60_1to4'
%upscan = 7;
%downscan = 9;

upinpath  = fullfile(basepath, samplename, ['S' num2str(upscan)]);
downinpath = fullfile(basepath, samplename, ['S' num2str(downscan)]);

upflist = dir(fullfile(upinpath, [samplename '*'] ));
downflist = dir(fullfile(downinpath, [samplename '*']));

upnum = size(upflist,1)
downnum = size(downflist,1)

numOK = 0;

if (upnum ~= downnum)
    sprintf('Error: These two scans have different lengths.')
else
    for jj = (0:upnum-1)  


    upname  =  [samplename '_s' num2str(upscan) '_' num2str(jj) '.tif'] ;
    downname = [samplename '_s' num2str(downscan) '_' num2str(jj) '.tif'] ;

    inname2 = upname;
    inname1 = downname;

    %   read in the data
    inframe2 = double(imread(fullfile(upinpath,upname)));
    inframe1 = double(imread(fullfile(downinpath,downname)));
    inheader1 = imfinfo(fullfile(downinpath,downname));
    
    % 1b) do the flatfield correction
    ffcorr = double(imread(fullfile(ffpath, ffname)));

    raw1 = inframe1; 
    raw2 = inframe2;

    inframe1 = raw1.*ffcorr;
    inframe2 = raw2.*ffcorr;

    % 2) determine which rows correspond to gaps
    copy1 = inframe1;
    %copy1(:,100) = 100;

    gaprows = [196 212 408 424 620 636 832 848]; % 4 gaps of 17 pixels


    % 3) subsitute rows from 2nd frame into 1st frame
    for gap =1:4
      for row = gaprows(2*gap-1)-extrapixel:gaprows(2*gap)+extrapixel
            copy1(row,:)=inframe2(row+nshift,:);
         % copy1(row,:)= 77;
      end
    end

    % 4) plot the data, if desired
    if graphicoutput
     scrsz = get(0,'ScreenSize');
     figure('Position',[100 scrsz(2)  scrsz(4)/3 scrsz(4)])
     daspect([1 1 1]);
     pbaspect([1 1 1]);
     subplot(3, 1, 1); 
     imagesc(raw1, clims)
     % Make sure underscore character appears correctly
     title(regexprep( inname1,'_', '\\_') )
     subplot(3,1,2); imagesc(raw2, clims)
     title(regexprep( inname2,'_', '\\_') )
     subplot(3,1,3); imagesc(copy1, clims)
     title([regexprep( inname1,'_', '\\_')  ' corrected'])
    
     figure
     daspect([1 1 1]);
     imagesc(copy1, clims)
     title([regexprep( inname1,'_', '\\_')  ' corrected;' ' extrapixel: ' num2str(extrapixel)])
  
    end % if graphicoutput

    % 5) write the result to a file

    % save current directory
    savepath = pwd;
    cd(basepath)
   
    if (exist('nogapdata')~=7)
     mkdir('nogapdata')
    end

    cd nogapdata
    outname = [inname1(1:end-4) '_ff_ep' num2str(extrapixel) ]
 
    imwrite2tif(copy1,inheader1, [ outname '.tif'], 'double' )
    % restore current path to original
    cd(savepath)
    numOK = numOK+1;
    end % jj loop
end % if (upnum ~= downnum)
end % function