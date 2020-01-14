% ***********************************************
% Copyright (c) 2020 UChicago Argonne, LLC
% See LICENSE file.
% ***********************************************
%
% GIXSDATA Class to visualize and process x-ray scattering data in 
%   transmission or reflection geometry.
%
%   OJB = GIXSDATA('IMAGE FILE') loads image from file into gixsdata and 
%   return handle OBJ.
%
%   OBJ = GIXSDATA initializes gixsdata with default parameters.
%
%   If handle OBJ exists, image DATA can be directly passed by OBJ.RawData
%   = DATA. 
%
%   CLEAR OBJ only removes the handle. Use DELETE(OBJ) to delete the
%   gixsdata object before clear handle OBJ.
%
%   See also COPYHOBJ, DELETE, IMAGESC, LINECUT, QMAPS, CQMAPS, GET_CMASK,
%   MASK_INITIALIZE, POST_IMAGESC, QMAPS_INITIALIZE, EFFICIENCY_CORRECTION,
%   POLARIZATION_CORRECTION, SOLIDANGLE_CORRECTION, RESHAPE_IMAGE,
%   FINDPEAK, ROICOUNT.
%
%   Type 'doc gixsdata' for reference page in Help browser

%   Zhang Jiang @8ID/APS/ANL
%   $Revision$  $Date: 2011/01/18 $
%   $Revision$  $Date: 2011/03/01 $
%       (1) Add property MaskedDataStats to store statistics of MaskedData.
%       (2) Add property FlatField for flat field correction. (ZJ and JS)
%       (3) Add property EfficiencyCorrection and method
%       efficiency_correction to correct images due to x-ray path 
%       absorption and detector sensor attenuation. (ZJ and JS)
%       (4) Add property PlotImageFlag to select data (MaskedData or 
%       SolidAngleCorrectedData) to be plotted. (ZJ)
%       (5) Add method reshape_image to convert image axis to Cartesian
%       coordinate. (ZJ)
%   $Revision$  $Date: 2011/09/16 $
%       (1) Add method findpeak to fit for 2D diffraction peak positions. 
%   $Revision$  $Date: 2011/10/07 $
%       (1) Add property PhiMode and set.PhiMode and get.PhiMap.
%   $Revision$  $Date: 2011/10/09 $
%       (1) Add property PolarizationMode and set.Polarization. Add
%       polarization correction to SolidAngleCorrectedData.
%   $Revision$  $Data: 2012/02/23 $
%       (1) Add PhiMode = 3 and 4 options
%   $Revision$  $Data: 2012/05/14 $
%       (1) To speed up reshape_image, add line space boundary check,
%       disable the collection of q map values of reshaped ROI, and flip x,
%       y variable if number of x is smaller than y.
%   $Revision$  $Data: 2012/05/14 $
%       (1) Add option 3 for PlotImageFlag to display RawData
%   $Revision$  $Data: 2012/07/27 $
%       (1) Calculate q maps using rotation matrices. Rotation method is
%       about 12% faster than old method.
%       (2) Change sign of qx to follow conventional definition of qx, ie.,
%       negative values above specular point.
%   $Revision$  $Data: 2012/08/06 $
%       (1) Add method mappixel to find pixel indices of given map values.
%   $Revision$  $Data: 2012/08/15 $
%       (1) Assign +/- sign to QrMap.
%   $Revision$  $Data: 2012/09/14 $
%       (1) Tick lable outside of plot.
%   $Revision$  $Data: 2013/01/25 $
%       (1) Include property ChiMap for polar angles.
%       (2) Include method CQMAPS to get the averaged map values within
%       defined region by mask and constraint.
%       (3) Add gixsdatacursor.
%   $Revision$  $Data: 2013/02/20 $
%       (1) Add an option output field - mapdata - to LINECUT, to calculate
%       corresponding averaged map values for each spot of the linecut.
%   $Revision$  $Data: 2014/06/11 $
%       (1) Add support for .edf image fomat (requested by CHZ)
%   $Revision$  $Data: 2014/12/1 $
%       (1) Add property for different types for Lorentz correction and
%       method lorentzfactor to calculate the Lorentz factor
%       (2) Add property for custom correction
%   $Revision$  $Data: 2015/03/13 $
%       (1) Fix the compatibility issue with 2014b when calling
%       gixsdatacursor.
%   $Revision$  $Data: 2015/12/30 $
%       (1) Fix xticklabel problem of off-centered lableling of 0 for the
%       x-axis.
%   $Revision$  $Data: 2016/04/19 $
%       (1) Fix the bug of AND/OR on linecut constraints. See get_cmask.
%   $Revision$  $Data: 2016/05/08 $
%       (1) Use pmedf_read.m by Petr Mikulik to read EDF files 
%   $Revision$  $Data: 2016/08/01 $
%       (1) Add NOFPIXLES output for linecut function.
%   $Revision$  $Data: 2017/09/05 $
%       (1) Also export the map name list as a field in the 4th output
%       argument in linecut method. 
%   $Revision$  $Data: 2017/11/09 $
%       (1) Allow inpainting to fill pixels with nan values after image
%       reshaping. If the reshaping is done for qr-qz, the forbidden region
%       is back filled with nans after inpainting.
%       (2) Limit the map span range to the extreme map values of the image
%       during reshaping.
%   $Revision$  $Data: 2019/08/23 $
%       (1) Enable Lambda and Eiger camera

classdef gixsdata < handle
    properties  % public
        % Type of camera. String must match 'Pilatus', 'MAR165' or 'Other'.
        % Default: 'Pilatus'
        Camera      = 'Pilatus';
        % Pixel sizes (1x2) are automatically set for 'Pilatus' and 'MAR165'.
        PixelSize = [0.172,0.172];
        % Beam zero (1x2) position in pixels.
        Beam0      = [NaN,NaN];
        % Sample to detector distance (mm).
        SDD         = NaN;
        % X-ray energy (keV). Default: 7.35
        XEnergy     = 10.86;
        % Scattering geometry. 1/[2] for transmission/[reflection].
        Geometry    = 2;
        % Phi mode. [1]/2/3/4 for (-180,180]/[0,360)/[-270,90)/[-90,270).
        PhiMode     = 1;        
        % Polarization mode. 1/[2]/3/4 for
        % none/[horizontal]/vertical/unpolarized
        PolarizationMode = 2; 
        % Horizontal polarization fraction. Value between [0,1]. Default: 1
        HorizontalPolarizationFraction = 1;
        % Incident angle for reflection geometry. It is not used for 
        % transmission geometry, and can set to NaN. Unit: degree.
        IncidentAngle = NaN;
        % Specular position (1x2) for reflection geometry. Can be NaN for
        % transmission geometry.
        Specular    = [NaN NaN];
        % Color limits type for plot. [1]/2 for [manual]/automatic.
        PlotCLimsType = 1;
        % Color limits (1x2) for plot. Values are set automatically based
        % current image intensity when PlotCLimsType=2.
        PlotCLims   = [1,2000];
        % Image display scale. 1/[2] for linear/[log].
        PlotScale   = 2;
        % Axis label style. [1]/2/3 for [pixel]/q/angle.
        PlotAxisLabel   = 1;
        % Plot flag. [1]/2/3 for MaskedData/SolidAngleCorrectedData/RawData
        PlotImageFlag = 1;        
        % Full image file name.
        ImFile      = '';
        % (mxn) Raw data.
        RawData     = single([]);
        % Type of Lorentz factor for correction. [1]/2/3/4 for no 
        % correction, in-plane (sample surface) randomly oriented 3D 
        % structures, in-plane randomly oriented 2D structures, and 
        % of 2D crystals, and and completely randomly oriented objects 
        % (such as powders; usually for SAXS and WAXS). Types 2 and 3 are 
        % for reflection geometry only.
        LorentzFactorType = 1; 
        % (mxn) Flat field to be multiplied by the image data. Default:
        % ones(m,n). The correction is automatically applied to both 
        % MaskedData and SolidAngleCorrectedData.        
        FlatField = [];
        % (mx) Custom correction matrix to be mulitiplied by the image
        % data; Default: ones(m,n). The correction is automatically applied
        % to SolidAngleCorrectedData.
        CustomCorrection = [];        
        % (2x3) Efficiency correction for path difference arising from path
        % medium (e.g. air) absorption and detector sensor attenuation.
        % These correction matrix is to be multiplied by the image data.
        % To apply the corrected data, obj.solidangle_correction(obj) or
        % obj.qmaps(obj) must be called be called after parameters in
        % EfficiencyCorrection is updated. Definition:
        %   {1,1}: path medium chemical formula. 'N78O21Ar1' for air;
        %   {1,2}: path medium mass density (g/cm^3);
        %   {1,2}: path medium length (mm);
        %   {2,1}: detector sensor chemical formula. 'Si' for Pilatus;
        %   {2,2}: detector sensor mass density (g/cm^3);
        %   {2,3}: detector sensor depth (mm);
        % Set chemical formula to 'N/A' if the corrrected is unwanted. For
        % path medium, alternatively {1,2} or {1,3} can be set to 0 to
        % disable the absorption correction.
        EfficiencyCorrection = {...
            'N78O21Ar',   1.1839e-3, 0;   ...
            'Si',         2.33,      0.32};
        % (mxn) Mask data (mxn).
        Mask        = [];
        % (mxn) Q map with respect to direct beam.
        QMap        = [];
        % (mxn) Qx map (only for reflection geometry).
        QxMap      = [];
        % (mxn) Qz map (only for reflection geometry).
        QzMap      = [];
        % (mxn) Qy map (only for reflection geometry).
        QyMap      = [];
        % (mxn) Qr map (only for reflection geometry).
        QrMap      = [];
        % (mxn) Phi (azimuthal angle) map with respect to direct beam.
        PhiMap     = [];
        % (mxn) 2Theta (out of plane angle) map (only for reflection geometry).
        TwoThetaMap   = [];
        % (mxn) Alpha_f (exit angle) (only for reflection geometry).
        AlphafMap  = [];
        % (mxn) Chi (polar angle) (only for reflection geometry).
        ChiMap = [];        
        % (mx1) Qx list for y axis label.
        Qx1DList   = [];
        % (1xn) Qy list for x axis label.
        Qy1DList   = [];
        % (mx1) Qz list for y axis label.
        Qz1DList   = [];
        % (1xn) 2Theta list for x axis label.
        TwoTheta1DList = [];
        % (mx1) Alpha_f list for y axis label.
        Alphaf1DList = [];
        % (mxn) Solid angle corrected data for linecut. Valid only when
        % QMap exists.
        SolidAngleCorrectedData = [];
        ImFileName      = '';   % Image file name.
    end
    properties (SetAccess = private)
        FigHandle       = [];   % Store the figure handle of the plotted image.
        ImFilePath      = '';   % Image file path.
%         ImFileName      = '';   % Image file name.
        ImFileExt       = '';   % Image file extension.
        ImFileInfo      = [];   % Image file information from image header.
    end
    properties (Dependent = true, SetAccess = private)
        % (mxn) Masked data after applying Mask to RawData. Values are
        % automatically updated when either Mask or RawData is changed.
        MaskedData  = [];
        % Statistics of MaskedData.
        MaskedDataStats = struct;
        % Dimension of image (1x2). Values are automatically detected.
        % Elements with 0 value in Mask are filled with NaN.
        ImDim;
        % X-ray wavelength (A). Value is automatically updated when energy
        % is changed.
        XWavelength;
    end
    % ===============
    methods % dynamic methods
        % --- interface constructor
        function obj = gixsdata(varargin)
            if nargin == 1
                if ischar(varargin{1}) && exist(varargin{1},'file')==2
                    obj.ImFile = varargin{1};
                elseif isnumeric(varargin{1})
                    obj.RawData = single(varargin{1});
                    obj.mask_initialize(obj);
                end
            end
        end
        % --- set Camera type
        function obj = set.Camera(obj,camera)
            if ~(strcmpi(camera,'Pilatus') || strcmpi(camera,'MAR165') || strcmpi(camera,'Lambda (Si)') ...
                    || strcmpi(camera,'Eiger (Si)') || strcmpi(camera,'Other'))
                error('Camera must be Pilatus, MAR165 or Other.');
            end
            obj.Camera = camera;
            switch obj.Camera
                case 'Pilatus'
                    obj.PixelSize = [0.172,0.172];
                    obj.EfficiencyCorrection{2,1} = 'Si';
                    obj.EfficiencyCorrection{2,2} = 2.33;
                    obj.EfficiencyCorrection{2,3} = 0.32;
                case 'MAR165'
                    obj.PixelSize = [0.079138,0.079138];
                    obj.EfficiencyCorrection{2,1} = 'N/A';
                    obj.EfficiencyCorrection{2,2} = NaN;
                    obj.EfficiencyCorrection{2,3} = NaN;
                case 'Lambda (Si)'
                    obj.PixelSize = [0.055,0.055];
                    obj.EfficiencyCorrection{2,1} = 'Si';
                    obj.EfficiencyCorrection{2,2} = 2.33;
                    obj.EfficiencyCorrection{2,3} = 0.3;                    
                case 'Eiger (Si)'
                    obj.PixelSize = [0.075,0.075];
                    obj.EfficiencyCorrection{2,1} = 'Si';
                    obj.EfficiencyCorrection{2,2} = 2.33;
                    obj.EfficiencyCorrection{2,3} = 0.45;  
            end
        end
        % --- set pixel size
        function obj = set.PixelSize(obj,pixelsize)
            if isempty(pixelsize) || ~isnumeric(pixelsize) || length(pixelsize)~=2 ...
                    || pixelsize(1)<0 || pixelsize(2)<0, return; end
            obj.PixelSize = pixelsize;
        end
        % --- set beam0
        function obj = set.Beam0(obj,beam0)
            if isempty(beam0) || ~isnumeric(beam0) || length(beam0)~=2,return; end
            obj.Beam0 = beam0;
        end
        % --- set SDD
        function obj = set.SDD(obj,sdd)
            if isempty(sdd) || ~isnumeric(sdd) || length(sdd)~=1 || sdd<=0, return; end
            obj.SDD = sdd;
        end
        % --- set energy
        function obj = set.XEnergy(obj,energy)
            if isempty(energy) || ~isnumeric(energy) || length(energy)~=1 || energy<=0, return; end
            obj.XEnergy = energy;
        end
        % --- set geometry
        function obj = set.Geometry(obj,geometry)
            if geometry == 1
                obj.PlotAxisLabel = 1; 
                if obj.LorentzFactorType == 2 || obj.LorentzFactorType == 3
                    obj.LorentzFactorType = 1;    
                end
            end
            obj.Geometry = geometry;
        end
        % --- set phi mode
        function obj = set.PhiMode(obj,phimode)
            if phimode ~= 1 && phimode ~=2 && ...
                    phimode ~=3 && phimode ~= 4, return; end
            obj.PhiMode = phimode;            
        end
        % --- set polarization mode
        function obj = set.PolarizationMode(obj,polarization)
            if polarization == 1 || polarization ==2 ...
            || polarization ==3 || polarization == 4
                obj.PolarizationMode = polarization;              
            end
        end        
        % --- set degree of polarization
        function obj = set.HorizontalPolarizationFraction(obj,hpf)
            if hpf<=1 && hpf>=0,
                obj.HorizontalPolarizationFraction = hpf;
            end
        end
        % --- set incident angle
        function obj = set.IncidentAngle(obj,angle)
            if isempty(angle) || ~isnumeric(angle) || length(angle)~=1, return; end
            obj.IncidentAngle = angle;
        end
        % --- specular position
        function obj=set.Specular(obj,specular)
            if isempty(specular) || ~isnumeric(specular) || length(specular)~=2, return;  end
            obj.Specular = specular;
        end
        % --- set plot axis label
        function obj = set.PlotAxisLabel(obj,plotaxislabel)
            if (obj.Geometry == 1 && plotaxislabel ~= 1) || isempty(obj.QMap)
                obj.PlotAxisLabel = 1;
                return;
            else
                obj.PlotAxisLabel = plotaxislabel;
            end
        end
        % --- set plot clims type
        function obj = set.PlotCLimsType(obj,climstype)
            if climstype~=1 && climstype~=2, return; end
            if climstype == 2
                data = obj.MaskedData(~isnan(obj.MaskedData));
                clims = [max(1,min(data)),  mean(data)*200];
                if ~isnan(sum(clims)) && clims(2)>=clims(1)
                    obj.PlotCLims = clims;
                end
            end
            obj.PlotCLimsType = climstype;
        end
        % --- set plot clims
        function obj = set.PlotCLims(obj,clims)
            if clims(2)<clims(1), return; end
            obj.PlotCLims = clims;
        end
        % --- set plot image flag
        function obj = set.PlotImageFlag(obj,flag)
            if isempty(obj.SolidAngleCorrectedData)
                obj.PlotImageFlag = 1;
            end
            if ~isnumeric(flag) || (flag~=1 && flag~=2 && flag~=3), return; end
            obj.PlotImageFlag = flag;
        end
        % --- set file name, automatically load image and initialize mask
        function obj = set.ImFile(obj,imfile)
            if exist(imfile,'file')~=2, return; end
            if ~isempty(which(imfile)), imfile = which(imfile); end
            [pathstr, name, ext] = fileparts(imfile);
            obj.ImFile = imfile;
            obj.ImFilePath = pathstr;
            obj.ImFileName = name;
            obj.ImFileExt = ext;
            [obj.RawData,obj.ImFileInfo] = loadimage(obj);
            if isequal(size(obj.RawData),size(obj.Mask))     % pass default mask;
                obj.Mask = (obj.RawData>=0) & obj.Mask;
            else
                obj.mask_initialize(obj);
            end
        end
        % --- set raw data
        function obj = set.RawData(obj,data)
            if isempty(data)
                obj.RawData = [];
                obj.SolidAngleCorrectedData = [];
                return;
            else
                obj.RawData = single(data);
            end
            % update FlatField
            if ~isequal(size(obj.FlatField),size(data))
                obj.FlatField = single(ones(size(data)));
            end
            % update CustomCorrection
            if ~isequal(size(obj.CustomCorrection),size(data))
                obj.CustomCorrection = single(ones(size(data)));
            end         
            % update q maps and solid angle corrected data
            if ~isequal(size(obj.Mask),size(data))
                obj.mask_initialize(obj);
                obj.qmaps_initialize(obj);
            end
            if isempty(obj.QMap)
                qmaps(obj);
            else
                obj.solidangle_correction(obj);
            end
        end
        % --- set mask
        function obj=set.Mask(obj,mask)
            obj.Mask = mask;
            if ~isempty(obj.RawData) && ~isequal(size(obj.Mask),size(obj.RawData))
                obj.mask_initialize(obj);
            end
        end
        % --- set efficiency correction
        function obj = set.EfficiencyCorrection(obj,effc)
            if ~iscell(effc) || ~isequal(size(effc),size(obj.EfficiencyCorrection)) || ...
                    ~ischar(effc{1,1}) || ~isnumeric(effc{1,2}) || length(effc{1,2})~=1 || ...
                    effc{1,2}<0 || ~isnumeric(effc{1,3}) || length(effc{1,3})~=1 || effc{1,3}<0 || ...
                    ~ischar(effc{2,1}) || ~isnumeric(effc{2,2}) || length(effc{2,2})~=1 || ...
                    effc{2,2}<0 || ~isnumeric(effc{2,3}) || length(effc{2,3})~=1 || ...
                    effc{2,3}<0
                return;
            end
            obj.EfficiencyCorrection = effc;
        end
        % --- set flat field
        function obj=set.FlatField(obj,ff)
            if ~isempty(obj.RawData) && ~isequal(size(ff),size(obj.RawData))
                obj.FlatField = single(ones(size(obj.RawData)));
            else
                obj.FlatField = single(ff);
            end
            obj.solidangle_correction(obj);
        end
        % --- set custom correction
        function obj=set.CustomCorrection(obj,cc)
            if ~isempty(obj.RawData) && ~isequal(size(cc),size(obj.RawData))
                obj.CustomCorrection = single(ones(size(obj.RawData)));
            else
                obj.CustomCorrection = single(cc);
            end
            obj.solidangle_correction(obj);
        end
        
        % --- set Lorentz factor
        function obj=set.LorentzFactorType(obj,lorentz)
            if lorentz == 1 || lorentz == 4 || ...
                    ((lorentz == 2 || lorentz == 3) && obj.Geometry == 2)
                obj.LorentzFactorType = lorentz;
            end
        end
        % --- get masked data
        function maskedData = get.MaskedData(obj)
            maskedData = obj.RawData;
            if isempty(maskedData), return; end
            maskedData(~obj.Mask) = NaN;
            maskedData = maskedData.*obj.FlatField;
        end
        % --- get statistics of masked data
        function stats = get.MaskedDataStats(obj)
            if isempty(obj.MaskedData), stats = struct; return; end
            % --- find min
            [v_min,~] = min(obj.MaskedData(:));
            ind_min = find(v_min==obj.MaskedData(:));
            [y_min,x_min] = ind2sub(size(obj.MaskedData),ind_min);
            p_min = [x_min,y_min];
            % --- find max
            [v_max,~] = max(obj.MaskedData(:));
            ind_max = find(v_max==obj.MaskedData(:));
            [y_max,x_max] = ind2sub(size(obj.MaskedData),ind_max);
            p_max = [x_max,y_max];
            % --- get total and averaged counts
            maskedData1D = obj.MaskedData(obj.Mask);
            v_tot = sum(maskedData1D);
            v_mean = mean(maskedData1D);
            v_std = std(maskedData1D);
            v_median = median(maskedData1D);
            n_nnz = nnz(maskedData1D);
            % --- assign values to stats
            stats.MinValue = v_min;
            stats.MinPixel = p_min;
            stats.MaxValue = v_max;
            stats.MaxPixel = p_max;
            stats.TotalValue    = v_tot;
            stats.MeanValue     = v_mean;
            stats.STDValue      = v_std;
            stats.MedianValue   = v_median;
            stats.NumOfNonZeros = n_nnz;
        end
        % --- get image size
        function imDim = get.ImDim(obj)
            imDim = fliplr(size(obj.RawData));
        end
        % --- get wavelength
        function xWavelength = get.XWavelength(obj)
            xWavelength =  12.3984171668278/obj.XEnergy;
        end
        % --- get phimap
        function phimap = get.PhiMap(obj)
            phimap = obj.PhiMap;
            % convert to (-180,180] scale
            phimap = angle(exp(1i*phimap*pi/180))*180/pi;
            switch obj.PhiMode
                case 1      % (-180,180]
                case 2      % [0,360)
                    index = (phimap<0);
                    phimap(index) = phimap(index) + 360;  
                case 3      % [-270,90)
                    index = (phimap<=180 & phimap>=90);
                    phimap(index) = phimap(index) - 360;  
                case 4      % [-90,270)
                    index = (phimap<=-90 & phimap>=-180);
                    phimap(index) = phimap(index) + 360;                      
            end
        end
        % --- duplicate gixsdata object
        function obj2 = copyhobj(obj1)
            % COPYHOBJ Duplicates gixsdata object.
            %   OBJ2 = COPYHOBJ(OBJ1) duplicate gixsdata object with handle
            %   OBJ1 to another object with a new handle OBJ2.
            %   FigHandle, ImFile, ImFilePath, ImFIleName, ImFileExt,
            %   ImFileInfo properties are not duplicated.
            obj2 = feval(class(obj1));
            m = metaclass(obj1);
            mp = findobj([m.Properties{:}],'SetAccess','public');
            p = {mp.Name};
            % --- for compatibility of new version of gixsdata with new
            % properties PhiMode and PolarizationMode, or ChiMap, LorentzFactorType
            if nnz(strcmpi(p,'PhiMode')) == 0
                copy_order = [1:8,32:-1:17,15:16,9:13];
            elseif nnz(strcmpi(p,'PolarizationMode')) == 0
                copy_order = [1:9,33:-1:18,16:17,10:14];
            elseif nnz(strcmpi(p,'HorizontalPolarizationFraction')) == 0                
                copy_order = [1:10,34:-1:19,17:18,11:15];
            elseif nnz(strcmpi(p,'ChiMap')) == 0                
                copy_order = [1:10,35:-1:19,17:18,11:15];
            elseif nnz(strcmpi(p,'LorentzFactorType')) == 0     % CustomCorrection field is added at the same time
                copy_order =  [1:11,36:-1:20,18:20,12:16]; 
            else
                %copy_order = [1:11,35:-1:20,18:19,12:16];
                copy_order =  [1:11,38:-1:21,18:20,12:16];                               
            end
            for ii=1:length(copy_order)
                obj2.(p{copy_order(ii)}) = obj1.(p{copy_order(ii)});
            end
        end
        % --- plot image
        function imagesc(obj)
            % IMAGESC Scales and displays image data.
            %   
            %   IMAGESC(OBJ) plot MaskedData or SolidAngleCorrectedData 
            %   (if exists) in current figure window using plot parameters
            %   PlotImageFlag, PlotCLimsType, PlotCLims, PlotScale, and 
            %   PlotAxisLabel. It overloads Matlab buildin function IMAGESC.
            obj.FigHandle = gcf;
            himg = imagesc([]);
            set(himg,'tag',['gixsdata:img:',obj.ImFileName]);
            set(get(himg,'Parent'),'Ydir','reverse');
            obj.post_imagesc(obj);
            axis image;
        end % method: imagesc
        function qmaps(obj)
            % QMAPS Calculates q and angle maps, 1D lists and 
            %   SolidAngleCorrectedData. It clears these properties if any 
            %   of the setup parameters are incomplete or invalid. QMAPS 
            %   needs to be manually called when setup parameters are 
            %   changed or updated.
            %
            %   QMAPS(OBJ) or OBJ.QMAPS where OBJ is the gixsdata handle. No
            %   output argument is necessary.
            nan_list = [obj.Beam0,obj.PixelSize,obj.SDD,obj.XEnergy];
            if obj.Geometry == 2
                nan_list = [nan_list,obj.IncidentAngle,obj.Specular];
            end
            if nnz(isnan(nan_list)) ~= 0
                obj.qmaps_initialize(obj);
                return;
            end
            qmaps_ker(obj);
        end % method: qmaps
        function varargout = linecut(obj,x,c,N,dataflag)
            % LINECUT Performs line cut on image data.
            %   [XDATA,YDATA] = LINECUT(OBJ,XFLAG,CONST,NOFPTS,DATAFLAG)
            %   perform line cut on data specified by DATAFLAG for gixsdata
            %   of handle OBJ with predictor parameter XFLAG, constraints
            %   CONST, total number of linecut points NOFPTS.
            %
            %   [XDATA,YDATA,NOFPIXELS] = LINECUT(OBJ,XFLAG,CONST,NOFPTS,DATAFLAG)
            %   returns the number pixels found and used to
            %   calculate the YDATA for each XDATA.
            %
            %   [XDATA,YDATA,NOFPT,MAPDATA] = LINECUT(OBJ,XFLAG,CONST,NOFPTS,DATAFLAG) 
            %   also returns map values corresponding to each XDATA.
            %   MAPDATA is a vector with 2 or 9 columns for transmission (q
            %   and phi) and reflection (all 9 maps) respectively.
            %
            %   XFLAG indicates the predictor parameter. It must be an 
            %   integer scalar. XFLAG chart:
            %       1 - q (total q)
            %       2 - phi (azimuthal angle with respect to direct beam)
            %       3 - qz (for reflection only)
            %       4 - qx (for reflection only)
            %       5 - qy (for reflection only)
            %       6 - qr (parallel q in surface plane; for reflection only)
            %       7 - 2Theta (scattering angle in surface plane; for
            %           reflection only)
            %       8 - Alpha_f (exit angle; for reflection only)
            %       9 - Chi (polar angle; for reflection only)
            %       10 - x pixel
            %       11 - y pixel
            %
            %   CONST specifies the constraints for linecut. It must be a
            %   mx4 array, where m is the number of constraints.
            %   Constraints create a second mask superimposed onto OBJ.Mask
            %   with logical AND operation. If CONST=[], only OBJ.Mask is 
            %   used. Each constraint (a row in CONST) has 4 elements:
            %       CONST = [OPERATOR,FLAG,LOWER,UPPER]
            %   where OPERATOR = 1 or 2 specifies either AND or OR logical
            %   operation is applied to the constraint. FLAG specifies the
            %   name of the constraint. It has the same chart as XFLAG. In
            %   addition, FLAG can also be 12, indicating that the
            %   constraint on that row is not in use. OPREATOR must be AND
            %   when FLAG = 12. LOWER and UPPER are lower and upper limits of
            %   the constraint. They are ignored if FLAG = 12. Constrains
            %   defined in CONST are applied in the following order.
            %   First, a mask array of all ones (each element equals 1) is
            %   created, which is then applied to the 1st constraint
            %   (on the 1st row) using specified logical OPERATOR on the
            %   1st row. This generates a 2nd mask which is then applied to
            %   the 2nd constraint on the 2nd row using specified OPERATOR.
            %   This generats a 3rd mask and so on. Last step is to
            %   apply the (m+1)th mask to main mask OBJ.Mask using AND
            %   operation, generating the final mask which is then used on 
            %   image data for the linecut.
            %
            %   NOFPTS is the total number of predictor points for the
            %   linecut.
            %
            %   DATAFLAG = 1 or 2 specifies either MaskedData or
            %   SolidAngleCorrectedData is used for the linecut. If
            %   SolidAngleCorrectedData does not exist, i.e., no maps,
            %   MaskedData will be used regardless DATAFLAG, and the 
            %   linecut can only be performed with XFLAG = 10 or 11.
            cmask = obj.get_cmask(obj,c);
            switch dataflag
                case 1
                    data = obj.MaskedData(cmask);
                case 2
                    data = obj.SolidAngleCorrectedData(cmask);
            end
            % --- generate xspan and xdata
            [xgrid,ygrid] = meshgrid(1:obj.ImDim(1),1:obj.ImDim(2));
            map_list = {'QMap','PhiMap','QzMap','QxMap','QyMap','QrMap','TwoThetaMap','AlphafMap','ChiMap'};
            switch x
                case num2cell(1:length(map_list))
                    x1DMap = obj.(map_list{x})(cmask);
                case 10
                    x1DMap = xgrid(cmask);
                case 11
                    x1DMap = ygrid(cmask);
            end
            xspan = linspace(min(x1DMap),max(x1DMap),N+1);
            % --- start collecting linecut data
            xdata = ones(N,1)*NaN;
            ydata = ones(N,1)*NaN;
            nofpixels = ones(N,1)*NaN;
            notnan_index = ~isnan(data);
            if nargout == 4             % calculate mapdata
                if obj.Geometry == 1
                    mapdata = ones(N,2)*NaN;
                    map_list_for_data = map_list(1:2);
                else
                    mapdata = ones(N,9)*NaN;
                    map_list_for_data = map_list(1:9);
                end
                mapdata_1DMAP = NaN*ones(length(x1DMap),length(map_list_for_data));
                for jj=1:length(map_list_for_data)
                    mapdata_1DMAP(:,jj) = obj.(map_list_for_data{jj})(cmask);
                end                
            end
            for ii=1:N
                tmp = find(x1DMap<xspan(ii+1) & x1DMap>=xspan(ii) & notnan_index);
                if ~isempty(tmp)
                    nofpixels(ii)  = numel(tmp);
                    ydata(ii)  = mean(data(tmp));
                    xdata(ii)  = mean(x1DMap(tmp));                    
                    if nargout == 4
                        mapdata(ii,:) = mean(mapdata_1DMAP(tmp,:),1);
                    end
                end
            end
            tmp = isnan(ydata);
            ydata(tmp) = [];
            xdata(tmp) = [];
            nofpixels(tmp) = [];
            varargout{1} = xdata;
            varargout{2} = ydata;
            if nargout >= 3
                varargout{3} = nofpixels;
            end
            if nargout == 4
                mapdata(tmp,:) = [];
                map.Data = mapdata;
                if obj.Geometry == 1
                    map.List = {'QMap','PhiMap'};
                else
                    map.List = {'QMap','PhiMap','QzMap','QxMap','QyMap','QrMap','TwoThetaMap','AlphafMap','ChiMap'};
                end
                varargout{4} = map;
            end
        end % method: linecut
        function [xdata,ydata,imgdata,countdata,imgdata_inpaint] = reshape_image(obj,params,dataFlag,varargin) %#ok<*INUSD>
            % RESHAPE_IMAGE Converts image into Cartesian coordinate system. 
            %   [XDATA,YDATA,IMGDATA,MAPDATA,COUNTDATA,IMGDATA_INPAINT] = RESHAPE_IMAGE(OBJ,PARAMS,DATAFLAG)
            %
            %   [XDATA,YDATA,IMGDATA,MAPDATA,COUNTDATA,IMGDATA_INPAINT] = RESHAPE_IMAGE(OBJ,PARAMS,DATAFLAG,INPAINT_METHOD)            
            %   converts image data specified by DATAFLAG for gixsdata of
            %   handle OBJ using parameters PARAMS to a new image IMGDATA
            %   with x, y axis XDATA and YDATA. MAPDATA is a structure
            %   composed of map fields that correspond to the averaged q
            %   and angle values of the pixles on the source image pixels 
            %   that are counted into pixel (x,y) of the new image. 
            %   COUNTDATA indicates the number of pixels on the source
            %   images used for the average for pixel (x,y) of the new
            %   image. IMGDATA_INPAINT is the inpainted data and is empty
            %   is inpaint method option is not specified.
            %
            %   PARAMS is a structure with fields X, Y, XNOfPts, YNOfPts, 
            %   XRange and YRange. Here, X and Y use the same chart as 
            %   XFLAG does in LINECUT method to specify the names of the 
            %   x and y axis. XNOfPts and YNOfPts are the numbers of points
            %   in XDATA and YDATA. XRange and YRange are the ranges of
            %   XDATA and YDATA.
            %
            %   DATAFLAG has the same usage as in LINECUT.      
            %
            %   INPAINT_METHOD specifies the method 0/1/2/3/4/5 for
            %   inpainting. method=0 is default. See function inpaint_nan
            %   for details.
            
            % construct xdata and ydata
            n = params.XNOfPts;
            m = params.YNOfPts;
            if n<m      % reshape an image of larger XNOfPts is faster;
                flip_flag = 1;          
            else
                flip_flag = 0;
            end
%             xdata = linspace(params.XRange(1),params.XRange(2),n);
%             ydata = linspace(params.YRange(1),params.YRange(2),m);
%             d_xdata = xdata(2)-xdata(1);
%             d_ydata = ydata(2)-ydata(1);
%             xdata_b = [xdata-d_xdata/2, xdata(end)+d_xdata/2];
%             ydata_b = [ydata-d_ydata/2, ydata(end)+d_ydata/2];
            % get source image data and map
            if dataFlag == 1
                I0 = obj.MaskedData;
            elseif dataFlag == 2
                I0 = obj.SolidAngleCorrectedData;
            end
            [xgrid,ygrid] = meshgrid(1:obj.ImDim(1),1:obj.ImDim(2));
            map_list = {'QMap','PhiMap','QzMap','QxMap','QyMap','QrMap','TwoThetaMap','AlphafMap','ChiMap'}; 
            switch params.X
                case num2cell(1:length(map_list))
                    xmap = obj.(map_list{params.X});
                case 10
                    xmap = xgrid;
                case 11
                    xmap = ygrid;
            end
            switch params.Y
                case num2cell(1:length(map_list))
                    ymap = obj.(map_list{params.Y});
                case 10
                    ymap = xgrid;
                case 11
                    ymap = ygrid;
            end
            % decide boundaries and span xdata and ydata
            
            xdata = linspace(max(params.XRange(1),min(xmap(:))),min(params.XRange(2),max(xmap(:))),n);
            ydata = linspace(max(params.YRange(1),min(ymap(:))),min(params.YRange(2),max(ymap(:))),m);
            d_xdata = xdata(2)-xdata(1);
            d_ydata = ydata(2)-ydata(1);
            xdata_b = [xdata-d_xdata/2, xdata(end)+d_xdata/2];
            ydata_b = [ydata-d_ydata/2, ydata(end)+d_ydata/2];            
            % flip xy if necessary
            if flip_flag
                foo = xmap; xmap = ymap; ymap = foo;
                foo = xdata_b; xdata_b = ydata_b; ydata_b = foo;
                foo = n; n=m; m=foo;
            end
            % fill unmasked region by NaN
%             xmap_dummy = xmap;  % dummy variables to label zero mask pixels
%             ymap_dummy = ymap;
%             xmap(~obj.Mask) = NaN;
%             ymap(~obj.Mask) = NaN;
%             
            I0(~obj.Mask) = NaN;
            I0_mask = nan(size(I0)) ;    % reshape mask file at the same time to identify the mask regions in the reshaped image
            I0_mask(~obj.Mask) = 1;
            % initialize
            imgdata         = NaN*ones(m,n);
            imgdata_mask    = NaN*ones(m,n);
%             imgdata_ismask0 = false(m,n);      % label empty pixels (those fall in the mask=0 and undefined regions)            
            countdata   = NaN*ones(m,n);
%            mapname = {'QMap','PhiMap','QzMap','QxMap','QyMap','QrMap','TwoThetaMap','AlphafMap'};
%             for mm=1:length(mapname)
%                 mapdata.(mapname{mm}) = NaN*ones(m,n);        % mapdata to store the mean q values of ROI                  
%             end
            ymap_min = min(ymap(:));
            ymap_max = max(ymap(:));
            for ii = 1:m
                if ydata_b(ii)>ymap_max || ydata_b(ii+1)<ymap_min
                    continue;
                end
                ind_y = (ymap>=ydata_b(ii) & ymap<ydata_b(ii+1));
                xmap_tmp = xmap(ind_y);
                if isempty(xmap_tmp), continue; end
                I0_tmp = I0(ind_y);
                I0_mask_tmp = I0_mask(ind_y);
%                 mask0_tmp = ~obj.Mask(ind_y);    % mask=0 pixles
%                 if any(mask0_tmp)
%                     keyboard
%                 end
%                 for mm = 1:length(mapname)
%                     if ~isempty(obj.(mapname{mm}))
%                         map0_tmp.(mapname{mm}) = obj.(mapname{mm})(ind_y) ;
%                     end
%                 end
                xmap_min = min(xmap_tmp(:));
                xmap_max = max(xmap_tmp(:));
                for jj=1:n
                    if xdata_b(jj)>xmap_max || xdata_b(jj+1)<xmap_min
                        continue;
                    end
                    ind_x = (xmap_tmp>=xdata_b(jj) & xmap_tmp<xdata_b(jj+1));
                    nnz_x = nnz(ind_x);
                    if nnz_x    % pixels has finite values
                        I_tmp = I0_tmp(ind_x);
                        I_mask_tmp = I0_mask_tmp(ind_x);
%                         mask_tmp = mask0_tmp(ind_x);
                        imgdata(ii,jj) = sum(I_tmp(:));
                        imgdata_mask(ii,jj) = sum(I_mask_tmp(:));
%                         if any(mask_tmp)
%                             imgdata_ismask0(ii,jj) = true;
%                         end
                        countdata(ii,jj) = nnz(ind_x);
%                         for mm=1:length(mapname)
%                             if ~isempty(obj.(mapname{mm}))
%                                 map_tmp = map0_tmp.(mapname{mm})(ind_x);
%                                 mapdata.(mapname{mm})(ii,jj) = sum(map_tmp(:));
%                             end
%                         end
                    end
                end
            end
            % normalize imgdata and mapdata
            imgdata = imgdata./countdata;
            imgdata_mask = imgdata_mask./countdata;
            countdata(isnan(countdata)) = 0;
            % flip xy back if necessary
            if flip_flag 
                imgdata = imgdata';
                imgdata_mask = imgdata_mask';                
                countdata = countdata';
%                 imgdata_ismask0 = imgdata_ismask0';
            end           
%             for mm=1:length(mapname)
%                 if isempty(obj.(mapname{mm})),
%                     mapdata.(mapname{mm}) = [];
%                 else
%                     mapdata.(mapname{mm}) = mapdata.(mapname{mm})./countdata;
%                 end
%             end
            if nargin == 4
                % inpaint in log scale
                imgdata_inpaint = inpaint_nans(imgdata,varargin{1});
                imgdata_inpaint(imgdata_inpaint<0) = 0;     % set negative inpainted values to 0
                if isequal([params.X,params.Y],[6,3])
                    qr = xdata;
                    qz = ydata;
                    imgdata_tmp = imgdata_inpaint;
                    qrqz_flag = 1;        % for [qr,qz];
                elseif isequal([params.X,params.Y],[3,6])
                    qr = ydata;
                    qz = xdata;
                    imgdata_tmp = imgdata_inpaint';
                    qrqz_flag = -1;       % for [qz,qr];
                else        % for no qr, qz reshaping
                    qrqz_flag = 0;
                end
                if qrqz_flag ~= 0 
                    k = 2*pi/obj.XWavelength;
                    alphai = obj.IncidentAngle;
                    qx=(sqrt(1-((qz(:)/k)-sind(alphai)).^2)-cosd(alphai))*k;   % qx boundary for forbidden region
                    qx_2D = repmat(qx,[1,length(qr)]);
                    [qr_2D,~] = meshgrid(qr,qz);
                    ind_nan = (abs(qr_2D - qx_2D) + abs(qr_2D + qx_2D)) <= (2*abs(qx_2D) + abs(qr(1)-qr(2))/100 );
                    imgdata_tmp(ind_nan) = nan;
                    if qrqz_flag == 1
                        imgdata_inpaint = imgdata_tmp;
                    elseif qrqz_flag == -1
                        imgdata_inpaint = imgdata_tmp';
                    end 
                end
                imgdata_inpaint(isfinite(imgdata_mask)) = NaN;
%                imgdata_inpaint(imgdata_ismask0) = NaN;
            else
                imgdata_inpaint = [];
            end
        end % method: reshape image
        function [I,roimap,mapname,pcounts,cmask] = roicount(obj,const,dataflag)
            % ROICOUNT Calculates the intensity of the pseduo point counter
            % defined by the constraints and mask.
            %
            %   [I,ROIMAP,MAPNAME,PCOUNTS] = ROICOUNT(OBJ,CONST,DATAFLAG) 
            %   returns the total intensity I (scalar), where OBJ the 
            %   gixsdata handle and CONST is the constraints. CONST=[] for
            %   no constraint. DATAFLAG has the same usage as in LINECUT. 
            %   ROIMAP is a list (4x1 or 11x1 depending on the scattering 
            %   geometry) of the averaged map values. The order of the list
            %   follows the XFLAG defined in method LINECUT. MAPNAME is the 
            %   list of corresponding map names. PCOUNTS is the number of 
            %   pixels in the ROI. CMASk is the logical matrix represnting
            %   the constrainted mask (mask && constraints).
            
            cmask = obj.get_cmask(obj,const);
            if obj.Geometry == 1
                mapname = {'QMap','PhiMap'};
            else
                mapname = {'QMap','PhiMap','QzMap','QxMap','QyMap','QrMap','TwoThetaMap','AlphafMap','ChiMap'};
            end
            roimap = NaN*ones(1,length(mapname)+2);
            for ii=1:length(mapname)
                roimap(ii) = mean(obj.(mapname{ii})(cmask));
            end
            [xgrid,ygrid] = meshgrid(1:obj.ImDim(1),1:obj.ImDim(2));
            roimap(end-1) = mean(xgrid(cmask));
            roimap(end) = mean(ygrid(cmask));
            mapname = [mapname,'X Pixel','Y Pixel'];
            if dataflag == 1
               I = obj.MaskedData(cmask) ;
            elseif dataflag == 2
                I = obj.SolidAngleCorrectedData(cmask);
            end
            I = sum(I(:));
            pcounts = nnz(cmask);
        end % method: roicount
        function [peak,fresult,vertex] = findpeak(obj,vertex,dataFlag,fitFlag,startpoint,opts)
            % FINDPEAK Finds 2D diffraction peak positions.
            %   [PEAK,FRESULT,VERTEX] = FINDPEAK(OBJ,VERTEX,DATAFLAG,
            %   FITFLAG,STARTPOINT,OPTS) finds peak position of a region of
            %   interest (ROI) defined in a polygon. The vertexes of the 
            %   polygon (x and y pixels) are defined by VERTEX (two column
            %   vector with the 1st and 2nd columns for x and y pixels).                         
            %
            %   DATAFLAG has the same usage as in method LINECUT.
            %
            %   FITFLAG = 1 or 2 specifies either COM or BVND model is used
            %   to determine the peak position.
            %
            %   STARTPOINT (1x9 vector) gives the initial values of the 
            %   BVND parameters for the fitting. These parameters are
            %   [A, Cx, Cy, C, Rho, Sigma_x, Sigam_y, X0, Y0]. STARTPOINT 
            %   can be [] (empty)for automatic assignment. It is not used 
            %   for the COM model.
            %
            %   OPTS gives the optimization parameters for the BVND 
            %   fittings. It can be defined using optimset function. Refer
            %   to fminsearch (Matlab built-in function) and fminsearchbnd 
            %   (by John D'Errico and included in the package) for usage.
            %   It can be [] (empty) for automatic assignment. It is not
            %   used for the COM model.
            %
            %   FINDPEAK returns the peak position (in x and y pixels) and 
            %   its map values in PEAK (structure). The fitting parameters 
            %   and results are given in FRESULT (structure). The vertexes
            %   for the ROI polygon are passed to VERTEX.
            
            % --- get roi for fitting
            x = vertex(:,1); y=vertex(:,2);
            x = round(x'); y = round(y');
            xpoly = [x, x(1)]; ypoly = [y,y(1)];
            [xgrid,ygrid] = meshgrid(1:obj.ImDim(1),1:obj.ImDim(2));
            in = inpolygon(xgrid(:),ygrid(:),xpoly,ypoly);
            BW = reshape(in,size(obj.RawData));
            roi = obj.Mask & BW;            
            % --- get roi data 
            if dataFlag == 1
                data = obj.MaskedData(roi);
            elseif dataFlag == 2
                data = obj.SolidAngleCorrectedData(roi);
            end
            data = double(data(:));
            xdata = xgrid(roi); xdata = xdata(:);
            ydata = ygrid(roi); ydata = ydata(:);
            % --- get rectangular roi
            xmin = min(x); ymin = min(y);
            xmax = max(x); ymax = max(y);
            [xdatarec,ydatarec] = meshgrid(xmin:xmax,ymin:ymax);
            if dataFlag == 1
                datarec = obj.MaskedData(ymin:ymax,xmin:xmax);
            elseif dataFlag == 2
                datarec = obj.SolidAngleCorrectedData(ymin:ymax,xmin:xmax);
            end
            datarec = double(datarec);
            % --- com or bvnd
            if fitFlag == 1         % use center of mass
                x0 = sum(xdata.*data)/sum(data);
                y0 = sum(ydata.*data)/sum(data);
            elseif fitFlag == 2    % use bvnd fitting
                % --- boundary
                lb =        [0     -Inf   -Inf  -Inf   -1   0   0    min(xdata)  min(ydata)];
                ub =        [Inf   Inf    Inf   Inf    1    Inf Inf  max(xdata)  max(ydata)];
                % --- start point
                if isempty(startpoint)
                    [~,x0_index] = max(data); x0_start = xdata(x0_index);
                    [~,y0_index] = max(data); y0_start = ydata(y0_index);                    
                    startpoint = [max(data)-min(data), 0,   0,  mean(data), 0, 1,1, x0_start, y0_start];
                end
                % --- options
                if isempty(opts)
                    opts = optimset('fminsearch');
                    opts.Display = 'off';
                    opts.TolX = 1e-6;
                    opts.TolFun = 1e-6;
                    opts.MaxFunEvals = 800;
                    opts.MaxIter = 600;
                end
                % --- start fitting
                [par,~,exitflag,output] = fminsearchbnd(@res_bvnd,startpoint,lb,ub,opts,xdata,ydata,data);                
                x0 = par(8);
                y0 = par(9);
                fval = bvnd(par,xdata,ydata);
                fvalrec = bvnd(par,xdatarec,ydatarec);
            else
                return;
            end
            % --- construct result
            peak.X = x0;
            peak.Y = y0;
            map_str_list = {'Q','Phi','Qz','Qx','Qy','Qr','TwoTheta','Alphaf'};
            for ii=1:length(map_str_list)
                map = map_str_list{ii};
                if ~isempty(obj.([map,'Map']))
                    peak.(map) = double(obj.([map,'Map'])(round(y0),round(x0)));
                else
                    peak.(map) = [];
                end
            end
            fresult = struct;
                fresult.imfilename = obj.ImFileName;
            if fitFlag == 1
                fresult.model = 'com';
            elseif fitFlag == 2
                fresult.model = 'bvnd';
                fresult.par = par;
                fresult.exitflag = exitflag;
                fresult.output = output;
                fresult.xdata = xdata;
                fresult.ydata = ydata;
                fresult.data = data;
                fresult.fval = fval;
                fresult.xdatarec = xdatarec;
                fresult.ydatarec = ydatarec;
                fresult.datarec = datarec;
                fresult.fvalrec = fvalrec;
            end
         end % method: findpeak
        function pixel = mappixel(obj,mapValue,mapFlag,directSearchFlag,opts)
            % MAPPIXEL Finds pixel indices of a point with given map values.
            %   PIXEL = MAPPIXEL(OBJ,MAPVALUE,MAPFLAG,DIRECTSEARCHFLAG,OPTS)
            %   maps out the pixel indices of a point given two map values.
            %
            %   MAPVALUE (1x2) provides two independent map (i.e.,
            %   [TwoTheta, Alphaf] or [Qz, Qy]) values of a point.
            %
            %   MAPFLAG (1x2) specifies the name of the two map values. It
            %   has the same usage as XFLAG in method linecut, except that
            %   the value can only be 1 through 8.
            %
            %   DIRECTSEARCHFLAG = 1 or 0 specifies the mapping method. 1
            %   is for a direct search, i.e., finding the nearest pixel 
            %   with integer pixel indices. 0 is for a search through
            %   optimization fitting for the pixel; its indices may not be
            %   integers. 
            %
            %   OPTS gives the optimization parameters for the optimization
            %   fitting. It can be defined using optimset function. Refer
            %   to fminsearch (Matlab built-in function) and fminsearchbnd 
            %   (by John D'Errico and included in the package) for usage.
            %   It can be [] (empty) for automatic assignment. It is not
            %   used when DIRECTSEARCHFLAG=1;
            %
            %   MAPPIXEL returns the PIXEL indices (1x2) if the mapping is
            %   valid and the pixel is within the dimension of the image.
            %   Otherwise, PIXEL=[].
            v1 = mapValue(1);
            v2 = mapValue(2);
            map_str_list = {'Q','Phi','Qz','Qx','Qy','Qr','TwoTheta','Alphaf'};
            map1 = obj.([map_str_list{mapFlag(1)},'Map']);
            map2 = obj.([map_str_list{mapFlag(2)},'Map']);
            %[xgrid,ygrid] = meshgrid(1:obj.ImDim(1),1:obj.ImDim(2));  
            % --- determine startpoint or direct pixel result automatically
            r2 = (map1-v1).^2+(map2-v2).^2;
            [r2_min,idx]=min(r2(:));
            [y0,x0]=ind2sub(size(r2),idx(1));
            % --- determine if [x0,y0] is good values by finding nearest
            % point in a smaller ROI (to save time)
            v1_near = map1(y0,x0);
            v2_near = map2(y0,x0);
            % define roi size
            dx = 10;        % roi x size is dx*2+1
            dy = 10;        % roi y size is dy*2+1
            roix = x0-dx:x0+dx;     % roi x index range
            roiy = y0-dy:y0+dy;     % roi y index range
            % correct roi box near image border
            if roix(1)<1        % roix offset
                roix_offset = 1-roix(1);
            elseif roix(end)>obj.ImDim(1)
                roix_offset = obj.ImDim(1)-roix(end);
            else
                roix_offset = 0;
            end
            roix = roix+roix_offset;
            if roiy(1)<1        % roiy offset
                roiy_offset = 1-roiy(1);
            elseif roiy(end)>obj.ImDim(2)
                roiy_offset = obj.ImDim(2)-roiy(end);
            else
                roiy_offset = 0;
            end
            roiy = roiy+roiy_offset;
            map1_roi = map1(roiy,roix);     % roi map
            map2_roi = map2(roiy,roix);     
            r2_near = (map1_roi-v1_near).^2+(map2_roi-v2_near).^2; 
            r2_near = sort(r2_near);
            r2_near_min = r2_near(2);
            if r2_min<=r2_near_min
                pixel = [x0,y0];
            else
                pixel = [];
            end
            if directSearchFlag == 1 || isempty(pixel)
                return;
            end
            % --- use square ROI around startpoint for fitting to save
            % computing time)
            [xgrid_roi,ygrid_roi] = meshgrid((-dx:dx)+roix_offset,(-dy:dy)+roiy_offset); % roi grid
            % --- boundary, new startpoint (roi), and fit options
            lb = [-dx -dy];         
            ub = [dx  dy];
            startpoint = [0,0];
            if isempty(opts)
                opts = optimset('fminsearch');
                opts.Display = 'off';
                opts.TolX = 1e-6;
                opts.TolFun = 1e-6;
                opts.MaxFunEvals = 800;
                opts.MaxIter = 600;
            end
            % --- start fitting
            [result,feval,exitflag,output] = fminsearchbnd(...
                @res_map,startpoint,lb,ub,opts,...
                map1_roi,map2_roi,xgrid_roi,ygrid_roi,v1,v2);
            if exitflag ~= 1
                pixel = [];
            else
                ptx = result(1)+x0;
                pty = result(2)+y0;
                pixel = [ptx,pty];
            end
        end % method: mappoint
        function c_lf = lorentzfactor(obj)
            % LORENTZFACTOR Calculate Lorentz factor with given type.
            switch obj.LorentzFactorType
                case 1  % no correction
                    c_lf = ones(size(obj.SolidAngleCorrectedData));
                case 2  % in-plane random 3D
                    c_lf = cosd(obj.IncidentAngle)*cosd(obj.AlphafMap).*sind(obj.TwoThetaMap);
                case 3  % in-plane random 2D
                    c_lf = sind(obj.TwoThetaMap);
                case 4  % completely random
                    k = 2*pi/obj.XWavelength;
                    transmission_th = asind(obj.QMap/(2*k));
                    c_lf = 4*(sind(transmission_th)).^2.*cosd(transmission_th);
            end
            c_lf = abs(c_lf); c_lf = max(c_lf(:))./c_lf; % inverse to get Lorentz factor
        end % method: lorentzfactor
    end % dynamic method
    methods (Static=true)
        function qmaps_initialize(obj)
            % QMAPS_INITIALIZE Clears all maps, 1D lists, and 
            %   SolidAngleCorrectedData.
            %
            %   OBJ.QMAPS_INITIALIZE(OBJ) where OBJ is the gixsdata handle.
            obj.QMap        = [];
            obj.QxMap      = [];
            obj.QzMap      = [];
            obj.QyMap      = [];
            obj.QrMap      = [];
            obj.PhiMap     = [];
            obj.TwoThetaMap   = [];
            obj.AlphafMap  = [];
            obj.ChiMap  = [];
            obj.Qx1DList   = [];
            obj.Qy1DList   = [];
            obj.Qz1DList   = [];
            obj.Alphaf1DList = [];
            obj.TwoTheta1DList = [];
            obj.SolidAngleCorrectedData = [];
            %obj.EfficiencyCorrection{1,3} = 0;
            %obj.EfficiencyCorrection{2,3} = Inf;
        end
        function solidangle_correction(obj)
            % SOLIDANGLE_CORRECTION Performs solid angle correction to 
            %   RawData, regardless the Mask. Flat field and efficiency 
            %   correction, if defined, are also applied. Given valid setup
            %   parameters, the correction is done automatically when 
            %   RawData is updated. If RawData is set to empty, 
            %   SolidAngleCorrectedData will be cleared.
            %
            %   OBJ.SOLIDANGLE_CORRECTION(OBJ), where OBJ is the gixsdata
            %   handle.
            nan_list = [obj.Beam0,obj.PixelSize,obj.SDD,obj.XEnergy];
            if nnz(isnan(nan_list)) ~= 0 || isempty(obj.QMap)
                obj.SolidAngleCorrectedData = [];
                return;
            end
            if isempty(obj.RawData)
                return;
            end            
            [x,y] = meshgrid(1:obj.ImDim(1),1:obj.ImDim(2));
            r_complex = obj.PixelSize(1)*(x-obj.Beam0(1)) + ...
                1i*obj.PixelSize(2)*(obj.Beam0(2)-y);
            r = abs(r_complex);
            R = sqrt(r.^2+obj.SDD^2);                               % pixel to sample distance
            theta = atan(r/obj.SDD);                                    % exit angle map
            dOmega = obj.PixelSize(1)*obj.PixelSize(2)*cos(theta)./R.^2;                        % solid angle map
            dOmega = dOmega/max(max(dOmega));                         % normalize solid angle
            obj.SolidAngleCorrectedData = single(obj.RawData./dOmega);
            % perform flat field correction
            if isequal(size(obj.SolidAngleCorrectedData),size(obj.FlatField))
                obj.SolidAngleCorrectedData = obj.SolidAngleCorrectedData.*obj.FlatField;
            end
            % perform custom correction
            if isequal(size(obj.SolidAngleCorrectedData),size(obj.CustomCorrection))
                obj.SolidAngleCorrectedData = obj.SolidAngleCorrectedData.*obj.CustomCorrection;
            end            
            % perform efficiency correction
            effcdata = obj.efficiency_correction(obj);
            if isequal(size(obj.SolidAngleCorrectedData),size(effcdata))
                obj.SolidAngleCorrectedData = obj.SolidAngleCorrectedData.*effcdata;
            end
            % perform polarization correction
            polardata = obj.polarization_correction(obj);
            if isequal(size(obj.SolidAngleCorrectedData),size(polardata))
                obj.SolidAngleCorrectedData = obj.SolidAngleCorrectedData./polardata;
            end
            % perform lorentz correction 
            CL = lorentzfactor(obj);
            obj.SolidAngleCorrectedData = obj.SolidAngleCorrectedData./CL;
        end % method: solid angle correction
        function c = efficiency_correction(obj)
            % EFFICIENCY_CORRECTION Calculates efficiency correction array
            %   from property EfficiencyCorreciton. It requires package 
            %   REFRAC to calculate the attenuation coefficients of the 
            %   path medium path and detector sensor materials. If the 
            %   package does not exist, an array of all ones is returned 
            %   without warnings.
            %
            %   C = OBJ.EFFICIENCY_CORRECTION(OBJ), where OBJ is the gixsdata
            %   handle and C is the correction array of the same dimension
            %   as the image data.
            
            % --- calculate diffraction theta
            k = 2*pi/obj.XWavelength;
            theta = 2*asin(obj.QMap/(2*k));
            % --- calculate corrections
            objeffc = obj.EfficiencyCorrection;
            % 1st correction coefficient for path aborption
            try
                r1 = refrac(objeffc{1,1},obj.XEnergy,objeffc{1,2});
                att1 = r1.attLength*10;     % convert from cm to mm
                c1 = 1./exp(-objeffc{1,3}./cos(theta)/att1);
            catch
                c1 = ones(size(obj.QMap));
            end
            % 2nd correction coefficient for chip attenuation
            try
                r2 = refrac(objeffc{2,1},obj.XEnergy,objeffc{2,2});
                att2 = r2.attLength*10;  % convert from cm to mm
                c2 = 1./(1-exp(-objeffc{2,3}./cos(theta)/att2));
            catch
                c2 = ones(size(obj.QMap));
            end
            c = c1.*c2;
            c = single(c/max(c(:)));  
        end % method: calculate efficiency correction coefficient
        function c = polarization_correction(obj)
            % POLARIZATION_CORRECTION Calculates polarization correction 
            %   array using polarization parameters from PolarizationMode 
            %   and HorizontalPolarizationFraction. 
            %
            %   C = OBJ.POLARIZATION_CORRECTION(OBJ), where OBJ is the gixsdata
            %   handle and C is the correction array of the same dimension
            %   as the image data.
                        
            k = 2*pi/obj.XWavelength;
            [x,y] = meshgrid(1:obj.ImDim(1),1:obj.ImDim(2));
            r_complex = obj.PixelSize(1)*(x-obj.Beam0(1)) + ...
                1i*obj.PixelSize(2)*(obj.Beam0(2)-y);              
            gamma = atan(real(r_complex)/obj.SDD);
            delta = atan(imag(r_complex)./sqrt(obj.SDD^2+real(r_complex).^2));
            uc = (1+cos(delta).^2.*cos(gamma).^2)/2; % for unpolarized
            hpf = obj.HorizontalPolarizationFraction;
            ch = (1-cos(delta).^2.*sin(gamma).^2);  % for fully horizontal
            cv = cos(delta).^2;                     % for fully vertical
            switch obj.PolarizationMode
                case 1
                    c = ones(obj.ImDim(2),obj.ImDim(1));
                case 2
                    c = hpf*ch + (1-hpf)*cv;
                case 3
                    c = hpf*cv + (1-hpf)*ch;
                case 4
                    c = uc;
            end
            c = single(c);
        end % method: calculate polarization correction coefficient
        function mask_initialize(obj)
            % MASK_INITIALIZE Initializes logical Mask from RawData. 
            %   Non-negative elements are denoted by 1 and the rest by 0. 
            %   User-defined Mask will be lost. 
            %
            %   OBJ.MASK_INITIALZE(OBJ), where OBJ is the gixsdata handle.
            obj.Mask = (obj.RawData>=0);
        end % method: mask_initialize
        function post_imagesc(obj)
            % POST_IMAGESC Updates the current figure with updated image 
            %   data (either MaskedData or SolidAngleCorrectedData) or 
            %   updated plot parameters.
            %
            %   OBJ.POST_IMAGESC(OBJ), where OBJ is the gixsdata handle.
            if obj.PlotImageFlag == 1 || isempty(obj.SolidAngleCorrectedData)
                cdata = obj.MaskedData;
            elseif obj.PlotImageFlag == 2
                cdata = obj.SolidAngleCorrectedData;
                cdata(~obj.Mask) = NaN;    
            elseif obj.PlotImageFlag == 3
                cdata = obj.RawData;
            end
            himg = findall(obj.FigHandle,'tag',['gixsdata:img:',obj.ImFileName]);
            haxes = get(himg,'parent');
            if obj.PlotScale==1
                %cdata = obj.MaskedData;
                cdata(isnan(cdata)) = -Inf;
                clims = obj.PlotCLims;
            elseif obj.PlotScale==2
                cdata = log10(cdata);
                cdata(~isreal(cdata)) = 0;      % for RawData<0
                cdata = real(cdata);                
                %cdata = log10(obj.MaskedData);
                cdata(isnan(cdata)) = -Inf;
                clims = log10(obj.PlotCLims);
                if ~isreal(clims(1)),clims(1) = min(eps,clims(2)-eps); end
            end
            set(himg,'CData',cdata);
            set(himg,'XData',1:obj.ImDim(1),'YData',1:obj.ImDim(2));
            set(haxes,'CLim',clims);
            set(haxes,'ydir','reverse');
            set(haxes,'xtickmode','auto','xticklabelmode','auto');
            set(haxes,'ytickmode','auto','yticklabelmode','auto');
            set(zoom(obj.FigHandle),'ActionPostCallback','');
            set(pan(obj.FigHandle),'ActionPostCallback','');
            switch obj.PlotAxisLabel
                case 1
                    xlabel_str = 'x pixel';
                    ylabel_str = 'y pixel';
                    axis_mode(haxes,xlabel_str,ylabel_str)
                case 2
                    xlabel_str = ['q_y (',char(197),'^{-1})'];
                    ylabel_str = ['q_z (',char(197),'^{-1})'];
                    axis_mode(haxes,xlabel_str,ylabel_str)
                    set(zoom(obj.FigHandle),'ActionPostCallback',{@zoompostcallback,haxes,obj.Qy1DList,obj.Qz1DList});
                    set(pan(obj.FigHandle),'ActionPostCallback',{@zoompostcallback,haxes,obj.Qy1DList,obj.Qz1DList});
                    feval(@zoompostcallback,[],[],haxes,obj.Qy1DList,obj.Qz1DList);
                case 3
                    xlabel_str = '2\theta (deg)';
                    ylabel_str = '\alpha_f (deg)';
                    axis_mode(haxes,xlabel_str,ylabel_str)
                    set(zoom(obj.FigHandle),'ActionPostCallback',{@zoompostcallback,haxes,obj.TwoTheta1DList,obj.Alphaf1DList});
                    set(pan(obj.FigHandle),'ActionPostCallback',{@zoompostcallback,haxes,obj.TwoTheta1DList,obj.Alphaf1DList});
                    feval(@zoompostcallback,[],[],haxes,obj.TwoTheta1DList,obj.Alphaf1DList);
            end
            str = [obj.ImFileName,obj.ImFileExt];
            title(haxes,titlestr(str));
            set(get(haxes,'parent'),'Name',str);
            % --- add data cursor
            hc = datacursormode(obj.FigHandle);
            set(hc,'UpdateFcn',{@gixsdatacursor,obj.FigHandle});
            set(obj.FigHandle,'UserData',obj);
        end % method: post_imagesc
        function cmask = get_cmask(obj,const)
            % GET_CMASK Calculates mask defined by linecut constraints and 
            %   Mask.
            %
            %   CMASK = OBJ.GET_CMASK(OBJ,CONST), where OBJ is the gixsdata
            %   handle and CONST is the linecut constraints. It returns a 
            %   logical array CMASK to be used for the linecut.
            %
            %   See also LINECUT.
            if isempty(const)
                cmask = obj.Mask;
                return;
            end
            nc = size(const,1);         % number of constraints
            cmask = nan;
            for ii=1:nc
                clower = min(const(ii,3:4));
                cupper = max(const(ii,3:4));
                switch const(ii,2)
                    case 1
                        tmp_cmask = obj.QMap>clower & obj.QMap<cupper;
                    case 2
                        tmp_cmask = obj.PhiMap>clower & obj.PhiMap<cupper;
                    case 3
                        tmp_cmask = obj.QzMap>clower & obj.QzMap<cupper;
                    case 4
                        tmp_cmask = obj.QxMap>clower & obj.QxMap<cupper & obj.AlphafMap>0;
                    case 5
                        tmp_cmask = obj.QyMap>clower & obj.QyMap<cupper;
                    case 6
                        tmp_cmask = obj.QrMap>clower & obj.QrMap<cupper;
                    case 7
                        tmp_cmask = obj.TwoThetaMap>clower & obj.TwoThetaMap<cupper;
                    case 8
                        tmp_cmask = obj.AlphafMap>clower & obj.AlphafMap<cupper;
                    case 9
                        tmp_cmask = obj.ChiMap>clower & obj.ChiMap<cupper;
                    case 10
                        tmp_cmask = false(size(obj.Mask));
                        tmp_cmask(:,clower:cupper) = 1;
                    case 11
                        tmp_cmask = false(size(obj.Mask));
                        tmp_cmask(clower:cupper,:) = 1;
                    case 12
                        continue;
                end
                if ~islogical(cmask) || ii==1
                    cmask = tmp_cmask;
                else
                    if const(ii,1) == 1
                        cmask = cmask & tmp_cmask;
                    elseif const(ii,1) == 2
                        cmask = cmask | tmp_cmask;
                    end
                end
            end
            if ~islogical(cmask)
                cmask = obj.Mask;                
            else
                cmask = cmask & obj.Mask;                            
            end
        end %method: get_cmask
    end % static method
end% classdef

% Internal functions
function qmaps_ker(obj)
% Called by QMAPS only.
k = 2*pi/obj.XWavelength;
[x,y] = meshgrid(1:obj.ImDim(1),1:obj.ImDim(2));
rx = (x-obj.Beam0(1))*obj.PixelSize(1);  
ry = (obj.Beam0(2)-y)*obj.PixelSize(2);
r_complex = rx+1i*ry;               % pixel position on camera wrt beam0
r = abs(r_complex);
R = sqrt(r.^2+obj.SDD^2);                               % pixel-to-sample distance
theta = atan(r/obj.SDD);                                    % exit angle map
obj.QMap = single(2*k*sin(theta/2));                                   % q map

% --- below 7 lines are commented out on 2014/12/1; solid angel correction
% will be done after qmaps is calculated
%dOmega = obj.PixelSize(1)*obj.PixelSize(2)*cos(theta)./R.^2;                        % solid angle map
%dOmega = dOmega/max(max(dOmega));                         % normalize solid angle
% obj.SolidAngleCorrectedData = single(obj.RawData./dOmega);
% obj.SolidAngleCorrectedData = obj.SolidAngleCorrectedData.*obj.FlatField;
% effcdata = obj.efficiency_correction(obj);
% polardata = obj.polarization_correction(obj);
% obj.SolidAngleCorrectedData = obj.SolidAngleCorrectedData.*effcdata./polardata;

% --- maps for specified geometry
if obj.Geometry == 1            % for transmission (no tilt is considered)
    obj.PhiMap = single(angle(r_complex)*180/pi); % phi map without tilt
    obj.TwoThetaMap = [];
    obj.AlphafMap = [];
    obj.ChiMap = [];    
    obj.QzMap = [];
    obj.QyMap = [];
    obj.QxMap = [];
    obj.QrMap = [];
    obj.Qx1DList   = [];
    obj.Qy1DList   = [];
    obj.Qz1DList   = [];
    obj.Alphaf1DList = [];
    obj.TwoTheta1DList = [];
elseif obj.Geometry == 2        % for reflection
    specular = (obj.Specular(1)-obj.Beam0(1)) + 1i*(obj.Beam0(2)-obj.Specular(2));
    chi = angle(specular)-pi/2;             % tilt angle
    r_complex = r_complex*exp(-1i*chi);      % rotate r to correct tilt
    obj.PhiMap = single(angle(r_complex)*180/pi);   % correct PhiMap with respect to sample horizon
    % assume detector is perpendicular to the beam
    gamma = atan(rx/obj.SDD);   
    delta = atan(ry./sqrt(obj.SDD^2+rx.^2));
    % incident wave vector in lab system
    kx_in = k;
    ky_in = 0;
    kz_in = 0;
    % scatterred wave vector in lab system
    kx_sc = k*cos(delta).*cos(gamma);
    ky_sc = k*cos(delta).*sin(gamma);
    kz_sc = k*sin(delta);
    % q maps in lab system
    qx = kx_sc - kx_in;
    qy = ky_sc - ky_in;
    qz = kz_sc - kz_in;
    % coordinate transformation matrix for the rotation of incident angle
    % the rotation angle is -alpha
    alpha = obj.IncidentAngle*pi/180;
        Ry_rot = [
        cos(alpha)   0   sin(-alpha);
        0            1   0;
        -sin(-alpha) 0   cos(alpha)];
    % coordinate transformation matrix for the rotation of chi
    % rotation angle is chi
    Rx_rot = [
        1   0           0;
        0   cos(chi)    -sin(chi);
        0   sin(chi)    cos(chi);];
    % total tranform matrix
    R_rot = Ry_rot*Rx_rot;
    % calculate q maps in sample system
    [qx,qy,qz] = frametransform(R_rot,qx,qy,qz); 
    obj.QxMap = single(qx);
    obj.QyMap = single(qy);
    obj.QzMap = single(qz);
    sign_idx = sign(obj.QyMap);
    sign_idx(sign_idx==0) = 1;    
    obj.QrMap = single(sign_idx.*sqrt(qx.^2+qy.^2));
    obj.AlphafMap   = single(asin(qz/k - sin(alpha))*180/pi);
    obj.TwoThetaMap = single(asin(qy/k./cos(obj.AlphafMap*pi/180))*180/pi);
    obj.ChiMap = single(180/pi*angle(obj.QzMap+1i*obj.QrMap));

%    % --- old geomeetric method to calcualte q maps; commented out on 2012/0727
%     r_horizon_spot = 1i*tan(obj.IncidentAngle*pi/180)*obj.SDD;         % sample horiztion on specular rod
%     r_horizon = (r_complex-r_horizon_spot);       % r with respect to horizon spot
%     r_tmp1 = real(r_complex)+r_horizon_spot;
%     R_tmp1 = sqrt(abs(r_tmp1).^2+obj.SDD^2);
%     R_horizon_spot = sqrt(abs(r_horizon_spot).^2+obj.SDD^2);
%     TwoThetaMap = single((180/pi)*sign(real(r_complex)).*acos((R_tmp1.^2+R_horizon_spot^2-real(r_complex).^2)./(2*R_tmp1*R_horizon_spot)));
%     AlphafMap = single((180/pi)*sign(imag(r_horizon)).*real(acos((R_tmp1.^2+R.^2-imag(r_horizon).^2)./(2*R_tmp1.*R))));
%     QzMap = single(k*(sin(AlphafMap*pi/180)+sin(obj.IncidentAngle*pi/180)));
%     QyMap = single(k*cos(AlphafMap*pi/180).*sin(TwoThetaMap*pi/180));
%     QxMap = -single(k*(cos(obj.IncidentAngle*pi/180)-cos(AlphafMap*pi/180).*cos(TwoThetaMap*pi/180)));
%     QrMap = single(sqrt(QyMap.^2 + QxMap.^2));

%     r_horizon_spot = 1i*tan(obj.IncidentAngle*pi/180)*obj.SDD;         % sample horiztion on specular rod
%     r_horizon = (r_complex-r_horizon_spot);       % r with respect to horizon spot
%     r_tmp1 = real(r_complex)+r_horizon_spot;
%     R_tmp1 = sqrt(abs(r_tmp1).^2+obj.SDD^2);
%     R_horizon_spot = sqrt(abs(r_horizon_spot).^2+obj.SDD^2);
%     obj.TwoThetaMap = single((180/pi)*sign(real(r_complex)).*acos((R_tmp1.^2+R_horizon_spot^2-real(r_complex).^2)./(2*R_tmp1*R_horizon_spot)));
%     obj.AlphafMap = single((180/pi)*sign(imag(r_horizon)).*real(acos((R_tmp1.^2+R.^2-imag(r_horizon).^2)./(2*R_tmp1.*R))));
%     obj.QzMap = single(k*(sin(obj.AlphafMap*pi/180)+sin(obj.IncidentAngle*pi/180)));
%     obj.QyMap = single(k*cos(obj.AlphafMap*pi/180).*sin(obj.TwoThetaMap*pi/180));
%     obj.QxMap = -single(k*(cos(obj.IncidentAngle*pi/180)-cos(obj.AlphafMap*pi/180).*cos(obj.TwoThetaMap*pi/180)));
%     obj.QrMap = single(sqrt(obj.QyMap.^2 + obj.QxMap.^2));
%     obj.ChiMap = single(180/pi*angle(obj.QzMap+1i*obj.QrMap));

  
    % 1D list for axis labeling
    obj.Alphaf1DList = 180/pi*(atan(  (obj.Beam0(2)-(1:obj.ImDim(2)))*obj.PixelSize(2)/obj.SDD))-obj.IncidentAngle;
    obj.TwoTheta1DList = 180/pi*atan( ((1:obj.ImDim(1))-obj.Beam0(1))*obj.PixelSize(1)/obj.SDD);
    obj.Qy1DList = k*sin(obj.TwoTheta1DList*pi/180);
    obj.Qz1DList = k*(sin(obj.Alphaf1DList*pi/180)+sin(obj.IncidentAngle*pi/180));
    obj.Qx1DList = -k*(cos(obj.IncidentAngle*pi/180)-cos(obj.Alphaf1DList*pi/180));
end
% solid angle correction
obj.solidangle_correction(obj);
end

function [im,info] = loadimage(obj)
if verLessThan('matlab','7.12')
    msg = 'MATLAB:tifftagsread:unknownTagPayloadState';
elseif verLessThan('matlab','8.5')
    msg = 'MATLAB:imagesci:tifftagsread:unknownTagPayloadState';
else
    msg = 'MATLAB:imagesci:tifftagsread:zeroComponentCount';
end
s = warning('query',msg);
warning('off',msg);
if strcmpi(obj.ImFileExt,'.tif') || strcmpi(obj.ImFileExt,'.tiff')
    im = single(imread(obj.ImFile));
    info = imfinfo(obj.ImFile);
elseif strcmpi(obj.ImFileExt,'.mat')
    im = load(obj.ImFile); tmp = fieldnames(im); im = single(im.(tmp{1}));
    info = [];
elseif strcmpi(obj.ImFileExt,'.cbf')
    tmp = cbfread(obj.ImFile);
    im = single(tmp.data');
    info = tmp.header;
elseif strcmpi(obj.ImFileExt,'.edf')
    %[tmp, info] = edfreader(obj.ImFile);
    [info, tmp] = pmedf_read(obj.ImFile);
    im = single(tmp);
elseif strcmpi(obj.ImFileExt,'.fits')
    im = fitsread(obj.ImFile,'image');
    info = fitsinfo(obj.ImFile);
else
    error('Invalid image type.');
end
warning(s.state,msg);
end

function b = titlestr(a)
if ~ischar(a)
    error('Invalid input argument.');
end
b = a;
c = strfind(a,'_');
for i = length(c):-1:1
    b = [b(1:c(i)-1),'\',b(c(i):end)];
end
end

function axis_mode(haxes,xlabel_str,ylabel_str)
% Only called by POST_IMAGESC
if ~strcmpi(get(get(haxes,'xlabel'),'string'),xlabel_str) || ...
        ~strcmpi(get(get(haxes,'ylabel'),'string'),ylabel_str)
    axis(haxes,'image');
    zoom(get(haxes,'parent'),'reset');
end
xlabel(haxes,xlabel_str);
ylabel(haxes,ylabel_str);
set(haxes,'tickdir','out');
end

function zoompostcallback(~,~,haxes,xdata,ydata)
% For axis label change
set(haxes,'xtickmode','auto','xticklabelmode','auto');
set(haxes,'ytickmode','auto','yticklabelmode','auto');
xlims0 = get(haxes,'xlim');
ylims0 = get(haxes,'ylim');
xlims = interp1(1:length(xdata),xdata,xlims0,'linear','extrap');
ylims = interp1(1:length(ydata),ydata,ylims0,'linear','extrap');
xtick0 = get(haxes,'XTick');
ytick0 = get(haxes,'YTick');
nxt = length(xtick0);
nyt = length(ytick0);
% calculate xtick
xdiff = diff(xlims);
xstep = sign(xdiff)*10^round(log10(abs(xdiff/nxt)));
xtick = ceil(xlims(1)/xstep)*xstep : xstep : floor(xlims(2)/xstep)*xstep;
if length(xtick)<nxt
    xtick = ceil(xlims(1)/xstep)*xstep : xstep/2 : floor(xlims(2)/xstep)*xstep;
end
while length(xtick)>10
    xstep = xstep*2;
    xtick = ceil(xlims(1)/xstep)*xstep : xstep : floor(xlims(2)/xstep)*xstep;
end
xtickpos = interp1(xdata,1:length(xdata),xtick,'linear','extrap');
% calculate ytick
ydiff = diff(ylims);
ystep = sign(ydiff)*10^round(log10(abs(ydiff/nxt)));
ytick = ceil(ylims(1)/ystep)*ystep : ystep : floor(ylims(2)/ystep)*ystep;
if length(ytick)<nyt
    ytick = ceil(ylims(1)/ystep)*ystep : ystep/2 : floor(ylims(2)/ystep)*ystep;
end
while length(ytick)>10
    ystep = ystep*2;
    ytick = ceil(ylims(1)/ystep)*ystep : ystep : floor(ylims(2)/ystep)*ystep;
end
ytickpos = interp1(ydata,1:length(ydata),ytick,'linear','extrap');
% set ticks
set(haxes,'xtickmode','manual','xticklabelmode','manual');
set(haxes,'ytickmode','manual','yticklabelmode','manual');
xtick = cellstr(strjust(num2str(xtick(:)),'left'));
ytick = cellstr(strjust(num2str(ytick(:)),'left'));
set(haxes,'xtick',xtickpos,'xticklabel',xtick);
set(haxes,'ytick',ytickpos,'yticklabel',ytick);
end

function z = bvnd(par,xi,yi)
% Bivariate normal distribution function
A = par(1);
Cx = par(2);
Cy = par(3);
C = par(4);
rho = par(5);
sigma_x = par(6);
sigma_y = par(7);
x0 = par(8);
y0 = par(9);
z = Cx*xi+Cy*yi+C+A*exp(-1/(2*(1-rho)^2)*( (xi-x0).^2/sigma_x^2 + (yi-y0).^2/sigma_y^2 - 2*rho*(xi-x0).*(yi-y0)/(sigma_x*sigma_y)));
end 

function res = res_bvnd(par,xi,yi,zi)
z = bvnd(par,xi,yi);
diff = ((z)-(zi)).^2; diff = diff(:);
diff(isnan(diff)) = [];
res = sum(diff);
end

function res =  res_map(x,map1_roi,map2_roi,xgrid_roi,ygrid_roi,v1,v2)
x0=x(1);  y0=x(2);
v1_tmp=interp2(xgrid_roi,ygrid_roi,map1_roi,x0,y0,'linear');
v2_tmp=interp2(xgrid_roi,ygrid_roi,map2_roi,x0,y0,'linear');
res = sqrt((v1_tmp-v1).^2+(v2_tmp-v2).^2);
end