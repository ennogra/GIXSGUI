function trueReflec = subtractreflectivitybackground(varargin)
% ***********************************************
% Copyright (c) 2020 UChicago Argonne, LLC
% See LICENSE file.
% ***********************************************
%
% SUBTRACTREFLECTIVITYBACKGROUND Subtract longitudial diffuse scattering
%   background from reflectivity.
%   Y = SUBTRACTREFLECTIVITYBACKGROUND(REF,DIFF1) if diffuse scattering is
%       measured only for one incident angle offset.
%   Y = SUBTRACTREFLECTIVITYBACKGROUND(REF,DIFF1,DIFF2) if diffuse
%       scattering is measured for both negative and positive offsets.
%
%   Input argument:
%       REF: Mx2 or Mx3 matrix with the 1st col for tth (twice of the
%           incident angle) or qz, 2nd column the intensity, and 3rd column
%           (optional) the absolute standard deviation.
%       DIFF1 and DIFF2: diffuse scattering matrix of the the same format
%           as REF.
%   
%   Output argument:
%       Y: background subtracted reflectivity (Mx3). The 3d col is zero if
%           standard deviation is not specified for the input.

%   Zhang Jiang @8ID/APS/ANL
%   $Revision: 1.0 $  $Date: 2015/12/19 $

if nargin == 2
    reflec = varargin{1};    
    posDiff = varargin{2};
    negDiff = [];
elseif nargin == 3
    reflec = varargin{1};        
    posDiff = varargin{2};    
    negDiff = varargin{3};
else
    error('Invalid input argument ...');
end

% --- if the 3rd column (error) does not exist, add zeros
if size(reflec,2) == 2
    reflec(:,3) = 0;
end
if size(posDiff,2) == 2
    posDiff(:,3) = 0;
end
if size(negDiff,2) == 2 && ~isempty(negDiff)
    negDiff(:,3) = 0;
end

% --- convert error to relative error
reflec(:,3) = reflec(:,3)./reflec(:,2);
posDiff(:,3) = posDiff(:,3)./posDiff(:,2);
if ~isempty(negDiff)
    negDiff(:,3) = negDiff(:,3)./negDiff(:,2);
end

% --- subtraction 
if isempty(negDiff)         % if only one diffuse exists 
    diffuse = posDiff;
    if reflec(1,1) < diffuse(1,1)
        diffuse  = [reflec(1,1),diffuse(1,2:3);diffuse];
    end
    if  reflec(end,1) > diffuse(end,1)
        diffuse = [diffuse;reflec(end,1),diffuse(end,2:3)];
    end
    interpDiff(:,1) = reflec(:,1);        % get angle of intercepted diffuse scattering
    %    interpDiff(:,2:3) = interp1(diffuse(:,1),diffuse(:,2:3),reflec(:,1),'linear');
    interpDiff(:,2:3) = interp1(diffuse(:,1),diffuse(:,2:3),reflec(:,1),'linear','extrap');
    trueReflec(:,1) = reflec(:,1);                            % true reflectivity angle
    trueReflec(:,2) = reflec(:,2)-interpDiff(:,2);            % true reflectivity intensity
    absReflecError = reflec(:,3).*reflec(:,2);                % absolute reflectivity error
    absInterpDiffError = interpDiff(:,3).*interpDiff(:,2);    % absolute interpolated diffuse error
    trueReflec(:,3) = sqrt(absReflecError.^2+absInterpDiffError.^2)./trueReflec(:,2); % releative true reflectivity error
else    % if both diffuse exist, calculate the average of the diffuse
    if reflec(1,1) < posDiff(1,1)
        posDiff  = [reflec(1,1),posDiff(1,2:3);posDiff];
    end
    if reflec(1,1) < negDiff(1,1)
        negDiff  = [reflec(1,1),negDiff(1,2:3);negDiff];
    end
    if  reflec(end,1) > posDiff(end,1)
        posDiff = [posDiff;reflec(end,1),posDiff(end,2:3)];
    end
    if  reflec(end,1) > negDiff(end,1)
        negDiff = [negDiff;reflec(end,1),negDiff(end,2:3)];
    end
    interpPosDiff(:,1) = reflec(:,1);       % set angle of intercepted positive diffuse scattering
    interpNegDiff(:,1) = reflec(:,1);       % set angle of intercepted negative diffuse scattering
    interpPosDiff(:,2:3) = interp1(posDiff(:,1),posDiff(:,2:3),reflec(:,1),'linear','extrap');    
    interpNegDiff(:,2:3) = interp1(negDiff(:,1),negDiff(:,2:3),reflec(:,1),'linear','extrap');            
    trueReflec(:,1) = reflec(:,1);                                                    % true reflectivity angle
    trueReflec(:,2) = reflec(:,2)-0.5*interpPosDiff(:,2)-0.5*interpNegDiff(:,2);      % true reflectivity intensity
    absReflecError = reflec(:,3).*reflec(:,2);                                        % absolute reflectivity error
    absInterpPosDiffError = interpPosDiff(:,3).*interpPosDiff(:,2);                   % absolute interpolated positive diffuse error    
    absInterpNegDiffError = interpNegDiff(:,3).*interpNegDiff(:,2);                   % absolute interpolated negtive diffuse error    
    trueReflec(:,3) = sqrt(absReflecError.^2+0.25*absInterpPosDiffError.^2+0.25*absInterpNegDiffError.^2)./trueReflec(:,2);
end

% --- remove negative data
trueReflec(trueReflec(:,2) < 0,:)=[];

% --- convert back to absolute error
trueReflec(:,3) = trueReflec(:,3).*trueReflec(:,2);
