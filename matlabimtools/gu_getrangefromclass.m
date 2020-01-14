function range = gu_getrangefromclass(I)
%GETRANGEFROMCLASS Get dynamic range of image based on its class.
%   RANGE = GETRANGEFROMCLASS(I) returns the dynamic range of the 
%   image I, based on its class type. 
%
%   Class Support
%   -------------
%   I can be numeric or logical. RANGE is a two-element vector of doubles.
%
%   Note
%   ----
%   For single and double data, GETRANGEFROMCLASS returns the range [0 1],
%   to be consistent with the way double and single images are interpreted
%   in MATLAB.  For integer data, GETRANGEFROMCLASS returns the range of
%   the class. For example, if the class is uint8, the dynamic range is
%   [0 255].
%
%   Example
%   -------
%       % Get the dynamic range of an int16 image.
%       CT = dicomread('CT-MONO2-16-ankle.dcm');
%       r = getrangefromclass(CT)
%
%   See also INTMIN, INTMAX.
  
%   Copyright 1993-2006 The MathWorks, Inc.
%   $Revision: 1.1.8.2 $  $Date: 2006/05/24 03:33:22 $
  
gu_iptchecknargin(1,1,nargin,mfilename);
% modified by ZJ
% iptchecknargin(1,1,nargin,mfilename);
% iptcheckinput(I,{'numeric','logical'}, {}, mfilename,'I',1);

if isinteger(I)
    classType = class(I);
    range = double([intmin(classType) intmax(classType)]);
else
    %double, single, or logical
    range = [0 1];
end 

