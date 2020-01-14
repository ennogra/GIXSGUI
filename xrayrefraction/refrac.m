function result = refrac(varargin)
% ***********************************************
% Copyright (c) 2020 UChicago Argonne, LLC
% See LICENSE file.
% ***********************************************
%
% REFRAC  Material property when interacted with x-rays.
%
%    Usage: result = refrac(Chemical Formula,Energy,MassDensity) 
%
%    Input: 
%       1) Chemical formula. Can be either cell array (for multiple 
%          formulas) or single string. Formula is case sensitive (e.g. CO 
%          for Carbon Monoxide vs Co for Cobalt);
%       2) Energy (0.03KeV~30KeV). Can be single or list.
%       3) List of mass densities (g/cm^3);
%    
%    Output: Can be either cell of structures (if input formular is cell 
%       array) or single structure (if input is single string) with the
%       following elements:
%       1) Chemical formula;
%       2) Molecular weight;
%       3) Number of electrons per molecule;
%       4) Mass density (g/cm^3);
%       5) Electron density (1/A^3);
%       6) Real part of the atomic scattering factor (f1)
%       7) Imaginary part of the atomic scattering factor (f2);
%       8) X-ray energy (KeV);
%       9) Corresponding X-ray wavelength (A);
%       10) Dispersion;
%       11) Absorption;
%       12) Critical angle (degree);
%       13) Attenuation length (cm);
%       14) Real part of scattering length density (SLD) (A^-2);
%       15) Imaginary part of SLD (A^-2).
%
%    Example 1: >> result=refrac({'H2O','Si3N4'},8.04778,[1,3.44])
%               Output: cell of structures 
%               result{2} = 
%                 chemFormula: 'Si3N4'
%             molecularWeight: 140.28
%           numberOfElectrons: 70
%                 massDensity: 3.44
%             electronDensity: 1.0337
%                          f1: 71.0209
%                          f2: 1.0482
%                      energy: 8.0478
%                  wavelength: 1.5406
%                  dispersion: 1.1164e-005
%                  absorption: 1.6477e-007
%               criticalAngle: 0.27074
%                   attLength: 0.0074403
%                       reSLD: 2.9554e-005
%                       imSLD: 4.362e-007
%    
%    Example 2: >>result=refrac('SiO2',8:0.5:10,2.33)
%               Output: structure 
%               result =
%                 chemFormula: 'SiO2'
%             molecularWeight: 60.0843
%           numberOfElectrons: 30
%                 massDensity: 2.33
%             electronDensity: 0.7006
%                          f1: [5x1 double]
%                          f2: [5x1 double]
%                      energy: [5x1 double]
%                  wavelength: [5x1 double]
%                  dispersion: [5x1 double]
%                  absorption: [5x1 double]
%               criticalAngle: [5x1 double]
%                   attLength: [5x1 double]
%                       reSLD: [5x1 double]
%                       imSLD: [5x1 double]
%
%    For more information about X-ray interactions with matter, go to
%       http://www.cxro.lbl.gov
%       http://www.nist.gov/
%   
%    Atomic scattering factor table is taken from the above two websites.
%
%    Zhang Jiang
%    $Revision: 1.0 $  $Date: 2004/10/10 $
%    $Revision: 1.1 $  $Date: 2005/05/25 $
%    $Revision: 1.2 $  $Date: 2008/08/28 $
%    $Revision: 1.3 $  $Date: 2008/10/21 $ Input chemical formula can be
%    either cell array or single char string.
%    $Revision: 1.4 $  $Date: 2015/07/29 $ Calculate summed f1 and f2.

if nargin ~= 3 
    error('Invalid number of input arguements.');
    return;
end
formulaStrCell = varargin{1};
energy = varargin{2};
massDensityList = varargin{3};
if ischar(formulaStrCell)
    formulaStrCell = {formulaStrCell};
elseif ~iscell(formulaStrCell);
    error('Need cell for chemical formula input argument.');
    return;
end
if ~isnumeric(energy) | isempty(energy)
    error('Invalid x-ray energy.');
    return;
end
if min(energy)<0.03 | max(energy)>30
    error('Energy is out of range 0.03KeV~30KeV.');
    return;
end
if ~isnumeric(massDensityList) 
    error('Invalid mass density.');
    return;
end
if length(formulaStrCell) ~= length(massDensityList)
    error('Input arguements do not match.');
    return;
end

% --- sort energy list and
energy = sort(energy);  % sort energy from min to max

for ii = 1:length(formulaStrCell)
    result{ii} = subrefrac(formulaStrCell{ii},energy,massDensityList(ii));
end
if ischar(varargin{1})
    result = result{ii};
end
return;


%===============================================================
% --- Sub function
%===============================================================
function xresult = subrefrac(formulaStr,energy,massDensity)
if ~ischar(formulaStr)
    error('Invalid chemical formula.');
    return;
end

% --------------------------------------------------------------
% some constants
% --------------------------------------------------------------
THOMPSON = 2.8179403227e-15;          % m
SPEEDOFLIGHT = 299792458;           % m/sec
PLANCK = 6.626068e-34;              % m^2*kg/sec
ELEMENTCHARGE = 1.60217646e-19;     % Coulombs
AVOGADRO = 6.02214199e23;           % mole^-1

% --------------------------------------------------------------
% convert energy to wavelength
% --------------------------------------------------------------
wavelength = (SPEEDOFLIGHT*PLANCK/ELEMENTCHARGE)./(energy'*1000.0);

% --------------------------------------------------------------
% determine elements and number of atoms in the formula
% --------------------------------------------------------------
nElements = 0;
formula.elements = {};
formula.nAtoms = {};
try
    for iFormulaStr = 1:length(formulaStr)
        formulaChar = formulaStr(iFormulaStr);
        if formulaChar <= 'Z' &  formulaChar >= 'A'
            nElements = nElements + 1;
            formula.elements{nElements} = formulaChar;
            formula.nAtoms{nElements} = '0';
        elseif formulaChar <= 'z' &  formulaChar >= 'a'...
                & ((formulaStr(iFormulaStr-1) <= 'Z'& formulaStr(iFormulaStr-1) >= 'A')...
                | (formulaStr(iFormulaStr-1) <= 'z'& formulaStr(iFormulaStr-1) >= 'a'))
            formula.elements{nElements} = ...
                [formula.elements{nElements},formulaChar];
        elseif (formulaChar <= '9' &  formulaChar >= '0') | formulaChar == '.'
            formula.nAtoms{nElements} = ...
                [formula.nAtoms{nElements},formulaChar];
        else
            error('Invalid chemical formula.');
            return;
        end
    end
catch
    error('Invalid chemical formula.');
    return;
end
for iElements = 1:nElements
    formula.nAtoms{iElements} = str2num(formula.nAtoms{iElements});
    if formula.nAtoms{iElements} == 0
        formula.nAtoms{iElements} = 1;
    end
end

% --------------------------------------------------------------
% read f1 and f2 from tables
% --------------------------------------------------------------
formula.f1f2Table = {};
for iElements = 1:nElements
%    file = fullfile(matlabroot,'toolbox','xraylabtool','xrayrefraction',...
%        'AtomicScatteringFactor',lower(formula.elements{iElements}));
    file = lower(formula.elements{iElements});
    file = [file,'.nff'];
    try 
        fid = fopen(file);
        fgetl(fid);
    catch
        error(['Element ''',formula.elements{iElements},...
            ''' is NOT in the table list.']);
        return;
    end
    table = cell2mat(textscan(fid,'%f %f %f'));
    table(find(table(:,1)<29),:) = [];
    formula.f1f2Table{iElements} = table;
    fclose(fid);
end

% --------------------------------------------------------------
% determine the atomic number and atomic weight
% --------------------------------------------------------------
formula.atomicNumber = {};
formula.atomicWeight = {};
%file = fullfile(matlabroot,'toolbox','xraylabtool','xrayrefraction','periodictable.mat');
file = 'periodictable.mat';
load(file);
for iElements = 1:nElements
    for iAtomicnum =1:length(elementAbbr)
        if strcmp(elementAbbr{iAtomicnum},formula.elements{iElements})
            formula.atomicNumber{iElements} = iAtomicnum;
            formula.atomicWeight{iElements} = atomicWeight(iAtomicnum);
            break;
        end
    end
end

% --------------------------------------------------------------
% determine molecular weight and number of electrons
% --------------------------------------------------------------
formula.molecularWeight = 0;
formula.numberOfElectrons = 0;
for iElements = 1:nElements
    formula.molecularWeight = formula.molecularWeight + ...
        formula.nAtoms{iElements}*formula.atomicWeight{iElements};
    formula.numberOfElectrons = formula.numberOfElectrons +...
        formula.atomicNumber{iElements}*formula.nAtoms{iElements};
end

% --------------------------------------------------------------
% interpolate to get f1 and f2 for given energies
% --------------------------------------------------------------
formula.interpf1 = {};
formula.interpf2 = {};
for iElements = 1:nElements
%     formula.interpf1{iElements} = exp(interp1(...
%         formula.f1f2Table{iElements}(:,1),...
%         log(formula.f1f2Table{iElements}(:,2)),...
%         energy'*1000,'cubic'));
%     formula.interpf2{iElements} = exp(interp1(...
%         formula.f1f2Table{iElements}(:,1),...
%         log(formula.f1f2Table{iElements}(:,3)),...
%         energy'*1000,'cubic'));
    formula.interpf1{iElements} = interp1(...
        formula.f1f2Table{iElements}(:,1),...
        formula.f1f2Table{iElements}(:,2),...
        energy'*1000,'pchip');
    formula.interpf2{iElements} = interp1(...
        formula.f1f2Table{iElements}(:,1),...
        formula.f1f2Table{iElements}(:,3),...
        energy'*1000,'pchip');
end

% --------------------------------------------------------------
% calculate dispersion, absorption, critical angle and attenuation length
% --------------------------------------------------------------
% initialize
xresult.chemFormula = formulaStr;
xresult.molecularWeight = formula.molecularWeight;
xresult.numberOfElectrons = formula.numberOfElectrons;
xresult.massDensity = massDensity;
xresult.electronDensity = 1e6*massDensity/formula.molecularWeight*AVOGADRO...
    *formula.numberOfElectrons/1e30;
xresult.f1 = zeros(length(energy),1);
xresult.f2 = zeros(length(energy),1);
xresult.energy      = energy';
xresult.wavelength  = wavelength*1e10;
xresult.dispersion  = zeros(length(energy),1);
xresult.absorption  = zeros(length(energy),1);
xresult.criticalAngle  = zeros(length(energy),1);
for iElements = 1:nElements
    xresult.dispersion = xresult.dispersion + ...
        wavelength.^2/(2*pi)*THOMPSON*AVOGADRO*massDensity*1e6 / ...
        formula.molecularWeight * ...
        formula.nAtoms{iElements} .* ...
        formula.interpf1{iElements};
    xresult.absorption = xresult.absorption + ...
        wavelength.^2/(2*pi)*THOMPSON*AVOGADRO*massDensity*1e6 / ...
        formula.molecularWeight * ...
        formula.nAtoms{iElements} .* ...
        formula.interpf2{iElements};
    xresult.f1 = xresult.f1 + formula.nAtoms{iElements} .* ...
        formula.interpf1{iElements};
    xresult.f2 = xresult.f2 + formula.nAtoms{iElements} .* ...
        formula.interpf2{iElements};
end
xresult.criticalAngle = sqrt((2*xresult.dispersion))*180/pi;
xresult.attLength = wavelength./xresult.absorption/(4*pi)*1e2;
xresult.reSLD            = xresult.dispersion*2*pi./wavelength.^2/1e20;
xresult.imSLD            = xresult.absorption*2*pi./wavelength.^2/1e20;
return;