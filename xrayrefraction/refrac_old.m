function xresult = refrac(varargin)
% ***********************************************
% Copyright (c) 2020 UChicago Argonne, LLC
% See LICENSE file.
% ***********************************************
%
% REFRAC  Material property function interacted with X-rays 
%    r = refrac('Formula String',Energy,MassDensity,'PlotStyle') returns
%    disperion absorption, critical angle and attenuation length
%    correponding to the enery input. Also plot if energy is a list.
%
%    Input: 
%       1) Chemical formula is case sensitive (e.g. CO for Carbon Monoxide
%           vs Co for Cobalt);
%       2) Energy (0.03KeV~30KeV). Can be a single energy or an energy list;
%       3) Mass density (g/cm^3);
%       4) PlotStyle can be 'logxy','linear','logx' or 'logy'. Default
%           is 'linear'.
%    
%    Output structure:
%       1) Chemical formula;
%       2) Molecular weight;
%       3) Number of electrons per molecule;
%       4) Mass density (g/cm^3);
%       5) Electron density (1/A^3);
%       6) X-ray energy (KeV);
%       7) Corresponding X-ray wavelength (A);
%       8) Dispersion;
%       9) Absorption;
%       10) Critical angle (degree);
%       11) Attenuation length (cm);
%       12) Real part of scattering length density (SLD) (A^-2);
%       13) Imaginary part of SLD (A^-2).
%       
%    Plot:
%       Only when energy is a list, plot dispersion, absorption, critical
%       angle and attenuation length.
%
%    Example 1: >> result=refrac('Si3N4',8.04778,3.44)
%           gives an output with structure below:
%               result = 
%                 chemFormula: 'Si3N4'
%             molecularWeight: 140.28
%           numberOfElectrons: 70
%                 massDensity: 3.44
%             electronDensity: 1.0337
%                      energy: 8.0478
%                  wavelength: 1.5406
%                  dispersion: 1.1164e-005
%                  absorption: 1.6477e-007
%               criticalAngle: 0.27074
%                   attLength: 0.0074403
%                       reSLD: 2.9554e-005
%                       imSLD: 4.362e-007
%    
%    Example 2: >>result=refrac('Si3N4',8:0.5:10,3.44)
%           gives an output with structure below together with a plot in
%           linear scale.
%               result =
%                 chemFormula: 'Si3N4'
%             molecularWeight: 140.28
%           numberOfElectrons: 70
%                 massDensity: 3.44
%             electronDensity: 1.0337
%                      energy: [5x1 double]
%                  wavelength: [5x1 double]
%                  dispersion: [5x1 double]
%                  absorption: [5x1 double]
%               criticalAngle: [5x1 double]
%                   attLength: [5x1 double]
%                       reSLD: [5x1 double]
%                       imSLD: [5x1 double]
%           result.energy, result.dispersion, etc. show the lists of this
%           values. >>result=refrac('Si3N4',8:0.5:10,3.44,'logxy') plots in
%           log-log scale.
%
%    For more information about X-ray interactions with matter, go to
%       http://www.cxro.lbl.gov
%       http://www.nist.gov/
%   
%    Atomic scattering factor table is taken from the above two websites.
%
%    Copyright 2004-2005 Zhang Jiang
%    $Revision: 1.0 $  $Date: 2004/10/10 18:52:19 $
%    $Revision: 1.1 $  $Date: 2005/05/25 17:04:11 $

if nargin ~= 3 & nargin ~= 4
    error('Incorrect inputment.');
    return;
end
formulaStr = varargin{1};
energy = varargin{2};
massDensity = varargin{3};
if ~ischar(formulaStr)
    error('Invalid chemical formula.');
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
if ~isnumeric(massDensity) | length(massDensity) > 1
    error('Invalid mass density.');
    return;
end
if nargin == 4
    plotStyle = varargin{4};
    if ~strcmp(plotStyle,'logxy') &...
            ~strcmp(plotStyle,'linear') &...
            ~strcmp(plotStyle,'logx') &...
            ~strcmp(plotStyle,'logy')
        error('Invalid PlotStyle.');
        return;
    end
else
    plotStyle = 'linear';
end

%===============================================================
% some constants
%===============================================================
THOMPSON = 2.81794092e-15;          % m
SPEEDOFLIGHT = 299792458;           % m/sec
PLANCK = 6.626068e-34;              % m^2*kg/sec
ELEMENTCHARGE = 1.60217646e-19;     % Coulombs
AVOGADRO = 6.02214199e23;           % mole^-1

%===============================================================
% sort energy list and convert to wavelength
%===============================================================
energy = sort(energy);  % sort energy from min to max
wavelength = (SPEEDOFLIGHT*PLANCK/ELEMENTCHARGE)./(energy'*1000.0);

%===============================================================
% determine elements and number of atoms in the formula
%===============================================================
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

%===============================================================
% read f1 and f2 from tables
%===============================================================
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

%===============================================================
% determine the atomic number and atomic weight
%===============================================================
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

%===============================================================
% determine molecular weight and number of electrons
%===============================================================
formula.molecularWeight = 0;
formula.numberOfElectrons = 0;
for iElements = 1:nElements
    formula.molecularWeight = formula.molecularWeight + ...
        formula.nAtoms{iElements}*formula.atomicWeight{iElements};
    formula.numberOfElectrons = formula.numberOfElectrons +...
        formula.atomicNumber{iElements}*formula.nAtoms{iElements};
end

%===============================================================
% interpolate to get f1 and f2 for given energies
%===============================================================
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
        energy'*1000,'cubic');
    formula.interpf2{iElements} = interp1(...
        formula.f1f2Table{iElements}(:,1),...
        formula.f1f2Table{iElements}(:,3),...
        energy'*1000,'cubic');
end

%===============================================================
% calculate dispersion, absorption, critical angle and attenuation length
%===============================================================
% initialize
xresult.chemFormula = formulaStr;
xresult.molecularWeight = formula.molecularWeight;
xresult.numberOfElectrons = formula.numberOfElectrons;
xresult.massDensity = massDensity;
xresult.electronDensity = 1e6*massDensity/formula.molecularWeight*AVOGADRO...
    *formula.numberOfElectrons/1e30;
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
end
xresult.criticalAngle = sqrt((2*xresult.dispersion))*180/pi;
xresult.attLength = wavelength./xresult.absorption/(4*pi)*1e2;
xresult.reSLD            = xresult.dispersion*2*pi./wavelength.^2/1e20;
xresult.imSLD            = xresult.absorption*2*pi./wavelength.^2/1e20;

%===============================================================
% plot
%===============================================================
% plot only when energy is a list
if length(energy) == 1
    return;
end
switch plotStyle
    case 'logxy'
        xscale = 'log';
        yscale = 'log';
    case 'logx'
        xscale = 'log';
        yscale = 'linear';
    case 'logy'
        xscale = 'linear';
        yscale = 'log';
    case 'linear'
        xscale = 'linear';
        yscale = 'linear';
end
xresultFig = figure(...
    'Name',['X-Ray Interactions with ',formulaStr],...
    'NumberTitle','off',...
    'PaperOrientation','landscape');
subplot(2,2,1);
% plot(xresult.energy,xresult.dispersion,'b-o',...
%     'MarkerEdge','r','MarkerSize',0.1);
plot(xresult.energy,xresult.dispersion);
xlabel('Energy (KeV)');
ylabel('Dispersion');
set(gca,'XScale',xscale,'YScale',yscale);
grid on;
box on;
title('Dispersion');
subplot(2,2,2);
% plot(xresult.energy,xresult.absorption,'b-o',...
%     'MarkerEdge','r','MarkerSize',0.1);
plot(xresult.energy,xresult.absorption);
xlabel('Energy (KeV)');
ylabel('Absorption');
set(gca,'XScale',xscale,'YScale',yscale);
grid on;
box on;
title('Absorption');
subplot(2,2,3);
% plot(xresult.energy,xresult.criticalAngle,'b-o',...
%     'MarkerEdge','r','MarkerSize',0.1);
plot(xresult.energy,xresult.criticalAngle);
xlabel('Energy (KeV)');
ylabel('Critical Angle (degree)');
set(gca,'XScale',xscale,'YScale',yscale);
grid on;
box on;
title('Critical Angle');
subplot(2,2,4);
% plot(xresult.energy,xresult.attLength,'b-o',...
%     'MarkerEdge','r','MarkerSize',0.1);
plot(xresult.energy,xresult.attLength);
xlabel('Energy (KeV)');
ylabel('Attenuation Length (cm)');
set(gca,'XScale',xscale,'YScale',yscale);
grid on;
box on;
title('Attenuation Length');
supertitle(['X-Ray Interactions with ',formulaStr,...
    ', \rho_{mass}=',num2str(xresult.massDensity),'g/cm^3']);

return;