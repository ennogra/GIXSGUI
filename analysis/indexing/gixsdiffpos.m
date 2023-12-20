function result = gixsdiffpos(lattice,sguvw,orientation,orientationmethod,hlist,klist,llist,k0,alpha_i,nfilm,qdeadband,qcutoff)
% ***********************************************
% Copyright (c) 2020 UChicago Argonne, LLC
% See LICENSE file.
% ***********************************************
%
% GIXSDIFFPOS Calculate the grazing-incidence x-ray diffraction positions
%   of polycrystalline films.
%   RESULT=GIXSDIFFPOS(LATTICE,SGUVW,ORIENTATION,ORIENTATIONMETHOD,HLIST,
%   KLIST,LLIST,K0,ALPHA_I,NFILM,QDEADBAND,QCUTOFF) calculates the GIXS
%   diffraction positions from a nanostructred film. Polycrystals in the film
%   have no preferred orientation in the plane of the substrate. Both Born
%   approximation (BA) and distorted wave Born approximation (DWBA) are 
%   used to calculate the angle and q positions of allowed diffractions. 
%   
%   LATTICE = [a,b,c,alpha,beta,gamma] defines the primitive lattice cell.
%   a,b,c are in A and angles are in degree. 
%
%   SGUVW defines the symmetry of the lattice. It is either an integer 
%   (1~230) for the space group, or an array (nx3) to define the lattice 
%   basis with n lattice "atoms" positioned at fractional coordinates [uvw] 
%   in the unit cell.
%
%   ORIENTATION (1x3) and ORIENTATIONMETHOD define the unit cell direction 
%   that is normal to the substrate surface. ORIENTATIONMETHOD=1 uses the 
%   Cartesian sample reference frame to define this direction. The initial
%   orientation of the unit cell is layer out with respect to the substrate
%   surface and the forward direction of the incident beam such that the 
%   a-b plane is on the substate plane and the lattice vector a and b are 
%   along and perpendicular to the forward direction of the incident beam, 
%   respectively. After the ORIENTATION is defined by an $[x,y,z] direction
%   in the Cartesian sample reference frame, the lattice will be 
%   automatically rotated so that the defined direction is along the 
%   substrate normal, and then the orientation of polycrystals is averaged 
%   over all in-plane angles for completely random in-plane orientations. 
%   A more convenient method to define the unit cell orientation is to do 
%   it in the lattice reference frame, i.e., ORIENTATIONMETHOD=2, where 
%   ORIENTATION is specified by [uvw], the indices of the crystallographic 
%   direction of a lattice.
%
%   HLIST, KLIST, LLIST are list of h,k,l values.
%
%   K0 is the wave vector of the incident beam.
%
%   ALPHG_I is the incident angle in degree.
%
%   NFILM is the refractive index of the film.
%
%   QDEADBAND specifies the q range within with the diffractions are
%   considered a multiplications, i.e., one diffraction corresponds to 
%   multiple (hkl) values. Typcal values are 1e-8 A^(-1).
%
%   QCUTOFF specifies the q cutoff below which the solved absolute q
%   (qx,qy,qz) values are set to zero. This is useful for large lattices
%   that may give various small qx, qy or qz values due to the precision of
%   equation solving algorithm. Typical values is 1e-6 A^-1.  
%
%   The output is a structure with fields of self-explanatory names.

%   Zhang Jiang @8ID/APS/ANL
%   $Revision: 1.0 $  $Date: 2012/08/08 $
%   $Revision: 1.1 $  $Date: 2017/01/10 $ Fix gixsdiffpos.m for correct
%       exit angles. When qz (BA) is less than minimum that can be achieve
%       for positive ext angles, the calculated exit angle should be
%       negative. Also convert the wave vector from vac to film in when
%       calcuation the q values for diffractions. Update (2018/02/07): (1)
%       Fix the typo in the Snell's law (Thansk to Ilja Gunkel). (2)
%       Include angles for BA in the output.

% use large number instead of inf for 2D lattices
lattice(isinf(lattice)) = 1/eps(eps); 

% lengths of hkl list
nh = length(hlist);
nk = length(klist);
nl = length(llist);

% convert k0 to k in the film
k_film = k0*real(nfilm);

% --- start calculate q in BA
q_ba = [];          % to store q in BA
miller_ba = [];     % to store miller
for ih=1:nh
    for ik=1:nk
        for il=1:nl
            miller = [hlist(ih),klist(ik),llist(il)];  % miller index
            [q,~]=diffposqba(...
                lattice,sguvw,orientation,orientationmethod,miller,k_film,alpha_i,qcutoff);
            if isempty(q), continue; end
            q_ba = [q_ba;q];
            miller_ba = [miller_ba;miller];
        end
    end
end

% --- return empty result f no diffration is found
if isempty(miller_ba)   
    result = [];
    return;
end

% --- find identical/degenerated diffractions
qq = [];      
mm = {};
while ~isempty(q_ba)
    qq = [qq;q_ba(1,:)];
    idx = rangesearch(q_ba(1,:),qdeadband,q_ba);
    if isempty(idx), idx = 1; end
    mm = [mm;miller_ba(idx,:)];
    q_ba(idx,:) = [];
    miller_ba(idx,:) = [];
end
[~,idx] = sort(sqrt(sum(qq.^2,2)));
q_ba = qq(idx,:);           % degenerated q (nx3) in BA
miller_ba = mm(idx);        % (nx1) cell with each element (mx3) stores m degenerated (hkl) values

% --- calculate angles and q in DWBA
qr_ba = sqrt(q_ba(:,1).^2+q_ba(:,2).^2);
af_ba = asin((q_ba(:,3)/k_film) - sind(alpha_i));       % BA angle ignoring the reflection effect % 2023/12/15
% af_ba = asin((q_ba(:,3)/k_film) - 0);       % BA angle ignoring the reflection effect % 2023/12/15

af1 = asin( real( sqrt( (q_ba(:,3)/k_film).^2+sin(alpha_i*pi/180)^2 - 2*q_ba(:,3)/k_film*sqrt(sin(alpha_i*pi/180)^2-1+nfilm^2) ) ) );
af2 = asin( real( sqrt( (q_ba(:,3)/k_film).^2+sin(alpha_i*pi/180)^2 + 2*q_ba(:,3)/k_film*sqrt(sin(alpha_i*pi/180)^2-1+nfilm^2))) );
% % check alpha_i and qz values in the film to correct the sign of af1 only
% $Date: 2017/01/10 $ & 2023/12/15 enable
alpha_i_film = acosd(cosd(alpha_i)/real(nfilm));         % snell's law
qz_film_min = k_film*nfilm*sind(alpha_i_film);          % minimum qz that can be reached in film, i.e. when alpha_f_film = 0
ind_qz_film_min = q_ba(:,3)<real(qz_film_min);          % here we ignore the GITSAXS peaks
af1(ind_qz_film_min) = -af1(ind_qz_film_min);

t_ba = acos(real((cos(af_ba).^2+cos(alpha_i*pi/180)^2-(qr_ba/k_film).^2)./(2*cos(af_ba)*cos(alpha_i*pi/180))));
t1 =   acos(real((cos(af1).^2  +cos(alpha_i*pi/180)^2-(qr_ba/k_film).^2)./(2*cos(af1)  *cos(alpha_i*pi/180))));
t2 =   acos(real((cos(af2).^2  +cos(alpha_i*pi/180)^2-(qr_ba/k_film).^2)./(2*cos(af2)  *cos(alpha_i*pi/180))));
af_ba_real = real(af_ba)*180/pi;
af1_real = real(af1)*180/pi;
af2_real = real(af2)*180/pi;
t_ba_real = real(t_ba)*180/pi;
t1_real = real(t1)*180/pi;
t2_real = real(t2)*180/pi;
% q(qx,qy,qz) in dwba for transmission
q_dwba_t = ones(size(q_ba));
q_dwba_t(:,1) = k0*cos(af1_real*pi/180).*cos(t1_real*pi/180) - k0*cos(alpha_i*pi/180);
q_dwba_t(:,2) = k0*cos(af1_real*pi/180).*sin(t1_real*pi/180);
q_dwba_t(:,3) = k0*(sin(af1_real*pi/180)+sin(alpha_i*pi/180));
% q(qx,qy,qz) in dwba for reflection
q_dwba_r = ones(size(q_ba));
q_dwba_r(:,1) = k0*cos(af2_real*pi/180).*cos(t2_real*pi/180) - k0*cos(alpha_i*pi/180);
q_dwba_r(:,2) = k0*cos(af2_real*pi/180).*sin(t2_real*pi/180);
q_dwba_r(:,3) = k0*(sin(af2_real*pi/180)+sin(alpha_i*pi/180));

% --- remove negative qz (in BA) values for GIXS; but keep those in af_ba
% and t_ba; 
index_q_ba = (q_ba(:,3)<0); 
%index_q_ba = [];
t1_real(index_q_ba) = [];
t2_real(index_q_ba) = [];
af1_real(index_q_ba) = [];
af2_real(index_q_ba) = [];

q_dwba_t(index_q_ba,:) = [];
q_dwba_r(index_q_ba,:) = [];
% miller_dwba = miller_ba(~index_q_ba);
miller_dwba = miller_ba;
miller_dwba(index_q_ba) = [];

% --- construct result structure

result.miller_ba_full = miller_ba;
result.angle_ba_full = [t_ba_real,af_ba_real];
result.q_ba_full = q_ba;
result.miller = miller_dwba;
% result.angle_ba = result.angle_ba_full(~index_q_ba,:);
result.angle_ba = result.angle_ba_full;
result.angle_ba(index_q_ba,:) = [];
result.angle_t = [t1_real,af1_real];
result.angle_r = [t2_real,af2_real];
result.q_dwba_t = q_dwba_t;
result.q_dwba_r = q_dwba_r;
% result.q_ba = q_ba(~index_q_ba,:);  
result.q_ba = q_ba;
result.q_ba(index_q_ba,:) = [];

%result.zeta = angle(result.q_ba(:,1)+1i*result.q_ba(:,2))*180/pi;    % angle in the substrate system that contribute to scattering 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- function to calculate q in ba. It returns q=[qx,qy,qz] values in
% substrate substrate system in BA approximation.)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [q,exitflag] = diffposqba(lattice,sguvw,orientation,orientationmethod,miller,k0,alpha_i,qcutoff)
% --- assign lattice parameters
a = lattice(1);        
b = lattice(2);
c = lattice(3);
alpha   = lattice(4);
beta    = lattice(5);
gamma   = lattice(6);

% --- direct structure matrix to convert vector from lattice system to
% substrate system
A=crystaldirectstructruematrix(a,b,c,alpha,beta,gamma);

% --- unit lattice vectors in lattice system
uvabc_lattice = [...
    1 0 0;
    0 1 0;
    0 0 1];     % 1/2/3 columns for a/b/c axes.

% --- lattice vectors in substrate system
vabc_substrate = A*uvabc_lattice;

% --- reflection conditions are determined either via space group or the
% lattice coordinates <uvw>.
if numel(sguvw) == 1 && sguvw>=1 && sguvw<=230 && floor(sguvw)==sguvw  
    issg = 1;
    sg = sguvw;
elseif size(sguvw,2) == 3     % use lattice points (nx3)
    issg = 0;
    uvw = sguvw;
else
    error('Invalude sg number or lattice point positions');
end

% --- unit axis of the substrate system. ex and ey are on the substrate
% plane, with ex along and qy normal to the beam direction. ez is normal
% to the substrate plane.
exyz_substrate =  [...
    1 0 0;
    0 1 0;
    0 0 1];     % 1/2/3 columns for ex/ey/ez axes.
ez_substrate = exyz_substrate(:,3);

% --- direction perpendicular to the substrate. There are 3 choices:
if numel(orientation)~=3 || nnz(orientation) == 0
    error('Invalid orientation');
end

switch orientationmethod
    case 1  % 1. substrate referecnce system [x'y'z']. 
        axis_norm = orientation(:);
    case 2  % 2. lattice coordinate <uvw> on lattice system. 
        axis_norm = orientation(:); % normal direction in lattice system
        axis_norm = A*axis_norm;    % normal direction in substrate
    case 3  % 3. plane index (hkl) parrallel to substrate. [hkl] direction 
        % may not perpendicular tot he plane (hkl).
        error('Under construction');
    otherwise
        error('Invalid lattice orientation');
end
axis_norm = axis_norm/sqrt(dot(axis_norm,axis_norm));

% --- rotate lattice vectors
u_rot = cross(axis_norm,ez_substrate);    % roation axis
if nnz(u_rot) ~= 0 % rotate only when aixs_norm is not along ez.
    u_rot = u_rot/sqrt(dot(u_rot,u_rot));  % convert to unit axis
    ux = u_rot(1);
    uy = u_rot(2);
    uz = u_rot(3);
    % calculate angles to be rotated wrt to u_rot (by righthand rule)
    theta = acos(dot(ez_substrate,axis_norm));
    % rotation matrix
    R = [...
        cos(theta)+ux^2*(1-cos(theta))      ux*uy*(1-cos(theta))-uz*sin(theta)  ux*uz*(1-cos(theta))+uy*sin(theta);...
        uy*ux*(1-cos(theta))+uz*sin(theta)  cos(theta)+uy^2*(1-cos(theta))      uy*uz*(1-cos(theta))-ux*sin(theta);...
        uz*ux*(1-cos(theta))-uy*sin(theta)  uz*uy*(1-cos(theta))+ux*sin(theta)  cos(theta)+uz^2*(1-cos(theta))];
    % new lattice vectors (in substrate system) after lattice rotation
    % 1/2/3 columns for a/b/c axes in substrate system
      vabc_substrate = R*vabc_substrate;    
end

% --- get miller index
miller = miller(:);

% --- check reflection condition
if issg     % use space group
    if sgreflection(miller,sg) == 0
        q = [];
        exitflag = [];
        return;
    end
else        % use manually built lattice
    SF = sum(exp(-2*pi*1i*uvw*miller));
    if abs(real(SF))<0.01
        q = [];
        exitflag = [];
        return;
    end
end

% --- calculate G vector in the substrate system
G_substrate = gvector(vabc_substrate,miller);


% --- calculate the q positions
options=optimset(...
    'Display','off',...
    'TolX',1e-8,...
    'TolFun',1e-8,...
    'MaxIter',1000,...
    'MaxFunEvals',2000); %,...
% q0 = rand(1,6);     % initial value
% [q,~] = fsolve(@ewald,q0,options,k0,alpha_i,G_substrate);  % Call solver
% q=q(:,1:3);           
exitflag = -2;
while exitflag <0
    q0 = rand(1,4);     % initial value    
    [q_tmp,~,exitflag] = fsolve(@ewald,q0,options,k0,alpha_i,G_substrate);  % Call solver
end

% --- construct result and apply qcutoff
q(1) = q_tmp(1);
q(2) = q_tmp(2);
q(3) = -sin(alpha_i*pi/180)*q_tmp(3)+cos(alpha_i*pi/180)*q_tmp(4);
% q(4) = q_tmp(3);    
% q(5) = q_tmp(2);
% q(6) = (q_tmp(1)-cos(alpha_i*pi/180)*q_tmp(3))/sin(alpha_i*pi/180);
q(abs(q)<=eps) = 0;         % set near-zero values to zero
q(abs(q)<=qcutoff) = 0;     % apply q tolerance


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- function to solve for q in ba 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F = ewald(x,k,alpha,G)

% --- incident angle
alpha = alpha*pi/180;

% --- in-plane and outplane values of G
SR_substrate = sqrt(G(1)^2+G(2)^2);
SZ_substrate = G(3);

% --- 6 equation method
% x(1), x(2), x(3) are qx, qy, qz in substrate frame.
% x(4), x(5), x(6) are qx, qy, qz in the lab frame.
% F = [...
%     x(1).^2 + x(2).^2 - SR_substrate.^2;...
%     x(3) - SZ_substrate;...
%     x(1) - cos(alpha)*x(4) - sin(alpha)*x(6);...
%     x(2) - x(5);...
%     x(3) + sin(alpha)*x(4) - cos(alpha)*x(6);...
%     (x(4) - k).^2 + x(5).^2 + x(6).^2 - k^2];

% --- 4 equation method (faster than 6 equation method)
% x(1), x(2) are for qx, qy in substrate frame
% x(3), x(4) are for qx, qz in the lab frame
F = [...
    x(1).^2 + x(2).^2 - SR_substrate.^2;...
    x(1) - cos(alpha)*x(3) - sin(alpha)*x(4);...
    SZ_substrate + sin(alpha)*x(3) - cos(alpha)*x(4);...
    (x(3) - k).^2 + x(2).^2 + x(4).^2 - k^2];