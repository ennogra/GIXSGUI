% ***********************************************
% Copyright (c) 2017 UChicago Argonne, LLC
% See LICENSE file.
% ***********************************************
%
% This example illustrates the usage of gixsdiffpos in the script mode to
% calcuate the diffraction postions

%% construct lattice and enter other paramters
a = 330; b=420; c=Inf;
alpha = 90; beta=alpha; gamma=alpha;
lattice = [a,b,c,alpha,beta,gamma];

%sguvw = [0 0 0; 1/2 1/2 0];    % use fractional basis coordinates 
sguvw = 35;             % use space group

orientationmethod = 2;      % use the unit cell frame
orientation = [0 1 0];      % [uvw]

k0 =  2*pi/1.6869;      % wave vector of incident beam
alpha_i = 0.23;         % incident angle
alpha_c = 0.19;         % average critical angle of the film

mu = 300;               % linear mass attenuation coefficient of the film (unit: 1/cm)                 
nfilm = 1-(alpha_c*pi/180)^2/2+1i*mu/1e8/(2*k0);    % index of refraction of the film

hlist = -4:4;
klist = -5:5;
llist = 0;

qdeadband = 1e-3;
qcutoff = 1e-3;        % q tolerance (Unit: A^-1)

%% calculate 
y = gixsdiffpos(lattice,sguvw,orientation,orientationmethod,hlist,klist,llist,k0,alpha_i,nfilm,qdeadband,qcutoff);

%% --- plot against angles
markersize = 6;
linewidth = 1;
figure
hold on;
plot(y.angle_t(:,1),y.angle_t(:,2),'ro','markersize',markersize,'linewidth',linewidth);
plot(y.angle_r(:,1),y.angle_r(:,2),'ks','markersize',markersize,'linewidth',linewidth);
plot(-y.angle_t(:,1),y.angle_t(:,2),'ro','markersize',markersize,'linewidth',linewidth);
plot(-y.angle_r(:,1),y.angle_r(:,2),'ks','markersize',markersize,'linewidth',linewidth);
hold off; box on;
legend('DWBA transmission','DWBA reflection');
xlabel('2\Theta (degree)');
ylabel('\alpha_f(degree)');
fontsize = 8;
% label the reflection channel
for ii=1:length(y.miller)
    hkl_str = cellfun(@(x)num2str(x,'%g%g%g'),y.miller(ii),'UniformOutput',0);
    text(-y.angle_t(ii,1)+0.1,y.angle_t(ii,2),hkl_str,'color','k','fontsize',fontsize);    
end
% label the transimission channel
for ii=1:length(y.miller)
    hkl_str = cellfun(@(x)num2str(x,'%g%g%g'),y.miller(ii),'UniformOutput',0);
    text(-y.angle_r(ii,1)+0.1,y.angle_r(ii,2),hkl_str,'color','r','fontsize',fontsize);    
end

%% plot against q
figure
hold on;
plot([-y.q_dwba_t(:,2);y.q_dwba_t(:,2)],[y.q_dwba_t(:,3);y.q_dwba_t(:,3)],'ro','markersize',markersize,'linewidth',linewidth);
plot([-y.q_dwba_r(:,2);y.q_dwba_r(:,2)],[y.q_dwba_r(:,3);y.q_dwba_r(:,3)],'ks','markersize',markersize,'linewidth',linewidth);
plot([-y.q_ba(:,2);y.q_ba(:,2)],[y.q_ba(:,3);y.q_ba(:,3)],'g^');
%plot(y.q_ba_full(:,2),y.q_ba_full(:,3),'bs');
hold off; box on;
set(gca,'ylim',[-0.01,max(y.q_ba(:,3))]);
legend('DWBA transmission','DWBA relfection','BA')
xlabel('q_y (A^{-1})')
ylabel('q_z (A^{-1})')
