clear
close all
fix(clock)
tic
addpath('reticolo_allege');
%% Parameter Define
% All variables are in SI unit (m)

n_top = 1.0;
n_org = 1.8;                             
n_Ag = 0.129807+3.09889i;
n_al= 0.83901+6.32423i;

% unit : m

param.wavelengths = 520*1e-9;     % Wavelength of light
param.periodx = 340*1e-9;          % Period of grating
param.periody = 340*1e-9;
param.period = [param.periodx,param.periody];


Grating_thickness = 20*1e-9;  % Thickness of grating
Ag_thickness = 15*1e-9;

param.nns = [10,10];                     % Fourier orders kept in RCWA (-nns:nns)
param.resolution = 10;             % # of points between 0th and 1st order in k-space
param.resolutiony = 10;

param.k0=2*pi/(param.wavelengths);

param.dx = param.wavelengths/param.periodx/n_org;
param.dy = param.wavelengths/param.periody/n_org;

param.uxsize = [-3*param.resolution:3*param.resolution]*param.dx/param.resolution;
param.uysize = [-3*param.resolutiony:3*param.resolutiony]*param.dy/param.resolutiony;

down.n = [n_org,[n_org,n_al],n_al];
down.thickness = [0,Grating_thickness,0];
down.grating = [0,1,0];
down.dutyx = [0,0.5,0];
down.dutyy = [0,0.5,0];

up.n = [n_org,n_Ag,n_top];
up.thickness = [0,Ag_thickness,0];

filename = ['RCWA\Org2D_nn' num2str(param.nns) '_res' num2str(param.resolution) '_period_' num2str(param.period*1e9) '_G',num2str(Grating_thickness*10^9), 'nm.mat'];

if isfile(filename)

load(filename)
timeLoad = toc

else

[rd,td,sd] = RCWA_DMM2D(param,down);
[ru,tu,su] = TMM_DMM2D(param,up);


timeRCWA = toc
save(filename)

end

%% AFTER PROCESSING PARAMETERS

param.x_dipole = 1/2*param.periodx; % x-position of dipole
param.y_dipole = 0*param.periody; % x-position of dipole

param.Org_thickness = 215*1e-9;
param.Dipole_pos = 55*1e-9;
param.uk = 1;

param.Dox = 1;
param.Doy = 1;
param.Doz = 1;

result = DMM_process2D(param,up,rd,ru,tu);
time = toc

toc