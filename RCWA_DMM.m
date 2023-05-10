function [rd,td,struc] = RCWA_DMM(param,down)

%% RCWA_DOWN PARAMETERS

uxsize = param.uxsize;
uysize = param.uysize;
nns = param.nns;
period = param.period;
wavelengths = param.wavelengths;

n_org = down.n(1);

rd = struct;
td = struct;
struc = zeros(round(sum(down.thickness*1e9)),100);
numlayer = length(down.grating);

%% Declare variables
R_amp_TMTM = zeros(length(uxsize), length(uysize),21);
R_amp_TETE = zeros(length(uxsize), length(uysize),21);
R_amp_TMTE = zeros(length(uxsize), length(uysize),21);
R_amp_TETM = zeros(length(uxsize), length(uysize),21);

E_I_TEu = zeros(length(uxsize), length(uysize));
E_I_TMu = zeros(length(uxsize), length(uysize));

E_R_TMTMu = zeros(length(uxsize), length(uysize),21);
E_R_TETEu = zeros(length(uxsize), length(uysize),21);
E_R_TMTEu = zeros(length(uxsize), length(uysize),21);
E_R_TETMu = zeros(length(uxsize), length(uysize),21);
% 
% T_amp_TMTM = zeros(length(uxsize), length(uysize),21);
% T_amp_TETE = zeros(length(uxsize), length(uysize),21);
% T_amp_TMTE = zeros(length(uxsize), length(uysize),21);
% T_amp_TETM = zeros(length(uxsize), length(uysize),21);
% 
% E_T_TMTMu = zeros(length(uxsize), length(uysize),21);
% E_T_TETEu = zeros(length(uxsize), length(uysize),21);
% E_T_TMTEu = zeros(length(uxsize), length(uysize),21);
% E_T_TETMu = zeros(length(uxsize), length(uysize),21);

if nns < 10
    R_amp_TMTM2 = zeros(length(uxsize), length(uysize),2*nns+1);
    R_amp_TETE2 = zeros(length(uxsize), length(uysize),2*nns+1);
    R_amp_TMTE2 = zeros(length(uxsize), length(uysize),2*nns+1);
    R_amp_TETM2 = zeros(length(uxsize), length(uysize),2*nns+1);
        
    E_R_TMTMu2 = zeros(length(uxsize), length(uysize),2*nns+1);
    E_R_TETEu2 = zeros(length(uxsize), length(uysize),2*nns+1);
    E_R_TMTEu2 = zeros(length(uxsize), length(uysize),2*nns+1);
    E_R_TETMu2 = zeros(length(uxsize), length(uysize),2*nns+1);

%     T_amp_TMTM2 = zeros(length(uxsize), length(uysize),2*nns+1);
%     T_amp_TETE2 = zeros(length(uxsize), length(uysize),2*nns+1);
%     T_amp_TMTE2 = zeros(length(uxsize), length(uysize),2*nns+1);
%     T_amp_TETM2 = zeros(length(uxsize), length(uysize),2*nns+1);
% 
%     E_T_TMTMu2 = zeros(length(uxsize), length(uysize),2*nns+1);
%     E_T_TETEu2 = zeros(length(uxsize), length(uysize),2*nns+1);
%     E_T_TMTEu2 = zeros(length(uxsize), length(uysize),2*nns+1);
%     E_T_TETMu2 = zeros(length(uxsize), length(uysize),2*nns+1);

end
%% Define Textures

textures = cell(1, numlayer);
j = 0;
k = 0; 
for i = 1:numlayer
    struc(k+1:k+down.thickness(i)*1e9,:) = down.n(i+j);
    if down.grating(i) == 0
        textures{i} = {down.n(i+j)};
    else
        textures{i} = {[0, down.duty(i)*period], [down.n(i+j), down.n(i+j+1)]};
        struc(k+1:k+down.thickness(i)*1e9,1:100*down.duty(i)) = down.n(i+j+1);
        j = j+1;
    end
    k = k+down.thickness(i)*1e9;
end
profile = {down.thickness, 1:numlayer};

parm = res0;
parm.res1.champ = 1;      % calculate precisely

parm.res2.tolh = inf;     % calculate evanescent fields
parm.res2.tolb = inf;     % calculate evanescent fields
%% RCWA_DOWN

if nns > 9


for uxindex = 1:length(uxsize)
parfor uyindex = 1:length(uysize)


[prv,vmax]=retio([],inf*1i);
retio;

ux = uxsize(uxindex);     % one point
uy = uysize(uyindex);     % one point


% To prvent division of 0 prob
if ux == 0
    ux = 1e-8;
end
if uy == 0
    uy = 1e-8;
end

% calculation of input theta and phi
theta = asind(sqrt(ux^2+uy^2));
phi = atand(uy/ux);

% k_parallel - added sign(ux) to distinguish -0.5 & 0.5
k_parallel = sign(ux) * n_org * sind(theta);
angle_delta = phi;

aa = res1(wavelengths,period,textures,nns,k_parallel,angle_delta,parm);
result = res2(aa, profile, parm);

% Calculate reflection coefficient of TE/TM plane wave input to TE/TM nth diffraction

% R_amp_TMTM0(uxindex,uyindex) = result.TMinc_top_reflected.amplitude_TM{0};
try
R_amp_TETE(uxindex,uyindex,:) = result.TEinc_top_reflected.amplitude_TE(nns+1-10:nns+1+10);
R_amp_TMTM(uxindex,uyindex,:) = result.TMinc_top_reflected.amplitude_TM(nns+1-10:nns+1+10);
R_amp_TETM(uxindex,uyindex,:) = result.TEinc_top_reflected.amplitude_TM(nns+1-10:nns+1+10);
R_amp_TMTE(uxindex,uyindex,:) = result.TMinc_top_reflected.amplitude_TE(nns+1-10:nns+1+10);

E_I_TEu(uxindex,uyindex) = result.TEinc_top.PlaneWave_TE_Eu(2);
E_I_TMu(uxindex,uyindex) = result.TMinc_top.PlaneWave_TM_Eu(1);

E_R_TETEu(uxindex,uyindex,:) = result.TEinc_top_reflected.PlaneWave_TE_Eu(nns+1-10:nns+1+10,2);
E_R_TMTMu(uxindex,uyindex,:) = result.TMinc_top_reflected.PlaneWave_TM_Eu(nns+1-10:nns+1+10,1);
E_R_TETMu(uxindex,uyindex,:) = result.TEinc_top_reflected.PlaneWave_TM_Eu(nns+1-10:nns+1+10,1);
E_R_TMTEu(uxindex,uyindex,:) = result.TMinc_top_reflected.PlaneWave_TE_Eu(nns+1-10:nns+1+10,2);

% T_amp_TETE(uxindex,uyindex,:) = result.TEinc_top_transmitted.amplitude_TE(nns+1-10:nns+1+10);
% T_amp_TMTM(uxindex,uyindex,:) = result.TMinc_top_transmitted.amplitude_TM(nns+1-10:nns+1+10);
% T_amp_TETM(uxindex,uyindex,:) = result.TEinc_top_transmitted.amplitude_TM(nns+1-10:nns+1+10);
% T_amp_TMTE(uxindex,uyindex,:) = result.TMinc_top_transmitted.amplitude_TE(nns+1-10:nns+1+10);
% 
% E_T_TETEu(uxindex,uyindex,:) = result.TEinc_top_transmitted.PlaneWave_TE_Eu(nns+1-10:nns+1+10,2);
% E_T_TMTMu(uxindex,uyindex,:) = result.TMinc_top_transmitted.PlaneWave_TM_Eu(nns+1-10:nns+1+10,1);
% E_T_TETMu(uxindex,uyindex,:) = result.TEinc_top_transmitted.PlaneWave_TM_Eu(nns+1-10:nns+1+10,1);
% E_T_TMTEu(uxindex,uyindex,:) = result.TMinc_top_transmitted.PlaneWave_TE_Eu(nns+1-10:nns+1+10,2);


catch
end

retio;
end

if uxindex == 10
   timefor1 = toc; 
   fprintf(['Expected RCWA time : ' num2str(timefor1*length(uxsize)/10) 's\n'])
end
end


else
   

for uxindex = 1:length(uxsize)
parfor uyindex = 1:length(uysize)

[prv,vmax]=retio([],inf*1i);
retio;
ux = uxsize(uxindex);     % one point
uy = uysize(uyindex);     % one point


% To prvent division of 0 prob
if ux == 0
    ux = 1e-8;
end
if uy == 0
    uy = 1e-8;
end

% calculation of input theta and phi
theta = asind(sqrt(ux^2+uy^2));
phi = atand(uy/ux);

% k_parallel - added sign(ux) to distinguish -0.5 & 0.5
k_parallel = sign(ux) * n_org * sind(theta);
angle_delta = phi;

aa = res1(wavelengths,period,textures,nns,k_parallel,angle_delta,parm);
result = res2(aa, profile, parm);

% Calculate reflection coefficient of TE/TM plane wave input to TE/TM nth diffraction

try
R_amp_TETE2(uxindex,uyindex,:) = result.TEinc_top_reflected.amplitude_TE(1:2*nns+1);
R_amp_TMTM2(uxindex,uyindex,:) = result.TMinc_top_reflected.amplitude_TM(1:2*nns+1);
R_amp_TETM2(uxindex,uyindex,:) = result.TEinc_top_reflected.amplitude_TM(1:2*nns+1);
R_amp_TMTE2(uxindex,uyindex,:) = result.TMinc_top_reflected.amplitude_TE(1:2*nns+1);

E_I_TEu(uxindex,uyindex) = result.TEinc_top.PlaneWave_TE_Eu(2);
E_I_TMu(uxindex,uyindex) = result.TMinc_top.PlaneWave_TM_Eu(1);

E_R_TETEu2(uxindex,uyindex,:) = result.TEinc_top_reflected.PlaneWave_TE_Eu(1:2*nns+1,2);
E_R_TMTMu2(uxindex,uyindex,:) = result.TMinc_top_reflected.PlaneWave_TM_Eu(1:2*nns+1,1);
E_R_TETMu2(uxindex,uyindex,:) = result.TEinc_top_reflected.PlaneWave_TM_Eu(1:2*nns+1,1);
E_R_TMTEu2(uxindex,uyindex,:) = result.TMinc_top_reflected.PlaneWave_TE_Eu(1:2*nns+1,2);
% 
% T_amp_TETE2(uxindex,uyindex,:) = result.TEinc_top_transmitted.amplitude_TE(nns+1-10:nns+1+10);
% T_amp_TMTM2(uxindex,uyindex,:) = result.TMinc_top_transmitted.amplitude_TM(nns+1-10:nns+1+10);
% T_amp_TETM2(uxindex,uyindex,:) = result.TEinc_top_transmitted.amplitude_TM(nns+1-10:nns+1+10);
% T_amp_TMTE2(uxindex,uyindex,:) = result.TMinc_top_transmitted.amplitude_TE(nns+1-10:nns+1+10);
% 
% E_T_TETEu2(uxindex,uyindex,:) = result.TEinc_top_transmitted.PlaneWave_TE_Eu(1:2*nns+1,2);
% E_T_TMTMu2(uxindex,uyindex,:) = result.TMinc_top_transmitted.PlaneWave_TM_Eu(1:2*nns+1,1);
% E_T_TETMu2(uxindex,uyindex,:) = result.TEinc_top_transmitted.PlaneWave_TM_Eu(1:2*nns+1,1);
% E_T_TMTEu2(uxindex,uyindex,:) = result.TMinc_top_transmitted.PlaneWave_TE_Eu(1:2*nns+1,2);



catch
end

retio;
end

% if uxindex == 10
%    timefor1 = toc; 
%    fprintf(['Expected time : ' num2str(timefor1*length(uxsize)/10) 's\n']) % give information of estimated calculation time
% end
end
    
    
end


if nns < 10
    R_amp_TETE(:,:,11-nns:11+nns) = R_amp_TETE2(:,:,:);
    R_amp_TETM(:,:,11-nns:11+nns) = R_amp_TETM2(:,:,:);
    R_amp_TMTE(:,:,11-nns:11+nns) = R_amp_TMTE2(:,:,:);
    R_amp_TMTM(:,:,11-nns:11+nns) = R_amp_TMTM2(:,:,:);
        
%     T_amp_TETE(:,:,11-nns:11+nns) = T_amp_TETE2(:,:,:);
%     T_amp_TETM(:,:,11-nns:11+nns) = T_amp_TETM2(:,:,:);
%     T_amp_TMTE(:,:,11-nns:11+nns) = T_amp_TMTE2(:,:,:);
%     T_amp_TMTM(:,:,11-nns:11+nns) = T_amp_TMTM2(:,:,:);

    E_R_TETEu(:,:,11-nns:11+nns) = E_R_TETEu2(:,:,:);
    E_R_TETMu(:,:,11-nns:11+nns) = E_R_TETMu2(:,:,:);
    E_R_TMTEu(:,:,11-nns:11+nns) = E_R_TMTEu2(:,:,:);
    E_R_TMTMu(:,:,11-nns:11+nns) = E_R_TMTMu2(:,:,:);

%     E_T_TETEu(:,:,11-nns:11+nns) = E_T_TETEu2(:,:,:);
%     E_T_TETMu(:,:,11-nns:11+nns) = E_T_TETMu2(:,:,:);
%     E_T_TMTEu(:,:,11-nns:11+nns) = E_T_TMTEu2(:,:,:);
%     E_T_TMTMu(:,:,11-nns:11+nns) = E_T_TMTMu2(:,:,:);
end

rd.TETE = R_amp_TETE .* E_R_TETEu ./ E_I_TEu;
rd.TMTM = R_amp_TMTM .* E_R_TMTMu ./ E_I_TMu;
rd.TETM = R_amp_TETM .* E_R_TETMu ./ E_I_TEu;
rd.TMTE = R_amp_TMTE .* E_R_TMTEu ./ E_I_TMu;

% td.TETE = T_amp_TETE .* E_T_TETEu ./ E_I_TEu;
% td.TMTM = T_amp_TMTM .* E_T_TMTMu ./ E_I_TMu;
% td.TETM = T_amp_TETM .* E_T_TETMu ./ E_I_TEu;
% td.TMTE = T_amp_TMTE .* E_T_TMTEu ./ E_I_TMu;

end
