function [rd,td,struc] = RCWA_DMM2D(param,down)

uxsize = param.uxsize;
uysize = param.uysize;
nns = param.nns;
period = param.period;
periodx = param.periodx;
periody = param.periody;
wavelengths = param.wavelengths;

n_org = down.n(1);

rd = struct;
td = struct;
struc = zeros(round(sum(down.thickness*1e9)),100,100);
numlayer = length(down.grating);

% RCWA_DOWN PARAMETERS

R_amp_TMTM = zeros(length(uxsize), length(uysize),13,13);
R_amp_TETE = zeros(length(uxsize), length(uysize),13,13);
R_amp_TMTE = zeros(length(uxsize), length(uysize),13,13);
R_amp_TETM = zeros(length(uxsize), length(uysize),13,13);

% T_amp_TMTM = zeros(length(uxsize), length(uysize),13,13);
% T_amp_TETE = zeros(length(uxsize), length(uysize),13,13);
% T_amp_TMTE = zeros(length(uxsize), length(uysize),13,13);
% T_amp_TETM = zeros(length(uxsize), length(uysize),13,13);

E_I_TEu = zeros(length(uxsize), length(uysize));
E_I_TMu = zeros(length(uxsize), length(uysize));

E_R_TMTMu = zeros(length(uxsize), length(uysize),13,13);
E_R_TETEu = zeros(length(uxsize), length(uysize),13,13);
E_R_TMTEu = zeros(length(uxsize), length(uysize),13,13);
E_R_TETMu = zeros(length(uxsize), length(uysize),13,13);

% E_T_TMTMu = zeros(length(uxsize), length(uysize),13,13);
% E_T_TETEu = zeros(length(uxsize), length(uysize),13,13);
% E_T_TMTEu = zeros(length(uxsize), length(uysize),13,13);
% E_T_TETMu = zeros(length(uxsize), length(uysize),13,13);
%% textures

textures = cell(1, numlayer);
j = 0;
k = 0; 
for i = 1:numlayer
    struc(k+1:k+down.thickness(i)*1e9,:,:) = down.n(i+j);
    if down.grating(i) == 0
        textures{i} = {down.n(i+j)};
    else
        textures{i} = {down.n(i+j),[0, 0, down.dutyx(i)*periodx, down.dutyy(i)*periody,down.n(i+j+1),1]};
        struc(k+1:k+down.thickness(i)*1e9,1:100*down.dutyx(i),1:100*down.dutyy(i)) = down.n(i+j+1);
        j = j+1;
    end
    k = k+down.thickness(i)*1e9;
end
profile = {down.thickness, 1:numlayer};

%% RCWA_DOWN

for uxindex = 1:length(uxsize)
parfor uyindex = 1:length(uysize)


[prv,vmax]=retio([],inf*1i);
retio;

parm = res0;
parm.res1.champ = 1;      % calculate precisely

parm.res2.tolh = inf;     % calculate evanescent fields
parm.res2.tolb = inf;     % calculate evanescent fields

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

for i=1:13
    for j = 1:13
        R_amp_TMTM(uxindex,uyindex,i,j) = sum(result.TMinc_top_reflected.amplitude_TM{i-7,j-7});
        R_amp_TETE(uxindex,uyindex,i,j) = sum(result.TEinc_top_reflected.amplitude_TE{i-7,j-7});
        R_amp_TMTE(uxindex,uyindex,i,j) = sum(result.TMinc_top_reflected.amplitude_TE{i-7,j-7});
        R_amp_TETM(uxindex,uyindex,i,j) = sum(result.TEinc_top_reflected.amplitude_TM{i-7,j-7});
        
        E_R_TETEu(uxindex,uyindex,i,j) = sum(sum(result.TEinc_top_reflected.PlaneWave_TE_Eu{i-7,j-7}));
        E_R_TETMu(uxindex,uyindex,i,j) = sum(sum(result.TEinc_top_reflected.PlaneWave_TM_Eu{i-7,j-7}));
        E_R_TMTEu(uxindex,uyindex,i,j) = sum(sum(result.TMinc_top_reflected.PlaneWave_TE_Eu{i-7,j-7}));
        E_R_TMTMu(uxindex,uyindex,i,j) = sum(sum(result.TMinc_top_reflected.PlaneWave_TM_Eu{i-7,j-7}));

%         T_amp_TMTM(uxindex,uyindex,i,j) = sum(result.TMinc_top_transmitted.amplitude_TM{i-7,j-7});
%         T_amp_TETE(uxindex,uyindex,i,j) = sum(result.TEinc_top_transmitted.amplitude_TE{i-7,j-7});
%         T_amp_TMTE(uxindex,uyindex,i,j) = sum(result.TMinc_top_transmitted.amplitude_TE{i-7,j-7});
%         T_amp_TETM(uxindex,uyindex,i,j) = sum(result.TEinc_top_transmitted.amplitude_TM{i-7,j-7});
%         
%         E_T_TETEu(uxindex,uyindex,i,j) = sum(sum(result.TEinc_top_transmitted.PlaneWave_TE_Eu{i-7,j-7}));
%         E_T_TETMu(uxindex,uyindex,i,j) = sum(sum(result.TEinc_top_transmitted.PlaneWave_TM_Eu{i-7,j-7}));
%         E_T_TMTEu(uxindex,uyindex,i,j) = sum(sum(result.TMinc_top_transmitted.PlaneWave_TE_Eu{i-7,j-7}));
%         E_T_TMTMu(uxindex,uyindex,i,j) = sum(sum(result.TMinc_top_transmitted.PlaneWave_TM_Eu{i-7,j-7}));
    end
end
E_I_TEu(uxindex,uyindex) = result.TEinc_top.PlaneWave_TE_Eu(2);
E_I_TMu(uxindex,uyindex) = result.TMinc_top.PlaneWave_TM_Eu(1);

retio;
end

if uxindex == 10
   timefor1 = toc; 
   fprintf(['Expected time : ' num2str(timefor1*length(uxsize)/10) 's\n'])
end
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
