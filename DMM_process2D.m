function [result] = DMM_process2D(param,up,rd,ru,tu)
%% AFTER PROCESSING PARAMETERS

small = 1e-12*1i;      % very small value to improve calculation speed
E_Coefficient = 1;
L_coefficient = 3/8;   % Factor at power inside
I_Coefficient = 3/16;  % Factor at power outside

loss = 0;
Dox = param.Dox;
Doy = param.Doy;
Doz = param.Doz;
x_dipole = param.x_dipole; % x-position of dipole
y_dipole = param.y_dipole;

uxsize = param.uxsize;
uysize = param.uysize;
resolution = param.resolution;
resolutiony = param.resolutiony;

k0 = param.k0;
n_org = up.n(1);
n_top = up.n(end);

Org_thickness = param.Org_thickness;
Dipole_pos = param.Dipole_pos;


%% Matrices

lenux = length(uxsize);
Diffsize = ceil(lenux/resolution);

lenuy = length(uysize);
Diffsizey = ceil(lenuy/resolutiony);

theta = zeros(2*Diffsize*Diffsizey,resolution,resolutiony);
bplus = zeros(2*Diffsize*Diffsizey,2*Diffsize*Diffsizey,resolution,resolutiony);
bminus = zeros(2*Diffsize*Diffsizey,2*Diffsize*Diffsizey,resolution,resolutiony);
e12 = zeros(2*Diffsize*Diffsizey,2*Diffsize*Diffsizey,resolution,resolutiony);
e21 = zeros(2*Diffsize*Diffsizey,2*Diffsize*Diffsizey,resolution,resolutiony);

r1 = zeros(2*Diffsize*Diffsizey,2*Diffsize*Diffsizey,resolution,resolutiony);
t1 = zeros(2*Diffsize*Diffsizey,2*Diffsize*Diffsizey,resolution,resolutiony);
r2 = zeros(2*Diffsize*Diffsizey,2*Diffsize*Diffsizey,resolution,resolutiony);
a2 = zeros(2*Diffsize*Diffsizey,2*Diffsize*Diffsizey,resolution,resolutiony);
a1a2 = zeros(2*Diffsize*Diffsizey,2*Diffsize*Diffsizey,resolution,resolutiony);
inva1a2 = zeros(2*Diffsize*Diffsizey,2*Diffsize*Diffsizey,resolution,resolutiony);

cplus = zeros(2*Diffsize*Diffsizey,2*Diffsize*Diffsizey,resolution,resolutiony);
cminus = zeros(2*Diffsize*Diffsizey,2*Diffsize*Diffsizey,resolution,resolutiony);
invcpm = zeros(2*Diffsize*Diffsizey,2*Diffsize*Diffsizey,resolution,resolutiony);
invcmp = zeros(2*Diffsize*Diffsizey,2*Diffsize*Diffsizey,resolution,resolutiony);

Eplusx = zeros(2*Diffsize*Diffsizey,resolution,resolutiony);
Eminusx = zeros(2*Diffsize*Diffsizey,resolution,resolutiony);
Eplusy = zeros(2*Diffsize*Diffsizey,resolution,resolutiony);
Eminusy = zeros(2*Diffsize*Diffsizey,resolution,resolutiony);
Eplusz = zeros(2*Diffsize*Diffsizey,resolution,resolutiony);
Eminusz = zeros(2*Diffsize*Diffsizey,resolution,resolutiony);

EinTEx = zeros(2*Diffsize*Diffsizey,resolution,resolutiony);
EinTMx = zeros(2*Diffsize*Diffsizey,resolution,resolutiony);
Linx = zeros(2*Diffsize*Diffsizey,resolution,resolutiony);

EinTEy = zeros(2*Diffsize*Diffsizey,resolution,resolutiony);
EinTMy = zeros(2*Diffsize*Diffsizey,resolution,resolutiony);
Liny = zeros(2*Diffsize*Diffsizey,resolution,resolutiony);

EinTEz = zeros(2*Diffsize*Diffsizey,resolution,resolutiony);
EinTMz = zeros(2*Diffsize*Diffsizey,resolution,resolutiony);
Linz = zeros(2*Diffsize*Diffsizey,resolution,resolutiony);

Eout1x = zeros(2*Diffsize*Diffsizey,resolution,resolutiony);
Eout2x = zeros(2*Diffsize*Diffsizey,resolution,resolutiony);
Eoutx = zeros(2*Diffsize*Diffsizey,resolution,resolutiony);
Ioutx = zeros(2*Diffsize*Diffsizey,resolution,resolutiony);
Koutx = zeros(2*Diffsize*Diffsizey,resolution,resolutiony);

Eout1y = zeros(2*Diffsize*Diffsizey,resolution,resolutiony);
Eout2y = zeros(2*Diffsize*Diffsizey,resolution,resolutiony);
Eouty = zeros(2*Diffsize*Diffsizey,resolution,resolutiony);
Iouty = zeros(2*Diffsize*Diffsizey,resolution,resolutiony);
Kouty = zeros(2*Diffsize*Diffsizey,resolution,resolutiony);

Eout1z = zeros(2*Diffsize*Diffsizey,resolution,resolutiony);
Eout2z = zeros(2*Diffsize*Diffsizey,resolution,resolutiony);
Eoutz = zeros(2*Diffsize*Diffsizey,resolution,resolutiony);
Ioutz = zeros(2*Diffsize*Diffsizey,resolution,resolutiony);
Koutz = zeros(2*Diffsize*Diffsizey,resolution,resolutiony);

Kouttempx = zeros(resolution*Diffsize*2,length(uysize));
Iouttempx = zeros(resolution*Diffsize*2,length(uysize));

Lxytempx = zeros(resolution*Diffsize*2,length(uysize));
Lxyx = zeros(length(uxsize),length(uysize));

Kouttempy = zeros(resolution*Diffsize*2,length(uysize));
Iouttempy = zeros(resolution*Diffsize*2,length(uysize));

Lxytempy = zeros(resolution*Diffsize*2,length(uysize));
Lxyy = zeros(length(uxsize),length(uysize));

Kouttempz = zeros(resolution*Diffsize*2,length(uysize));
Iouttempz = zeros(resolution*Diffsize*2,length(uysize));

Lxytempz = zeros(resolution*Diffsize*2,length(uysize));
Lxyz = zeros(length(uxsize),length(uysize));

%% r1 r2
for z = 1:resolution
for zy = 1:resolutiony
for uyindex = 1:Diffsizey
for uxindex = 1:Diffsize
if (uxindex-1)*resolution+z<=lenux
if (uyindex-1)*resolutiony+zy<=lenuy
    
    ux = uxsize((uxindex-1)*resolution+z);
    uy = uysize((uyindex-1)*resolutiony+zy);

    if ux == 0
        ux = 1e-8;
    end
    if uy == 0
        uy = 1e-8;
    end
    ull = sqrt(ux^2+uy^2);
    uz = sqrt(1-ull^2);
    if abs(uz) < 1e-4
        uz = 1e-4;
    end

    % reflectance by RCWA
    % r1 = r_up
    % t1 = t_up

    for px = -10 : 10
    for py = -10 : 10
        D_change = px * Diffsize + py;

        if uxindex + px < Diffsize + 1
        if uyindex + py < Diffsize + 1
        if uxindex + px > 0
        if uyindex + py > 0
            r1((uxindex-1)*Diffsize+uyindex,(uxindex-1)*Diffsize+uyindex+D_change,z,zy) = ru.TETE((uxindex-1)*resolution+z,(uyindex-1)*resolutiony+zy,7+px,7+py);
            r1((uxindex-1)*Diffsize+uyindex+Diffsize*Diffsizey,(uxindex-1)*Diffsize+uyindex+D_change,z,zy) = ru.TMTE((uxindex-1)*resolution+z,(uyindex-1)*resolutiony+zy,7+px,7+py);
            r1((uxindex-1)*Diffsize+uyindex,(uxindex-1)*Diffsize+uyindex+Diffsize*Diffsizey+D_change,z,zy) = ru.TETM((uxindex-1)*resolution+z,(uyindex-1)*resolutiony+zy,7+px,7+py);
            r1((uxindex-1)*Diffsize+uyindex+Diffsize*Diffsizey,(uxindex-1)*Diffsize+uyindex+Diffsize*Diffsizey+D_change,z,zy) = ru.TMTM((uxindex-1)*resolution+z,(uyindex-1)*resolutiony+zy,7+px,7+py);
        end
        end
        end
        end
        
    end
    end

    for px = -10 : 10
    for py = -10 : 10
        D_change = px * Diffsize + py;

        if uxindex + px < Diffsize + 1
        if uyindex + py < Diffsize + 1
        if uxindex + px > 0
        if uyindex + py > 0
            t1((uxindex-1)*Diffsize+uyindex,(uxindex-1)*Diffsize+uyindex+D_change,z,zy) = tu.TETE((uxindex-1)*resolution+z,(uyindex-1)*resolutiony+zy,7+px,7+py);
            t1((uxindex-1)*Diffsize+uyindex+Diffsize*Diffsizey,(uxindex-1)*Diffsize+uyindex+D_change,z,zy) = tu.TMTE((uxindex-1)*resolution+z,(uyindex-1)*resolutiony+zy,7+px,7+py);
            t1((uxindex-1)*Diffsize+uyindex,(uxindex-1)*Diffsize+uyindex+Diffsize*Diffsizey+D_change,z,zy) = tu.TETM((uxindex-1)*resolution+z,(uyindex-1)*resolutiony+zy,7+px,7+py);
            t1((uxindex-1)*Diffsize+uyindex+Diffsize*Diffsizey,(uxindex-1)*Diffsize+uyindex+Diffsize*Diffsizey+D_change,z,zy) = tu.TMTM((uxindex-1)*resolution+z,(uyindex-1)*resolutiony+zy,7+px,7+py);
        end
        end
        end
        end
        
    end
    end
    % r2 = diffraction at bottom
    for px = -10 : 10
    for py = -10 : 10
        D_change = px * Diffsize + py;

        if uxindex + px < Diffsize + 1
        if uyindex + py < Diffsize + 1
        if uxindex + px > 0
        if uyindex + py > 0
            r2((uxindex-1)*Diffsize+uyindex,(uxindex-1)*Diffsize+uyindex+D_change,z,zy) = rd.TETE((uxindex-1)*resolution+z,(uyindex-1)*resolutiony+zy,7+px,7+py);
            r2((uxindex-1)*Diffsize+uyindex+Diffsize*Diffsizey,(uxindex-1)*Diffsize+uyindex+D_change,z,zy) = rd.TMTE((uxindex-1)*resolution+z,(uyindex-1)*resolutiony+zy,7+px,7+py);
            r2((uxindex-1)*Diffsize+uyindex,(uxindex-1)*Diffsize+uyindex+Diffsize*Diffsizey+D_change,z,zy) = rd.TETM((uxindex-1)*resolution+z,(uyindex-1)*resolutiony+zy,7+px,7+py);
            r2((uxindex-1)*Diffsize+uyindex+Diffsize*Diffsizey,(uxindex-1)*Diffsize+uyindex+Diffsize*Diffsizey+D_change,z,zy) = rd.TMTM((uxindex-1)*resolution+z,(uyindex-1)*resolutiony+zy,7+px,7+py);
        end
        end
        end
        end
        
    end
    end
end
end
end
end
end
end

%% b+ b- e12 e21 theta theta_air

for z = 1:resolution
for zy = 1:resolutiony
for uyindex = 1:Diffsizey
for uxindex = 1:Diffsize

if (uxindex-1)*resolution+z<=lenux
if (uyindex-1)*resolutiony+zy<=lenuy

    ux = uxsize((uxindex-1)*resolution+z);
    uy = uysize((uyindex-1)*resolutiony+zy);

    if ux == 0
        ux = 1e-8;
    end
    if uy == 0
        uy = 1e-8;
    end
    ull = sqrt(ux^2+uy^2);
    uz = sqrt(1-ull^2);
    if abs(uz) < 1e-4
        uz = 1e-4;
    end

    theta((uxindex-1)*Diffsize+uyindex,z,zy) = asind(sqrt(ux^2+uy^2));
    theta((uxindex-1)*Diffsize+uyindex+Diffsize*Diffsizey,z,zy) = asind(sqrt(ux^2+uy^2));
end
end
end
end
        % propagation matrices
        % bplus : from dipole to top layer
        % bmins : from dipole to bottom layer
        % e12 = e21 : from top to bottom or bottom to top
        bplus(:,:,z,zy) = diag(exp(1i*k0*n_org*cosd(theta(:,z,zy))*(Org_thickness-Dipole_pos)));
        bminus(:,:,z,zy) = diag(exp(1i*k0*n_org*cosd(theta(:,z,zy))*(Dipole_pos)));
        
        e12(:,:,z,zy) = diag(exp(1i*k0*n_org*cosd(theta(:,z,zy))*(Org_thickness)));
        e21(:,:,z,zy) = diag(exp(1i*k0*n_org*cosd(theta(:,z,zy))*(Org_thickness)));
end
end

theta_air = asind(sind(theta)*n_org/n_top);

%% E and I, a1a2, c+, c-


for z = 1:resolution
for zy = 1:resolutiony
for uyindex = 1:Diffsizey
for uxindex = 1:Diffsize
if (uxindex-1)*resolution+z<=lenux
if (uyindex-1)*resolutiony+zy<=lenuy

    ux = uxsize((uxindex-1)*resolution+z);
    uy = uysize((uyindex-1)*resolutiony+zy);
    if ux == 0
        ux = 1e-8;
    end
    if uy == 0
        uy = 1e-8;
    end
    ull = sqrt(ux^2+uy^2);
    uz = sqrt(1-ull^2);
    if abs(uz) < 1e-4
        uz = 1e-4;
    end

    Expe = exp(-1i*k0*n_org*(x_dipole*ux+y_dipole*uy));
   
    % E field ( x-dipole )
    % See pdf slide
    % Expe : position factor e^-ikr

    if Dox == 1
        Eplusx((uxindex-1)*Diffsize+uyindex,z,zy) = E_Coefficient * uy/uz/ull * Expe; %% uy / uz / ull
        Eplusx((uxindex-1)*Diffsize+uyindex+Diffsize*Diffsizey,z,zy) = -E_Coefficient * ux/ull * Expe; %% ux / ull
        Eminusx((uxindex-1)*Diffsize+uyindex,z,zy) = E_Coefficient * uy/uz/ull * Expe;
        Eminusx((uxindex-1)*Diffsize+uyindex+Diffsize*Diffsizey,z,zy) = E_Coefficient * ux/ull * Expe;
    end
    % E field ( y-dipole )
    if Doy == 1
        Eplusy((uxindex-1)*Diffsize+uyindex,z,zy) = E_Coefficient * -ux/uz/ull * Expe;
        Eplusy((uxindex-1)*Diffsize+uyindex+Diffsize*Diffsizey,z,zy) = -E_Coefficient * uy/ull * Expe;
        Eminusy((uxindex-1)*Diffsize+uyindex,z,zy) = E_Coefficient * -ux/uz/ull * Expe;
        Eminusy((uxindex-1)*Diffsize+uyindex+Diffsize*Diffsizey,z,zy) = E_Coefficient * uy/ull * Expe;
    end
    % E field ( z-dipole )
    if Doz ==1
        Eplusz((uxindex-1)*Diffsize+uyindex,z,zy) = E_Coefficient * 0 * Expe;
        Eplusz((uxindex-1)*Diffsize+uyindex+Diffsize*Diffsizey,z,zy) = E_Coefficient * ull/uz * Expe;
        Eminusz((uxindex-1)*Diffsize+uyindex,z,zy) = E_Coefficient * 0 * Expe;
        Eminusz((uxindex-1)*Diffsize+uyindex+Diffsize*Diffsizey,z,zy) = E_Coefficient * ull/uz * Expe;
    end
    end
end
end
end
% a1a2,a2 -> for calculation of transmittance (transmitted out)
a2(:,:,z,zy) = (r2(:,:,z,zy))*(e21(:,:,z,zy));
a1a2(:,:,z,zy) = (r1(:,:,z,zy))*(e12(:,:,z,zy))*(a2(:,:,z,zy));
inva1a2(:,:,z,zy) = inv(eye(2*Diffsize*Diffsizey)-(a1a2(:,:,z,zy)));

% c's -> for calculation of electric field at dipole
cplus(:,:,z,zy) = (bplus(:,:,z,zy))*(r1(:,:,z,zy))*(bplus(:,:,z,zy));
cminus(:,:,z,zy) = (bminus(:,:,z,zy))*(r2(:,:,z,zy))*(bminus(:,:,z,zy));
invcpm(:,:,z,zy) = inv(eye(2*Diffsize*Diffsizey)-(cplus(:,:,z,zy))*(cminus(:,:,z,zy)));
invcmp(:,:,z,zy) = inv(eye(2*Diffsize*Diffsizey)-(cminus(:,:,z,zy))*(cplus(:,:,z,zy)));


% Eout
if Dox == 1
Eout1x(:,z,zy) = (Eplusx(:,z,zy)).'*(bplus(:,:,z,zy))*(inva1a2(:,:,z,zy))*(t1(:,:,z,zy));
Eout2x(:,z,zy) = (Eminusx(:,z,zy)).'*(bminus(:,:,z,zy))*(a2(:,:,z,zy))*(inva1a2(:,:,z,zy))*(t1(:,:,z,zy));
Eoutx(:,z,zy) = Eout1x(:,z,zy) + Eout2x(:,z,zy);

% Ioutx(:,z,zy) = abs(Eoutx(:,z,zy)).^2 .* I_Coefficient .* max(0,cosd(theta_air(:,z,zy)).^2); 
Ioutx(:,z,zy) = abs(Eoutx(:,z,zy)).^2 .* n_top/n_org .* cosd(theta_air(:,z,zy)) ./ cosd(theta(:,z,zy))*  I_Coefficient .* max(0,cosd(theta_air(:,z,zy)).^2);     
Koutx(:,z,zy) = Ioutx(:,z,zy)./ max(0,cosd(theta_air(:,z,zy)).^2);

% Dipole position E field (TE/TM divided)
EinTEx(:,z,zy) = (Eplusx(:,z,zy)).' + ( (Eminusx(:,z,zy)).'*(cminus(:,:,z,zy)) + (Eplusx(:,z,zy)).'*(cplus(:,:,z,zy))*(cminus(:,:,z,zy))  )*(invcpm(:,:,z,zy));
EinTEx(:,z,zy) = EinTEx(:,z,zy).' + ( (Eplusx(:,z,zy)).'*(cplus(:,:,z,zy)) + (Eminusx(:,z,zy)).'*(cminus(:,:,z,zy))*(cplus(:,:,z,zy))  )*(invcmp(:,:,z,zy));
EinTMx(:,z,zy) = -(Eplusx(:,z,zy)).' - ( (Eminusx(:,z,zy)).'*(cminus(:,:,z,zy)) + (Eplusx(:,z,zy)).'*(cplus(:,:,z,zy))*(cminus(:,:,z,zy))  )*(invcpm(:,:,z,zy));
EinTMx(:,z,zy) = EinTMx(:,z,zy).' + ( (Eplusx(:,z,zy)).'*(cplus(:,:,z,zy)) + (Eminusx(:,z,zy)).'*(cminus(:,:,z,zy))*(cplus(:,:,z,zy))  )*(invcmp(:,:,z,zy));

end

if Doy == 1
Eout1y(:,z,zy) = (Eplusy(:,z,zy)).'*(bplus(:,:,z,zy))*(inva1a2(:,:,z,zy))*(t1(:,:,z,zy));
Eout2y(:,z,zy) = (Eminusy(:,z,zy)).'*(bminus(:,:,z,zy))*(a2(:,:,z,zy))*(inva1a2(:,:,z,zy))*(t1(:,:,z,zy));
Eouty(:,z,zy) = Eout1y(:,z,zy) + Eout2y(:,z,zy);

Iouty(:,z,zy) = abs(Eouty(:,z,zy)).^2 .* n_top/n_org .* cosd(theta_air(:,z,zy)) ./ cosd(theta(:,z,zy))*  I_Coefficient .* max(0,cosd(theta_air(:,z,zy)).^2);     
Kouty(:,z,zy) = Iouty(:,z,zy)./ max(0,cosd(theta_air(:,z,zy)).^2);

EinTEy(:,z,zy) = (Eplusy(:,z,zy)).' + ( (Eminusy(:,z,zy)).'*(cminus(:,:,z,zy)) + (Eplusy(:,z,zy)).'*(cplus(:,:,z,zy))*(cminus(:,:,z,zy))  )*(invcpm(:,:,z,zy));
EinTEy(:,z,zy) = EinTEy(:,z,zy).' + ( (Eplusy(:,z,zy)).'*(cplus(:,:,z,zy)) + (Eminusy(:,z,zy)).'*(cminus(:,:,z,zy))*(cplus(:,:,z,zy))  )*(invcmp(:,:,z,zy));
EinTMy(:,z,zy) = -(Eplusy(:,z,zy)).' - ( (Eminusy(:,z,zy)).'*(cminus(:,:,z,zy)) + (Eplusy(:,z,zy)).'*(cplus(:,:,z,zy))*(cminus(:,:,z,zy))  )*(invcpm(:,:,z,zy));
EinTMy(:,z,zy) = EinTMy(:,z,zy).' + ( (Eplusy(:,z,zy)).'*(cplus(:,:,z,zy)) + (Eminusy(:,z,zy)).'*(cminus(:,:,z,zy))*(cplus(:,:,z,zy))  )*(invcmp(:,:,z,zy));

end

if Doz == 1
Eout1z(:,z,zy) = (Eplusz(:,z,zy)).'*(bplus(:,:,z,zy))*(inva1a2(:,:,z,zy))*(t1(:,:,z,zy));
Eout2z(:,z,zy) = (Eminusz(:,z,zy)).'*(bminus(:,:,z,zy))*(a2(:,:,z,zy))*(inva1a2(:,:,z,zy))*(t1(:,:,z,zy));
Eoutz(:,z,zy) = Eout1z(:,z,zy) + Eout2z(:,z,zy);

Ioutz(:,z,zy) = abs(Eoutz(:,z,zy)).^2 .* n_top/n_org .* cosd(theta_air(:,z,zy)) ./ cosd(theta(:,z,zy))*  I_Coefficient .* max(0,cosd(theta_air(:,z,zy)).^2);     
Koutz(:,z,zy) = Ioutz(:,z,zy)./ max(0,cosd(theta_air(:,z,zy)).^2);

EinTEz(:,z,zy) = (Eplusz(:,z,zy)).' + ( (Eminusz(:,z,zy)).'*(cminus(:,:,z,zy)) + (Eplusz(:,z,zy)).'*(cplus(:,:,z,zy))*(cminus(:,:,z,zy))  )*(invcpm(:,:,z,zy));
EinTEz(:,z,zy) = EinTEz(:,z,zy).' + ( (Eplusz(:,z,zy)).'*(cplus(:,:,z,zy)) + (Eminusz(:,z,zy)).'*(cminus(:,:,z,zy))*(cplus(:,:,z,zy))  )*(invcmp(:,:,z,zy));
EinTMz(:,z,zy) = (Eplusz(:,z,zy)).' + ( (Eminusz(:,z,zy)).'*(cminus(:,:,z,zy)) + (Eplusz(:,z,zy)).'*(cplus(:,:,z,zy))*(cminus(:,:,z,zy))  )*(invcpm(:,:,z,zy));
EinTMz(:,z,zy) = EinTMz(:,z,zy).' + ( (Eplusz(:,z,zy)).'*(cplus(:,:,z,zy)) + (Eminusz(:,z,zy)).'*(cminus(:,:,z,zy))*(cplus(:,:,z,zy))  )*(invcmp(:,:,z,zy));

end


end
end

%% Lin generator


% Power at the dipole position
% Ein * some u factors
for z = 1:resolution
for zy = 1:resolutiony
for uyindex = 1:Diffsizey
for uxindex = 1:Diffsize
if (uxindex-1)*resolution+z<=lenux
if (uyindex-1)*resolutiony+zy<=lenuy
    
    ux = uxsize((uxindex-1)*resolution+z);
    uy = uysize((uyindex-1)*resolutiony+zy);
    Expe2 = exp(+1i*k0*n_org*(x_dipole*ux+y_dipole*uy));
    
    if ux == 0
        ux = 1e-8;
    end
    if uy == 0
        uy = 1e-8;
    end
    
    ull = sqrt(ux^2+uy^2);
    uz = sqrt(1-ull^2);
    
    if Dox == 1
        Linx((uxindex-1)*Diffsize+uyindex,z,zy) = real(EinTEx((uxindex-1)*Diffsize+uyindex,z,zy) * uy/ull * Expe2); % x_dipole TE % uy/ull
        Linx((uxindex-1)*Diffsize+uyindex+Diffsize*Diffsizey,z,zy) = real(EinTMx((uxindex-1)*Diffsize+uyindex+Diffsize*Diffsizey,z,zy) * (uz*ux/ull) * Expe2); % x_dipole TM uz*ux/ull
    end
    
    if Doy ==1
        Liny((uxindex-1)*Diffsize+uyindex,z,zy) = real(EinTEy((uxindex-1)*Diffsize+uyindex,z,zy) * -ux/ull * Expe2); % y_dipole TE
        Liny((uxindex-1)*Diffsize+uyindex+Diffsize*Diffsizey,z,zy) = real(EinTMy((uxindex-1)*Diffsize+uyindex+Diffsize*Diffsizey,z,zy) * (uz*uy/ull) * Expe2); % y_dipole TM
    end
    
    if Doz ==1
        Linz((uxindex-1)*Diffsize+uyindex+Diffsize*Diffsizey,z,zy) = real(EinTMz((uxindex-1)*Diffsize+uyindex+Diffsize*Diffsizey,z,zy) * ull * Expe2); % z_dipole TM
    end
end
end
end
end
end
end


%% Make L,K matrix
Koutx(isnan(Koutx))=0;
Kouty(isnan(Kouty))=0;
Koutz(isnan(Koutz))=0;


% grouped ux -> real ux
% K * uz

for z = 1:resolution
for zy = 1:resolutiony
for uyindex = 1:Diffsizey
for uxindex = 1:Diffsize
if (uxindex-1)*resolution+z<=lenux
if (uyindex-1)*resolutiony+zy<=lenuy
    
    ux = uxsize((uxindex-1)*resolution+z);
    uy = uysize((uyindex-1)*resolutiony+zy);
    
    if ux == 0
        ux = 1e-8;
    end
    if uy == 0
        uy = 1e-8;
    end
    ull = sqrt(ux^2+uy^2);
    uz = sqrt(1-ull^2);
    uz_air = sqrt(1- (n_org/n_top)^2 * ull^2);
    if abs(uz) < 1e-6
        uz = 0;
    end
    
    
    if Dox ==1
    Iouttempx((uxindex-1)*resolution+z,(uyindex-1)*resolutiony+zy) = Ioutx((uxindex-1)*Diffsize+uyindex,z,zy)*uz/uz_air;
    Iouttempx((uxindex-1)*resolution+z+Diffsize*resolution,(uyindex-1)*resolutiony+zy) = Ioutx((uxindex-1)*Diffsize+uyindex+Diffsize*Diffsizey,z,zy)*uz/uz_air;
     
    Kouttempx((uxindex-1)*resolution+z,(uyindex-1)*resolutiony+zy) = Koutx((uxindex-1)*Diffsize+uyindex,z,zy)*uz;
    Kouttempx((uxindex-1)*resolution+z+Diffsize*resolution,(uyindex-1)*resolutiony+zy) = Koutx((uxindex-1)*Diffsize+uyindex+Diffsize*Diffsizey,z,zy)*uz;
    
    Lxytempx((uxindex-1)*resolution+z,(uyindex-1)*resolutiony+zy) = Linx((uxindex-1)*Diffsize+uyindex,z,zy);
    Lxytempx((uxindex-1)*resolution+z+Diffsize*resolution,(uyindex-1)*resolutiony+zy) = Linx((uxindex-1)*Diffsize+uyindex+Diffsize*Diffsizey,z,zy);
    if uz == 0
       Lxytempx((uxindex-1)*resolution+z,(uyindex-1)*resolutiony+zy) = 0;
       Lxytempx((uxindex-1)*resolution+z+Diffsize*resolution,(uyindex-1)*resolutiony+zy) = 0;
     
    end
    end
    
    if Doy ==1
    Iouttempy((uxindex-1)*resolution+z,(uyindex-1)*resolutiony+zy) = Iouty((uxindex-1)*Diffsize+uyindex,z,zy);
    Iouttempy((uxindex-1)*resolution+z+Diffsize*resolution,(uyindex-1)*resolutiony+zy) = Iouty((uxindex-1)*Diffsize+uyindex+Diffsize*Diffsizey,z,zy);
    
    Kouttempy((uxindex-1)*resolution+z,(uyindex-1)*resolutiony+zy) = Kouty((uxindex-1)*Diffsize+uyindex,z,zy)*uz;
    Kouttempy((uxindex-1)*resolution+z+Diffsize*resolution,(uyindex-1)*resolutiony+zy) = Kouty((uxindex-1)*Diffsize+uyindex+Diffsize*Diffsizey,z,zy)*uz;
    
    Lxytempy((uxindex-1)*resolution+z,(uyindex-1)*resolutiony+zy) = Liny((uxindex-1)*Diffsize+uyindex,z,zy);
    Lxytempy((uxindex-1)*resolution+z+Diffsize*resolution,(uyindex-1)*resolutiony+zy) = Liny((uxindex-1)*Diffsize+uyindex+Diffsize*Diffsizey,z,zy);
    end
    
    if Doz ==1
    Iouttempz((uxindex-1)*resolution+z,(uyindex-1)*resolutiony+zy) = Ioutz((uxindex-1)*Diffsize+uyindex,z,zy);
    Iouttempz((uxindex-1)*resolution+z+Diffsize*resolution,(uyindex-1)*resolutiony+zy) = Ioutz((uxindex-1)*Diffsize+uyindex+Diffsize*Diffsizey,z,zy);
    
    Kouttempz((uxindex-1)*resolution+z,(uyindex-1)*resolutiony+zy) = Koutz((uxindex-1)*Diffsize+uyindex,z,zy)*uz;
    Kouttempz((uxindex-1)*resolution+z+Diffsize*resolution,(uyindex-1)*resolutiony+zy) = Koutz((uxindex-1)*Diffsize+uyindex+Diffsize*Diffsizey,z,zy)*uz;
    
    Lxytempz((uxindex-1)*resolution+z,(uyindex-1)*resolutiony+zy) = Linz((uxindex-1)*Diffsize+uyindex,z,zy);
    Lxytempz((uxindex-1)*resolution+z+Diffsize*resolution,(uyindex-1)*resolutiony+zy) = Linz((uxindex-1)*Diffsize+uyindex+Diffsize*Diffsizey,z,zy);
    end
    
    end
end
end
end
end
end

if Dox == 1
Lxyx = ( Lxytempx(1:length(uxsize),1:length(uysize)) +Lxytempx(Diffsize*resolution+1:Diffsize*resolution+length(uxsize),1:length(uysize)) )* L_coefficient ;

Ioutuxuyx = Iouttempx(1:length(uxsize),1:length(uysize));
Ioutuxuy2x = Iouttempx(Diffsize*resolution+1:Diffsize*resolution+length(uxsize),1:length(uysize));
Ioutnormx = (Ioutuxuyx'+Ioutuxuy2x');

Koutuxuyx = Kouttempx(1:length(uxsize),1:length(uysize));
Koutuxuy2x = Kouttempx(Diffsize*resolution+1:Diffsize*resolution+length(uxsize),1:length(uysize));
Koutnormx = (Koutuxuyx'+Koutuxuy2x');
end

if Doy == 1
Lxyy = ( Lxytempy(1:length(uxsize),1:length(uysize)) +Lxytempy(Diffsize*resolution+1:Diffsize*resolution+length(uxsize),1:length(uysize)) )* L_coefficient ;

Ioutuxuyy = Iouttempy(1:length(uxsize),1:length(uysize));
Ioutuxuy2y = Iouttempy(Diffsize*resolution+1:Diffsize*resolution+length(uxsize),1:length(uysize));
Ioutnormy = (Ioutuxuyy'+Ioutuxuy2y');

Koutuxuyy = Kouttempy(1:length(uxsize),1:length(uysize));
Koutuxuy2y = Kouttempy(Diffsize*resolution+1:Diffsize*resolution+length(uxsize),1:length(uysize));
Koutnormy = (Koutuxuyy'+Koutuxuy2y');
end

if Doz == 1
Lxyz = ( Lxytempz(1:length(uxsize),1:length(uysize)) +Lxytempz(Diffsize*resolution+1:Diffsize*resolution+length(uxsize),1:length(uysize)) )* L_coefficient ;

Ioutuxuyz = Iouttempz(1:length(uxsize),1:length(uysize));
Ioutuxuy2z = Iouttempz(Diffsize*resolution+1:Diffsize*resolution+length(uxsize),1:length(uysize));
Ioutnormz = (Ioutuxuyz'+Ioutuxuy2z');

Koutuxuyz = Kouttempz(1:length(uxsize),1:length(uysize));
Koutuxuy2z = Kouttempz(Diffsize*resolution+1:Diffsize*resolution+length(uxsize),1:length(uysize));
Koutnormz = (Koutuxuyz'+Koutuxuy2z');
end
%% Error 

Lxyx(isnan(Lxyx)) = 0;

% Lxyx(Lxyx<-1000) = -1000;
% Lxyx(Lxyx>1000) = 1000;


%% Figure
% close all

% figure()
% imagesc(Ioutnormx)
% caxis([0 1])
% title('I far-field')
% 
% figure()
% imagesc(Koutnormx)
% caxis([0 1])
% title('K far-field')
if Dox == 1
LEEx = sum(sum(Koutnormx))/sum(sum(Lxyx));
Purcellx = sum(sum(Lxyx))*(uxsize(2)-uxsize(1))*(uysize(2)-uysize(1))*2/pi;
Transx = LEEx * Purcellx;

[LEEx,Purcellx,Transx]
end

if Doy == 1
LEEy = sum(sum(Koutnormy))/sum(sum(Lxyy));
Purcelly = sum(sum(Lxyy))*(uxsize(2)-uxsize(1))*(uysize(2)-uysize(1))*2/pi;
Transy = LEEy * Purcelly;

[LEEy,Purcelly,Transy]
end

if Doz == 1
LEEz = sum(sum(Koutnormz))/sum(sum(Lxyz));
Purcellz = sum(sum(Lxyz))*(uxsize(2)-uxsize(1))*(uysize(2)-uysize(1))*2/pi;
Transz = LEEz * Purcellz;

[LEEz,Purcellz,Transz]
end
% figure
% plot(uysize*2,Ioutnormx(:,4))
% xlim([-1 1])

Trans = (Transx+Transy+Transz)/3;
PurcellA = (Purcellx+Purcelly+Purcellz)/3;
LEE = Trans / PurcellA;

Transh = (Transx+Transy)/2;
Purcellh = (Purcellx+Purcelly)/2;
LEEh = Transh / Purcellh;

[LEE,PurcellA,Trans]

[LEEx,LEEy,LEEz,LEE]
[Purcellx,Purcelly,Purcellz,PurcellA,resolution]
% LEEy = sum(sum(Koutnormy))/sum(sum(Lxyy))
% Purcelly = sum(sum(Lxyy))*(uxsize(2)-uxsize(1))*(uysize(2)-uysize(1))*2/pi
[a,b] = min(abs(uxsize-n_top/n_org));
[c,d] = min(abs(uxsize+n_top/n_org));
% figure
% imagesc(Ioutnormx(d:b,d:b))
% imagesc(Ioutnormx(d:b,d:b)+Ioutnormy(d:b,d:b)+Ioutnormz(d:b,d:b))
% xlabel(['LEE = ' num2str(LEEx) ', Purcell = ' num2str(Purcellx)] )
% title([num2str(x_dipole/periodx) ', ' num2str(y_dipole/periody)])
% % figure
% % imagesc(Lxyx')
% % caxis([0 1.5])
% tgprintf_PCH('NoAir JM Done');
% toc



Ksum = (Koutnormx+Koutnormy+Koutnormz)/3;
Isum = (Ioutnormx+Ioutnormy+Ioutnormz)/3;
c_DMM = sum(sum(Ksum))*n_org^2*(uxsize(2)-uxsize(1))*(uysize(2)-uysize(1));

Ifin = Isum * LEE / c_DMM;

figure(1)
imagesc(Ifin(d:b,d:b))
colorbar

if Dox == 1
LEEx = sum(sum(Koutnormx))/sum(sum(Lxyx));
Purcellx = sum(sum(Lxyx))*(uxsize(2)-uxsize(1))*(uysize(2)-uysize(1))*2/pi;
Transx = LEEx * Purcellx;

[LEEx,Purcellx,Transx];
end

if Doy == 1
LEEy = sum(sum(Koutnormy))/sum(sum(Lxyy));
Purcelly = sum(sum(Lxyy))*(uxsize(2)-uxsize(1))*(uysize(2)-uysize(1))*2/pi;
Transy = LEEy * Purcelly;

[LEEy,Purcelly,Transy];
end

if Doz == 1
LEEz = sum(sum(Koutnormz))/sum(sum(Lxyz));
Purcellz = sum(sum(Lxyz))*(uxsize(2)-uxsize(1))*(uysize(2)-uysize(1))*2/pi;
Transz = LEEz * Purcellz;

[LEEz,Purcellz,Transz];
end

Trans = (Transx+Transy+Transz)/3;
PurcellA = (Purcellx+Purcelly+Purcellz)/3;
LEE = Trans / PurcellA;

[a,b] = min(abs(uxsize-n_top/n_org));
[c,d] = min(abs(uxsize+n_top/n_org));

Ksum = (Koutnormx+Koutnormy+Koutnormz)/3;
Isum = (Ioutnormx+Ioutnormy+Ioutnormz)/3;
Lsum = (Lxyx+Lxyy+Lxyz)/3;
c_DMM = sum(sum(Ksum))*n_org^2*(uxsize(2)-uxsize(1))*(uysize(2)-uysize(1));

Ifin = Isum * LEE / c_DMM;

result.xdip.LEE = LEEx;
result.ydip.LEE = LEEy;
result.zdip.LEE = LEEz;
result.tdip.LEE = LEE;

result.xdip.Purcell = Purcellx;
result.ydip.Purcell = Purcelly;
result.zdip.Purcell = Purcellz;
result.tdip.Purcell = PurcellA;

result.xdip.T = Transx;
result.ydip.T = Transy;
result.zdip.T = Transz;
result.tdip.T = Trans;

result.xdip.K = Koutnormx;
result.ydip.K = Koutnormy;
result.zdip.K = Koutnormz;
result.tdip.K = Ksum;

result.xdip.I = Ioutnormx;
result.ydip.I = Ioutnormy;
result.zdip.I = Ioutnormz;
result.tdip.I = Isum;

result.xdip.L = Lxyx;
result.ydip.L = Lxyy;
result.zdip.L = Lxyz;
result.tdip.L = Lsum;

result.Ifin = Ifin;

result.d = d;
result.b = b;
end