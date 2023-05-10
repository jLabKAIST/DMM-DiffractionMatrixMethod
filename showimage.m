function [struc] = showimage(param,up,down,su,sd)
% show structure of periodically corrugated OLEDs
% (show real part of refractive indices)
struc = zeros((sum(up.thickness)+sum(down.thickness)+param.Org_thickness)*1e9+20,100);
struc(:,:) = up.n(1);
struc(1:10,:) = up.n(end);
struc(11:10+round(sum(up.thickness)*1e9),:) = flip(su);
struc(end-9-round(sum(down.thickness)*1e9):end-10,:) = sd;

struc(end-9:end,:) = down.n(end);
figure
imagesc(real(struc))
hold on;
if param.singledip == 1
    text(param.x_dipole/param.period*100,round(1e9*(param.Org_thickness-param.Dipole_pos+sum(up.thickness))+10), 'O','Color', 'r')
else
    yline(round(1e9*(param.Org_thickness-param.Dipole_pos+sum(up.thickness))+10), 'Color', 'r');
end
end