function [u,uK,K_prime,uK_u,uK_uK,uK_K_prime] = Get_uK(param,up,down,result)

n_org = up.n(1);
n_glass = up.n(end);

duk = 0.003;

uK_u = 0:duk:2;

uK_Kh = zeros(length(uK_u),1);
uK_Kh_prime = zeros(length(uK_u),1);

uK_K = zeros(length(uK_u),1);
uK_K_prime = zeros(length(uK_u),1);

%% interp
newrange =-2.5:0.0005:2.5;
[X,Y] = meshgrid(param.uxsize);
[Xq,Yq] = meshgrid(newrange);
Lq = interp2(X,Y,result.tdip.L,Xq,Yq);
Kq = interp2(X,Y,result.tdip.K,Xq,Yq);

for i = 1:length(newrange)
for j = 1:length(newrange)

    ux = newrange(i);
    uy = newrange(j);
    
    if ux == 0
        ux = 1e-8;
    end
    if uy == 0
        uy = 1e-8;
    end
    ull = sqrt(ux^2+uy^2);
    
    uK_num = ceil(ull/duk);
    
    if uK_num <= length(uK_u)
       uK_K(uK_num) = uK_K(uK_num)+Lq(i,j)*(newrange(2)-newrange(1))*(newrange(2)-newrange(1))*2/pi; 
       uK_K_prime(uK_num) = uK_K_prime(uK_num)+Kq(i,j)*(newrange(2)-newrange(1))*(newrange(2)-newrange(1))*2/pi; 
    end
    
end
end 

%%

uK_u = uK_u + 0.5*duk;

num_out = ceil(n_glass/n_org/duk);
num_wg = ceil(1/duk);

Purcellh2 = sum(uK_Kh);

LEEh2 = sum(uK_Kh_prime(1:num_out))/Purcellh2;
Effh_out = sum(uK_Kh(1:num_out))/Purcellh2;
Effh_abs = Effh_out - LEEh2;
Effh_wg = sum(uK_Kh(num_out+1:num_wg))/Purcellh2;
Effh_sp = sum(uK_Kh(num_wg+1:end))/Purcellh2;

Purcell2 = sum(uK_K);

LEE2 = sum(uK_K_prime(1:num_out))/Purcell2;
Eff_out = sum(uK_K(1:num_out))/Purcell2;
Eff_abs = Eff_out - LEE2;
Eff_wg = sum(uK_K(num_out+1:num_wg))/Purcell2;
Eff_sp = sum(uK_K(num_wg+1:end))/Purcell2;


uK_Kh = uK_Kh * 1/2/duk;
uK_uKh = uK_Kh;
uK_Kh = uK_Kh ./ uK_u';

uK_Kh_prime = uK_Kh_prime ./ uK_u'* 1/2/duk;

uK_K = uK_K * 1/2/duk;
uK_uK = uK_K;
uK_K = uK_K ./ uK_u';

uK_K_prime = uK_K_prime ./ uK_u'* 1/2/duk;

% uK_u = uK_u + 0.5*duk;
% close all

CPS_test

%% Figure2

figure();
hold on
plot(uK_u,(max(2*uK_u.*uK_K',1e-14)),'LineWidth',2)
hold off
hAx=gca;                    % create an axes
hAx.LineWidth=1.1;
hAx.FontSize = 14;
ylim([10^-4 10^2])

hold on
plot(u,(2*uK),'Linewidth',2)
lgd = legend('Corr','Planar');
xlim([0 2])
hold off

title('K')
    
set(gca, 'YScale', 'log')
%%
figure();
polarplot(real(asin(n_org*uK_u)),real(max(uK_K_prime',1e-14)),'LineWidth',2,'Color','b')
hold on
polarplot(real(asin(n_org*u)),real(K_prime),'Linewidth',2,'Color','r')
polarplot(-real(asin(n_org*uK_u)),real(max(uK_K_prime',1e-14)),'LineWidth',2,'Color','b')
polarplot(-real(asin(n_org*u)),real(K_prime),'Linewidth',2,'Color','r')
lgd = legend('Corr','Planar');

set(gca,'ThetaZeroLocation','top')
thetalim([-90 90])

end
