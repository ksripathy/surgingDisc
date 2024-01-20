% % This is the script to calculate the force field using the Momentum integral of control-volume contour method


close all 
clear all 



%% set plot properties 
m = {'p' 'o' '*' 'x' 's' 'd' '^' 'v' '>' '<' '+' 'h' '.'};
c = colormap(lines);
L = repmat({'-' '--' ':' '-.'},1,3);
set(0,'DefaultLinelineWidth',2)
set(0,'DefaultAxesFontSize',15)
set(0,'DefaultTextFontSize',15)
set(0,'defaultTextInterpreter','LaTex');
close gcf


%% Data Directory
datadir.main  = './';
%% steady path 
datadir.sub1  = '/Steady';
datadir.sub2 = './';



datadir.file = 'F*.*';
datafiles = dir(fullfile(datadir.main,datadir.sub1,datadir.sub2,datadir.file));

[COLORMAP_u,COLORMAP_v,COLORMAP_w,COLORMAP_p]=self_colormap();

% %% % %  initial conditions, using the average value from all the cases 
% Rho_inf =1.2065; % % kg/m^3
% % V_inf=5.89;   % % m/s 
% T=19.52+273.15; % % plotting 
%  % % K 
%  R_const= 8.314;  % %  % Pa*m^3*mol^-1*K^-1
%  m_mol=29/1000; % % the mol mass of air kg/mol
% 
% P_inf=Rho_inf*R_const*T/m_mol; % % % Pa  ideal gas law 




%% Physical properties
load('initial_cond.mat')
P_inf = p_m800(2);
V_inf = V_m800(2);
Rho_inf = rho_m800(2);
D = 60/100; %[m]


Q_inf =(1/2)*Rho_inf*(V_inf.^2);
Area_disc= pi*(D/2)^2;
% return
% % % % should the force on the nacelle be substracted? 
% cd_cube = 1.05; % % 1.05 or 2.05 I'm not so sure 
% % % % The area is 6.5cm x5cm 
% Area_nacelle = 6.5*5/100/100;   % % m^2 
% Fd=cd_cube*(Q_inf*Area_nacelle);
% 
% Area_ratio=Area_nacelle/Area_disc;  % % % the ratio of nacelle area and disc area 


% % % % general observations 
%%%%% The CMI method is very sensitive to the location of boudary, the
%%%%% shorter the wake, the lower the CT_CMI
%%%%% The AMI method is less sensitive to the location where is treated as
%%%%% the far wake. The shorter the wake, the slightly lower CT_AMI

% return 
for jj= 1%:4
    
    datadir.file = datafiles(jj).name;  
    load(fullfile(datadir.main,datadir.sub1,datadir.sub2,datadir.file))

    
% return
%% % % change the normalized variables to the actual unit , remove the negetive part, which influence the calculation of Poisson solver.
    Vx=Vx*V_inf;
    Vr=Vr*V_inf;
    X=X*D;
    R=R*D;
    
    
    
    % % % % % Test the smoothing, doesn't influence much the methods of
    % AMI, and CMI 
%     Vx=smooth2a(Vx,1,1);
%     Vr=smooth2a(Vr,1,1);
    
    
%     %% stiching for the last FOV  
%     delR=R(1,2:end)-R(1,1:end-1);
%     vals=delR(delR~=0);
%     [cid rid]=find(delR~=0);
%     if jj<4
%      R(:,rid(end)+1:end)=R(:,rid(end)+1:end)-0.5*vals(end); 
%     else
%      R(:,rid(2)+1:end)=R(:,rid(2)+1:end)-0.5*vals(2);       
%     end 
    
    
%     return 

    
% %% % to set Mask 
    Mask=ones(size(X));
%     
    % Using the mask 
%       Mask = 1-imdilate(imclose(Vx==0,strel('disk',20)),strel('disk',15)); % % creating a mask on the disk and reflection region, from Karanbir 
  % Wei's new mask 
   Mask(:,1:1200) = 1-imdilate(imclose(Vx(:,1:1200)==0,strel('disk',20)),strel('disk',10)); % creating a mask on the disk and reflection region

   if jj==3
      Mask(694:741,528:587)=0;
  else
  end 
  
% % % % % to check the if the mask only cover the disk, 
% A=imclose(Vx==0,strel('disk',20));
% B=strel('disk',10);
% C=imdilate(A,B);
% figure
% imshow(A)
% figure
% imshow(C)
% % figure
% imshow(1-Mask)
%     return 
    
    %%  Navier Stokes equations

dx=X(1,2)-X(1,1);
dr=R(2,1)-R(1,1);
% return 
dUdt=0;
dVdt=0;

  
[dUdx,dUdr]=gradient(Vx,dx,dr);
[dVdx,dVdr]=gradient(Vr,dx,dr);


dPdx=-Rho_inf*(dUdt+Vx.*dUdx+Vr.*dUdr);
dPdr=-Rho_inf*(dVdt+Vx.*dVdx+Vr.*dVdr);

[dP2dx2,~]=gradient(dPdx,dx,dr);
[~,dP2dr2]=gradient(dPdr,dx,dr);
    
% % check the size of mask, make the new mask depend on dp2dx2+dp2dr2.  don't know why this mask doesn't work 

 F_int= dP2dx2+dP2dr2;

%  Fig1=figure('units','centimeters','papersize',[24 12],'position',[3 7 24 12])
%  pcolor(X,R,F_int);
%   axis equal
%  axis([-0.52 0.95 -0.2 0.7])
%  shading flat
%  colormap(COLORMAP_w)
%  hcb= colorbar
%  title(hcb,'$\bar{w}\frac{D}{U_{\infty}}$','interpreter','latex')
%  caxis([-20 20])
%  box on
%  xlabel('x/D [-]')
%  ylabel('y/D [-]')
%  
 
%  return 
 
%  % % % adjust the mask 
%  delta_r=R(1,2:end)-R(1,1:end-1);
%  idd0=find(delta_r~=0, 3, 'last')+1;
%  
%  Idx1 = find(X(1,:)<-0.04454, 1, 'last');
%  Idx2 = find(X(1,:)<-0.009904, 1, 'last');
%  Idx3 = find(X(1,:)<0.006179, 1, 'last');
%  Idx4 = find(X(1,:)<0.1506, 1, 'last');
%  
%  
%  Idr1 = find(R(:,Idx1)>-0.05629, 1, 'last');
%  Idr2 = find(R(:,Idx1)>0.2981, 1, 'last');
%  Idr3 = find(R(:,Idx1)>0.03588, 1, 'last');
%  Idr4 = find(R(:,idd0(2))>0.03588, 1, 'last');
%   
%  
%   Mask(Idr1:end,Idx1:Idx2)=0;
%   Mask(Idr2:end,Idx2:Idx3)=0;
%   Mask(Idr3:end,Idx3:idd0(2)-1)=0;
%   Mask(Idr4:end,idd0(2):Idx4)=0;
% % % %   return 
%   figure
% imshow(1-Mask)
% % 
%  return 

%% Pressure from Bernoulli equation
PBern=P_inf+(1/2)*Rho_inf*(V_inf.^2-(Vx.^2+Vr.^2)); % Pressure using Bernoulli Equation
CPBern = (PBern-P_inf).*Mask/Q_inf;
PBern=PBern.*Mask;



% % calculated pressure field from the Poisson equation, left and top of the real boundaries in are
% % Dirichlet boundary (left is the inlet of the flow, right is the outlet of the flow, top is the outer boundary of the CV, down is the center axis of the slice of the AD)
% % in the equation below, the order are for left,right,top,bottom of the data matrix.
% % % % 
%% by applying the mask 
P=Poisson_plane_nw(Mask,dPdx,dP2dx2,dPdr,dP2dr2,PBern,dx,dr,'Dir','Nx','Dir','Nx');



% % %% % % 
% % % % % by split the domain into an upstream and a downstream domain. 
% P=zeros(size(X));
% % % % find the index of the column without 0 values 
% ind=find(sum(Vx(:,1:1000)==0,1)>0);
% ind_up=min(ind)-1;
% ind_down=max(ind)+1;
% 
% % ind = find(X(1,:)<0, 1, 'last');
% % % offset=[0:10:50];
% % % for ii=1:length(offset)
% % ind_up=ind-20;
% % ind_down=ind+20;
% P(:,1:ind_up)=Poisson_plane_nw(Mask(:,1:ind_up),dPdx(:,1:ind_up),dP2dx2(:,1:ind_up),dPdr(:,1:ind_up),dP2dr2(:,1:ind_up),PBern(:,1:ind_up),dx,dr,'Dir','Nx','Dir','Nx');
% P(:,ind_down:end)=Poisson_plane_nw(Mask(:,ind_down:end),dPdx(:,ind_down:end),dP2dx2(:,ind_down:end),dPdr(:,ind_down:end),dP2dr2(:,ind_down:end),PBern(:,ind_down:end),dx,dr,'Nx','Nx','Dir','Nx');    

CP = (P-P_inf).*Mask/Q_inf;

%% % Force determination, momentum integration

%% % % % generate the contour for momentum integration  

% x_sta=[-0.45:0.05:-0.1];
% x_end=[0.05:0.05:0.75]; 
% r_bound=[0.5:0.02:0.7];
% 
if jj==1
    r_bound=0.66;
elseif jj==2
    r_bound=0.64;
elseif jj==3   
    r_bound=0.62;
else 
    r_bound=0.58;
end 
    
% for ii=1:length(x_sta)

x1=[-0.1:0.001:0.7]*D;
y1=[0:.001:r_bound]*D;

xline= [zeros(size(y1))+min(x1) x1 zeros(size(y1))+max(x1)]; 
yline= [ y1 zeros(size(x1))+max(y1)  fliplr(y1)]; 
zline = yline.*0;
dyx=sqrt((yline(2:end) - yline(1:end-1)).^2+(xline(2:end) - xline(1:end-1)).^2);

AREA_MIDPOINT = 2*pi*( yline(2:end)+yline(1:end-1))/2.*dyx;
TOTALAREA=sum(AREA_MIDPOINT);

ny=-gradient(xline);
nx=gradient(yline);
normvec=sqrt(nx.^2+ny.^2);
NORVEC=-[nx;ny;nx.*0]./normvec;

%% velocity and pressure at the contour 
UC_VEC(1,:)=griddata(X,R,Vx,xline,yline);
UC_VEC(2,:)=griddata(X,R,Vr,xline,yline);
UC_VEC(3,:)=0;
PC=griddata(X,R,P,xline,yline);
% return 
F_CMI=F_cal_CMI([xline;yline;zline],UC_VEC,PC,NORVEC);
F_CMI_MID=(F_CMI(:,1:end-1)+F_CMI(:,2:end))/2;


FAREA=F_CMI_MID.*AREA_MIDPOINT;
LOAD=sum(FAREA,2);
CT=sum(FAREA,2)/(0.5*pi*Rho_inf*(D/2)^2*V_inf^2); % % velocity is already normalized 
% CT_CMI(jj,ii)=CT(1);
CT_CMI(jj)=CT(1);







% return 

% % % % plotting 
% % %  
% 
 Fig1=figure('units','centimeters','papersize',[24 12],'position',[3 5 24 12])
 hold on
 pcolor(X/D,R/D,CP);
 axis equal
 axis([-0.52 0.95 -0.2 0.7])
 shading flat
 colormap(COLORMAP_p)
 hcb= colorbar
 title(hcb,'$C_p$','interpreter','latex')
 caxis([-1.0 1.0])
 box on 
 xlabel('x/D [-]')
 ylabel('y/D [-]')
%  
%  
%  Fig1=figure('units','centimeters','papersize',[24 12],'position',[3 5 24 12])
%  hold on
%  pcolor(X(:,1:ind_up)/D,R(:,1:ind_up)/D,P(:,1:ind_up));
%  axis equal
% %  axis([-0.52 0.95 -0.2 0.7])
%  shading flat
%  colormap(COLORMAP_p)
%  hcb= colorbar
%  title(hcb,'$C_p$','interpreter','latex')
% %  caxis([-1.0 1.0])
%  box on 
%  xlabel('x/D [-]')
%  ylabel('y/D [-]')
%  % % %  
%  Fig2=figure('units','centimeters','papersize',[24 12],'position',[3 5 24 12])
%  hold on
%  pcolor(X(:,ind_down:end)/D,R(:,ind_down:end)/D,P(:,ind_down:end));
%  axis equal
% %  axis([-0.52 0.95 -0.2 0.7])
%  shading flat
%  colormap(COLORMAP_p)
%  hcb= colorbar
%  title(hcb,'$C_p$','interpreter','latex')
% %  caxis([-1.0 1.0])
%  box on 
%  xlabel('x/D [-]')
%  ylabel('y/D [-]')
% % %  
% %  
%  
% Fig2=figure
% plot(R(1:500,443),CP(1:500,443))
% hold on
% plot(R(1:500,584),CP(1:500,584))
% xlabel('R/D [-]')
% ylabel('Cp [-]')
% legend('x/D = -0.04','x/D = 0.105')
% axis([0 0.6 -1.2 1])


Mask(Mask == 0) = NaN;
 Fig2=figure('units','centimeters','papersize',[24 12],'position',[3 5 24 12])
 hold on
 pcolor(X/D,R/D,Vx.*Mask/V_inf);
%  plot(Xw,Rw)
 axis equal
axis([-0.52 0.95 -0.2 0.7])
 shading flat
 colormap(COLORMAP_u)
hcb= colorbar
title(hcb,'$\bar{u}/U_{\infty}$','interpreter','latex')
 caxis([0 1.2])
 box on 
 xlabel('x/D [-]')
 ylabel('y/D [-]')

%   h = vline(x1,'k-') 
%   h = vline(x2,'k--')
%   h = vline(x3,'k--')
%   h = vline(x4,'k--') 
%   h = hline(r1,'r-')
%   h = hline(r2,'r-')
%   h = hline(r3,'r-')
%   
%   return 

 Fig3=figure('units','centimeters','papersize',[24 12],'position',[3 7 24 12])
 pcolor(X/D,R/D,Vr.*Mask/V_inf);
  axis equal
 axis([-0.52 0.95 -0.2 0.7])
 shading flat
 colormap(COLORMAP_v)
hcb= colorbar
title(hcb,'$\bar{v}/U_{\infty}$','interpreter','latex')
 caxis([-0.2 0.2])

 box on
 xlabel('x/D [-]')
 ylabel('y/D [-]')
 
 
% clear UC_VEC PC NORVEC 
% clear UC_VEC2 PC2  NORVEC2

% end 
end 




return 
figure
R=R./D;
hold on 
plot(R(1:680,end),ct(:,1:680))
% plot(R(1:600,end),Fx1(:,1:600))
grid on 
legend



for k= 1:4
Fig1=figure
hold on 
plot(0.3+0.05*[1:10],Ct_AMI_linear(k,:),'-','color',c(k,:))
plot(0.3+0.05*[1:10],Ct_AMI_annular(k,:),'--','color',c(k,:))

title(strcat('case ',num2str(k)))
legend('AMI linear','AMI annulus','Location','Best')
xlabel('x/D')
ylabel('$C_t$')
grid on 
box on 
printpdf(Fig1, ['Ct_AMI',num2str(k)])
end 

 
for kk= 1:4
Fig2=figure
hold on
plot(0.3+0.05*[1:10],CT_CMI_linear(kk,:),'-','color',c(kk,:))
plot(0.3+0.05*[1:10],CT_CMI_annular(kk,:),'--','color',c(kk,:))

legend('CMI linear','CMI annulus','Location','Best')
title(strcat('case ',num2str(kk)))
xlabel('x/D')
ylabel('$C_t$')
grid on 
box on 
printpdf(Fig2, ['Ct_CMI',num2str(kk)])
end 

return 

printpdf(Fig1, ['Ct_AMI',num2str(k)])
printpdf(Fig2, ['Ct_CMI',num2str(kk)])
