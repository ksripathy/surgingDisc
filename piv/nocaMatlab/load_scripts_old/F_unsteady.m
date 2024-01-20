% % This is the script to calculate the force field 

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
%% unsteady path 
datadir.sub1  = '/Unsteady';
datadir.sub2 = '/1_5Hz';

% % % determine the period time 
if datadir.sub2 == '/5_0Hz'
    T=1/5.0;
    N_max=2348;
elseif datadir.sub2 == '/3_0Hz'
    T=1/3.0;
    N_max=2358;
elseif datadir.sub2 == '/1_5Hz'
    T=1/1.5;
    N_max=2348;
end 
delta_t=T/10; 
DT=2*delta_t;


datadir.file = 'F*.*';
datafiles = dir(fullfile(datadir.main,datadir.sub1,datadir.sub2,datadir.file));


[COLORMAP_u,COLORMAP_v,COLORMAP_w,COLORMAP_p]=self_colormap();

% %% % %  initial conditions, using the average value from all the cases 
% Rho_inf =1.2065; % % kg/m^3
% % V_inf=5.89;   % % m/s 
% T=19.52+273.15;  % % K 
%  R_const= 8.314;  % %  % Pa*m^3*mol^-1*K^-1
%  m_mol=29/1000; % % the mol mass of air kg/mol
% P_inf=Rho_inf*R_const*T/m_mol; % % % Pa  ideal gas law 

%% Physical properties
load('initial_cond.mat')
P_inf = p_m800(2);   % % This absolute value of P_inf doesn't change the solution, only the pressure gradient matters, which is determined by the vlocity field. 
V_inf = V_m800(2);
Rho_inf = rho_m800(2);
D = 60/100; %[m]


Q_inf =(1/2)*Rho_inf*(V_inf.^2);
Area_disc= pi*(D/2)^2;

% % % % should the force on the nacelle be substracted? 
% cd_cube = 1.05; % % 1.05 or 2.05 I'm not so sure 
% % % % The area is 6.5cm x5cm 
% Area_nacelle = 6.5*5/100/100;   % % m^2 
% Fd=cd_cube*(Q_inf*Area_nacelle);
% 
% Area_ratio=Area_nacelle/Area_disc;  % % % the ratio of nacelle area and disc area 
% 



for K=   1%:11  
 

if K>1
    DATA1=load(fullfile(datadir.main,datadir.sub1,datadir.sub2,datafiles(K-1).name));
    DATA1.U=DATA1.Vx(:,1:N_max)*V_inf;
    DATA1.V=DATA1.Vr(:,1:N_max)*V_inf;
    DATA1.X=DATA1.X(:,1:N_max)*D;
    DATA1.R=DATA1.R(:,1:N_max)*D;
end 



if K<11
    DATA2=load(fullfile(datadir.main,datadir.sub1,datadir.sub2,datafiles(K+1).name)); 
    DATA2.U=DATA2.Vx(:,1:N_max)*V_inf;
    DATA2.V=DATA2.Vr(:,1:N_max)*V_inf;
    DATA2.X=DATA2.X(:,1:N_max)*D;
    DATA2.R=DATA2.R(:,1:N_max)*D; 
end 
 



DATA=load(fullfile(datadir.main,datadir.sub1,datadir.sub2,datafiles(K).name));
DATA.U=DATA.Vx(:,1:N_max)*V_inf;
DATA.V=DATA.Vr(:,1:N_max)*V_inf;
DATA.X=DATA.X(:,1:N_max)*D;
DATA.R=DATA.R(:,1:N_max)*D;
    
    
Mask=ones(size(DATA.X));

%% Wei's new mask 
% if K==1  %%  if (datadir.sub2 == '/3_0Hz') &&  (K==1) 
% Mask(:,1:1200) = 1-imdilate(imclose(DATA.U(:,1:1200)<0,strel('disk',20)),strel('disk',10)); % creating a mask on the disk and reflection region
% else 
Mask(:,1:1200) = 1-imdilate(imclose(DATA.U(:,1:1200)==0,strel('disk',20)),strel('disk',10)); % creating a mask on the disk and reflection region
% end 





% % % to ch%eck the if the mask only cover the disk 
% A=imclose(Vx==0,strel('disk',20));
% B=strel('disk',10);
% C=imdilate(A,B);
% figure
% imshow(A)
% figure
% imshow(C)
% figure
% imshow(1-Mask)
%     
%     
   

    %%  Navier Stokes equations

dx=DATA.X(1,2)-DATA.X(1,1);
dr=DATA.R(2,1)-DATA.R(1,1);


%% % % % For the unsteady case, dUdt != 0   calculate dUdt= (U(phase+1)-U(phase-1))/(2*Delta_T )
if K==1
DUDT=(DATA2.U-DATA.U)/DT/2;
DVDT=(DATA2.V-DATA.V)/DT/2;
elseif K==11
DUDT=(DATA.U-DATA1.U)/DT/2;
DVDT=(DATA.V-DATA1.V)/DT/2;
else 
DUDT=(DATA2.U-DATA1.U)/DT;
DVDT=(DATA2.V-DATA1.V)/DT;
end 



[dUdx,dUdr]=gradient(DATA.U,dx,dr);
[dVdx,dVdr]=gradient(DATA.V,dx,dr);


dPdx=-Rho_inf*(DUDT+DATA.U.*dUdx+DATA.V.*dUdr);
dPdr=-Rho_inf*(DVDT+DATA.U.*dVdx+DATA.V.*dVdr);

[dP2dx2,~]=gradient(dPdx,dx,dr);
[~,dP2dr2]=gradient(dPdr,dx,dr);
   
%% Pressure from Bernoulli equation
PBern=P_inf+(1/2)*Rho_inf*(V_inf.^2-(DATA.U.^2+DATA.V.^2)); % Pressure using Bernoulli Equation
PBern=PBern.*Mask;
CPBern = (PBern-P_inf).*Mask/Q_inf;


%% by applying the mask 
P=Poisson_plane_nw(Mask,dPdx,dP2dx2,dPdr,dP2dr2,PBern,dx,dr,'Dir','Nx','Dir','Nx');

%% split the domain 
% % % calculated pressure field from the Poisson equation, left and top are
% % % Dirichlet boundary 
% P=zeros(size(DATA.X));
% ind = find(DATA.X(1,:)<0, 1, 'last');
% % ind1=find(sum(DATA.U(:,1:1000)==0,1)>0);
% ind_up=ind-15;
% ind_down=ind+30;
% % num_skip=15; % % % the minimal num_skip=5
% P(:,1:ind_up)=Poisson_plane_nw(Mask(:,1:ind_up),dPdx(:,1:ind_up),dP2dx2(:,1:ind_up),dPdr(:,1:ind_up),dP2dr2(:,1:ind_up),PBern(:,1:ind_up),dx,dr,'Dir','Nx','Dir','Nx');
% P(:,ind_down:end)=Poisson_plane_nw(Mask(:,ind_down:end),dPdx(:,ind_down:end),dP2dx2(:,ind_down:end),dPdr(:,ind_down:end),dP2dr2(:,ind_down:end),PBern(:,ind_down:end),dx,dr,'Nx','Nx','Dir','Nx');   

CP = (P-P_inf).*Mask/Q_inf;


%  Fig1=figure('units','centimeters','papersize',[24 12],'position',[3 5 24 12])
%  hold on
%  pcolor(X./D,R./D,CP);
%  axis equal
%  axis([-0.52 2.0 -0.2 0.7])
%  shading flat
%  colormap(COLORMAP_p)
%  hcb= colorbar
%  title(hcb,'$C_p$','interpreter','latex')
%  caxis([-1.0 1.0])
%  box on 
%  xlabel('x/D [-]')
%  ylabel('y/D [-]')

%% % % % generate the contour for momentum integration  

% x_sta=[-0.1];
x_end=[1.0];%[0.3:0.1:1.5];
for ii=1:length(x_end)
x1=[-0.1:0.001:x_end(ii)]*D;
y1=[0:.001:0.65]*D;

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

% %  % velocity and pressure at the contour 
UC_VEC(1,:)=griddata(DATA.X,DATA.R,DATA.U,xline,yline);
UC_VEC(2,:)=griddata(DATA.X,DATA.R,DATA.V,xline,yline);
UC_VEC(3,:)=0;
PC=griddata(DATA.X,DATA.R,P,xline,yline);

%% % % unsteady term, surface integration 
if K==1    
UC_VEC_tp1(1,:)=griddata(DATA2.X,DATA2.R,DATA2.U,xline,yline);
UC_VEC_tp1(2,:)=griddata(DATA2.X,DATA2.R,DATA2.V,xline,yline);
UC_VEC_tp1(3,:)=0;
UC_VEC_tm1(1,:)=griddata(DATA.X,DATA.R,DATA.U,xline,yline);
UC_VEC_tm1(2,:)=griddata(DATA.X,DATA.R,DATA.V,xline,yline);
UC_VEC_tm1(3,:)=0;
F_CMI=F_cal_CMI_unsteady([xline;yline;zline],UC_VEC,PC,NORVEC,UC_VEC_tp1,UC_VEC_tm1,DT/2);
elseif K==11
UC_VEC_tp1(1,:)=griddata(DATA.X,DATA.R,DATA.U,xline,yline);
UC_VEC_tp1(2,:)=griddata(DATA.X,DATA.R,DATA.V,xline,yline);
UC_VEC_tp1(3,:)=0;
UC_VEC_tm1(1,:)=griddata(DATA1.X,DATA1.R,DATA1.U,xline,yline);
UC_VEC_tm1(2,:)=griddata(DATA1.X,DATA1.R,DATA1.V,xline,yline);
UC_VEC_tm1(3,:)=0;
F_CMI=F_cal_CMI_unsteady([xline;yline;zline],UC_VEC,PC,NORVEC,UC_VEC_tp1,UC_VEC_tm1,DT/2);    
else 
UC_VEC_tp1(1,:)=griddata(DATA2.X,DATA2.R,DATA2.U,xline,yline);
UC_VEC_tp1(2,:)=griddata(DATA2.X,DATA2.R,DATA2.V,xline,yline);
UC_VEC_tp1(3,:)=0;
UC_VEC_tm1(1,:)=griddata(DATA1.X,DATA1.R,DATA1.U,xline,yline);
UC_VEC_tm1(2,:)=griddata(DATA1.X,DATA1.R,DATA1.V,xline,yline);
UC_VEC_tm1(3,:)=0;
F_CMI=F_cal_CMI_unsteady([xline;yline;zline],UC_VEC,PC,NORVEC,UC_VEC_tp1,UC_VEC_tm1,DT);       
end 
        

F_CMI_MID=(F_CMI(:,1:end-1)+F_CMI(:,2:end))/2;


FAREA=F_CMI_MID.*AREA_MIDPOINT;
LOAD=sum(FAREA,2);

% % %% % % % the volume integration of the unsteady term. 
% % [XL,RL] = meshgrid(x1, y1);
% % U2 = griddata(DATA2.X,DATA2.R,DATA2.U,XL,RL);
% % V2 = griddata(DATA2.X,DATA2.R,DATA2.V,XL,RL);
% % W2=zeros(size(U2));
% % U1 = griddata(DATA1.X,DATA1.R,DATA1.U,XL,RL);
% % V1 = griddata(DATA1.X,DATA1.R,DATA1.V,XL,RL);
% % W1=zeros(size(U1));
% % F_vol=CMI_unsteady_term(U2,V2,W2,U1,V1,W1,XL,RL,DT);
% % CT=(F_vol+sum(FAREA,2))/(0.5*Rho_inf*pi*(D/2)^2*V_inf^2);


CT=sum(FAREA,2)/(0.5*Rho_inf*pi*(D/2)^2*V_inf^2);
% CT_CMI(K)=CT(1);

CT_CMI(K,ii)=CT(1);
% clear UC_VEC PC NORVEC 
clear UC_VEC PC NORVEC UC_VEC_tp1  UC_VEC_tm1
end 

%% % plotting 
 
%  Fig1=figure('units','centimeters','papersize',[24 12],'position',[3 5 24 12])
%  hold on
%  pcolor(DATA.X/D,DATA.R/D,CP);
%  axis equal
%  axis([-0.52 1.8 -0.2 0.7])
%  shading flat
%  colormap(COLORMAP_p)
%  hcb= colorbar
%  title(hcb,'$C_p$','interpreter','latex')
%  caxis([-1.0 1.0])
%  box on 
%  xlabel('x/D [-]')
%  ylabel('y/D [-]')
%  
%  
%  
% Fig2=figure
% plot(R(1:500,443),CP(1:500,443))
% hold on
% plot(R(1:500,584),CP(1:500,584))
% xlabel('R/D [-]')
% ylabel('Cp [-]')
% legend('x/D = -0.04','x/D = 0.105')
% axis([0 0.6 -1.2 1])
% printpdf(Fig1, ['Cp',datadir.sub2(2:end),num2str(K)])
% clear UC_VEC PC NORVEC UC_VEC_tp1  UC_VEC_tm1
% clear UC_VEC2 PC2  NORVEC2 UC_VEC2_tp1 UC_VEC2_tm1
% end 
end 


CT_CMI=CT_CMI';


