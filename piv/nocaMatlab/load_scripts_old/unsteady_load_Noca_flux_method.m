%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This script determin the loads from velocity fields and their derivatives using the3D Noca method to determine 
 % % % from 
% @article{NOCA1999,
% abstract = {We continue with the 1997 work of Noca et al. and offer some additional closed-form expressions (and their derivations) for the evaluation of time-dependent forces on a body in an incompressible, viscous, and rotational flow, which require only the knowledge of the velocity field (and its derivatives) in a finite and arbitrarily chosen region enclosing the body. In particular, we offer an expression for the force which only depends on the velocity field (and its derivatives) on thesurface of an arbitrary control volume. These expressions are particularly useful for experimental techniques like Digital Particle Image Velocimetry (DPIV) which provide time sequences of 2-D velocity fields but not pressure fields. For some common flow situations (freely moving objects, flexible bodies, flying and swimming animals, low Reynolds number flows, soap film tunnels), these techniques may be more viable than traditional methods (strain gages). The formulations can also be of some interest to the Computational Fluid Dynamics (CFD) community, especially when pressure is not evaluated explicitly, such as in vorticity-based algorithms. From a theoretical point of view, they provide an explicit relation between loading and flow structure. In the present work, the formulations are tested on a numerical flow simulation using a high-resolution vortex method and experimentally with DPIV on a circular cylinder flow. {\textcopyright} 1999 Academic Press.},
% author = {Noca, F. and Shiels, D. and Jeon, D.},
%  doi = {10.1006/jfls.1999.0219},
%  journal = {Journal of Fluids and Structures},
% month = {jul},
% number = {5},
%  pages = {551--578},
% title = {A COMPARISON OF METHODS FOR EVALUATING TIME-DEPENDENT FLUID DYNAMIC FORCES ON BODIES, USING ONLY VELOCITY FIELDS AND THEIR DERIVATIVES},
% url = {http://linkinghub.elsevier.com/retrieve/pii/S0889974699902190},
% volume = {13},
% year = {1999}
%  }
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



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
datadir.sub2 = '/5_0Hz';


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

%% Physical properties
load('initial_cond.mat')
P_inf = p_m800(2);
V_inf = V_m800(2);
Rho_inf = rho_m800(2);
D = 60/100; %[m]


Q_inf =(1/2)*Rho_inf*(V_inf.^2);
Area_disc= pi*(D/2)^2;

% % % should the force on the nacelle be substracted? 
cd_cube = 1.05; % % 1.05 or 2.05 I'm not so sure 
% % % The area is 6.5cm x5cm 
Area_nacelle = 6.5*5/100/100;   % % m^2 
Fd=cd_cube*(Q_inf*Area_nacelle);

Area_ratio=Area_nacelle/Area_disc;  % % % the ratio of nacelle area and disc area 

% % % % % matrixIn using a mean filter over a rectangle of size (2*Nr+1)-by-(2*Nc+1) 
Nr=5;
Nc=5;



%%% read data
for K=1:11





if K==11
    DATA2=load(fullfile(datadir.main,datadir.sub1,datadir.sub2,datafiles(1).name))
else 
    DATA2=load(fullfile(datadir.main,datadir.sub1,datadir.sub2,datafiles(K+1).name))
end 
    DATA2.Vx=smooth2a(DATA2.Vx,Nr,Nc);
    DATA2.Vr=smooth2a(DATA2.Vr,Nr,Nc);
    
    DATA2.U=DATA2.Vx(:,1:N_max);
    DATA2.V=DATA2.Vr(:,1:N_max);
    DATA2.X=DATA2.X(:,1:N_max);
    DATA2.R=DATA2.R(:,1:N_max);




if K==1
    DATA1=load(fullfile(datadir.main,datadir.sub1,datadir.sub2,datafiles(11).name))
else 
    DATA1=load(fullfile(datadir.main,datadir.sub1,datadir.sub2,datafiles(K-1).name))
end 

    DATA1.Vx=smooth2a(DATA1.Vx,Nr,Nc);
    DATA1.Vr=smooth2a(DATA1.Vr,Nr,Nc);
    
    DATA1.U=DATA1.Vx(:,1:N_max);
    DATA1.V=DATA1.Vr(:,1:N_max);
    DATA1.X=DATA1.X(:,1:N_max);
    DATA1.R=DATA1.R(:,1:N_max);
%     DATA1.R=flip(DATA1.R,1);

DUDT=(DATA2.U-DATA1.U)/DT;
DVDT=(DATA2.V-DATA1.V)/DT;

% 
% % % % %  for steady case, time derivative is 0. 
% DUDT=DATA.Vx*0.0;
% DVDT=DATA.Vx*0.0;


DATA=load(fullfile(datadir.main,datadir.sub1,datadir.sub2,datafiles(K).name))
    
    DATA.Vx=smooth2a(DATA.Vx,Nr,Nc);
    DATA.Vr=smooth2a(DATA.Vr,Nr,Nc);   

    DATA.U=DATA.Vx(:,1:N_max);
    DATA.V=DATA.Vr(:,1:N_max);
    DATA.X=DATA.X(:,1:N_max);
    DATA.R=DATA.R(:,1:N_max);

    
% outfilename=['GRID6_CT077_k1_axi0' '_subset'];
% load(outfilename)
% 
% DUDT=(DATA.Vx-DATA.Vx)/DT;
% DVDT=(DATA.V-DATA.V)/DT;



y1=0:.001:0.65;
x1=[-0.36:0.001:0.1];  % % % % % tests show that the load is very sensitive to
% the right boundary.

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



DATA.U(1,:)=DATA.U(2,:);
DATA.V(1,:)=DATA.V(2,:).*0;


DUDT(1,:)=DATA.U(2,:)*0;
DVDT(1,:)=DATA.V(2,:).*0;




[OMEGAZ,cav]= curl(DATA.X,DATA.R,DATA.U,DATA.V);
VORTICITY(3,:) =griddata(DATA.X,DATA.R,OMEGAZ,xline,yline);
% VORTICITY(3,:) =interp2(DATA.X,DATA.R,OMEGAZ,xline,yline);
VORTICITY(2,:)=VORTICITY(3,:).*0; 
VORTICITY(1,:)=VORTICITY(2,:);

[ddx_u,ddy_u] = gradient(DATA.U,0.01,0.01);
[ddx_v,ddy_v] = gradient(DATA.V,0.01,0.01);

[ddxddx_u,ddxddy_u] = gradient(ddx_u,0.01,0.01);
[ddxddy_u,ddyddy_u] = gradient(ddy_u,0.01,0.01);

[ddxddx_v,ddxddy_v] = gradient(ddx_v,0.01,0.01);
[ddxddy_v,ddyddy_v] = gradient(ddy_v,0.01,0.01);

UC_VEC(1,:)=griddata(DATA.X,DATA.R,DATA.U,xline,yline);
UC_VEC(2,:)=griddata(DATA.X,DATA.R,DATA.V,xline,yline);
UC_VEC(3,:)=UC_VEC(2,:).*0;

% DUDT_VEC2(1,:)=interp2(DATA.X,DATA.R,DUDT,xline,yline);
% DUDT_VEC2(2,:)=interp2(DATA.X,DATA.R,DVDT,xline,yline);
DUDT_VEC2(1,:)=griddata(DATA.X,DATA.R,DUDT,xline,yline);
DUDT_VEC2(2,:)=griddata(DATA.X,DATA.R,DVDT,xline,yline);
DUDT_VEC2(3,:)=DUDT_VEC2(2,:).*0;



PARTIAL1st(1,:)=griddata(DATA.X,DATA.R,ddx_u,xline,yline);
PARTIAL1st(2,:)=griddata(DATA.X,DATA.R,ddx_v,xline,yline);
PARTIAL1st(3,:)=griddata(DATA.X,DATA.R,ddy_u,xline,yline);
PARTIAL1st(4,:)=griddata(DATA.X,DATA.R,ddy_v,xline,yline);

PARTIAL2nd(1,:)=griddata(DATA.X,DATA.R,ddxddx_u,xline,yline);
PARTIAL2nd(2,:)=griddata(DATA.X,DATA.R,ddxddy_v,xline,yline);
PARTIAL2nd(3,:)=griddata(DATA.X,DATA.R,ddyddy_u,xline,yline);
PARTIAL2nd(4,:)=griddata(DATA.X,DATA.R,ddxddx_v,xline,yline);
PARTIAL2nd(5,:)=griddata(DATA.X,DATA.R,ddxddy_u,xline,yline);
PARTIAL2nd(6,:)=griddata(DATA.X,DATA.R,ddyddy_v,xline,yline);

% return 
[UVEC, PARTIALVEC] = convert_from_cylindrical_to_cartesian(UC_VEC,PARTIAL1st,PARTIAL2nd);
[DUDTVEC2, PARTIALVEC] = convert_from_cylindrical_to_cartesian(DUDT_VEC2,PARTIAL1st,PARTIAL2nd);

%DUDT=(UVEC2-UVEC)/DT;


[FLUX, F1, F2, F3, F4, F5, F6] = calculate_flux_3D(...
    [xline;yline;zline],UVEC,DUDTVEC2,VORTICITY,PARTIALVEC,1.8e-5,NORVEC);


FLUXMID= (FLUX(:,1:end-1)+FLUX(:,2:end))/2;

FAREA=FLUXMID.*AREA_MIDPOINT;
LOAD=sum(FAREA,2)
CT=sum(FAREA,2)/(0.5*pi*0.5^2)
LOADF1=sum((F1(:,1:end-1)+F1(:,2:end))/2.*AREA_MIDPOINT,2)/(0.5*pi*0.5^2);
LOADF2=sum((F2(:,1:end-1)+F2(:,2:end))/2.*AREA_MIDPOINT,2)/(0.5*pi*0.5^2);
LOADF3=sum((F3(:,1:end-1)+F3(:,2:end))/2.*AREA_MIDPOINT,2)/(0.5*pi*0.5^2);
LOADF4=sum((F4(:,1:end-1)+F4(:,2:end))/2.*AREA_MIDPOINT,2)/(0.5*pi*0.5^2);
LOADF5=sum((F5(:,1:end-1)+F5(:,2:end))/2.*AREA_MIDPOINT,2)/(0.5*pi*0.5^2);
LOADF6=sum((F6(:,1:end-1)+F6(:,2:end))/2.*AREA_MIDPOINT,2)/(0.5*pi*0.5^2);
LDS= [LOADF1 LOADF2 LOADF3 LOADF4 LOADF5 LOADF6] 
% 
% figure
% hold on
% pcolor(DATA.X,DATA.R,OMEGAZ)
% shading flat
% colorbar
% quiver(DATA.X,DATA.R,DATA.U,DATA.V)
% plot(xline,yline,'g--')
% daspect([1 1 1])



CTARRAY(K)=CT(1);
end



% CTARRAY=CTARRAY(2:8);

% CTNUM= 7/9 + 1/9*sin(pi*1*2*([1:1:7]* 0.04+22)) ;
figure
hold on
plot(CTARRAY, 'bo -' )
% plot(CTNUM, 'r+-' )
legend ('Noca', 'CFD' )




% 
% figure
% hold on
% plot(yline,UVEC(1,:))
% 
% 
% 
% figure
% title('vorticity')
% hold on
% plot(yline,VORTICITY(3,:),'+-')









