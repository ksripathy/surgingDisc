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

caseDir = "Z:\home\ksripathy\phd\surgingDisc\piv\data\matlab\p70Case7.mat";

%staticCt(caseDir,165,165,200,15)

minHOffsetArr = 200:300;
CtArr = zeros(3,length(minHOffsetArr));

for i=1:length(minHOffsetArr)

    [CtArr(:,i)] = staticCt(caseDir,165,165,minHOffsetArr(1,i),15);

end

disp([mean(CtArr(1,:)),std(CtArr(1,:))])


%%Functions
function [CT] = staticCt(dir, minVOffset, maxVOffset, minHOffset, maxHoffset)
    dataDir = dir;
    
    [COLORMAP_u,COLORMAP_v,COLORMAP_w,COLORMAP_p]=self_colormap();
    
    P_inf = 1.017e5;
    D = 0.2;
    
    Area_disc= pi*(D/2)^2;
    
    dataParent = load(dataDir);
    Rho_inf = dataParent.rhoInf;
    deltaX = dataParent.deltaX;
    deltaY = dataParent.deltaY;
    
    DATA = dataParent.phase1;
    
    V_inf = DATA.Vinf;
    %V_inf = 3.096507;
    Q_inf =(1/2)*Rho_inf*(V_inf.^2);
    
    DATA.U=DATA.Vx;
    DATA.V=DATA.Vr;
    DATA.X=DATA.X;
    DATA.R=DATA.Y;
    
    % % % % %  for steady case, time derivative is 0. 
    DUDT=DATA.U*0.0;
    DVDT=DATA.V*0.0;
    
    %Index corresponding to disc centre
    rc = dataParent.rcIndex;
    xc = dataParent.xcIndex;
    
    %Control volume bound in indices
    rMax = rc-maxVOffset; %Lowest index corresponds to highest position for radial direction
    rMin = rc+minVOffset;
    xMin = xc-minHOffset;
    xMax = xc+maxHoffset;
    xSize = size(xMin:xMax);
    rSize = size(rMax:rMin);
    
    %Indices spanning the boundary of control volume
    rIndices = [rMin:-1:rMax ones(xSize,"int64")*rMax rMax:rMin ones(xSize,"int64")*rMin];
    xIndices = [ones(rSize,"int64")*xMin xMin:xMax ones(rSize,"int64")*xMax xMax:-1:xMin];
    
    %Coordinates corresponding to control volume boundaries
    xline = DATA.X(rc,xIndices);
    yline = DATA.Y(rIndices,xc).';
    zline = zeros(size(xline));
    
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
    DVDT(1,:)=DATA.V(2,:)*0;
    
    [OMEGAZ,cav]= curl(DATA.X,DATA.R,DATA.U,DATA.V);
    
    VORTICITY = zeros(3,length(xline));
    
    for i=1:length(xline)
        
        VORTICITY(3,i) = OMEGAZ(rIndices(1,i),xIndices(1,i));
    
    end
    
    [ddx_u,ddy_u] = gradient(DATA.U,deltaX,deltaY);
    [ddx_v,ddy_v] = gradient(DATA.V,deltaX,deltaY);
    
    [ddxddx_u,ddxddy_u] = gradient(ddx_u,deltaX,deltaY);
    [ddxddy_u,ddyddy_u] = gradient(ddy_u,deltaX,deltaY);
    
    [ddxddx_v,ddxddy_v] = gradient(ddx_v,deltaX,deltaY);
    [ddxddy_v,ddyddy_v] = gradient(ddy_v,deltaX,deltaY);
    
    UC_VEC = zeros(3,length(xline));
    DUDT_VEC2 = zeros(3,length(xline));
    PARTIAL1st = zeros(4,length(xline));
    PARTIAL2nd = zeros(6,length(xline));
    
    for i=1:length(xline)
    
        UC_VEC(1,i) = DATA.U(rIndices(1,i),xIndices(1,i));
        UC_VEC(2,i) = DATA.V(rIndices(1,i),xIndices(1,i));
    
        DUDT_VEC2(1,i) = DUDT(rIndices(1,i),xIndices(1,i));
        DUDT_VEC2(2,i) = DVDT(rIndices(1,i),xIndices(1,i));
    
        PARTIAL1st(1,i) = ddx_u(rIndices(1,i),xIndices(1,i));
        PARTIAL1st(2,i) = ddx_v(rIndices(1,i),xIndices(1,i));
        PARTIAL1st(3,i) = ddy_u(rIndices(1,i),xIndices(1,i));
        PARTIAL1st(4,i) = ddy_v(rIndices(1,i),xIndices(1,i));
    
        PARTIAL2nd(1,i) = ddxddx_u(rIndices(1,i),xIndices(1,i));
        PARTIAL2nd(2,i) = ddxddy_v(rIndices(1,i),xIndices(1,i));
        PARTIAL2nd(3,i) = ddyddy_u(rIndices(1,i),xIndices(1,i));
        PARTIAL2nd(4,i) = ddxddx_v(rIndices(1,i),xIndices(1,i));
        PARTIAL2nd(5,i) = ddxddy_u(rIndices(1,i),xIndices(1,i));
        PARTIAL2nd(6,i) = ddyddy_u(rIndices(1,i),xIndices(1,i));
    
    end
    
    [UVEC, PARTIALVEC] = convert_from_cylindrical_to_cartesian(UC_VEC,PARTIAL1st,PARTIAL2nd);
    [DUDTVEC2, PARTIALVEC] = convert_from_cylindrical_to_cartesian(DUDT_VEC2,PARTIAL1st,PARTIAL2nd);
    
    [FLUX, F1, F2, F3, F4, F5, F6] = calculate_flux_3D(...
        [xline;yline;zline],UVEC,DUDTVEC2,VORTICITY,PARTIALVEC,1.8e-5,NORVEC);
    
    
    FLUXMID= (FLUX(:,1:end-1)+FLUX(:,2:end))/2;

    FAREA=FLUXMID.*AREA_MIDPOINT;
    
    LOAD = sum(FAREA,2);
    CT = LOAD/(0.5*Rho_inf*V_inf^2*Area_disc);
    
    % fig = figure;
    % ax = gca;
    % hold on
    % pcolor(DATA.X,DATA.R,OMEGAZ)
    % shading flat
    % colorbar
    % ax.CLim = [-125,125];
    % %quiver(DATA.X,DATA.R,DATA.U,DATA.V)
    % plot(xline,yline,'g--')
    % daspect([1 1 1])

end





