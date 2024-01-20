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

caseDir = "Z:\home\ksripathy\phd\surgingDisc\piv\data\matlab\p70Case6.mat";

%staticCt1(caseDir,0,165,150,75) + staticCt2(caseDir,165,0,150,100)

minHOffsetArr = 50:5:350;
maxHOffsetArr = 15:5:315;
CtArr1 = zeros(3,length(minHOffsetArr));
CtArr2 = zeros(3,length(minHOffsetArr));

for i=1:length(minHOffsetArr)

    [CtArr1(:,i)] = staticCt1(caseDir,165,165,minHOffsetArr(1,i),75) + staticCt2(caseDir,165,165,minHOffsetArr(1,i),75);

end

for i=1:length(maxHOffsetArr)

    %[CtArr2(:,i)] = staticCt1(caseDir,165,165,100,maxHOffsetArr(1,i)) + staticCt2(caseDir,165,165,100,maxHOffsetArr(1,i));

end

disp([mean(CtArr1(1,:)),std(CtArr1(1,:))])
disp([mean(CtArr2(1,:)),std(CtArr2(1,:))])


%%Functions
function [CT] = staticCt1(dir, minVOffset, maxVOffset, minHOffset, maxHoffset)
    dataDir = dir;
    
    [COLORMAP_u,COLORMAP_v,COLORMAP_w,COLORMAP_p]=self_colormap();
    
    P_inf = 1.017e5;
    D = 0.2;
    
    Area_disc= pi*(D/2)^2;
    
    dataParent = load(dataDir);
    %Rho_inf = dataParent.rhoInf;
    Rho_inf = 1;
    deltaX = dataParent.deltaX;
    deltaY = dataParent.deltaY;
    
    DATA = dataParent.phase1;
    
    V_inf = DATA.Vinf;
    %V_inf = 3.096507
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
    xc = dataParent.xcIndex(1);
    %disp(xc)
    %disp(DATA.X(rc,xc))
    
    %Control volume bound in indices
    rMax = rc-maxVOffset; %Lowest index corresponds to highest position for radial direction
    rMin = rc+minVOffset;
    xMin = xc-minHOffset;
    xMax = xc+maxHoffset;
    %disp(xMin)
    %disp(xMax)
    xSize = size(xMin:xMax);
    rSize = size(rMax:rMin);
    %disp(DATA.X(rc,xMin))
    %disp(DATA.X(rc,xMax))

    x1=DATA.X(rc,xMin):0.001:DATA.X(rc,xMax);
    y1=DATA.Y(rMin,xc):0.001:DATA.Y(rMax,xc);

    xline= [zeros(size(y1))+min(x1) x1 zeros(size(y1))+max(x1) ];  
    yline= [ y1 zeros(size(x1))+max(y1)  fliplr(y1) ]; 
    zline = yline.*0;
    %disp(length(xline))
    AREA_MIDPOINT = [pi*(y1(2:end).^2 - y1(1:end-1).^2) 0 2*pi*(max(y1) - min(y1))*(x1(2:end) - x1(1:end-1)) 0 fliplr(pi*(y1(2:end).^2 - y1(1:end-1).^2))];
    %disp(AREA_MIDPOINT);
    
    dyx=sqrt((yline(2:end) - yline(1:end-1)).^2+(xline(2:end) - xline(1:end-1)).^2);
    
    %AREA_MIDPOINT = 2*pi*( yline(2:end)+yline(1:end-1))/2.*dyx;
    %disp(AREA_MIDPOINT);
    %disp(length(AREA_MIDPOINT))
    SURFACEAREA=sum(AREA_MIDPOINT);
    CROSSAREA = 0.25 * pi * (DATA.Y(rMax,xc) - DATA.Y(rMin,xc))^2;
    
    ny=-gradient(xline);
    nx=gradient(yline);
    normvec=sqrt(nx.^2+ny.^2);
    NORVEC=-[nx;ny;nx.*0]./normvec;
    
    DATA.U(1,:)=DATA.U(2,:);
    DATA.V(1,:)=DATA.V(2,:).*0;
    
    DUDT(1,:)=DATA.U(2,:)*0;
    DVDT(1,:)=DATA.V(2,:)*0;
    
    [OMEGAZ,cav]= curl(DATA.X,DATA.R,DATA.U,DATA.V);
    %OMEGAZ = DATA.Vortz;
    
    VORTICITY = zeros(3,length(xline));
    
    VORTICITY(3,:) =griddata(DATA.X,DATA.R,OMEGAZ,xline,yline);
    
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

    UC_VEC(1,:)=griddata(DATA.X,DATA.R,DATA.U,xline,yline);
    UC_VEC(2,:)=griddata(DATA.X,DATA.R,DATA.V,xline,yline);

    DUDT_VEC2(1,:)=griddata(DATA.X,DATA.R,DUDT,xline,yline);
    DUDT_VEC2(2,:)=griddata(DATA.X,DATA.R,DVDT,xline,yline);

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
    
    [UVEC, PARTIALVEC] = convert_from_cylindrical_to_cartesian(UC_VEC,PARTIAL1st,PARTIAL2nd);
    [DUDTVEC2, PARTIALVEC] = convert_from_cylindrical_to_cartesian(DUDT_VEC2,PARTIAL1st,PARTIAL2nd);
    
    [FLUX, F1, F2, F3, F4, F5, F6] = calculate_flux_3D(...
        [xline;yline;zline],UVEC,DUDTVEC2,VORTICITY,PARTIALVEC,1.8e-5,NORVEC);
    
    
    FLUXMID= (FLUX(:,1:end-1)+FLUX(:,2:end))/2;
    %disp(length(FLUXMID))
    FAREA = FLUXMID.*abs(AREA_MIDPOINT);
    
    %LOAD = sum(FLUXMID.*dyx,2)*pi*D/2;
    LOAD = sum(FAREA,2);
    %CT = LOAD/(0.5*Rho_inf*V_inf^2*Area_disc);
    CT = LOAD/(0.5*Rho_inf*V_inf^2*Area_disc);
    
    % fig = figure
    % ax = gca;
    % hold on
    % pcolor(DATA.X,DATA.R,OMEGAZ)
    % shading flat
    % colorbar
    % ax.CLim = [-125,125];
    % %quiver(DATA.X,DATA.R,DATA.U,DATA.V)
    % quiver(DATA.X(1:25:end,1:25:end),DATA.R(1:25:end,1:25:end),DATA.U(1:25:end,1:25:end),DATA.V(1:25:end,1:25:end))
    % plot(xline,yline,'g--')
    % plot(xline,yline,'g--')
    % daspect([1 1 1])

end

function [CT] = staticCt2(dir, minVOffset, maxVOffset, minHOffset, maxHoffset)
    dataDir = dir;
    
    [COLORMAP_u,COLORMAP_v,COLORMAP_w,COLORMAP_p]=self_colormap();
    
    P_inf = 1.017e5;
    D = 0.2;
    
    Area_disc= pi*(D/2)^2;
    
    dataParent = load(dataDir);
    %Rho_inf = dataParent.rhoInf;
    Rho_inf = 1;
    deltaX = dataParent.deltaX;
    deltaY = dataParent.deltaY;
    
    DATA = dataParent.phase1;
    
    V_inf = DATA.Vinf;
    %V_inf = 3.096507
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
    xc = dataParent.xcIndex(1);
    
    %Control volume bound in indices
    rMax = rc-maxVOffset; %Lowest index corresponds to highest position for radial direction
    rMin = rc+minVOffset;
    xMin = xc-minHOffset;
    xMax = xc+maxHoffset;
    xSize = size(xMin:xMax);
    rSize = size(rMax:rMin);

    x1=DATA.X(rc,xMin):0.001:DATA.X(rc,xMax);
    y1=DATA.Y(rMin,xc):0.001:DATA.Y(rMax,xc);

    xline= [zeros(size(y1))+max(x1) fliplr(x1) zeros(size(y1))+min(x1) ];  
    yline= [ fliplr(y1) zeros(size(x1))+min(y1)  y1 ]; 
    zline = yline.*0;
    %disp(length(xline))
    AREA_MIDPOINT = [pi*(y1(2:end).^2 - y1(1:end-1).^2) 0 2*pi*(max(y1) - min(y1))*(x1(2:end) - x1(1:end-1)) 0 fliplr(pi*(y1(2:end).^2 - y1(1:end-1).^2))];
    %disp(AREA_MIDPOINT);
    
    dyx=sqrt((yline(2:end) - yline(1:end-1)).^2+(xline(2:end) - xline(1:end-1)).^2);
    
    %AREA_MIDPOINT = 2*pi*( yline(2:end)+yline(1:end-1))/2.*dyx;
    %disp(AREA_MIDPOINT);
    %disp(length(AREA_MIDPOINT))
    SURFACEAREA=sum(AREA_MIDPOINT);
    CROSSAREA = 0.25 * pi * (DATA.Y(rMax,xc) - DATA.Y(rMin,xc))^2;
    
    ny=-gradient(xline);
    nx=gradient(yline);
    normvec=sqrt(nx.^2+ny.^2);
    NORVEC=-[nx;ny;nx.*0]./normvec;
    
    DATA.U(1,:)=DATA.U(2,:);
    DATA.V(1,:)=DATA.V(2,:).*0;
    
    DUDT(1,:)=DATA.U(2,:)*0;
    DVDT(1,:)=DATA.V(2,:)*0;
    
    [OMEGAZ,cav]= curl(DATA.X,DATA.R,DATA.U,DATA.V);
    %OMEGAZ = DATA.Vortz;
    
    VORTICITY = zeros(3,length(xline));
    
    VORTICITY(3,:) =griddata(DATA.X,DATA.R,OMEGAZ,xline,yline);
    
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

    UC_VEC(1,:)=griddata(DATA.X,DATA.R,DATA.U,xline,yline);
    UC_VEC(2,:)=griddata(DATA.X,DATA.R,DATA.V,xline,yline);

    DUDT_VEC2(1,:)=griddata(DATA.X,DATA.R,DUDT,xline,yline);
    DUDT_VEC2(2,:)=griddata(DATA.X,DATA.R,DVDT,xline,yline);

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
    
    [UVEC, PARTIALVEC] = convert_from_cylindrical_to_cartesian(UC_VEC,PARTIAL1st,PARTIAL2nd);
    [DUDTVEC2, PARTIALVEC] = convert_from_cylindrical_to_cartesian(DUDT_VEC2,PARTIAL1st,PARTIAL2nd);
    
    [FLUX, F1, F2, F3, F4, F5, F6] = calculate_flux_3D(...
        [xline;yline;zline],UVEC,DUDTVEC2,VORTICITY,PARTIALVEC,1.5e-5,NORVEC);
    
    
    FLUXMID= (FLUX(:,1:end-1)+FLUX(:,2:end))/2;
    %disp(length(FLUXMID))
    FAREA = FLUXMID.*abs(AREA_MIDPOINT);
    
    %LOAD = sum(FLUXMID.*dyx,2)*pi*D/2;
    LOAD = sum(FAREA,2);
    %CT = LOAD/(0.5*Rho_inf*V_inf^2*Area_disc);
    CT = LOAD/(0.5*Rho_inf*V_inf^2*Area_disc);
    
    % fig = figure
    % ax = gca;
    % hold on
    % pcolor(DATA.X,DATA.R,OMEGAZ)
    % shading flat
    % colorbar
    % ax.CLim = [-125,125];
    % %quiver(DATA.X,DATA.R,DATA.U,DATA.V)
    % quiver(DATA.X(1:25:end,1:25:end),DATA.R(1:25:end,1:25:end),DATA.U(1:25:end,1:25:end),DATA.V(1:25:end,1:25:end))
    % plot(xline,yline,'g--')
    % plot(xline,yline,'g--')
    % daspect([1 1 1])

end





