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

caseDir = "Z:\home\ksripathy\phd\surgingDisc\piv\data\matlab\p70Case5.mat";

dynCT = dynamicCt(caseDir,165,165,150,150);
save("results2\p70Case5.mat","dynCT","-v7")

% minHOffsetArr = 50:5:350;
% CtArr = zeros(3,12,length(minHOffsetArr));
% 
% for i=1:length(minHOffsetArr)
% 
%     [CtArr(:,:,i)] = dynamicCt(caseDir,165,165,minHOffsetArr(1,i),25);
% 
% end
% 
% dynCt = zeros(2,12);
% 
% 
% for i=1:12
% 
%     dynCt(1,i) = mean(CtArr(1,i,:));
%     dynCt(2,i) = std(CtArr(1,i,:));
% 
% end
% 
% disp(dynCt(1,:));
% disp(dynCt(2,:));


%disp([mean(CtArr(1,:)),std(CtArr(1,:))])

%dynCT = dynamicCt(caseDir, 165, 165, 200, 25);
%disp(dynCT(:,4))


%%Functions
function [CTArr] = dynamicCt(dir, minVOffset, maxVOffset, minHOffset, maxHoffset)
    dataDir = dir;
    
    [COLORMAP_u,COLORMAP_v,COLORMAP_w,COLORMAP_p]=self_colormap();
    
    P_inf = 1.017e5;
    D = 0.2;
    
    Area_disc= pi*(D/2)^2;
    
    dataParent = load(dataDir);
    Rho_inf = 1;
    dispVol1 = dataParent.por * Area_disc * 3e-3;
    %disp(dispVol1)
    dispVol2 = Area_disc * 3e-3;
    dispVol2
    freq = dataParent.freq;
    T = 1/freq;
    dt = T/10;
    DT = 2 * dt;
    deltaX = dataParent.deltaX;
    deltaY = dataParent.deltaY;

    % % % % % matrixIn using a mean filter over a rectangle of size (2*Nr+1)-by-(2*Nc+1) 
    Nr=19;
    Nc=10;

    CTArr = zeros(3,12);
    CTFArr = zeros(3,12);
    CTBArr = zeros(3,12);

    for j=1:12
   
        if j==12
            DATA2 = eval("dataParent.phase"+num2str(3));
            DATA2.surgeVel=dataParent.surgeVel(3);

        else
            DATA2 = eval("dataParent.phase"+num2str(j+1));
            DATA2.surgeVel=dataParent.surgeVel(j+1);

        end

        DATA2.Vx=smoothdata2(DATA2.Vx,"movmean",Nr);
        DATA2.Vr=smoothdata2(DATA2.Vr,"movmean",Nr);
    
        DATA2.U=DATA2.Vx;
        DATA2.V=DATA2.Vr;
        DATA2.X=DATA2.X;
        DATA2.R=DATA2.Y;
        

        if j==1
            DATA1 = eval("dataParent.phase"+num2str(10));
            DATA1.surgeVel=dataParent.surgeVel(10);

        else
            DATA1 = eval("dataParent.phase"+num2str(j-1));
            DATA1.surgeVel=dataParent.surgeVel(j-1);

        end

        DATA1.Vx=smoothdata2(DATA1.Vx,"movmean",Nr);
        DATA1.Vr=smoothdata2(DATA1.Vr,"movmean",Nr);

        DATA1.U=DATA1.Vx;
        DATA1.V=DATA1.Vr;
        DATA1.X=DATA1.X;
        DATA1.R=DATA1.Y;

        DATA = eval("dataParent.phase"+num2str(j));
    
        V_inf = DATA.Vinf;
        Q_inf =(1/2)*Rho_inf*(V_inf.^2);

        DATA.Vx=smoothdata2(DATA.Vx,"movmean",Nr);
        DATA.Vr=smoothdata2(DATA.Vr,"movmean",Nr);   
    
        DATA.U=DATA.Vx;
        DATA.V=DATA.Vr;
        DATA.X=DATA.X;
        DATA.R=DATA.Y;
        DATA.surgeAccl=dataParent.surgeAccl(j);

        % % % % %  for steady case, time derivative is 0. 
        % DUDT=DATA.U*0.0;
        % DVDT=DATA.V*0.0;
        % 
        DUDT=(DATA2.U-DATA1.U)/DT;
        DVDT=(DATA2.V-DATA1.V)/DT;
        DUDTF=(DATA2.U-DATA.U)/dt;
        DVDTF=(DATA2.V-DATA.V)/dt;
        DUDTB=(DATA.U-DATA1.U)/dt;
        DVDTB=(DATA.V-DATA1.V)/dt;
        DSurgeVelDT = DATA.surgeAccl;
        DSurgeVelDT1 = 0.2 * (DATA2.surgeVel - DATA1.surgeVel)/DT;%0.2 axialInd for P70
        %disp(DSurgeVelDT1)
        DSurgeVelDT2 = (DATA2.surgeVel - DATA1.surgeVel)/DT;
        
    
        %Index corresponding to disc centre
        rc = dataParent.rcIndex;
        xc = dataParent.xcIndex(j);
        
        %Control volume bound in indices
        rMax = rc-maxVOffset; %Lowest index corresponds to highest position for radial direction
        rMin = rc+minVOffset;
        xMin = xc-minHOffset;
        xMax = xc+maxHoffset;
        xSize = size(xMin:xMax);
        rSize = size(rMax:rMin);

        %disp([j,xMin, xMax]);
        
        %Indices spanning the boundary of control volume
        rIndices = [rMin:-1:rMax ones(xSize,"int64")*rMax rMax:rMin ones(xSize,"int64")*rMin];
        xIndices = [ones(rSize,"int64")*xMin xMin:xMax ones(rSize,"int64")*xMax xMax:-1:xMin];
        
        %Coordinates corresponding to control volume boundaries
        xline = DATA.X(rc,xIndices);
        yline = DATA.Y(rIndices,xc).';
        zline = zeros(size(xline));
        
        dyx=sqrt((yline(2:end) - yline(1:end-1)).^2+(xline(2:end) - xline(1:end-1)).^2);
        
        %AREA_MIDPOINT = abs(2*pi*( yline(2:end)+yline(1:end-1))/2.*dyx);
        %TOTALAREA=sum(AREA_MIDPOINT);

        x1 = DATA.X(rc,xMin:xMax);
        y1 = DATA.Y(rMin:-1:rMax,xc).';
        AREA_MIDPOINT = abs([pi*(y1(2:end).^2 - y1(1:end-1).^2) 0 2*pi*(max(y1) - DATA.Y(rc,xc))*(x1(2:end) - x1(1:end-1)) 0 pi*(y1(end:-1:2).^2 - y1(end-1:-1:1).^2) 0 2*pi*(DATA.Y(rc,xc) - min(y1))*(x1(end:-1:2) - x1(end-1:-1:1))  ]);

        
        ny=-gradient(xline);
        nx=gradient(yline);
        normvec=sqrt(nx.^2+ny.^2);
        NORVEC=-[nx;ny;nx.*0]./normvec;
    
        DATA.U(1,:)=DATA.U(2,:);
        DATA.V(1,:)=DATA.V(2,:).*0;
        
        DUDT(1,:)=DATA.U(2,:)*0;
        DVDT(1,:)=DATA.V(2,:)*0;
        DUDTF(1,:)=DATA.U(2,:)*0;
        DVDTF(1,:)=DATA.V(2,:)*0;
        DUDTB(1,:)=DATA.U(2,:)*0;
        DVDTB(1,:)=DATA.V(2,:)*0;
    
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
        DUDT_VEC = zeros(3,length(xline));
        DUDTF_VEC = zeros(3,length(xline));
        DUDTB_VEC = zeros(3,length(xline));
        PARTIAL1st = zeros(4,length(xline));
        PARTIAL2nd = zeros(6,length(xline));
    
        for i=1:length(xline)
        
            UC_VEC(1,i) = DATA.U(rIndices(1,i),xIndices(1,i));
            UC_VEC(2,i) = DATA.V(rIndices(1,i),xIndices(1,i));
        
            DUDT_VEC(1,i) = DUDT(rIndices(1,i),xIndices(1,i));
            DUDT_VEC(2,i) = DVDT(rIndices(1,i),xIndices(1,i));

            DUDTF_VEC(1,i) = DUDTF(rIndices(1,i),xIndices(1,i));
            DUDTF_VEC(2,i) = DVDTF(rIndices(1,i),xIndices(1,i));

            DUDTB_VEC(1,i) = DUDTB(rIndices(1,i),xIndices(1,i));
            DUDTB_VEC(2,i) = DVDTB(rIndices(1,i),xIndices(1,i));
        
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
        [DUDTVEC, PARTIALVEC] = convert_from_cylindrical_to_cartesian(DUDT_VEC,PARTIAL1st,PARTIAL2nd);
        [DUDTFVEC, PARTIALVEC] = convert_from_cylindrical_to_cartesian(DUDTF_VEC,PARTIAL1st,PARTIAL2nd);
        [DUDTBVEC, PARTIALVEC] = convert_from_cylindrical_to_cartesian(DUDTB_VEC,PARTIAL1st,PARTIAL2nd);
        
        [FLUX, F1, F2, F3, F4, F5, F6] = calculate_flux_3D(...
            [xline;yline;zline],UVEC,DUDTVEC,VORTICITY,PARTIALVEC,1.5e-5,NORVEC);

        [FLUXF, F1, F2, F3, F4, F5, F6] = calculate_flux_3D(...
            [xline;yline;zline],UVEC,DUDTFVEC,VORTICITY,PARTIALVEC,1.5e-5,NORVEC);
        
        [FLUXB, F1, F2, F3, F4, F5, F6] = calculate_flux_3D(...
            [xline;yline;zline],UVEC,DUDTBVEC,VORTICITY,PARTIALVEC,1.5e-5,NORVEC);

        FLUXMID= (FLUX(:,1:end-1)+FLUX(:,2:end))/2;
        FLUXFMID= (FLUXF(:,1:end-1)+FLUXF(:,2:end))/2;
        FLUXBMID= (FLUXB(:,1:end-1)+FLUXB(:,2:end))/2;
        FAREA=0.5*FLUXMID.*AREA_MIDPOINT;%Averaging the force contribution from both halves of the disc
        FAREAF=0.5*FLUXFMID.*AREA_MIDPOINT;
        FAREAB=0.5*FLUXBMID.*AREA_MIDPOINT;
        %disp(sum(FAREA,2));
        FVOL = - Rho_inf * dispVol2 * DSurgeVelDT;
        %disp(FVOL);
        %bodyAcclForce = [bodyX; 0; 0];
        bodyAcclForce=0;
        
        LOAD=sum(FAREA,2) + bodyAcclForce;
        LOADF=sum(FAREAF,2);
        LOADB=sum(FAREAB,2);
        %LOAD = sum(FLUXMID.*dyx,2)*pi*D/2 + bodyAcclForce;
        CT = LOAD/(0.5*Rho_inf*V_inf^2*Area_disc);
        CTF=LOADF/(0.5*Rho_inf*V_inf^2*Area_disc);
        CTB=LOADB/(0.5*Rho_inf*V_inf^2*Area_disc);
        %CT(1,:) = 0.5 * CT(1,:); %Averaging the force contribution from both halves of the disc
        CTArr(:,j) = CT;
        CTFArr(:,j) = CTF;
        CTBArr(:,j) = CTB;

        fig = figure;
        ax = gca;
        hold on
        pcolor(DATA.X,DATA.R,OMEGAZ)
        shading flat
        colorbar
        ax.CLim = [-125,125];
        %quiver(DATA.X,DATA.R,DATA.U,DATA.V)
        plot(xline,yline,'g--')
        daspect([1 1 1])

    end
    
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

    %CTArr
    %CTFArr
    %CTBArr

end





