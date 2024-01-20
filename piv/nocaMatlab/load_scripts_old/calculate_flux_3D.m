%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function calculates the flux term in the FLUX EQUATION method of
%reference : A comparison of methods for evaluating time-dependent fluid
%dynamic forces on bodies, using only velocity fields and their derivatives
%F.Noca D. Shiels and D. Jeon
%Journal of Fluid and Structures 1999, 13, 551-578
%
%The method is apllied to the 3D flow case
%
% Developed by: Carlos Simao Ferreira, DUWIND TUDelft
% version: 2018-10-18
% 
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [FLUX, F1, F2, F3, F4, F5, F6] = calculate_flux_3D(coordinates,velocity,time_derivatives,vorticity,viscous_stress_tensor_terms,mu,normaldirection)
%   inputs:
%           U,V,W are the velocities in X,Y,Z direction at the point on the
%           surface.
%           N1,N2,N3 are components of the vector normal to the surfacein X,Y,Z 
%           direction, at the point on the surface.
%           

N=3; % dimensions of space

% normal direction
n1 = normaldirection(1,:);
n2 = normaldirection(2,:);
n3 = normaldirection(3,:);

% coordinates
x = coordinates(1,:);
y = coordinates(2,:);
z = coordinates(3,:);

% velocity terms
u = velocity(1,:);
v = velocity(2,:);
w = velocity(3,:);


% time dependent terms
dudt = time_derivatives(1,:);
dvdt = time_derivatives(2,:);
dwdt = time_derivatives(3,:);

% vorticity
omega_x = vorticity(1,:);
omega_y = vorticity(2,:);
omega_z = vorticity(3,:);

% allocate stress tensor terms
ddx_u = viscous_stress_tensor_terms(1,:);
ddx_v = viscous_stress_tensor_terms(2,:);
ddx_w = viscous_stress_tensor_terms(3,:);
ddy_u = viscous_stress_tensor_terms(4,:);
ddy_v = viscous_stress_tensor_terms(5,:);
ddy_w = viscous_stress_tensor_terms(6,:);
ddz_u = viscous_stress_tensor_terms(7,:);
ddz_v = viscous_stress_tensor_terms(8,:);
ddz_w = viscous_stress_tensor_terms(9,:);

ddxddx_u = viscous_stress_tensor_terms(10,:);
ddyddx_v = viscous_stress_tensor_terms(11,:);
ddyddy_u = viscous_stress_tensor_terms(12,:);
ddzddx_w = viscous_stress_tensor_terms(13,:);
ddzddz_u = viscous_stress_tensor_terms(14,:);
ddxddx_v = viscous_stress_tensor_terms(15,:);
ddxddy_u = viscous_stress_tensor_terms(16,:);
ddyddy_v = viscous_stress_tensor_terms(17,:);
ddzddy_w = viscous_stress_tensor_terms(18,:);
ddzddz_v = viscous_stress_tensor_terms(19,:);
ddxddx_w = viscous_stress_tensor_terms(20,:);
ddxddz_u = viscous_stress_tensor_terms(21,:);
ddyddy_w = viscous_stress_tensor_terms(22,:);
ddyddz_v = viscous_stress_tensor_terms(23,:);
ddzddz_w = viscous_stress_tensor_terms(24,:);






% we calculate the first term, F1

F1 = [ (u.^2+v.^2+w.^2)/2 .* n1  ;...
       (u.^2+v.^2+w.^2)/2 .* n2   ;...
       (u.^2+v.^2+w.^2)/2 .* n3  ]; 

% we calculate the second term, F2

F2 = [ -n1.*(u.^2)-n2.*u.*v-n3.*u.*w  ;...
       -n1.*u.*v-n2.*(v.^2)-n3.*v.*w ;...
         -n1.*u.*w-n2.*v.*w-n3.*(w.^2)];  

% vorticity dependent terms

F3 = [ -(y.*omega_z-z.*omega_y).*(u.*n1+v.*n2+n3.*w).*(1/(N-1))  ;...
        (x.*omega_z-z.*omega_x).*(u.*n1+v.*n2+n3.*w).*(1/(N-1))  ;...
       -(x.*omega_y-y.*omega_x).*(u.*n1+v.*n2+n3.*w).*(1/(N-1))];  



F4 = [ -(v.*z-w.*y).*(n1.*omega_x+omega_y.*n2+n3.*omega_z).*(1/(N-1)) ;...
        (u.*z-w.*x).*(n1.*omega_x+omega_y.*n2+n3.*omega_z).*(1/(N-1));...
       -(u.*y-v.*x).*(n1.*omega_x+omega_y.*n2+n3.*omega_z).*(1/(N-1))];  


% time dependent terms   
   
F5 =  [  (-(n1.*dudt+n2.*dvdt+n3.*dwdt).*(N-1).*x+(y.*n2+n3.*z).*dudt-n1.*(y.*dvdt+z.*dwdt)).*(1/(N-1)) ;...
        (-(n1.*dudt+n2.*dvdt+n3.*dwdt).*(N-1).*y+(x.*n1+n3.*z).*dvdt-n2.*(x.*dudt+z.*dwdt)).*(1/(N-1));...
       (-(n1.*dudt+n2.*dvdt+n3.*dwdt).*(N-1).*z+(x.*n1+y.*n2).*dwdt-n3.*(x.*dudt+y.*dvdt)).*(1/(N-1))];     



% viscous terms

F6a = [ n1.*(x.*(2.*ddxddx_u+ddyddx_v+ddyddy_u+ddzddx_w+ddzddz_u)+y.*(ddxddx_v+ddxddy_u+2.*ddyddy_v+ddzddy_w+ddzddz_v)+z.*(ddxddx_w+ddxddz_u+ddyddy_w+ddyddz_v+2.*ddzddz_w)).*(1/(N-1))  ;...
        n2.*(x.*(2.*ddxddx_u+ddyddx_v+ddyddy_u+ddzddx_w+ddzddz_u)+y.*(ddxddx_v+ddxddy_u+2.*ddyddy_v+ddzddy_w+ddzddz_v)+z.*(ddxddx_w+ddxddz_u+ddyddy_w+ddyddz_v+2.*ddzddz_w)).*(1/(N-1));...
        n3.*(x.*(2.*ddxddx_u+ddyddx_v+ddyddy_u+ddzddx_w+ddzddz_u)+y.*(ddxddx_v+ddxddy_u+2.*ddyddy_v+ddzddy_w+ddzddz_v)+z.*(ddxddx_w+ddxddz_u+ddyddy_w+ddyddz_v+2.*ddzddz_w)).*(1/(N-1))];     


F6b = [ -(2.*ddxddx_u+ddyddx_v+ddyddy_u+ddzddx_w+ddzddz_u).*(x.*n1+y.*n2+n3.*z).*(1/(N-1))  ;...
        -(ddxddx_v+ddxddy_u+2.*ddyddy_v+ddzddy_w+ddzddz_v).*(x.*n1+y.*n2+n3.*z).*(1/(N-1));...
        -(ddxddx_w+ddxddz_u+ddyddy_w+ddyddz_v+2.*ddzddz_w).*(x.*n1+y.*n2+n3.*z).*(1/(N-1))];     


F6c = [2.*n1.*ddx_u+(ddx_v+ddy_u).*n2+n3.*(ddx_w+ddz_u);...
       (ddx_v+ddy_u).*n1+2.*n2.*ddy_v+n3.*(ddy_w+ddz_v);...
       (ddx_w+ddz_u).*n1+(ddy_w+ddz_v).*n2+2.*n3.*ddz_w];


F6= mu*(F6a + F6b + F6c);



% total flux element

FLUX=F1+F2+F3+F4+F5+F6;












