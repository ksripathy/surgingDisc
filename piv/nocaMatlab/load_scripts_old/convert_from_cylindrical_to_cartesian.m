%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function converts the units and partial derivatives of the flow in an
%axisymmetric flow from cylindrical coordinates to cartesian
%
%The method is apllied to the 3D flow case
%
% Developed by: Carlos Simao Ferreira, DUWIND TUDelft
% version: 2018-11-12
% 
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






function [UVEC, PARTIALVEC] = convert_from_cylindrical_to_cartesian(UC_VEC,PARTIAL1st,PARTIAL2nd)

%%% zeros
ZEROVEC = zeros(size(UC_VEC(2,:)));

%% from a 2D space to a 3D space, convert velocity field
UVEC = [UC_VEC(1,:); UC_VEC(2,:); ZEROVEC];

%% first order partial derivatives

% allocate stress tensor terms
ddx_u = PARTIAL1st(1,:);
ddx_v = PARTIAL1st(2,:);
ddx_w = ZEROVEC;
ddy_u = PARTIAL1st(3,:);
ddy_v = PARTIAL1st(4,:);
ddy_w = ZEROVEC;
ddz_u = ZEROVEC;
ddz_v = ZEROVEC;
ddz_w = ZEROVEC;

ddxddx_u = PARTIAL2nd(1,:);
ddyddx_v = PARTIAL2nd(2,:);
ddyddy_u = PARTIAL2nd(3,:);
ddzddx_w = ZEROVEC;
ddzddz_u = ZEROVEC;
ddxddx_v = PARTIAL2nd(4,:);
ddxddy_u = PARTIAL2nd(5,:);
ddyddy_v = PARTIAL2nd(6,:);
ddzddy_w = ZEROVEC;
ddzddz_v = ZEROVEC;
ddxddx_w = ZEROVEC;
ddxddz_u = ZEROVEC;
ddyddy_w = ZEROVEC;
ddyddz_v = ZEROVEC;
ddzddz_w = ZEROVEC;

PARTIALVEC=[...
ddx_u;...
ddx_v ;...
ddx_w ;...
ddy_u ;...
ddy_v;...
ddy_w ;...
ddz_u ;...
ddz_v ;...
ddz_w ;...
ddxddx_u ;...
ddyddx_v ;...
ddyddy_u ;...
ddzddx_w ;...
ddzddz_u ;...
ddxddx_v ;...
ddxddy_u ;...
ddyddy_v ;...
ddzddy_w ;...
ddzddz_v ;...
ddxddx_w ;...
ddxddz_u ;...
ddyddy_w ;...
ddyddz_v;...
ddzddz_w];






