close all
clear all

% datadir.main  = 'C:\Users\kiransripathy\AppData\Local\Packages\CanonicalGroupLimited.Ubuntu_79rhkp1fndgsc\LocalState\ext4.vhdx\home\ksripathy\phd\surgingDisc\piv\data';
% datadir.sub1  = 'matlab/';
% datadir.sub2 = './';
% 
% datadir.file = 'P45Case0.mat';
% datafiles = dir(fullfile(datadir.main,datadir.sub1,datadir.sub2,datadir.file));

datadir = "\\wsl.localhost\Ubuntu\home\ksripathy\phd\surgingDisc\piv\data\matlab\p45Case6.mat";
data = load(datadir);

var = 1; 
temp = eval("data.phase"+ num2str(var));



%DATA=load(fullfile(datadir.main,datadir.sub1,datadir.sub2,datadir.file));