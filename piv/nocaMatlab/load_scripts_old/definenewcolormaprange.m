function [COLOR2] = definenewcolormaprange(COLOR1, LEVELS,TOP,MID,BOT)


%%% COLOR1 - original colormap
%%% LEVELS - number of levls of color
%%% TOP - top bound
%%% MID - mid value
%%% BOT - bottom value


%%% size of the original colormap
[M,N]=size(COLOR1);

%%%
COLOR2=zeros(LEVELS,3);


N1=round(((MID-BOT)/(TOP-BOT))*LEVELS);
N2=LEVELS-N1;

IND1=(0:1:round(M/2)-1)';
IND2=(0:1:(N1-1))*max(IND1)/(N1-1);

for L=1:1:3
COLOR2(1:N1,L)=interp1(IND1,COLOR1(1:length(IND1),L),IND2);
end

IND3=(0:1:(M-round(M/2))-1)';
IND2=(0:1:(N2-1))*max(IND3)/(N2-1);

for L=1:1:3
COLOR2(N1+1:LEVELS,L)=interp1(IND1,COLOR1(end-length(IND3)+1:end,L),IND2);
end

























