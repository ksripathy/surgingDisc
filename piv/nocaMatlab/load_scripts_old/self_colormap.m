
function[COLORMAP_u,COLORMAP_v,COLORMAP_w,COLORMAP_p]=self_colormap()


load('COLORMAP_HOT_COLD.mat')
RED_MAX=1.2;
WHITE=1;
BLUE_MIN=0;
NUMBER_LEVELS=50;
COLORMAP_u =definenewcolormaprange(COLORMAP_HOT_COLD,NUMBER_LEVELS,RED_MAX,WHITE,BLUE_MIN);


RED_MAX=0.2;
WHITE=0;
BLUE_MIN=-0.2;
NUMBER_LEVELS=50;
COLORMAP_v =definenewcolormaprange(COLORMAP_HOT_COLD,NUMBER_LEVELS,RED_MAX,WHITE,BLUE_MIN);


RED_MAX=20;
WHITE=0;
BLUE_MIN=-20;
NUMBER_LEVELS=50000;
COLORMAP_w =definenewcolormaprange(COLORMAP_HOT_COLD,NUMBER_LEVELS,RED_MAX,WHITE,BLUE_MIN);

RED_MAX=1.0;
WHITE=0;
BLUE_MIN=-1.0;
NUMBER_LEVELS=50;
COLORMAP_p =definenewcolormaprange(COLORMAP_HOT_COLD,NUMBER_LEVELS,RED_MAX,WHITE,BLUE_MIN);