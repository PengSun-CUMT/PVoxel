load("Pmodel.mat");
load("Tmodel.mat");
load("ZWDmodel.mat");
load("ud2008.mat");
B=30;%latitude in degree;
L=117;%longitude in degree;
H=100;%height in m;
doy=123;%day of year;
HType=0;%type of height; 0:orthometric height ; 1:ellipsolidal height
%P:atmospheric pressure (hPa)
%T:atmospheric temperature(K)
%ZWD:zenith wet delay(m) 
[P,T, ZWD]=PVoxel(B,L,H,HType,doy,Pmodel, Tmodel,ZWDmodel,ud2008);