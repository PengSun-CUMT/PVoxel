function [P,T, ZWD]=PVoxel2(B,L,H,HType,doy, Pmodel, Tmodel,ZWDmodel, ud2008)
%input:
%B: latitude in degree [-90 90]
%L: longitude in degree [-180 180] or[0 360]
%H：height of station in m
%Htype: type of height 0 or 1
      %0:orthometric height
      %1:ellipsolidal height
%doy:day of year
%output:
%P:atmospheric pressure (hPa)
%T:atmospheric temperature(K)
%ZWD:zenith wet delay(m) 

%%
%Contact Peng Sun (sunpcxs@gmail.com) if you find any bugs or if you need any helps
%reference 

P00model=Pmodel(:,1:5);
P04model=Pmodel(:,6:10);
P10model=Pmodel(:,11:15);
P15model=Pmodel(:,16:20);
T00model=Tmodel(:,1:5);
T04model=Tmodel(:,6:10);
T10model=Tmodel(:,11:15);
T15model=Tmodel(:,16:20);
ZWD00model=ZWDmodel(:,1:5);
ZWD04model=ZWDmodel(:,11:15);
gama00model=ZWDmodel(:,6:10);%gama=γ+1
gama04model=ZWDmodel(:,16:20);
% mean gravity in m/s^2
gm = 9.80665;
% molar mass of dry air in kg/mol
dMtr = 28.965*10^-3;
% universal gas constant in J/K/mol
Rg = 8.31432;
ct=gm*dMtr/Rg;
ppod=90-B;
plon=L;
if(L<0)
    plon=L+360;
end
% find the indx (line in the grid file) of the nearest point
% changed for the 1 degree grid
ipod = floor(ppod+1);
ilon = floor(plon+1);

% normalized (to one) differences, can be positive or negative
% changed for the 1 degree grid
diffpod = (ppod - (ipod - 0.5));
difflon = (plon - (ilon - 0.5));
% changed for the 1 degree grid
if ipod == 181
    ipod = 180;
end
if ilon == 361
    ilon = 1;
end
if ilon == 0
    ilon = 360;
end

% get the number of the corresponding line
% changed for the 1 degree grid
indx(1) = (ipod - 1)*360 + ilon;

% near the poles: nearest neighbour interpolation, otherwise: bilinear
% with the 1 degree grid the limits are lower and upper
bilinear = 0;
if ppod > 0.5 && ppod < 179.5
    bilinear = 1;
    ipod1 = ipod + sign(diffpod);
    ilon1 = ilon + sign(difflon);
    % changed for the 1 degree grid
    if ilon1 == 361
        ilon1 = 1;
    end
    if ilon1 == 0
        ilon1 = 360;
    end
    % changed for the 1 degree grid
    indx(2) = (ipod1 - 1)*360 + ilon;  % along same longitude
    indx(3) = (ipod  - 1)*360 + ilon1; % along same polar distance
    indx(4) = (ipod1 - 1)*360 + ilon1; % diagonal
end

dnpod1 = abs(diffpod); % distance nearer point
dnpod2 = 1 - dnpod1;   % distance to distant point
dnlon1 = abs(difflon);
dnlon2 = 1 - dnlon1;
if(HType==1)
    unduIdx=ud2008(indx);
    R1 = dnpod2*unduIdx(1)+dnpod1*unduIdx(2);
    R2 = dnpod2*unduIdx(3)+dnpod1*unduIdx(4);
    undu = dnlon2*R1+dnlon1*R2;
    H=H-undu;
end
cosyP=cos((doy-P00model(indx,3))*2*pi/365.25);
coshyP=cos((doy-P00model(indx,5))*4*pi/365.25);
P00=P00model(indx,1)+P00model(indx,2).*cosyP+P00model(indx,4).*coshyP;
cosyT=cos((doy-T00model(indx,3))*2*pi/365.25);
coshyT=cos((doy-T00model(indx,5))*4*pi/365.25);
T00=T00model(indx,1)+T00model(indx,2).*cosyT+T00model(indx,4).*coshyT;
cosyP=cos((doy-P04model(indx,3))*2*pi/365.25);
coshyP=cos((doy-P04model(indx,5))*4*pi/365.25);
P04=P04model(indx,1)+P04model(indx,2).*cosyP+P04model(indx,4).*coshyP;
cosyT=cos((doy-T04model(indx,3))*2*pi/365.25);
coshyT=cos((doy-T04model(indx,5))*4*pi/365.25);
T04=T04model(indx,1)+T04model(indx,2).*cosyT+T04model(indx,4).*coshyT;
cosyP=cos((doy-P10model(indx,3))*2*pi/365.25);
coshyP=cos((doy-P10model(indx,5))*4*pi/365.25);
P10=P10model(indx,1)+P10model(indx,2).*cosyP+P10model(indx,4).*coshyP;
cosyT=cos((doy-T10model(indx,3))*2*pi/365.25);
coshyT=cos((doy-T10model(indx,5))*4*pi/365.25);
T10=T10model(indx,1)+T10model(indx,2).*cosyT+T10model(indx,4).*coshyT;
cosyP=cos((doy-P15model(indx,3))*2*pi/365.25);
coshyP=cos((doy-P15model(indx,5))*4*pi/365.25);
P15=P15model(indx,1)+P15model(indx,2).*cosyP+P15model(indx,4).*coshyP;
cosyT=cos((doy-T15model(indx,3))*2*pi/365.25);
coshyT=cos((doy-T15model(indx,5))*4*pi/365.25);
T15=T15model(indx,1)+T15model(indx,2).*cosyT+T15model(indx,4).*coshyT;
cosy=cos((doy-ZWD00model(indx,3))*2*pi/365.25);
coshy=cos((doy-ZWD00model(indx,5))*4*pi/365.25);
ZWD00=ZWD00model(indx,1)+ZWD00model(indx,2).*cosy+ZWD00model(indx,4).*coshy;
cosy=cos((doy-gama00model(indx,3))*2*pi/365.25);
coshy=cos((doy-gama00model(indx,5))*4*pi/365.25);
gama00=gama00model(indx,1)+gama00model(indx,2).*cosy+gama00model(indx,4).*coshy;
cosy=cos((doy-ZWD04model(indx,3))*2*pi/365.25);
coshy=cos((doy-ZWD04model(indx,5))*4*pi/365.25);
ZWD04=ZWD04model(indx,1)+ZWD04model(indx,2).*cosy+ZWD04model(indx,4).*coshy;
cosy=cos((doy-gama04model(indx,3))*2*pi/365.25);
coshy=cos((doy-gama04model(indx,5))*4*pi/365.25);
gama04=gama04model(indx,1)+gama04model(indx,2).*cosy+gama04model(indx,4).*coshy;

btvoxel(1,:)=(T10-T15)/5000;
btvoxel(2,:)=(T04-T10)/6000;
btvoxel(3,:)=(T00-T04)/4000;
% case of nearest neighbourhood
if bilinear == 0
    if(H==4000)
        P=P04;
        T=T04;
    elseif(H==10000)
        P=P10;
        T=T10;
    elseif(H<=0)
        bt=btvoxel(3,1);
        dh=H-0;
        P=P00*(1-bt*dh/T00)^(ct/bt);
        T=T00-bt*dh;
    elseif(H>=15000)
        bt=btvoxel(1,1);
        dh=H-15000;
        P=P15*exp(-ct*dh/T15);
        T=T15-bt*dh;
    else
        if(H>0 && H<4000)
            dh1=H;
            dh2=H-4000;
            bt=btvoxel(3,1);
            P1=P00*(1-bt*dh1/T00)^(ct/bt);
            P2=P04*(1-bt*dh2/T04)^(ct/bt);
            T=T00-bt*dh1;
        elseif(H>4000 && H<10000)
            dh1=H-4000;
            dh2=H-10000;
            bt=btvoxel(2,1);
            P1=P04*(1-bt*dh1/T04)^(ct/bt);
            P2=P10*(1-bt*dh2/T10)^(ct/bt);
            T=T04-bt*dh1;
        elseif(H>10000 && H<15000)
            dh1=H-10000;
            dh2=H-15000;
            P1=P10*exp(-ct*dh1/T10);%Tv~=T at 10 km;
            P2=P15*exp(-ct*dh2/T15);%Tv~=T at 15 km;
            bt=btvoxel(1,1);
            T=T10-bt*dh1;
        end
        Q1=1/dh1/(1/dh1+1/abs(dh2));
        Q2=1/abs(dh2)/(1/dh1+1/abs(dh2));
        P=P1*Q1+P2*Q2;
    end
    if(H<=0)
        ZWD=ZWD00*(P/P00)^gama00;
    elseif(H>=4000)
        ZWD=ZWD04*(P/P04)^gama04;
    else
        ZWD1=ZWD00*(P/P00)^gama00;
        ZWD2=ZWD04*(P/P04)^gama04;
        dh1=H;
        dh2=abs(H-4000);
        Q1=1/dh1/(1/dh1+1/abs(dh2));
        Q2=1/abs(dh2)/(1/dh1+1/abs(dh2));
        ZWD=ZWD1*Q1+ZWD2*Q2;
    end
else
    Ptemp=[0 0 0 0];
    Ttemp=[0 0 0 0];
    ZWDtemp=[0 0 0 0];
    for i=1:4
        if H==4000
            Ptemp(i)=P04(i);
            Ttemp(i)=T04(i);
        elseif (H==10000)
            Ptemp(i)=P10(i);
            Ttemp(i)=T10(i);
        elseif(H<=0)
            dh=H-0;
            bt=btvoxel(3,i);
            Ptemp(i)=P00(i)*(1-bt*dh/T00(i))^(ct/bt);
            Ttemp(i)=T00(i)-bt*dh;
        elseif(H>=15000)
            dh=H-15000;
            bt=btvoxel(1,i);
            Ptemp(i)=P15(i)*exp(-ct*dh/T15(i));
            Ttemp(i)=T15(i)-bt*dh;
        else
            if(H>0 && H<4000)
                dh1=H;
                dh2=H-4000;
                bt=btvoxel(3,i);
                p1=P00(i)*(1-bt*dh1/T00(i))^(ct/bt);
                p2=P04(i)*(1-bt*dh2/T04(i))^(ct/bt);
                Ttemp(i)=T00(i)-bt*dh1;
            elseif(H>4000 && H<10000)
                dh1=H-4000;
                dh2=H-10000;
                bt=btvoxel(2,i);
                p1=P04(i)*(1-bt*dh1/T04(i))^(ct/bt);
                p2=P10(i)*(1-bt*dh2/T10(i))^(ct/bt);
                Ttemp(i)=T04(i)-bt*dh1;
            elseif(H>10000 && H<15000)
                bt=btvoxel(1,i);
                dh1=H-10000;
                dh2=H-15000;
                p1=P10(i)*exp(-ct*dh1/T10(i));%Tv~=T at 10 km;
                p2=P15(i)*exp(-ct*dh2/T15(i));%Tv~=T at 15 km;
                Ttemp(i)=T10(i)-bt*dh1;
            end
            Q1=1/dh1/(1/dh1+1/abs(dh2));
            Q2=1/abs(dh2)/(1/dh1+1/abs(dh2));
            Ptemp(i)=p1*Q1+p2*Q2;

        end
        if(H<=0)
            ZWDtemp(i)=ZWD00(i)*(Ptemp(i)/P00(i))^gama00(i);
        elseif(H>=4000)
            ZWDtemp(i)=ZWD04(i)*(Ptemp(i)/P04(i))^gama04(i);
        else
            ZWD1=ZWD00(i)*(Ptemp(i)/P00(i))^gama00(i);
            ZWD2=ZWD04(i)*(Ptemp(i)/P04(i))^gama04(i);
            dh1=H;
            dh2=abs(H-4000);
            Q1=1/dh1/(1/dh1+1/abs(dh2));
            Q2=1/abs(dh2)/(1/dh1+1/abs(dh2));
            ZWDtemp(i)=ZWD1*Q1+ZWD2*Q2;
        end
    end

    % pressure
    R1 = dnpod2*Ptemp(1)+dnpod1*Ptemp(2);
    R2 = dnpod2*Ptemp(3)+dnpod1*Ptemp(4);
    P = dnlon2*R1+dnlon1*R2;
    % temperature
    R1 = dnpod2*Ttemp(1)+dnpod1*Ttemp(2);
    R2 = dnpod2*Ttemp(3)+dnpod1*Ttemp(4);
    T = dnlon2*R1+dnlon1*R2;
    % temperature
    R1 = dnpod2*ZWDtemp(1)+dnpod1*ZWDtemp(2);
    R2 = dnpod2*ZWDtemp(3)+dnpod1*ZWDtemp(4);
    ZWD = dnlon2*R1+dnlon1*R2;
end