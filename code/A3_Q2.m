clear
m0 = 9.1093837015E-31;
mn = 0.26*m0;
T = 300;
kB = 1.38064852E-23;
tmn = 0.2E-12;
q = 1.60217662E-19;

vth = ((2*kB*T)/mn)^0.5; %thermal velocity

mfp = vth*tmn; %mean free path

v = vth;

Vx = 0.5;

xmax = 200E-9; %max positions
L =xmax;
ymax = 100E-9;
W = ymax;

Np = 10000; % # particles, want 1000-10000
Nplot = 10;

Px = xmax*rand(Np,1);
vx = v*(randn(Np,1)-0.5);
% vacc = sqrt((2*q*V)/mn);

Py = ymax*rand(Np,1);
vy = v*(randn(Np,1)-0.5); %initial velocities

n = 1E19; % electron concentration in m^−2
A = xmax*ymax; % cross-sectional area in m^2
% I = n*A*vx*q; % drift current calculation

dt = 0.01*(ymax/v); %time step

tstop = 1000; %simulation time

I = zeros(tstop,1);

Ppx = Px; %previous postions
Ppy = Py;

box1 = 1; %turning on box 1
box2 = 1; %turning on box 2

figure(1)
xlabel('X (m)')
ylabel('Y (m)')
hold on
axis([0 xmax 0 ymax]);

if box1 == 1
b1 = rectangle('Position',[0.8E-7,0,0.4E-7,0.4E-7]); %box1 position (bottom)
end

if box2 == 1
b2 = rectangle('Position',[0.8E-7,0.6E-7,0.4E-7,0.4E-7]); %box2 position (top)
end

inbox = Px > 0.8E-7 & Px < 1.2E-7 & (Py < 0.4E-7 | Py > 0.6E-7);
while sum(inbox) > 0
    Px(inbox) = rand(sum(inbox),1)*xmax;
    Py(inbox) = rand(sum(inbox),1)*ymax;
    inbox = Px > 0.8E-7 & Px < 1.2E-7 & (Py < 0.4E-7 | Py > 0.6E-7);
end


nx = 20;
ny = 20;

V = zeros(nx,ny);
G = sparse(ny*nx,ny*nx);
F = zeros(ny*nx,1);
Cond = 1; %conductivity outside the boxes


CondBox = 10^-2; %conductivity inside the boxes



for i = 1:nx
    for j = 1:ny
        
        l = (L/nx)*i;
        w = (W/ny)*j;
        
        if l > 0.8E-7 && l < 1.2E-7 && (w > 0.6E-7 || w < 0.4E-7)
            C(i,j) = CondBox;
        else
            C(i,j) = Cond;
        end
    end
end

case1 = 1; %keep equal to 1 for pt. 2

for i = 1:nx
    for j = 1:ny
        n = j + (i-1)*ny; %mapping equation

        if i == 1
                F(n) = Vx;
                G(n,n) = C(i,j);
        elseif i == nx
            if case1 == 1
                F(n) = 0;
                G(n,n) = C(i,j);
            else
                F(n) = Vx;
                G(n,n) = C(i,j);
            end
        elseif j == 1
            if case1 == 1
                F(n) = 0;
                nxm = j + ((i-1)-1)*ny; %(i-1,j)
                nxp = j + ((i+1)-1)*ny; %(i+1,j)
                nyp = (j+1) + (i-1)*ny; %(i,j+1)
                G(n,n) = -(C(i-1,j) + C(i+1,j) + C(i,j+1));
                G(n,nxm) = C(i-1,j);
                G(n,nxp) = C(i+1,j);
                G(n,nyp) = C(i,j+1);
            else
                F(n) = 0;
                G(n,n) = C(i,j);
            end
        elseif j == ny
            if case1 == 1
                F(n) = 0;
                nxm = j + ((i-1)-1)*ny; %(i-1,j)
                nxp = j + ((i+1)-1)*ny; %(i+1,j)
                nym = (j-1) + (i-1)*ny; %(i,j-1)
                G(n,n) = -(C(i-1,j) + C(i+1,j) + C(i,j-1));
                G(n,nxm) = C(i-1,j);
                G(n,nxp) = C(i+1,j);
                G(n,nym) = C(i,j-1);
            else
                F(n) = 0;
                G(n,n) = C(i,j);
            end
            
        else
            nxm = j + ((i-1)-1)*ny; %(i-1,j)
            nxp = j + ((i+1)-1)*ny; %(i+1,j)
            nym = (j-1) + (i-1)*ny; %(i,j-1)
            nyp = (j+1) + (i-1)*ny; %(i,j+1)
            
            G(n,n) = -(C(i-1,j) + C(i+1,j) + C(i,j-1) + C(i,j+1));
            G(n,nxm) = C(i-1,j);
            G(n,nxp) = C(i+1,j);
            G(n,nym) = C(i,j-1);
            G(n,nyp) = C(i,j+1);
            
        end
    end
end

axis([0 L 0 W]);

M = G\F;
for i = 1:nx
    for j = 1:ny
        n = j + (i-1)*ny;
        V(i,j) = M(n);
    end
end


a = ymax;
b = xmax/2;
x = linspace(-b,b,nx);
y = linspace(0,a,ny);

figure(2)
[Gy,Gx] = meshgrid(x,y);
surf(Gy,Gx,V)
title('Surface Potential Plot')
xlabel('Y')
ylabel('X')
zlabel('V(x,y)')


[Ey, Ex] = gradient(V); %plotting gradient
   
dx = xmax/nx;
dy = ymax/ny;

Ex = -Ex/dx;
Ey = -Ey/dy;

figure(3)
quiver(-Ex',-Ey',10) %vector arrows
title('Electric Field Gradient')
xlabel('X')
ylabel('Y')



lx = linspace(0,xmax,nx);
ly = linspace(0,ymax,ny);

[LX,LY] = meshgrid(lx,ly);


c = hsv(Nplot);



for i = 1:tstop
    Ppx = Px;
    Ppy = Py;
    
    Px = Px + vx*dt;
    Py = Py + vy*dt;
    

    xBC = 1; %x boundary conditions ON
    %xBC = 0; %x boundary conditions OFF
   
    if xBC == 1
        ix1 = Px < 0;
        Px(ix1) = Px(ix1) + xmax;
        Ppx(ix1) = Ppx(ix1) + xmax;
        NumPartLeft = sum(ix1);
       
        ix2 = Px > xmax;
        Px(ix2) = Px(ix2) - xmax;
        Ppx(ix2) = Ppx(ix2) - xmax;
        NumPartRight = sum(ix2);
    elseif xBC == 0
        ix1 = Px < 0;
        Px(ix1) = Px(ix1);
        Ppx(ix1) = Ppx(ix1);
        
        ix2 = Px > xmax;
        Px(ix2) = Px(ix2);
        Ppx(ix2) = Ppx(ix2);
    end
    
    iy1 = Py < 0 | Py > ymax;
    vy(iy1) = -vy(iy1);

    
    scatter = 1; %scattering is ON
    %scatter = 0; %scattering is OFF
    
    if scatter == 1
        Psc = 1 - exp(-(dt/tmn));
    elseif scatter == 0
        Psc = 0;
    end
    
    std = sqrt((kB*T)/mn);
    
    Ex_p= interp2(LX,LY,Ex.',Ppx, Ppy);
    Ey_p= interp2(LX,LY,Ey.',Ppx, Ppy);
    Ex_p(isnan(Ex_p)) = 0;
    Ey_p(isnan(Ey_p)) = 0;
    
    Fx = Ex_p*q;
    Fy = Ey_p*q;
    ax = Fx/mn;
    ay = Fy/mn;
    
    isc = Psc > rand(Np,1);
    vx = vx + ax*dt;
    vy = vy + ay*dt;
    Px = Px + vx*dt + 0.5*ax*(dt)^2;
    Py = Py + vy*dt + 0.5*ay*(dt)^2;
    vx(isc) = randn(sum(isc),1)*std;
    vy(isc) = randn(sum(isc),1)*std;
    
    
    if box1 == 1
        inbox1 = Px > 0.8E-7 & Px < 1.2E-7 & Py < 0.4E-7; %in box1
        LorR = Ppx < 0.8E-7 | Ppx > 1.2E-7; %coming from L or R
        Px(inbox1 & LorR) = Ppx(inbox1 & LorR);
        vx(inbox1 & LorR) = -vx(inbox1 & LorR);
        Py(inbox1 & ~LorR) = Ppy(inbox1 & ~LorR);
        vy(inbox1 & ~LorR) = -vy(inbox1 & ~LorR);
    end
    
    if box2 == 1
        inbox2 = Px > 0.8E-7 & Px < 1.2E-7 & Py > 0.6E-7; %in box2
        LorR = Ppx < 0.8E-7 | Ppx > 1.2E-7; %coming from L or R
        Px(inbox2 & LorR) = Ppx(inbox2 & LorR);
        vx(inbox2 & LorR) = -vx(inbox2 & LorR);
        Py(inbox2 & ~LorR) = Ppy(inbox2 & ~LorR);
        vy(inbox2 & ~LorR) = -vy(inbox2 & ~LorR);
    end
    
    
    figure(1)
    xlabel('X (m)')
    ylabel('Y (m)')
    hold on
    axis([0 xmax 0 ymax]);
    
    for j = 1:Nplot
    plot([Ppx(j),Px(j)]',[Ppy(j),Py(j)]','color',c(j,:));
    end
    
    
    vavg = mean(sqrt(vx.^2 + vy.^2)); %average velocity
    
    TSi = ((vavg.^2)*mn)/(2*kB); %temperature of the Si semiconductor
    
    DirectionOverTime = (NumPartRight - NumPartLeft)/dt;
    
    %I(i) = n*A*vx(i)*q; % drift current calculation
    I(i) = n*q*DirectionOverTime; % updated drift current calculation
    
    pause(0.001)
end

hold off




figure(4)
plot(1:tstop,I)
title('Drift Current in the X Direction Over Time')
xlabel('Time (timesteps)')
ylabel('X Current (A)')



Xbins = discretize(Px, 30); %bins for the temp map
Ybins = discretize(Py, 30);

num_bin = zeros(30,30);
Tmap = zeros(30,30); %for loops checking the temp of the particles in each of the x & y bins
for x = 1:30
    for y= 1:30
        num_bin(x,y) = sum(Xbins==x & Ybins==y);
        Vmean = mean(sqrt(vx(Xbins == x & Ybins == y).^2 + vy(Xbins == x & Ybins == y).^2));
        Tmap(x,y) = ((Vmean).^2*mn)/(2*kB);
    end
end


figure(5)
surf(num_bin/Np)
title('Electron Density Map')
xlabel('Y (m)')
ylabel('X (m)')
zlabel('Number of Electrons')

figure(6)
surf(Tmap)
title('Temperature Map')
xlabel('Y (m)')
ylabel('X (m)')
zlabel('Temperature (K)')
