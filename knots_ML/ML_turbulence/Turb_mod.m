function [turb] = Turb_mod(H,V,SR,Pixel)
%%H:  Number of pixels in the horizontal direction
%%V:  Number of pixels in the vertical direction
%%SR: Strehl ratio with values in the inverval [0,1]
%%w0: Beam size at the SLM in  mm 
%%Pixel size in um

Size = min(H,V);

% Number of points for square area
Delta = 1/Pixel/Size; % increment size for x and y
r0 = Delta*(SR^(5/6)/(1-SR^(5/6)))^(3/5); % Fried's Parameter

% put zero (origin) between samples to avoid singularity
[nx,ny] = meshgrid((1:Size)-Size/2-1/2);
Modgrid = real(exp(-1i*pi*(nx+ny)));
rr = (nx.^2 + ny.^2)*Delta^2;

% Square root of the Kolmogorov spectrum:
qKol = 0.1516*Delta/r0^(5/6)*rr.^(-11/12);

f0 = (randn(Size)+1i*randn(Size)).*qKol/sqrt(2);
f1 = fft2(f0).*Modgrid;

ary = [-0.25,-0.25,-0.25,-0.125,-0.125,-0.125,0,0,0,0,0.125,0.125,0.125,0.25,0.25,0.25];
bry = [-0.25,0,0.25,-0.125,0,0.125,-0.25,-0.125,0.125,0.25,-0.125,0,0.125,-0.25,0,0.25];
dary = [0.25,0.25,0.25,0.125,0.125,0.125,0.25,0.125,0.125,0.25,0.125,0.125,0.125,0.25,0.25,0.25];
dbry = [0.25,0.25,0.25,0.125,0.125,0.125,0.25,0.125,0.125,0.25,0.125,0.125,0.125,0.25,0.25,0.25];

ss = (ary.*ary+bry.*bry)*Delta^2;
qsKol = 0.1516*Delta/r0^(5/6)*ss.^(-11/12);
f0 = (randn(1,16)+1i*randn(1,16)).*qsKol/sqrt(2);
fn = f1; % zeros(Size);
for pp = 1:16
  eks = exp(1i*2*pi*(nx*ary(pp)+ny*bry(pp))/Size);
  fn = fn + f0(pp)*eks*dary(pp)*dbry(pp);
end

ff = zeros(Size,Size);
ff((Size/2-Size/2+1):(Size/2+Size/2),(Size/2-Size/2+1):(Size/2+Size/2)) = real(fn);
yo=(H-Size)/2+1;
yf=yo+Size-1;
xo=(V-Size)/2+1;
xf=xo+Size-1;
T=zeros(V,H);
turb0 = ff(:,1:Size);
T(xo:xf,yo:yf)=turb0;
turb=T;
end

