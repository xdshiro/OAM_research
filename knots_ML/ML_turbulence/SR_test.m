clearvars;

ft2 = @(x) fftshift(fft2(fftshift(x)));
ift2 = @(x) fftshift(ifft2(fftshift(x)));

%% Parameters [um]
% N = 512;
N = 128;
L = 20; %Length of background in 
Ld = 532e-3; %Incident wavelength of light
wo = 1.3; %Use this
zR = pi*wo^2/Ld; %Rayleigh length
n0 = 1.0;
params = [wo,Ld,zR,n0];

%% Axis: 
dx = L/N;
x = -L/2:dx:L/2-dx; y = x;
[X,Y] = meshgrid(x,y); 
[theta,r] = cart2pol(X,Y);

%% Optical beam

LG =@(p,m,z) sqrt(factorial(p)/(pi*factorial(abs(m)+p))).*r.^abs(m)./(wo^(abs(m)+1)).*exp(1i*m*theta) ...
          .*(1-1i*z/zR)^p/(1+1i*z/zR)^(p+abs(m)+1).*exp(-r.^2./(2*wo^2*(1+1i*z/zR))).*LaguerrePoly([p,abs(m)],r.^2./(wo^2*(1+z^2/zR^2)));

Uo = LG(0,0,0);

%% Propagation: 

pad = 2*N; crop = 10;
U = ift2(padarray(Uo,[pad pad],0,'both'));
U = U(size(U,1)/2-N/2:size(U,1)/2+N/2-1 , size(U,1)/2-N/2:size(U,1)/2+N/2-1);

P = dx*1.3e-1; SR = .95;

figure(1); clf(1);
img1 = imagesc(zeros(N));
axis image; axis xy; colormap parula(1024); colorbar

figure(2); clf(2);
img2 = imagesc(zeros(N));
axis image; axis xy; colormap jet(1024); colorbar

figure(3); clf(3);
img3 = imagesc(zeros(N));
% imagesc(abs(U).^2)
axis image; axis xy; colormap jet(1024); colorbar

figure(4); clf(4);
% imagesc(abs(Uo).^2);
img4 = imagesc(zeros(N));
axis image; axis xy; colormap jet(1024); colorbar
ii = 0;

n_screens = 100;
for iz = 1:1e5
    [turb] = Turb_mod(N,N,SR,P);
    turb = mod(turb,2*pi);

    U0 = ft2(U);
    Uz = ft2(U.*exp(1i*turb));

    U0 = U0(size(U0,1)/2-crop+1:size(U0,1)/2+crop , size(U0,1)/2-crop+1:size(U0,1)/2+crop);
    Uz = Uz(size(Uz,1)/2-crop+1:size(Uz,1)/2+crop , size(Uz,1)/2-crop+1:size(Uz,1)/2+crop);

    xo = linspace(-L/2,L/2,size(Uz,1)); yo = xo;
    [Xo,Yo] = meshgrid(xo,yo);

    U0 = interp2(Xo,Yo,U0,X,Y,'cubic');
    Uz = interp2(Xo,Yo,Uz,X,Y,'cubic');

    Iz = abs(Uz(N/2+4,N/2+4)).^2;
    I0 = abs(U0(N/2+4,N/2+4)).^2;
    ratio = Iz/I0;
    img1.CData = turb; drawnow;
    if (ratio >= SR-SR*0.01) && (ratio <= SR+SR*0.01)
        ii = ii+1; disp(ii);
%         img1.CData = turb;
        img2.CData = angle(Uz);
        img3.CData = abs(Uz).^2;
        img4.CData = abs(U0).^2;
        drawnow
        SRlist(ii) = ratio;
        turb_set(:,:,ii) = turb;
        if ii == n_screens
            break
        end
    else
        continue
    end
    

    pause(1)
end

SRm = mean(SRlist)

% name = sprintf('turb_set_SRm_%d.mat',SR);
% save(name,'turb_set')

D_r0 = ((1-SRm^(5/6))/SRm^(5/6))^(3/5);



%% LG Functions
function out = LaguerreGaussian(params,p,l,r,theta,z)
    w0 = params(1);
    lambda = params(2);
    zR = params(3);
    n = params(4);
    k = 2*pi/lambda;

    Npl = sqrt((2*factorial(p))./(pi*factorial(p + abs(l))));
    Gouy = -(2*p + abs(l) + 1).*atan2(z,zR);
    if z ~= 0
        wz = w0*sqrt(1+(z./zR).^2);
        Rz = z.*(1+(zR./z).^2);
        Phi_C = (k*(r.^2)./(2*Rz));
    else
        wz = w0;
        Phi_C = 0;
    end

    out = (Npl./wz).*((sqrt(2)*r)./wz).^(abs(l))...
        .*LaguerrePoly([p,abs(l)],2*r.^2./wz.^2)...
        .*exp(-r.^2./wz.^2).*exp(1i*(l*theta + Phi_C + Gouy));

end

function y = LaguerrePoly(params,x)

    n=params(1)-1;
    k=params(2);

    L0 = 1;
    L1 = 1 + k - x;
    L2 = zeros(size(x));
    switch (n+1)
        case 0
            y = L0;
        case 1
            y = L1;
        otherwise
            for p = 1:n
                %             disp(p)
                L2 = ((2*p + 1 + k - x).*L1 - (p+k)*L0)./(p+1);
                L0 = L1;
                L1 = L2;
            end
            y = L2;
    end
end