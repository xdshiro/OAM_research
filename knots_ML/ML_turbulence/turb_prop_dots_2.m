clearvars;
% N = 512; %sampling number
N = 128;
tic
num = 100;
num = 1;

% for w = 1.3:0.1:1.6
for ii = 1:num
    
    ft2 = @(x) fftshift(fft2(fftshift(x)));
    ift2 = @(x) fftshift(ifft2(fftshift(x)));
    
    %% Parameters [um]
    
    L = 20; %Length of background in
    Ld = 532e-3; %Incident wavelength of light
    ko = (2*pi)/Ld; %wavenumber
    wo = 1.3;
    zR = pi*wo^2/Ld; %Rayleigh length
    Z = 0; Nz = 1; dz = Z/Nz; %destination of z and step size
    n0 = 1.0;
    params = [wo,Ld,zR,n0];
    
    %% Axis:
    
    dx = L/N;
    x = -L/2:dx:L/2-dx; y = x;
    [X,Y] = meshgrid(x,y);
    [theta,r] = cart2pol(X,Y);
    
    %% Optical Knots
    Gaussian = @(XO,YO) exp(-((X-XO).^2+(Y-YO).^2)./(0.05*wo)^2);
    
    z_min = -Z;
    LG =@(p,m,z) sqrt(factorial(p)/(pi*factorial(abs(m)+p))).*r.^abs(m)./(wo^(abs(m)+1)).*exp(1i*m*theta) ...
        .*(1-1i*z/zR)^p/(1+1i*z/zR)^(p+abs(m)+1).*exp(-r.^2./(2*wo^2*(1+1i*z/zR))).*LaguerrePoly([p,abs(m)],r.^2./(wo^2*(1+z^2/zR^2)));
    
    
    % Uo = @(z) 0.61.*LG(0,0,z) - 2.56*LG(1,0,z) +6.15*LG(2,0,z) -6.35*LG(3,0,z) +2.92*LG(4,0,z) -0.61*LG(5,0,z) -2.45*LG(0,5,z); %Cinquefoil
    % Uo = @(z) 1.71.*LG(0,0,z) - 5.66*LG(1,0,z) +6.38*LG(2,0,z) -2.30*LG(3,0,z) -4.36*LG(0,3,z); % Trefoil
    % Uo = @(z) 2.63.*LG(0,0,z) -6.32.*LG(1,0,z) +4.21.*LG(2,0,z) -5.95.*LG(0,2,z); %Hopf
    
    w = 1.376;
    
    a00 = 1 - 2*w^2 + 2*w^4;
    a01 = 2*w^2 - 4*w^4;
    a02 = 2*w^4;
    a20 = 4*sqrt(2)*w^2;
    
    atot2 = a00^2 + a01^2 + a02^2 + a20^2;
    
    an_00 = a00/sqrt(atot2)*10;
    an_01 = a01/sqrt(atot2)*10;
    an_02 = a02/sqrt(atot2)*10;
    an_20 = a20/sqrt(atot2)*10;
    
    Uo = @(z) an_00*LG(0,0,z) +an_01*LG(1,0,z) +an_02*LG(2,0,z) +an_20*LG(0,2,z); % Hopf
    
    
    
    %% 3foil
    
    %         w = 1.05;
    %
    %         a00 = 1 - w^2 - 2*w^4 + 6*w^6;
    %         a01 = w^2*(1 + 4*w^2 - 18*w^4);
    %         a02 = -2*w^4*(1 - 9*w^2);
    %         a03 = -6*w^6;
    %         a30 = -8*sqrt(6)*w^3;
    %
    %         atot2 = a00^2 + a01^2 + a02^2 + a03^2 + a30^2;
    %
    %         an_00 = 1*a00/sqrt(atot2)*10;
    %         an_01 = 1*a01/sqrt(atot2)*10;
    %         an_02 = 1*a02/sqrt(atot2)*10;
    %         an_03 = 1*a03/sqrt(atot2)*10;
    %         an_30 = a30/sqrt(atot2)*10;
    %
    %         Uo = @(z) an_00*LG(0,0,z) +an_01*LG(1,0,z) +an_02*LG(2,0,z) +an_03*LG(3,0,z) +an_30*LG(0,3,z); % Trefoil (w)
    
    %                     Uo = @(z) 1.29*LG(0,0,z) -3.95*LG(1,0,z) +7.49*LG(2,0,z) -3.28*LG(3,0,z) -3.98*LG(0,3,z); % Trefoil our optm
    
    %             Uo = @(z) 1.51.*LG(0,0,z) - 5.06*LG(1,0,z) +7.23*LG(2,0,z) -2.03*LG(3,0,z) -3.37*LG(0,3,z); % Trefoil
    
    %         Uo = @(z) 0.61*LG(0,0,z) -2.56*LG(1,0,z) +6.15*LG(2,0,z) -6.35*LG(3,0,z) +2.92*LG(4,0,z) -0.61*LG(5,0,z) -2.45*LG(0,3,z); %Trefoil single shot
    
    %     Uo = @(z) 3.59.*LG(0,0,z) -6.31*LG(1,0,z) +5.47*LG(2,0,z) +5.0*LG(0,2,z); % Hopf
    
    %         Uo = @(z) LG(0,0,z);
    
    
    %% Propagation:
    
    pad = 2*N; SR = 0.95; crop = 10;
    
    load(sprintf('turb_set_SRm_%d.mat',SR))
    
%     figure(1); clf(1);
    %         img3 = imagesc(zeros(N));
    imagesc(abs(Uo(0)))
%     axis image; axis xy; colormap parula(1024); colorbar
    
%     figure(2); clf(2);
%     img2 = imagesc(zeros(N));
%     axis image; axis xy; colormap jet(1024);
    
%     figure(3); clf(3);
%     [C_real, h_real]  = contour(x,y,ones(N,N), [0 0], '--', 'LineWidth', 1);
%     hold on
%     [C_imag, h_imag]  = contour(x,y,ones(N,N), [0 0], '-', 'LineWidth', 1);
%     positive = plot([0],[0],'r.','markersize',22);
%     hold off
%     axis image; axis square; colormap gray;
    
    xx = L/2;
    
    turb = turb_set(:,:,ii);
    disp(ii)
    img3.CData = turb;
    drawnow
    z = z_min;
    
    for iz = 1:2*Nz
        
        z = z + dz;
        Uz = ift2(padarray(Uo(z),[pad pad],0,'both'));
        Uz = Uz(size(Uz,1)/2-N/2:size(Uz,1)/2+N/2-1 , size(Uz,1)/2-N/2:size(Uz,1)/2+N/2-1);
        Uz = ft2(Uz.*exp(1i*turb));
        Uz = Uz(size(Uz,1)/2-crop+1:size(Uz,1)/2+crop , size(Uz,1)/2-crop+1:size(Uz,1)/2+crop);
        
        xo = linspace(-L/2,L/2,size(Uz,1)); yo = xo;
        [Xo,Yo] = meshgrid(xo,yo);
        Uz = interp2(Xo,Yo,Uz,X,Y,'cubic');
        
%         phase = (angle(Uz) + pi)./(2*pi);
        
%         xi = real(Uz); eta = imag(Uz);
        
%         C1 = contourcs(x,y,xi, [0 0]);    C2 = contourcs(x,y,eta, [0 0]);
%         
%         x1 = []; y1 = [];    x2 = []; y2 = [];
%         
%         for i = 1:length(C1)
%             x1 = cat(2, x1, NaN, C1(i).X);        y1 = cat(2, y1, NaN, C1(i).Y);
%         end
%         for i = 1:length(C2)
%             x2 = cat(2, x2, NaN, C2(i).X);        y2 = cat(2, y2, NaN, C2(i).Y);
%         end
%         
%         [xout,yout] = intersections(x1, y1, x2, y2 ,1);
%         
%         num_cargas = length(xout);
%         
%         if (num_cargas ~= 0)
%             for kk = 1:num_cargas
%                 if (isnan(yout(kk)))
%                     break;
%                 end
%                 p_charges{iz}(kk,:) = [xout(kk),yout(kk)];
%             end
%         end
%         
%         h_real.ZData = xi;    h_imag.ZData = eta;
%         
%         if num_cargas ~= 0
%             posx = p_charges{iz}(:,1); posy = p_charges{iz}(:,2);
%             set(positive, 'XData', posx, 'YData', posy);
%         else
%             posx = xx; posy = xx;
%             set(positive, 'XData', posx, 'YData', posy);
%         end
%         
%         knot = zeros(N);
%         if num_cargas ~= 0
%             for mm = 1:size(posx,1)
%                 xp = posx(mm,1); yp = posy(mm,1);
%                 knot = knot + Gaussian(xp,yp);
%             end
%         else
%             knot = zeros(N);
%         end
%         
%         %     img1.CData = abs(knot).^2;
%         img2.CData = abs(Uz).^2;
%         
%         %     Efield(:,:,iz) = Uz;
%         isoI(:,:,iz) = abs(knot).^2;
%         drawnow
%         %     pause(eps)
    end
    
%     name = sprintf('hopfnoturb_w_%d_iso.mat',w);
%     save(name,'isoI','-v7');
    
    %         name = sprintf('3foil_w_%d_SR_%d_num_%d.mat',w,SR,ii);
    %         save(name,'p_charges')
    
    name = sprintf('..\\data\\test\\Efield_%d_SR_%d.mat',ii,SR);
    save(name,'Uz')
    
    clearvars -except N w ii num
    close all;
end
% end
toc

%% LG Functions
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