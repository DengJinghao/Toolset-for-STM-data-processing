clear variables;
close all

a_ref1 = 0.369;%lattice constant of VCl3 on NbSe2, in terms of nm
k1 = 2*pi/(a_ref1);

atomicposition = [];
PathWay = 'E:\onedrive\OneDrive - whu.edu.cn\STM data\LHE\LHE\20211227crcl3\';
BaseName = 'Vcl300361';
%choose the pathway and basename by replacing the variant
Rawdata = readmatrix(strcat(PathWay,BaseName,'.txt'));
%change files into variables in forms of 'double'

[L_Y,L_X] = size(Rawdata);
LenX = 11.07;
LenY = 3.29;
%nm, the actual size of the image
X = 0:LenX/(L_X-1):LenX;
Y = 0:LenY/(L_Y-1):LenY;

fftim = fftshift(fft2(Rawdata));%计算2d fft并将其零频移动至中心
fftimmodulus = abs(fftim);
qxmax = 0.5/(LenY/(L_Y-1));
qymax = 0.5/(LenX/(L_X-1));
qy = -qymax:1/LenY:qymax;
qx = -qxmax:1/LenX:qxmax;
%define k-space or q-space

figure(1);
imagesc(Rawdata)
title('click to choose the starting point')
axis equal
[xzero,yzero] = ginput(1);
waitforbuttonpress %pause the program and wait until you press the buttom
close(1)

Rawdatashift = Rawdata(yzero:L_Y,xzero:L_X);
[L_Ycut,L_Xcut] = size(Rawdatashift);
Xcut = X(1:L_Xcut);
Ycut = Y(1:L_Ycut);
%re-define the actual size of the cutted image

figure(1);
imagesc(Xcut,Ycut,Rawdatashift)
title('Rawdata');
axis equal

%theta = 90;%twist angle betweem reference lattice and horizonal aligned lattice
Zamp = max(max(Rawdata));
for ii = 1:L_Xcut
    for jj = 1:L_Ycut
        T1(jj,ii) = Zamp*exp(-sqrt(-1).*k1.*(Xcut(ii)+Ycut(jj).*0));
    end
end
T = T1;%+ Tinner1 + Tinner2 + Tinner3;%reference lattice matrix
P_ref = double(imregionalmax(real(T)));
P_real = double(imregionalmax(Rawdatashift));
%find local maxima, which indicates the atomic positions
%figure(5);
%hold on
%imagesc(P_ref)
%scatter(P_)
%hold off

figure(2);
imagesc(Xcut,Ycut,real(T))
title('reference lattice');
axis equal

tmulti =imag(T).*Rawdatashift;%multiply reference lattice and real lattice
tmulti2 = real(T).*Rawdatashift;

%r = 240;
%coarsening_length = 60;
%for kk = 1:L_Ycut
%    data = tmulti(kk,:);
%    udisp_field1(kk,:) = Gaussianfilter(data,r,coarsening_length);
%    udisp_field1(kk,:) = lowpass(data,0.5,2);
%    data = tmulti2(kk,:);
%    udisp_field2(kk,:) = lowpass(data,0.5,2);   
%end
%udisp_field1 = imgaussfilt(tmulti,60);
%udisp_field2 = imgaussfilt(tmulti2,60);
tmulti_fft = (fft(tmulti,[],2));
tmulti_fft(:,10:852) = 0;
udisp_field1 = abs(ifft(tmulti_fft,862,2));

tmulti_fft2 = (fft(tmulti2,[],2));
tmulti_fft2(:,10:852) = 0;
udisp_field2 = abs(ifft(tmulti_fft2,862,2));
%windowed fft to filter the atomic fluctuations

phase = atan(udisp_field1./udisp_field2);
[strainrawx,strainrawy] = gradient(phase,(LenX/(L_X-1)).*k1,(LenY/(L_Y-1)).*k1);
figure(3);
imagesc(Xcut,Ycut,tmulti)
title('image*raw');
axis equal

figure(4);
imagesc(Xcut,Ycut,udisp_field1)
title('filtered image*raw')
axis equal

figure(5);
imagesc(Xcut,Ycut,phase)
title('phase')
axis equal

figure(6);
imagesc(Xcut,Ycut,strainrawx)
title('strainy')
axis equal

figure(7);
imagesc(Xcut,Ycut,udisp_field1)
title('strain')
axis equal

save ycrcl2driftfield.txt -ascii strainrawx
