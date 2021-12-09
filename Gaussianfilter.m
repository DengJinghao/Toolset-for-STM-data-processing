%针对一维数组的高斯滤波, 两端的扩充矩阵元取值为数组的start与end值
%Gaussianfilter(data, r, coarsening_length)
%key parameters:
%length: coarsening length in terms of pixel
%r: integrating scale, r = 4*length, usually
%   for r = 4*length, the templet contains 95.45% of the total area
%   for r = 6*length, the templet contains 99.73% of the total area
%data: experimental data for filtering process
function y_filted = Gaussianfilter(data, r, coarsening_length)

%creating templet of 1D gaussian function
for x = 1:2*r-1
    gaussianfunc(x) = (1./(sqrt(2*pi).*coarsening_length)).*exp(-(x-r).^2./(2.*coarsening_length.^2));
end

%filtering process
y_filted = data
for x = 1:r-1
    zer = data(1).*ones(1,r-x);
    data_cut = data(1:x+r-1);
    A = [zer data_cut];
    y_filted(x) = A*gaussianfunc.'
end
for x = r:length(data)-r+1
    y_filted(x) = data(x-r+1 : x+r-1)*gaussianfunc.'
end
for x = length(data)-r+2:length(data)
    zer = data(length(data)).*ones(1,r-1-length(data)+x);
    data_cut = data(x-r+1 : length(data));
    A = [data_cut zer];
    y_filted(x) = A*gaussianfunc.'
end
end
