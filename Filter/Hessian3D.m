
function [Dxx, Dyy, Dzz, Dxy, Dxz, Dyz] = Hessian3D(Volume,Sigma)
%  This function Hessian3D filters the image with an Gaussian kernel
%  followed by calculation of 2nd order gradients, which aprroximates the
%  2nd order derivatives of the image.
%  该函数用高斯核对图像进行滤波，然后计算二阶梯度，近似地逼近图像的二阶导数。
% 
% [Dxx, Dyy, Dzz, Dxy, Dxz, Dyz] = Hessian3D(I,Sigma)
% 
% inputs,
%   I : The image volume, class preferable double or single 
%   %double:8byte 双精度浮点数 single:4byte 单精度浮点数
%   Sigma : The sigma of the gaussian kernel used. If sigma is zero
%           no gaussian filtering.  高斯核中的sigma
%
% outputs,
%   Dxx, Dyy, Dzz, Dxy, Dxz, Dyz: The 2nd derivatives 二阶导数
%
% Function is written by D.Kroon University of Twente (June 2009)
% defaults
if nargin < 2, Sigma = 1; end  %nargin - 函数输入参数数目

if(Sigma>0)
    F=imgaussian(Volume,Sigma);
%     F=Volume;
else
    F=Volume;
end

% Create first and second order diferentiations 求一阶导数和二阶导数
Dz=gradient3(F,'z');
Dzz=(gradient3(Dz,'z'));
clear Dz;

Dy=gradient3(F,'y');
Dyy=(gradient3(Dy,'y'));
Dyz=(gradient3(Dy,'z'));
clear Dy;

Dx=gradient3(F,'x');
Dxx=(gradient3(Dx,'x'));
Dxy=(gradient3(Dx,'y'));
Dxz=(gradient3(Dx,'z'));
clear Dx;

function D = gradient3(F,option)
% This function does the same as the default matlab "gradient" function
% but with one direction at the time, less cpu and less memory usage.
% 和matlab自带函数gradient一样，但是使用cpu和内存较少
% Example:
%
% Fx = gradient3(F,'x');

[k,l,m] = size(F);
D  = zeros(size(F),class(F)); 

switch lower(option)
case 'x'
    % Take forward differences on left and right edges
    D(1,:,:) = (F(2,:,:) - F(1,:,:));
    D(k,:,:) = (F(k,:,:) - F(k-1,:,:));
    % Take centered differences on interior points
    D(2:k-1,:,:) = (F(3:k,:,:)-F(1:k-2,:,:))/2;
case 'y'
    D(:,1,:) = (F(:,2,:) - F(:,1,:));
    D(:,l,:) = (F(:,l,:) - F(:,l-1,:));
    D(:,2:l-1,:) = (F(:,3:l,:)-F(:,1:l-2,:))/2;
case 'z'
    D(:,:,1) = (F(:,:,2) - F(:,:,1));
    D(:,:,m) = (F(:,:,m) - F(:,:,m-1));
    D(:,:,2:m-1) = (F(:,:,3:m)-F(:,:,1:m-2))/2;
otherwise
    disp('Unknown option')
end
        