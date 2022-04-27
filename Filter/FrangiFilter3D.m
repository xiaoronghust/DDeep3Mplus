function [Iout,whatScale,Voutx,Vouty,Voutz]=FrangiFilter3D(I,options)
%
% This function FRANGIFILTER3D uses the eigenvectors of the Hessian to
% compute the likeliness of an image region to vessels, according
% to the method described by Frangi
% 这个函数FRANGIFILTER3D使用Hessian的特征向量，按照Frangi描述的方法计算图像区域对血管的概率

% [J,Scale,Vx,Vy,Vz] = FrangiFilter3D(I, Options)
%
% inputs,
%   I : The input image volume (vessel volume)
%   Options : Struct with input options,
%       .FrangiScaleRange : The range of sigmas used, default [1 8]
%       .FrangiScaleRatio : Step size between sigmas, default 2
%       .FrangiAlpha : Frangi vesselness constant, treshold on Lambda2/Lambda3
%					   determines if its a line(vessel) or a plane like structure
%					   default .5;第二个比率的阈值，Lambda2/Lambda3，决定了是线状斑点还是板状斑点
%       .FrangiBeta  : Frangi vesselness constant, which determines the deviation
%					   from a blob like structure, default .5;决定了从水滴状结构的误差，默认值0.5;
%       .FrangiC     : Frangi vesselness constant which gives
%					   the threshold between eigenvalues of noise and 
%					   vessel structure. A thumb rule is dividing the 
%					   the greyvalues of the vessels by 4 till 6, default 500;
%                      噪声特征值与血管结构之间的阈值。经验法则是将血管的灰度值除以4到6，默认为500;
%       .BlackWhite : Detect black ridges (default) set to true, for
%                       white ridges set to false.检测黑色信号(默认)设置为真，白色信号设置为假。
%       .verbose : Show debug information, default true %显示调试信息，默认为真
%
% outputs,
%   J : The vessel enhanced image (pixel is the maximum found in all scales)血管增强图像(像素是所有尺度中发现的最大值)
%   Scale : Matrix with the scales on which the maximum intensity 
%           of every pixel is found具有每个像素的最大强度的尺度的矩阵
%   Vx,Vy,Vz: Matrices with the direction of the smallest eigenvector, pointing
%				in the direction of the line/vessel.矩阵具有最小特征向量的方向，指向线/血管的方向。
%
% Literature, 
%	Manniesing et al. "Multiscale Vessel Enhancing Diffusion in 
%		CT Angiography Noise Filtering"
%
% Example,
%   % compile needed mex file
%   mex eig3volume.c
%
%   load('ExampleVolumeStent');
%   
%   % Frangi Filter the stent volume
%   options.BlackWhite=false;
%   options.FrangiScaleRange=[1 1];
%   Vfiltered=FrangiFilter3D(V,options);
%
%   % Show maximum intensity plots of input and result
%   figure, 
%   subplot(2,2,1), imshow(squeeze(max(V,[],2)),[])
%   subplot(2,2,2), imshow(squeeze(max(Vfiltered,[],2)),[])
%   subplot(2,2,3), imshow(V(:,:,100),[])
%   subplot(2,2,4), imshow(Vfiltered(:,:,100),[])
%
% Written by D.Kroon University of Twente (May 2009)

% Constants vesselness function

% 设置默认参数
defaultoptions = struct('FrangiScaleRange', [1 10], 'FrangiScaleRatio', 2, 'FrangiAlpha', 0.3, 'FrangiBeta', 0.5, 'FrangiC', 1000, 'verbose',true,'BlackWhite',true);

% Process inputs
if(~exist('options','var')) % 每一项都不存在
    options=defaultoptions; 
else %其中存在某些项
    tags = fieldnames(defaultoptions); %返回字符矢量元胞数组，这些字符矢量包含结构体s中的字段名称
    for i=1:length(tags)
         if(~isfield(options,tags{i})),  options.(tags{i})=defaultoptions.(tags{i}); end %如果不包含其中某些参数，就用默认参数替代
    end
    if(length(tags)~=length(fieldnames(options)))  %出现未定义的参数
        warning('FrangiFilter3D:unknownoption','unknown options found');
    end
end

% Use single or double for calculations
if(~isa(I,'double')), I=single(I); end %如果a不是double类型 

sigmas=options.FrangiScaleRange(1):options.FrangiScaleRatio:options.FrangiScaleRange(2); %在(1)-(2)之间，隔ratio取多少
sigmas = sort(sigmas, 'ascend'); %沿其大小不为 1 的第一个数组维度以升序对 A 元素进行排序。

% Frangi filter for all sigmas
for i = 1:length(sigmas)
    % Show progress
    if(options.verbose)
        disp(['Current Frangi Filter Sigma: ' num2str(sigmas(i)) ]);
    end
    
    % Calculate 3D hessian
    [Dxx, Dyy, Dzz, Dxy, Dxz, Dyz] = Hessian3D(I,sigmas(i));

    if(sigmas(i)>0)
        % Correct for scaling
        c=(sigmas(i)^2);
        Dxx = c*Dxx; Dxy = c*Dxy;
        Dxz = c*Dxz; Dyy = c*Dyy;
        Dyz = c*Dyz; Dzz = c*Dzz;
    end
    
    % Calculate eigen values 计算特征值
    if(nargout>2)
        [Lambda1,Lambda2,Lambda3,Vx,Vy,Vz]=eig3volume(Dxx,Dxy,Dxz,Dyy,Dyz,Dzz);
    else
        [Lambda1,Lambda2,Lambda3]=eig3volume(Dxx,Dxy,Dxz,Dyy,Dyz,Dzz);
    end
    
    % Free memory
    clear Dxx Dyy  Dzz Dxy  Dxz Dyz;

    % Calculate absolute values of eigen values  计算特征值的绝对值
    LambdaAbs1=abs(Lambda1);
    LambdaAbs2=abs(Lambda2);
    LambdaAbs3=abs(Lambda3);

    % The Vesselness Features
    Ra=LambdaAbs2./LambdaAbs3;
    Rb=LambdaAbs1./sqrt(LambdaAbs2.*LambdaAbs3);

    % Second order structureness. S = sqrt(sum(L^2[i])) met i =< D
    S = sqrt(LambdaAbs1.^2+LambdaAbs2.^2+LambdaAbs3.^2);
    A = 2*options.FrangiAlpha^2; B = 2*options.FrangiBeta^2;  C = 2*options.FrangiC^2;
	
    % Free memory
    clear LambdaAbs1 LambdaAbs2 LambdaAbs3

    %Compute Vesselness function
    expRa = (1-exp(-(Ra.^2./A)));
    expRb =    exp(-(Rb.^2./B));
    expS  = (1-exp(-S.^2./(2*options.FrangiC^2)));
    % keyboard
    % Free memory
    clear S A B C Ra Rb

    %Compute Vesselness function
    % Voxel_data = expRa.* expRb.* expS;
    Voxel_data = expRa.* expS;
    
    % Free memory
    clear expRa expRb expRc;
    
    if(options.BlackWhite)
        Voxel_data(Lambda2 < 0)=0; Voxel_data(Lambda3 < 0)=0;
    else
        Voxel_data(Lambda2 > 0)=0; Voxel_data(Lambda3 > 0)=0;
    end
        
    % Remove NaN values
    Voxel_data(~isfinite(Voxel_data))=0;

    % Add result of this scale to output
    if(i==1)
        Iout=Voxel_data;
        if(nargout>1)
            %函数输出参数的数目。
            whatScale = ones(size(I),class(Iout));
        end
        if(nargout>2)
            Voutx=Vx; Vouty=Vy; Voutz=Vz;
        end
    else
        if(nargout>1)
            whatScale(Voxel_data>Iout)=i;
        end
        if(nargout>2)
            Voutx(Voxel_data>Iout)=Vx(Voxel_data>Iout);
            Vouty(Voxel_data>Iout)=Vy(Voxel_data>Iout);
            Voutz(Voxel_data>Iout)=Vz(Voxel_data>Iout);
        end
        % Keep maximum filter response
        Iout=max(Iout,Voxel_data);
    end   
end

