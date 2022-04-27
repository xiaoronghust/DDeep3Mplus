clc
clear

% I1=double(imread ('vessel.png'));
% Ivessel=FrangiFilter2D(I1);
% figure,
% subplot(1,2,1), imshow(I1,[]);
% subplot(1,2,2), imshow(Ivessel,[0 0.25]);

% compile needed mex file
% mex eig3volume.c

%  load('ExampleVolumeStent');
%   
%  % Frangi Filter the stent volume
%  options.BlackWhite=false;
%  options.FrangiScaleRange=[1 1];
%  Vfiltered=FrangiFilter3D(V,options);
% 
%  % Show maximum intensity plots of input and result
%  figure(1), 
%  subplot(2,2,1), imshow(squeeze(max(V,[],2)),[])
%  subplot(2,2,2), imshow(squeeze(max(Vfiltered,[],2)),[])
%  subplot(2,2,3), imshow(V(:,:,100),[])
%  subplot(2,2,4), imshow(Vfiltered(:,:,100),[])

% I=double(imread ('XHQ1um224.tif'));
% subplot(1,2,1), imshow(I,[]);

% fname='18459\VED\15-ved-20200710.tif';
fname='G:\PAPER DATA\Synthetic(ground truth)\p_checked6_mouse_RGC_uw\091201c1-downsample\091201c1.tif.v3dpbd-downsample-x5y5times-x435y429.tif';
info = imfinfo(fname); %����һ���ṹ�壬�ýṹ����ֶΰ����й�ͼ���ļ� filename �е�ͼ�����Ϣ��
num_images = numel(info); %�������� A �е�Ԫ����Ŀ n
for i=1:num_images
    a= imread(fname,i);
    I(:,:,i)=a;
end


%% Filter the stent volume.
options.BlackWhite=false;
options.FrangiScaleRange=[1 1];

% options.BlackWhite=false;%����ɫ
% options.EnhancementScaleRange=[1 5];
% options.EnhancementType = 'zhou';

%  Ifiltered=EnhancementFilter3D(I,options);
% t1 = clock;
Ifiltered=FrangiFilter3D(I,options);
% t2 = clock;
% etime(t1,t2)
% Ifiltered=vesselness3D(I, 1:5,  0.75, true);
% figure(1), 
% subplot(2,2,1), imshow(squeeze(max(I,[],1)),[])
% subplot(2,2,2), imshow(squeeze(max(Ifiltered,[],1)),[])
% subplot(2,2,3), imshow(I(:,:,150),[])
% subplot(2,2,4),imshow(Ifiltered(:,:,140),[])
% I140 = Ifiltered(:,:,140);
% 
% Pro_of_max_of_origin = squeeze(max(I,[],2));
% Pro_of_max_of_origin_adjust = imadjust(Pro_of_max_of_origin)*255000;
% Pro_of_max_of_result = double(squeeze(max(Ifiltered,[],2))*255000);
% Pro_of_max_of_result_adjust = imadjust(Pro_of_max_of_result)*255000;
% % imwrite(Pro_of_max_of_origin,'18459-z11-y22-15_22_11-Pro_of_max_of_origin.tif')
% % imwrite(Pro_of_max_of_result,'18459-z11-y22-15_22_11-Pro_of_max_of_result.tif')
% figure(1),
% subplot(2,2,1), imshow(Pro_of_max_of_origin,[]),title('origin max')
% subplot(2,2,2), imshow(Pro_of_max_of_result,[]),title('result max')
% subplot(2,2,3), imshow(Pro_of_max_of_origin_adjust,[]),title('origin max adjust')
% subplot(2,2,4), imshow(Pro_of_max_of_result_adjust,[]),title('result max adjust')

% 
% %%������
% D_Threshold = 50;
% %������������������������������������������������������������������������������������������������������������������������
% Ifiltered = Ifiltered*25500;
% Ifiltered(D_high>D_Threshold) = Volume(D_high>D_Threshold);


fname='*.tif';
% img_name = 'D:\Liaojiangshan\Code\frangiFilter2D\18459\18459-z11-y22-15_22_11-Frangi-sigma12.tif';
img_name = 'G:\D trans\Data-Analysis\done_neuron4\neuron4-Frangi-sigma11-alpha0.3-c1000-multiply2550000-20210401.tif';
num_images = size(Ifiltered,3);
for i=1:num_images
% % ������˰���Ŀ׶�   
%     J = Ifiltered(:,:,i);
%  %   J = imfill(J);
%     J = J *25500;
    J = Ifiltered(:,:,i)*2550000;
    imwrite(uint8(J),img_name,'WriteMode','append');
end

