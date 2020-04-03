%AUTHOR:
    %DANIEL TOVBIS (2019)
%DESCRIPTION:
    % This function will register images using four consecutive registrations.
    %"The Grand Image Registration Script" has the fields of this function
    %populated by default to process the entire dataset.
%INPUTS:
    %directory (string): directory where image files are stored. By default
    %these are stored in Nerve Images\directory
    %Image files must be numbered sequentially (1.tif,2.tif, etc)
    %Ground truths should be stored in the folder *directory\Ground Truths, with
    %the filename *directorygts.mat
    %initialimage (integer): The initial image file to start registering from
    %numimages (integer): Total number of images to register (including initial image)
    %isihc (binary) : Whether IHC images are being registered or not (isihc=1 for IHC, =0 for HE)
%Outputs are:
    %regim: array of registered images
    %SSIMorig: Array of SSIMs for original, unregistered images
    %SSIMreg: Array of SSIMs for registered images
    %MSEorig: Array of MSEs for original, unregistered images
    %MSEreg: Array of MSEs for registered images
    %gtsreg: Registered ground truth array (used in subsequent processing steps)
function [regim,SSIMorig,SSIMreg,MSEorig,MSEreg,gtsreg]=regist4x75(directory,initialimage,numimages,isihc)
tic
currentdir=pwd;
%% Registration Parameters
%Set up registration parameters
regim=cell(1,numimages);
[optimizer,metric] = imregconfig('monomodal'); %Images are monomodal, that is, both taken with the same modaliy (histology)
optimizer.MinimumStepLength = 1e-6;
optimizer.MaximumStepLength = 0.00625;
optimizer.RelaxationFactor = 0.6;
optimizer.MaximumIterations=250; %Number of iterations for registration
%% Registration Parameters (affine)
%Different parameters for affine registration (a bit more conservative to reduce distortion)
[optimizeraffine,metricaffine] = imregconfig('monomodal'); %Images are monomodal, that is, both taken with the same modaliy (histology)
optimizeraffine.MinimumStepLength = 1e-7;
optimizeraffine.MaximumStepLength = 0.001675;
optimizeraffine.RelaxationFactor = 0.7;
optimizeraffine.MaximumIterations=200;
%% Generate origim and regim array
alltruths=importdata(strcat(currentdir,'\Nerve Images\', directory, '\Ground Truths\', directory, 'gts.mat')); %Get all ground truths for current nerve
for i=1:numimages
    origim{i}=imread(strcat(currentdir,'\Nerve Images\', directory, '\',num2str(i+(initialimage-1)), '.tif')); %An array of unregistered images
    regim{i}=imread(strcat(currentdir,'\Nerve Images\', directory, '\',num2str(i+(initialimage-1)), '.tif'));  %Also an array of unregistered images (but this will be changed later)
    %% Generate ground truth array
    % Get array of only ground truths that are relevant
    gtsreg{i}=alltruths{i+(initialimage-1)};
end
%% Do registration
if isihc==1
    for i=1:numimages-1 %Minus 1 So we don't repeat the loop for the last image as a "fixed" image
        disp(['Registering Image ' num2str(i+1)])
        contour=improcessihc(regim{i}); %Process IHC image
        fixed=255.*uint8(contour); %convert to uint8
        contour2=improcessihc(regim{i+1}); %Process next image
        moving =255.*uint8(contour2); %Process the moving image (initially unregistered)
        groundtruth=gtsreg{i+1};
        % First registration (translation only)
        fixed25=imresize(fixed,0.25); %Downscale 25%
        moving25=imresize(moving,0.25);
        tform25= imregtform(moving25, fixed25, 'translation', optimizer, metric,'DisplayOptimization',false,'PyramidLevels',4); %Get the transformation function
        moving_t=imwarp(moving,tform25,'OutputView',imref2d(size(fixed)),'FillValues',234); %Warp the previously registered image to create a new rgb image
        groundtruth_t=imwarp(groundtruth,tform25,'OutputView',imref2d(size(fixed))); %Warp the ground truth
        original_t= imwarp(regim{i+1},tform25,'OutputView',imref2d(size(fixed)),'FillValues',[221;219;220]);
        % Second registration (translation + rotation)
        fixed50=imresize(fixed,0.5); %Downscale 50%
        moving50=imresize(moving_t,0.5);
        tform50 = imregtform(moving50, fixed50, 'rigid', optimizer, metric,'DisplayOptimization',false,'PyramidLevels',4); %Get the transformation function
        moving_tr=imwarp(moving_t,tform50,'OutputView',imref2d(size(fixed)),'FillValues',234); %Warp the previously registered image to create a new rgb image
        groundtruth_tr= imwarp(groundtruth_t,tform50,'OutputView',imref2d(size(fixed)));
        original_tr= imwarp(original_t,tform50,'OutputView',imref2d(size(fixed)),'FillValues',[221;219;220]);
        % Third registration (translation + rotation + scale)
        fixed75=fixed; %No downscaling for IHC images (variable name is the same for convenience)
        moving75=moving_tr;
        tform75 = imregtform(moving75, fixed75, 'similarity', optimizer, metric,'DisplayOptimization',false,'PyramidLevels',4); %Get the transformation function
        moving_trs=imwarp(moving_tr,tform75,'OutputView',imref2d(size(fixed)),'FillValues',234);
        groundtruth_trs= imwarp(groundtruth_tr,tform75,'OutputView',imref2d(size(fixed)));
        original_trs= imwarp(original_tr,tform75,'OutputView',imref2d(size(fixed)),'FillValues',[221;219;220]);
        % Fourth registration (translation + rotation + scale +shear)
        tform100 = imregtform(moving_trs, fixed, 'affine', optimizeraffine, metricaffine,'DisplayOptimization',false,'PyramidLevels',4); %Get the transformation function
        ifreg=imwarp(original_trs,tform100,'OutputView',imref2d(size(fixed)),'FillValues',[221;219;220]);
        gtsreg{i+1}=imwarp(groundtruth_trs,tform100,'OutputView',imref2d(size(fixed))); %Warp the ground truth
        regim{i+1}=ifreg; %Replace the image
    end
else %For H&E images
    for i=1:numimages-1
        disp(['Registering Image ' num2str(i+1)])
        fixed =improcess(regim{i},7,7); % Process the fixed image (already registered)
        moving =improcess(regim{i+1},7,7); %Process the moving image (initially unregistered)
        groundtruth=gtsreg{i+1};
        % First registration (translation only)
        fixed25=imresize(fixed,0.25); %Downscale 25%
        moving25=imresize(moving,0.25);
        tform25= imregtform(moving25, fixed25, 'translation', optimizer, metric,'DisplayOptimization',false,'PyramidLevels',4); %Get the transformation function
        moving_t=imwarp(moving,tform25,'OutputView',imref2d(size(fixed)),'FillValues',234); %Warp the previously registered image to create a new rgb image
        groundtruth_t=imwarp(groundtruth,tform25,'OutputView',imref2d(size(fixed)));
        original_t= imwarp(regim{i+1},tform25,'OutputView',imref2d(size(fixed)),'FillValues',[221;219;220]);
        % Second registration (translation + rotation)
        fixed50=imresize(fixed,0.5); %Downscale 50%
        moving50=imresize(moving_t,0.5);
        tform50 = imregtform(moving50, fixed50, 'rigid', optimizer, metric,'DisplayOptimization',false,'PyramidLevels',4); %Get the transformation function
        moving_tr=imwarp(moving_t,tform50,'OutputView',imref2d(size(fixed)),'FillValues',234); %Warp the previously registered image to create a new rgb image
        groundtruth_tr= imwarp(groundtruth_t,tform50,'OutputView',imref2d(size(fixed)));
        original_tr= imwarp(original_t,tform50,'OutputView',imref2d(size(fixed)),'FillValues',[221;219;220]);
        % Third registration (translation + rotation + scale)
        fixed75=imresize(fixed,0.75); %Downscale 75%
        moving75=imresize(moving_tr,0.75);
        tform75 = imregtform(moving75, fixed75, 'similarity', optimizer, metric,'DisplayOptimization',false,'PyramidLevels',4); %Get the transformation function
        moving_trs=imwarp(moving_tr,tform75,'OutputView',imref2d(size(fixed)),'FillValues',234);
        groundtruth_trs= imwarp(groundtruth_tr,tform75,'OutputView',imref2d(size(fixed)));
        original_trs= imwarp(original_tr,tform75,'OutputView',imref2d(size(fixed)),'FillValues',[221;219;220]);
        % Fourth registration (translation + rotation + scale +shear
        tform100 = imregtform(moving_trs, fixed, 'affine', optimizeraffine, metricaffine,'DisplayOptimization',false,'PyramidLevels',4); %Get the transformation function
        ifreg=imwarp(original_trs,tform100,'OutputView',imref2d(size(fixed)),'FillValues',[221;219;220]);
        gtsreg{i+1}=imwarp(groundtruth_trs,tform100,'OutputView',imref2d(size(fixed))); %Warp the ground truth
        regim{i+1}=ifreg; %Replace the image
    end
end
%% Check SSIM Values
%Calculate statistics
for i=2:numimages
    if isihc==1
        ref=origim{i-1};
        test=origim{i};
        refreg=regim{i-1};
        testreg=regim{i};
    else
        ref=improcess(origim{i-1},7,7);
        test=improcess(origim{i},7,7);
        refreg=improcess(regim{i-1},7,7);
        testreg=improcess(regim{i},7,7);
    end
    SSIMorig(i)=ssim(test,ref);
    SSIMreg(i)=ssim(testreg,refreg);
end
mean(SSIMorig(2:numimages))
std(SSIMorig(2:numimages))
mean(SSIMreg(2:numimages))
std(SSIMreg(2:numimages))
%% Check MSE values
for i=2:numimages
    if isihc==1
        ref=origim{i-1};
        test=origim{i};
        refreg=regim{i-1};
        testreg=regim{i};
    else
        ref=improcess(origim{i-1},7,7);
        test=improcess(origim{i},7,7);
        refreg=improcess(regim{i-1},7,7);
        testreg=improcess(regim{i},7,7);
    end
    MSEorig(i)=immse(test,ref);
    MSEreg(i)=immse(testreg,refreg);
end
mean(MSEorig(2:numimages))
std(MSEorig(2:numimages))
mean(MSEreg(2:numimages))
std(MSEreg(2:numimages))
toc
end