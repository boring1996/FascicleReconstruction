%AUTHOR:
    %DANIEL TOVBIS (2019)
%DESCRIPTION:
    %This function takes in the registered images and their ground truths (from %regist4x75)
    %and the desired neural network and does segmentation
    %When segmenting, do not use a network that was trained on the same nerve slice
    %e.g.,if segmenting prox1, use the prox2dist3 network
    %The CreateMaskforallimages_allmethods script should populate all these
    %values by default.
%INPUTS:
    %images: Array of registered images (for a single segment)
    %network: Network used for detetection
    %gtsreg: Array of registered images (for a single segment)
%OUTPUTS:
    %maskedImages: An intermediate result; Array of images with connective tissue masked out
    %calccontour: Contour after segmentation and detection
    %meanIOU: Mean IOU for the segment
    %stdIOU: Standard Deviation of IOU for the segment
function [maskedImages,calccontour,meanIOU,stdIOU] = NN_ACoutsidein(images,network,gtsreg)
tic
%% Detection, Mask Generation
numimages=size(images,2);
annotation = sprintf('fascicle'); %Annotation label (not strictly necessary)
mask=cell(1,numimages); %Initialize Mask Matrix
colourmat=[221;219;220]; %Colour of background
for j=1:numimages
    testImage = images{j};
    [bboxes,~,~] = detect(network,testImage); %do detection
    outputImage = insertObjectAnnotation(testImage, 'rectangle', bboxes, annotation,'LineWidth',3,'Color','Black','TextColor','White');
    figure
    imshow(outputImage)
    mask{j}=false(size(outputImage(:,:,1))); %initialize masks
    for i=1:size(bboxes)
        bwmaskfin{j}=false(size(outputImage(:,:,1)));
        strelsize=round(0.4*sqrt((bboxes(i,3))^2+(bboxes(i,4))^2));
        bwrect{i}=imrect(gca,bboxes(i,:));
        bwmask{i}=createMask(bwrect{i});
        props=regionprops(bwmask{i},'Centroid');
        centre=floor(props.Centroid);
        bwmaskfin{j}(centre(2),centre(1))=1;
        mask{j}=or(mask{j},imdilate(bwmaskfin{j},strel('disk',strelsize,0))); %generate circles
        maskedImages{j} = bsxfun(@times, images{j}, cast(mask{j},class(images{j}))); %remove connective tissue
        for k=1:3 %Change black background to greyish
            temp=maskedImages{j}(:,:,k);
            temp(temp==0)=colourmat(k);
            maskedImages{j}(:,:,k)=temp;
        end
    end
    close
end
%% Segmentation
maxIterations = 500; %Iterations for Segmentation
for k=1:numimages
    disp(['Processing Image ' num2str(k)])
    Iobrcbr2=improcess(maskedImages{k},7,7,0); %Process the image
    contractionbias(k)=0.5; %set contraction bias
    calccontour{k}=activecontour(Iobrcbr2,mask{k},maxIterations,'Chan-Vese','smoothfactor',0.6,'contractionbias',contractionbias(k));
    calccontour{k}=imfill(calccontour{k},'holes'); %fill holes
    calccontour{k}=bwareaopen(calccontour{k},500); %delete small objects
    CC_diameters{k}=bwconncomp(calccontour{k});
    diameters{k}=regionprops(CC_diameters{k},'MinorAxisLength');
    smallestdiameter=min([diameters{k}.MinorAxisLength]);
    largeststrelsize=floor(0.5*smallestdiameter); %get smallest diameter in image
    calccontour{k}=imopen(calccontour{k},strel('disk',largeststrelsize,0)); %image opening to split fascicles
end
%% IOU calculation
for k=1:numimages
    int=and(calccontour{k},gtsreg{k});
    uni=or(calccontour{k},gtsreg{k});
    intoveruni(k)=sum(int(:))/sum(uni(:));
end
meanIOU=mean(intoveruni)
stdIOU=std(intoveruni)
%%
elapsed=toc
end