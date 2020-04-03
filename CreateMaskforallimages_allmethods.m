%AUTHOR:
    %DANIEL TOVBIS (2019)
%DESCRIPTION:
    %This script contains a variety of functions for detecting/segmenting
    %fascicles in H&E and IHC images, as well as otsu's method and k-means
    %detection (used for comparison purposes).
    %You can run the appropriate section after generating the required sets (by using regist4x75) and
    %ensuring the proper variables are loaded (the neural networks)
    %All loops are set to run assuming the full data set has been registered.
    %Only run one section at a time, then save the variables (as they use the same variable names and
    %some information will be overwritten).
%% Processing for H&E slices
%INPUTS:
    %network: A trained neural network (.mat) (from trainFascicleRCNN_Noaugment_NSG)
    %regim: Cell array of registered images (from regist4x75)
    %gtsreg: Cell array of registered ground truths (from regist4x75)
%OUTPUTS:
    %calccontour_oi: Cell array of segmented images
    %meanIOU_oi: Array of mean IOU per segment.
    %stdIOU_oi: Array of standard deviation of IOU per segment.
    %detectiontime: Total processing time
currentdir=pwd;
for j=1:11
    if j==1
        network=fasterrcnnvgg16noaugment_prox2dist3;
    elseif j==6
        clear network
        network=fasterrcnnvgg16noaugment_prox1dist3;
    elseif j==9
        clear network
        network=fasterrcnnvgg16noaugment_prox1prox2;
    end
    tic
    disp(['Creating Mask For Segment ' num2str(j)])
    [maskedImages_oi{j},calccontour_oi{j},meanIOU_oi(j),stdIOU_oi(j)] = NN_ACoutsidein(regim{j},network,gtsreg{j});
    detectiontime(j)=toc
end
%% Processing for IHC slices
%INPUTS:
    %regim: Cell array of registered images (from regist4x75)
    %gtsreg: Cell array of registered ground truths (from regist4x75)
%OUTPUTS:
    %calccontour_ihc: Cell array of segmented images
    %meanIOU_ihc: Array of mean IOU per segment.
    %stdIOU_ihc: Array of standard deviation of IOU per segment.
    %detectiontime_ihc: Total processing time
for j=1:10
    tic
    disp(['Creating Mask For Segment ' num2str(j)])
    currentimages=regim{j};
    currentgts=gtsreg{j};
    numimages=length(currentimages);
    for k=1:numimages
        disp(['Processing Image ' num2str(k)])
        currentcontour{k}=improcessihc(currentimages{k});
        int=and(currentcontour{k},currentgts{k});
        uni=or(currentcontour{k},currentgts{k});
        intoveruni(j,k)=sum(int(:))/sum(uni(:));
    end
    detectiontime_ihc(j)=toc;
    meanIOU_ihc(j)=mean(intoveruni(j,1:numimages));
    stdIOU_ihc(j)=std(intoveruni(j,1:numimages));
    calccontour_ihc{j}=currentcontour;
    clear currentcontour
end
%% Otsu's Method
%INPUTS:
    %gtsreg: Cell array of registered ground truths (from regist4x75)
    %maskedImages (Otsu's method/Kmeans): cell array of images with connective tissue removed (from NN_ACoutsidein)
%OUTPUTS:
    %calccontour_otsu: Cell array of segmented images
    %meanIOU_otsu: Array of mean IOU per segment.
    %stdIOU_otsu: Array of standard deviation of IOU per segment.
    %detectiontime_otsu: Total processing time
currentdir=pwd;
for j=1:11
    tic
    disp(['Creating Mask For Segment ' num2str(j)])
    currentimages=maskedImages{j};
    currentgts=gtsreg{j};
    numimages=length(currentimages);
    for k=1:numimages
        disp(['Processing Image ' num2str(k)])
        Iobrcbr2=improcess(currentimages{k},7,7,0);
        initbin=imcomplement(imbinarize(Iobrcbr2));
        binarea=bwareaopen(initbin,1000,4);
        CC_diameters=bwconncomp(binarea);
        diameters=regionprops(CC_diameters,'EquivDiameter');smallestdiameter=min([diameters.EquivDiameter]);
        largeststrelsize=floor(0.35*smallestdiameter);
        binareafill=imfill(binarea,'holes');
        calccontour_otsu_thissegment{k}=imopen(binareafill,strel('disk',largeststrelsize,0));
        int=and(calccontour_otsu_thissegment{k},currentgts{k});
        uni=or(calccontour_otsu_thissegment{k},currentgts{k});
        intoveruni(j,k)=sum(int(:))/sum(uni(:));
    end
    meanIOU_otsu(j)=mean(intoveruni(j,:))
    stdIOU_otsu(j)=std(intoveruni(j,:))
    detectiontime_otsu(j)=toc
    calccontour_otsu{j}=calccontour_otsu_thissegment;
end
%% K means
%INPUTS:
    %gtsreg: Cell array of registered ground truths (from regist4x75)
    %maskedImages : cell array of images with connective tissue removed (from NN_ACoutsidein)
    %numberofimages: Number of images per segment (should not be changed if
    %using default data set)
%OUTPUTS:
    %calccontour_kmeans: Cell array of segmented images
    %meanIOU_kmeans: Array of mean IOU per segment.
    %stdIOU_kmeans: Array of standard deviation of IOU per segment.
    %detectiontime_kmeans: Total processing time
currentdir=pwd;
numberofimages=[9,13,13,13,12,12,12,12,12,12,12];
for j=1:11
    tic
    disp(['Creating Mask For Segment ' num2str(j)])
    currentsegment=maskedImages{j};
    currentgts=gtsreg{j};
    for k=1:numberofimages(j)
        disp(['Processing Image ' num2str(k)])
        rgb = currentsegment{k};
        % Image preprocessing
        Iobrcbr2=improcess(rgb,7,7,0);
        % Kmeans segmentation
        J=imadjust(Iobrcbr2); % Enhance contrast
        Id = im2double(J);
        c = kmeans(Id(:), 3,'MaxIter',5000,'Replicates',10);
        p = reshape(c, size(Id));
        iback=mat2gray(p);
        fasc=cell(1,5,1);
        fasc(1:5)={zeros(size(iback))};
        truefasc=zeros(size(iback));
        fascprops(1:5)={[]};
        objects=zeros(1,5);
        for i=1:5
            fasc{i}(iback==(-0.25)+i*0.25)=1;
            objects(i)=fascconncomps{i}.NumObjects;
        end
        objects=sort(objects);
        maxobj=objects(1,end-1);
        for q=1:5
            if fascconncomps{q}.NumObjects==maxobj
                truefasc(iback==(-0.25)+q*0.25)=1;
            end
        end
        %Create structuring element
        se2 = strel('disk', 8);
        %Image opening
        Io2 = imopen(truefasc, se2);
        % Image eroding
        Ie2 = imerode(truefasc, se2);
        % Reconstruct image using erosion to enhance contrast
        Iobr2 = imreconstruct(Ie2, truefasc);
        Iobrd2 = imdilate(Iobr2, se2);
        Iobrcbrfin = imreconstruct(imcomplement(Iobrd2), imcomplement(Iobr2));
        Iobrcbrfin = imcomplement(Iobrcbrfin);
        % Final filtering
        CC=bwconncomp(Iobrcbrfin);
        props=regionprops(CC,'Centroid','Image','Eccentricity','PixelIdxList');
        base=false(size(Iobrcbrfin));
        for i=1:length(props)
            if props(i).Eccentricity<0.868
                base(props(i).PixelIdxList)=1;
            end
        end
        calccontour_kmeans_currentsegment{k}=base;
        int=and(calccontour_kmeans_currentsegment{k},currentgts{k});
        uni=or(calccontour_kmeans_currentsegment{k},currentgts{k});
        intoveruni(j,k)=sum(int(:))/sum(uni(:));
    end
    meanIOU_kmeans(j)=mean(intoveruni(j,1:numberofimages(j)))
    stdIOU_kmeans(j)=std(intoveruni(j,:))
    detectiontime_kmeans(j)=toc
    calccontour_kmeans{j}=calccontour_kmeans_currentsegment;
end