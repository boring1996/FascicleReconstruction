%AUTHOR: 
    %DANIEL TOVBIS (2019)
%DESCRIPTION:
    %This function will generate fascicle models using interpmask. It may
    %also fix holes. Fascicles are split apart after reconstruction by
    %using the watershed-erode-watershed method.
%INPUTS:
    %contour: Cell array of segmented fascicle contours (from NN_ACoutsidein or IHC
    %segmentation)
    %dofixing (logical): 1 to do hole fixing, 0 to ignore it.
%OUTPUTS:
    %model: The reconstructed fascicle models, stored as a cell array.
    %reconstructiontime: Array of reconstruction times for each segment.
function [model,reconstructiontime]=fascrecon_modelonly(contour,dofixing)
numsegments=length(contour);
for j=1:numsegments
    disp(['Reconstructing Segment ' num2str(j)])
    if isempty(contour{j})==1
        model{j}=[];
        continue
    end
    tic
    currentsegment=contour{j};
    %% Create Model
    numimages=size(currentsegment,2);
    for v=1:numimages
        BWorig(:,:,v)= masksmooth(imresize(currentsegment{v},0.3), 5, 0.7); %resize the model to save space/time
    end
    %% Fix Holes
    Num_gts=size(BWorig,3);
    BW=BWorig;
    if dofixing==1
        for t=(1:Num_gts-1)
            CC_prerecon{t}=bwconncomp(BW(:,:,(t))); %Acquire info about connected components
            %CC contains connectivity, image size, number of objects, and the Pixel
            %Index List for each object
            props_prerecon{t}=regionprops(CC_prerecon{t},'Centroid','EquivDiameter','Area','Image','PixelIdxList'); %Acquire properties of connected components
            Num_objects_on_slice=length(props_prerecon{t});
            currslice=BW(:,:,t); %Temporary object to store the current slice
            for u=1:Num_objects_on_slice %For all the objects on the current slice
                overlapfound=0;
                x=t; %t represents the variable of the current slice, x is the test slice
                while overlapfound==0 %As long we haven't found an overlap
                    if x+1<=Num_gts %As long as we can find an overlap within the slices we have
                        testslice=BW(:,:,x+1); %Look at the next slice
                        idxs=props_prerecon{t}(u).PixelIdxList; %Find the pixel indices of the current object
                        values_at_next_slice=testslice(idxs); %Check the values of those indices at the next slice
                        sumvalues=sum(values_at_next_slice);
                        if sumvalues==0 %If it is totally empty
                            x=x+1; %Set nextslice to nextslice+1 and repeat
                        else %If there is some overlap
                            overlapfound=1; %We exit the while loop
                            break
                        end
                    else %If we get to the last slice and there is no overlap
                        break %Just do nothing and exit the while loop
                    end
                end
                if x~=t && x+1<=Num_gts  %If there wasn't an overlap found on the next slice, but there was one found eventually
                    currsliceatidx=logical(zeros(size(BW(:,:,t))));
                    currsliceatidx(idxs)=1;
                    testsliceatidx=logical(zeros(size(BW(:,:,t))));
                    testsliceatidx(idxs)=values_at_next_slice;
                    %Values of the current slice at the indicies
                    BWtempinterp=cat(3,currsliceatidx,testsliceatidx); %Concatenate the matrices
                    layersbetween=(x-t)*100; %Number of layers between object on current slice and object on
                    tempinterp=interpmask(1:2,BWtempinterp,linspace(1,2,2*layersbetween)); %Interpolate between them
                    tempinterp100=tempinterp(:,:,100);
                    newobj=find(tempinterp100==1);
                    nextslice=BW(:,:,t+1); %Get the slice that follows the current one
                    nextslice(newobj)=1; %Take the 100th slice from the interpolated array, apply it to the following slice
                    BW(:,:,t+1)=nextslice; %Apply the change to the master array
                else %If there was an overlap found on the next slice, or no overlap found don't do anything
                end
                clear testslice nextslice currsliceatidx BWtempinterp tempinterp values_at_next_slice idxs
            end
        end
    end
    %% Interpolate
    finlayers=50*numimages-50; %The number of layers we want to have in the end.
    BWout = interpmask(1:numimages, BW, linspace(1,numimages,finlayers),'linear'); %Do the interpolation
    se = strel('disk', 1);
    BWout_morphopen=BWout;
    for i=1:finlayers %Small Morphological Opening to break apart errant mergers
        BWout_morphopen(:,:,i)=imopen(BWout(:,:,i),se);
    end
    %% Watershed splitting
    %Adapted from: https://blogs.mathworks.com/steve/2013/11/19/watershed-transform-question-from-tech-support/
    for i=1:finlayers
        D = -bwdist(~BWout_morphopen(:,:,i)); %Distance to border of fascicle
        Ld = watershed(D); %Get watershed of distance transform.
        mask = imextendedmin(D,2); %Get extended minima
        D2 = imimposemin(D,mask); %Impose the minima
        Ld2 = watershed(D2); %Watershed the modified distance transform
        bw3 = BWout_morphopen(:,:,i);
        bw3(Ld2 == 0) = 0;
        BWout_morphopen_watershed(:,:,i)=bw3;
    end
    %% Erode
    % We erode the watershedded fascicles maximally (without removing them) to separate them as much as possible
    % using the diameter of the smallest fascicle in the image
    % This value is unique to each slice and is stored so they can later be
    % dilated again
    for i=1:finlayers
        CC_diameters=bwconncomp(BWout_morphopen_watershed(:,:,i));
        diameters=regionprops(BWout_morphopen_watershed(:,:,i),'MinorAxisLength');
        smallestdiameter=min([diameters.MinorAxisLength]);
        largeststrelsize(i)=floor(0.3*smallestdiameter);
        BWout_morphopen_watershed_eroded(:,:,i)=imerode(BWout_morphopen_watershed(:,:,i),strel('disk',largeststrelsize(i),0));
        clear diameters CC_diameters smallestdiameter
    end
    %% Watershed again
    % Another watershedding on the now eroded fascicles helps separate them as much as possible
    for i=1:finlayers
        D_2 = -bwdist(~BWout_morphopen_watershed_eroded(:,:,i));
        Ld_2 = watershed(D_2);
        mask_2 = imextendedmin(D_2,2);
        D2_2 = imimposemin(D_2,mask_2);
        Ld2_2 = watershed(D2_2);
        bw3_2 = BWout_morphopen_watershed_eroded(:,:,i);
        bw3_2(Ld2_2 == 0) = 0;
        BWout_morphopen_watershed_eroded_watershed(:,:,i)=bw3_2;
    end
    %% Dilate
    %Using the stored strel sizes, every layer is dilated again.
    for i=1:finlayers
        currentmodel(:,:,i)=imdilate(BWout_morphopen_watershed_eroded_watershed(:,:,i),strel('disk',largeststrelsize(i),0));
    end
    %% Save everything
    model{j}=currentmodel;
    reconstructiontime(j)=toc;
    clearvars BWorig  BWout_morphopen_watershed BW BWout_morphopen_watershed_eroded BWout_morphopen_watershed_eroded_watershed BWout_morphopen currentmodel
end
end