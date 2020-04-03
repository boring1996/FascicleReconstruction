%AUTHOR:
%DANIEL TOVBIS (2019)
%DESCRIPTION:
    %This script will train the VGG-16 RCNN to detect fascicles. You only train
%with two sets of training data at a time (for testing on the third- you
%obviously don't want to test a network on the data it was trained on!)
%This was done on the NeuroScience Gateway but can also work on any machine
%assuming the file structures are set up correctly.
%INPUTS:
    %Training Data: The list of images and bounding boxes, stored as a .mat file
    %Pretrained Neural Network: A pretrained neural network. Currently using vgg16.
    %Training Options: Type help trainingOptions for full description.
    %Positive/Negative Overlap Range: When training determines what percent overlap to consider a positive
    %vs a negative result. see help trainFasterRCNNObjectDetector
    %for more details.
%OUTPUTS:
    %A trained neural network (.mat), named by the user.
%% Input: Load training data and pretrained neural net
currentdir=pwd;
load([currentdir,'/Nerve Images/Proximal1HE/trainingdataProximal1HE.mat']);
load([currentdir,'/Nerve Images/Proximal2HE/trainingdataProximal2HE.mat']);
load([currentdir,'/Nerve Images/Distal3HE/trainingdataDistal3HE.mat']);
load([currentdir, '/vgg16.mat']);
%% Change data source to fit user's file structure
%You can do this as part of this code, or you can do it seperately on your
%own and save the files so you don't have to do this every time.
% N=size(gTruthProximal1HE.DataSource,1);
% for q=1:N
%     gTruthProximal1HE.DataSource{q}=strrep(gTruthProximal1HE.DataSource{q},'M:\Peripheral Nerve Studies\Daniel\Nerve Processing',currentdir);
% end
% clearvars N q
% N=size(gTruthProximal2HE.DataSource,1);
% for q=1:N
%     gTruthProximal2HE.DataSource{q}=strrep(gTruthProximal2HE.DataSource{q},'M:\Peripheral Nerve Studies\Daniel\Nerve Processing',currentdir);
% end
% clearvars N q
% N=size(gTruthDistal3HE.DataSource,1);
% for q=1:N
%     gTruthDistal3HE.DataSource{q}=strrep(gTruthDistal3HE.DataSource{q},'M:\Peripheral Nerve Studies\Daniel\Nerve Processing',currentdir);
% end
% clearvars N q
%% set up the training data
%Comment one out so training is only done on two
trainingData1 = objectDetectorTrainingData(gTruthProximal1HE);
trainingData2= objectDetectorTrainingData(gTruthProximal2HE);
%trainingData3= objectDetectorTrainingData(gTruthDistal3HE);
trainingDatafull=vertcat(trainingData1,trainingData2); %Make sure you change this to fit the training data you're using)
%% Reformat data to fit NSG syntax
%Reformat the training data to fit with NSG's syntax. NSG uses a linux
%machine, but this code was written on Windows, so backslashes need to be changed to forward slashes.
%The directory is also different on their server so we need to change any
%references to the user directory (M drive in the training data)
%This is optional depending on how you set the data up.
N=height(trainingDatafull);
for q=1:N
    trainingDatafull.imageFilename{q}=strrep(trainingDatafull.imageFilename{q},'M:\Peripheral Nerve Studies\Daniel\Nerve Processing',currentdir);
    trainingDatafull.imageFilename{q}=strrep(trainingDatafull.imageFilename{q},'\','/');
end
%% Input: Training Options
options = trainingOptions('sgdm', ...
    'MiniBatchSize', 32, ...
    'InitialLearnRate', 1e-3, ...
    'LearnRateSchedule', 'piecewise', ...
    'LearnRateDropFactor', 0.1, ...
    'LearnRateDropPeriod', 5, ...
    'MaxEpochs', 10, ...
    'ExecutionEnvironment','auto',...
    'Verbose', true,...
    'VerboseFrequency',100);
%% Train an R-CNN object detector. This will take several hours
%Make sure to rename the variable to whatever you want it to be.
fasterrcnnvgg16noaugment_p1p2 = trainFasterRCNNObjectDetector(trainingDatafull, net, options, ...
    'NegativeOverlapRange', [0 0.5], 'PositiveOverlapRange',[0.5 1])
save('fasterrcnnvgg16noaugment_p1p2','fasterrcnnvgg16noaugment_p1p2')