%AUTHOR: DANIEL TOVBIS (2019)
%DESCRIPTION: This script contains options for registering H&E and IHC
%images using regist4x75. Please only run one subsection, appropriate to
% the desired image type, then save or rename the output files as the
% output names are the same.
%INPUTS: 
%See regist4x75 documentation for description of inputs.
%OUTPUTS:
%See regist4x75 documentation for description of outputs.
%mean(SSIM/MSE)(orig/reg): Mean SSIM or MSE for original or registered
%images, stored as a cell array.
%% IHC
tic
[regim{1},SSIMorig{1},SSIMreg{1},MSEorig{1},MSEreg{1},gtsreg{1}]=regist4x75('Distal1IHC',1,9,1); %10-12 excluded due to damage
[regim{2},SSIMorig{2},SSIMreg{2},MSEorig{2},MSEreg{2},gtsreg{2}]=regist4x75('Distal1IHC',14,11,1);%13 excluded due to damage
[regim{3},SSIMorig{3},SSIMreg{3},MSEorig{3},MSEreg{3},gtsreg{3}]=regist4x75('Distal1IHC',26,11,1);%25 excluded due to damage
[regim{4},SSIMorig{4},SSIMreg{4},MSEorig{4},MSEreg{4},gtsreg{4}]=regist4x75('Distal1IHC',37,12,1);%42 is damaged, may affect outcome
[regim{5},SSIMorig{5},SSIMreg{5},MSEorig{5},MSEreg{5},gtsreg{5}]=regist4x75('Distal1IHC',51,10,1);%49 and 50 excluded due to damage
[regim{6},SSIMorig{6},SSIMreg{6},MSEorig{6},MSEreg{6},gtsreg{6}]=regist4x75('Distal2IHC',1,12,1);
[regim{7},SSIMorig{7},SSIMreg{7},MSEorig{7},MSEreg{7},gtsreg{7}]=regist4x75('Distal2IHC',13,12,1);
[regim{8},SSIMorig{8},SSIMreg{8},MSEorig{8},MSEreg{8},gtsreg{8}]=regist4x75('Distal2IHC',25,12,1);
[regim{9},SSIMorig{9},SSIMreg{9},MSEorig{9},MSEreg{9},gtsreg{9}]=regist4x75('Distal4IHC',1,16,1); %17 and 18 excluded due to damage
[regim{10},SSIMorig{10},SSIMreg{10},MSEorig{10},MSEreg{10},gtsreg{10}]=regist4x75('Distal4IHC',20,17,1); %19 excluded due to damage
toc
%% HE
tic
[regim{1},SSIMorig{1},SSIMreg{1},MSEorig{1},MSEreg{1},gtsreg{1}]=regist4x75('Proximal1HE',2,9,0); %1,11,12,and 13 excluded for damage
[regim{2},SSIMorig{2},SSIMreg{2},MSEorig{2},MSEreg{2},gtsreg{2}]=regist4x75('Proximal1HE',14,13,0);
[regim{3},SSIMorig{3},SSIMreg{3},MSEorig{3},MSEreg{3},gtsreg{3}]=regist4x75('Proximal1HE',27,13,0);
[regim{4},SSIMorig{4},SSIMreg{4},MSEorig{4},MSEreg{4},gtsreg{4}]=regist4x75('Proximal1HE',40,13,0);
[regim{5},SSIMorig{5},SSIMreg{5},MSEorig{5},MSEreg{5},gtsreg{5}]=regist4x75('Proximal1HE',54,12,0); %53 excluded for damage. There are a lot of shit slices in this segment, and we should really excluded more of them but it would make a contiguous model very hard.
[regim{6},SSIMorig{6},SSIMreg{6},MSEorig{6},MSEreg{6},gtsreg{6}]=regist4x75('Proximal2HE',1,12,0);
[regim{7},SSIMorig{7},SSIMreg{7},MSEorig{7},MSEreg{7},gtsreg{7}]=regist4x75('Proximal2HE',13,12,0);
[regim{8},SSIMorig{8},SSIMreg{8},MSEorig{8},MSEreg{8},gtsreg{8}]=regist4x75('Proximal2HE',25,12,0);
[regim{9},SSIMorig{9},SSIMreg{9},MSEorig{9},MSEreg{9},gtsreg{9}]=regist4x75('Distal3HE',1,12,0);
[regim{10},SSIMorig{10},SSIMreg{10},MSEorig{10},MSEreg{10},gtsreg{10}]=regist4x75('Distal3HE',13,12,0);
[regim{11},SSIMorig{11},SSIMreg{11},MSEorig{11},MSEreg{11},gtsreg{11}]=regist4x75('Distal3HE',25,12,0);
toc
%% Stats
for k=1:11 %use 1:10 for IHC
meanSSIMorig(k)=mean(SSIMorig{k});
meanSSIMreg(k)=mean(SSIMreg{k});
meanMSEorig(k)=mean(MSEorig{k});
meanMSEreg(k)=mean(MSEreg{k});
end