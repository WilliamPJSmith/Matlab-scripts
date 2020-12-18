%%% ScriptToGenerateSIFigure_Ravm
% Created by WPJS, Summer 2018

% To Run:
% i) Open script in a new Matlab window
% ii) Ensure dirRoot (path to raw image data) is correctly specified
% iii) Under the "Editor" tab, hit "Run"

%% 1) Loading images

% Location of raw image files
dirRoot = '/Volumes/EHD_MAC_2/Experimental_data/SI_processing/';

% Load raw images
[A1,A2,A3,A4,A5,...
    B1,B2,B3,B4,B5,...
    C1,C2,C3,C4,C5,...
    D1,D2,D3,D4,D5,...
    E1,E2,E3,E4,E5] = LoadAllImages(dirRoot);

%% 2) Segmenting images (Otsu thresholding)
[M_A1] = SegmentImage(A1);
[M_A2] = SegmentImage(A2);
[M_A3] = SegmentImage(A3);
[M_A4] = SegmentImage(A4);
[M_A5] = SegmentImage(A5);

[M_B1] = SegmentImage(B1);
[M_B2] = SegmentImage(B2);
[M_B3] = SegmentImage(B3);
[M_B4] = SegmentImage(B4);
[M_B5] = SegmentImage(B5);

[M_C1] = SegmentImage(C1);
[M_C2] = SegmentImage(C2);
[M_C3] = SegmentImage(C3);
[M_C4] = SegmentImage(C4);
[M_C5] = SegmentImage(C5);

[M_D1] = SegmentImage(D1);
[M_D2] = SegmentImage(D2);
[M_D3] = SegmentImage(D3);
[M_D4] = SegmentImage(D4);
[M_D5] = SegmentImage(D5);

[M_E1] = SegmentImage(E1);
[M_E2] = SegmentImage(E2);
[M_E3] = SegmentImage(E3);
[M_E4] = SegmentImage(E4);
[M_E5] = SegmentImage(E5);

%% 3) Measuring average segregation index across segmented images
% /!\ This step is time-consuming. Uncomment only to re-compute SI data
% from scratch (otherwise, load precomputed SI data)
load('./SI_Data')

% set span in pixels
% span = 5;
% [~, ~, FullData_5]   = GetSITraceForSpan(M_A1,M_A2,M_A3,M_A4,M_A5,M_B1,M_B2,M_B3,M_B4,M_B5,M_C1,M_C2,M_C3,M_C4,M_C5,M_D1,M_D2,M_D3,M_D4,M_D5,M_E1,M_E2,M_E3,M_E4,M_E5,span);
% 
% span = 10;
% [~, ~, FullData_10]   = GetSITraceForSpan(M_A1,M_A2,M_A3,M_A4,M_A5,M_B1,M_B2,M_B3,M_B4,M_B5,M_C1,M_C2,M_C3,M_C4,M_C5,M_D1,M_D2,M_D3,M_D4,M_D5,M_E1,M_E2,M_E3,M_E4,M_E5,span);
% 
% span = 20;
% [~, ~, FullData_20]   = GetSITraceForSpan(M_A1,M_A2,M_A3,M_A4,M_A5,M_B1,M_B2,M_B3,M_B4,M_B5,M_C1,M_C2,M_C3,M_C4,M_C5,M_D1,M_D2,M_D3,M_D4,M_D5,M_E1,M_E2,M_E3,M_E4,M_E5,span);
% 
% span = 50;
% [~, ~, FullData_50]   = GetSITraceForSpan(M_A1,M_A2,M_A3,M_A4,M_A5,M_B1,M_B2,M_B3,M_B4,M_B5,M_C1,M_C2,M_C3,M_C4,M_C5,M_D1,M_D2,M_D3,M_D4,M_D5,M_E1,M_E2,M_E3,M_E4,M_E5,span);
% 
% span = 100;
% [~, ~, FullData_100] = GetSITraceForSpan(M_A1,M_A2,M_A3,M_A4,M_A5,M_B1,M_B2,M_B3,M_B4,M_B5,M_C1,M_C2,M_C3,M_C4,M_C5,M_D1,M_D2,M_D3,M_D4,M_D5,M_E1,M_E2,M_E3,M_E4,M_E5,span);

%% 4) Plot SI values for each grid replicate, as a function of grid motif and SI span

% plotting options
close all, figure('color','w')
LW = {'LineWidth',2.5};
LW2 = {'LineWidth',3.5};
FS = {'FontSize',16};
FS2 = {'FontSize',20};
LO = {'Location','SouthEast'};

% set plot colorscheme
colors = [[31,119,180] ./ 255;...
[255,127,14] ./ 255;...
[44,160,44] ./ 255;...
[214,39,40] ./ 255;...
[148,103,189] ./ 255];

% scatterplotting measured SI values
xvec = repmat([1,2,3,4,5]',[1,5]);
sz = 50;
h1 = scatter(xvec(:), FullData_5(:), sz, 'o', 'MarkerEdgeColor', colors(1,:)); hold on
h2 = scatter(xvec(:), FullData_10(:), sz, 'o', 'MarkerEdgeColor', colors(2,:));
h3 = scatter(xvec(:), FullData_20(:), sz, 'o', 'MarkerEdgeColor', colors(3,:));
h4 = scatter(xvec(:), FullData_50(:), sz, 'o', 'MarkerEdgeColor', colors(4,:));
h5 = scatter(xvec(:), FullData_100(:), sz, 'o', 'MarkerEdgeColor', colors(5,:));

% plot formatting
hax = gca; set(hax, LW{:}, FS{:})
xlabel('Motif',FS{:})
ylabel('S.I.',FS{:})
ylim([0.45,1.0])
yticks(0.5:0.1:1.0)
pbaspect([5,3,1])
xticks([1 2 3 4 5])
xlim([0.5,5.5])
SIRefVec = fliplr([0.967, 0.936, 0.900, 0.820, 0.545]);
% plot(SIRefVec,'k--',LW{:})

%% 5) Repeat for reference grids

% generate reference grids
[GridA,GridB,GridC,GridD,GridE] = GetReferenceGrids();

% compute SI values (same method as above, but 1 replicate only)
% [MeansVec_5_REF] = GetSITraceForSpan_REF(GridA,GridB,GridC,GridD,GridE,5);
% [MeansVec_10_REF] = GetSITraceForSpan_REF(GridA,GridB,GridC,GridD,GridE,10);
% [MeansVec_20_REF] = GetSITraceForSpan_REF(GridA,GridB,GridC,GridD,GridE,20);
% [MeansVec_50_REF] = GetSITraceForSpan_REF(GridA,GridB,GridC,GridD,GridE,50);
% [MeansVec_100_REF] = GetSITraceForSpan_REF(GridA,GridB,GridC,GridD,GridE,100);

load('./SI_Data_Ref')

% plotting reference SIs as dotted lines
plot(fliplr([1,2,3,4,5]),MeansVec_5_REF,'o--',LW{:},'color',colors(1,:)), hold on
plot(fliplr([1,2,3,4,5]),MeansVec_10_REF,'o--',LW{:},'color',colors(2,:))
plot(fliplr([1,2,3,4,5]),MeansVec_20_REF,'o--',LW{:},'color',colors(3,:))
plot(fliplr([1,2,3,4,5]),MeansVec_50_REF,'o--',LW{:},'color',colors(4,:))
plot(fliplr([1,2,3,4,5]),MeansVec_100_REF,'o--',LW{:},'color',colors(5,:))

% format figure legend
hl = legend('7.6','15.2','30.3','75.8','151.5');
set(hl,FS{:},LW{:},'Location','SouthEast'), legend boxoff
title(hl,'Span [\mum]')
box on

% save figure
saveas(gcf,'SI_Figure.fig')






%%% SUPPORTING FUNCTIONS %%%

function [SI_matrix] = GetDiscreteSIMatrix(M, span)
%% GETDISCRETESIMATRIX computes a discrete analogue to Nadel's SI index, for a given segmented confocal image.
[im_height, im_width] = size(M);
SI_matrix = zeros(im_height-(2*span),im_width-(2*span)); % preallocate matrix of local SI indices
N = ((2*span)+1)^2 - 1;                                  % number of neighbours in a square of edge length 2*span+1

% go through each pixel in the image, except for outer border of width <Span>
for i=(1+span):(im_height-span)
    for j=(1+span):(im_width-span)
        
        % the focal pixel
        focal_type = M(i,j);
        
        % check if focal pixel is empty (local SI = NaN)
        if focal_type==0
            SI_matrix(i-span,j-span) = NaN;
            continue
        end
        
        % if focal pixel not empty, measure local SI around this point,
        % ignoring empty / ambiguous pixels
        M_sub = M(i-span:i+span,j-span:j+span);
        num_off = sum(M_sub(:)==0);
        num_focal = sum(M_sub(:)==focal_type) - 1;
        SI_matrix(i-span,j-span) = num_focal / (N-num_off);
    end
end

end

function [MeansVec, StDevVec, FullData] = GetSITraceForSpan(M_A1,M_A2,M_A3,M_A4,M_A5,M_B1,M_B2,M_B3,M_B4,M_B5,M_C1,M_C2,M_C3,M_C4,M_C5,M_D1,M_D2,M_D3,M_D4,M_D5,M_E1,M_E2,M_E3,M_E4,M_E5,span)
%% GETSITRACEFORSPAN creates a local segregation index (SI) map for each image, reporting the average SI value over each individual image
% SI is computed for a given span value (in pixels); e.g. 10 pixels ->
% report the degree of segregation for a square area measuring 2*10+1 = 21
% pixels along a side.

% local segregation index (SI) map for each image
[SI_A1] = GetDiscreteSIMatrix(M_A1, span);
[SI_A2] = GetDiscreteSIMatrix(M_A2, span);
[SI_A3] = GetDiscreteSIMatrix(M_A3, span);
[SI_A4] = GetDiscreteSIMatrix(M_A4, span);
[SI_A5] = GetDiscreteSIMatrix(M_A5, span);

[SI_B1] = GetDiscreteSIMatrix(M_B1, span);
[SI_B2] = GetDiscreteSIMatrix(M_B2, span);
[SI_B3] = GetDiscreteSIMatrix(M_B3, span);
[SI_B4] = GetDiscreteSIMatrix(M_B4, span);
[SI_B5] = GetDiscreteSIMatrix(M_B5, span);

[SI_C1] = GetDiscreteSIMatrix(M_C1, span);
[SI_C2] = GetDiscreteSIMatrix(M_C2, span);
[SI_C3] = GetDiscreteSIMatrix(M_C3, span);
[SI_C4] = GetDiscreteSIMatrix(M_C4, span);
[SI_C5] = GetDiscreteSIMatrix(M_C5, span);

[SI_D1] = GetDiscreteSIMatrix(M_D1, span);
[SI_D2] = GetDiscreteSIMatrix(M_D2, span);
[SI_D3] = GetDiscreteSIMatrix(M_D3, span);
[SI_D4] = GetDiscreteSIMatrix(M_D4, span);
[SI_D5] = GetDiscreteSIMatrix(M_D5, span);

[SI_E1] = GetDiscreteSIMatrix(M_E1, span);
[SI_E2] = GetDiscreteSIMatrix(M_E2, span);
[SI_E3] = GetDiscreteSIMatrix(M_E3, span);
[SI_E4] = GetDiscreteSIMatrix(M_E4, span);
[SI_E5] = GetDiscreteSIMatrix(M_E5, span);

MeanSIVec_A = [nanmean(nanmean(SI_A1))...
               nanmean(nanmean(SI_A2))...
               nanmean(nanmean(SI_A3))...
               nanmean(nanmean(SI_A4))...
               nanmean(nanmean(SI_A5))];...
               
MeanSIVec_B = [nanmean(nanmean(SI_B1))...
               nanmean(nanmean(SI_B2))...
               nanmean(nanmean(SI_B3))...
               nanmean(nanmean(SI_B4))...
               nanmean(nanmean(SI_B5))];...
               
MeanSIVec_C = [nanmean(nanmean(SI_C1))...
               nanmean(nanmean(SI_C2))...
               nanmean(nanmean(SI_C3))...
               nanmean(nanmean(SI_C4))...
               nanmean(nanmean(SI_C5))];...
               
MeanSIVec_D = [nanmean(nanmean(SI_D1))...
               nanmean(nanmean(SI_D2))...
               nanmean(nanmean(SI_D3))...
               nanmean(nanmean(SI_D4))...
               nanmean(nanmean(SI_D5))];...

MeanSIVec_E = [nanmean(nanmean(SI_E1))...
               nanmean(nanmean(SI_E2))...
               nanmean(nanmean(SI_E3))...
               nanmean(nanmean(SI_E4))...
               nanmean(nanmean(SI_E5))];...

MeansVec = mean([MeanSIVec_A; MeanSIVec_B; MeanSIVec_C; MeanSIVec_D; MeanSIVec_E],2);
StDevVec =  std([MeanSIVec_A; MeanSIVec_B; MeanSIVec_C; MeanSIVec_D; MeanSIVec_E],[],2);
FullData = [MeanSIVec_A; MeanSIVec_B; MeanSIVec_C; MeanSIVec_D; MeanSIVec_E];
end

function [MeansVec] = GetSITraceForSpan_REF(GridA,GridB,GridC,GridD,GridE,span)

[SI_A_REF] = GetDiscreteSIMatrix(GridA, span);
[SI_B_REF] = GetDiscreteSIMatrix(GridB, span);
[SI_C_REF] = GetDiscreteSIMatrix(GridC, span);
[SI_D_REF] = GetDiscreteSIMatrix(GridD, span);
[SI_E_REF] = GetDiscreteSIMatrix(GridE, span);

MeansVec = [nanmean(nanmean(SI_A_REF)), nanmean(nanmean(SI_B_REF)), nanmean(nanmean(SI_C_REF)), nanmean(nanmean(SI_D_REF)), nanmean(nanmean(SI_E_REF))];
end

function [M] = SegmentImage(A)
%% SEGMENTIMAGE uses Otsu's thresholding method to segment raw images, assigning pixels according to thresholded red and green channels

% Compute thresholds for red and green channels separately
[level_R,~] = graythresh(A(:,:,1));
[level_G,~] = graythresh(A(:,:,2));

% Binary images showing pixels > respective threshold
BW_R = im2bw(A(:,:,1),level_R);
BW_G = im2bw(A(:,:,2),level_G);

% Composite image with the following encoding:
% 0: EMPTY (no signal, or ambigious signal)
% 1: RED (only red channel above threshold)
% 2: GREEN (only green channel above threshold))

M = BW_R + 2*BW_G;
M(M>2) = 0;
end

function [A1,A2,A3,A4,A5,...
          B1,B2,B3,B4,B5,...
          C1,C2,C3,C4,C5,...
          D1,D2,D3,D4,D5,...
          E1,E2,E3,E4,E5] = LoadAllImages(dirRoot)
%% LOADALLIMAGES imports raw images of printed communities for processing

% File labels
dirA = '0_545/'; dirB = '0_820/'; dirC = '0_900/'; dirD = '0_936/'; dirE = '0_967/';
fileRootA = '20180430_bn52';
fileRootB = '20180525_bn61';
fileRootC = '20180512_bn58';
fileRootD = '20180513_bn59';
fileRootE = '20180503_bn54';

% Loading images
A1 = imread(strcat(dirRoot, dirA, fileRootA, '_1.tif'));
A2 = imread(strcat(dirRoot, dirA, fileRootA, '_2.tif'));
A3 = imread(strcat(dirRoot, dirA, fileRootA, '_3.tif'));
A4 = imread(strcat(dirRoot, dirA, fileRootA, '_4.tif'));
A5 = imread(strcat(dirRoot, dirA, fileRootA, '_5.tif'));

B1 = imread(strcat(dirRoot, dirB, fileRootB, '_1.tif'));
B2 = imread(strcat(dirRoot, dirB, fileRootB, '_2.tif'));
B3 = imread(strcat(dirRoot, dirB, fileRootB, '_3.tif'));
B4 = imread(strcat(dirRoot, dirB, fileRootB, '_4.tif'));
B5 = imread(strcat(dirRoot, dirB, fileRootB, '_5.tif'));

C1 = imread(strcat(dirRoot, dirC, fileRootC, '_1.tif'));
C2 = imread(strcat(dirRoot, dirC, fileRootC, '_2.tif'));
C3 = imread(strcat(dirRoot, dirC, fileRootC, '_3.tif'));
C4 = imread(strcat(dirRoot, dirC, fileRootC, '_4.tif'));
C5 = imread(strcat(dirRoot, dirC, fileRootC, '_5.tif'));

D1 = imread(strcat(dirRoot, dirD, fileRootD, '_1.tif'));
D2 = imread(strcat(dirRoot, dirD, fileRootD, '_2.tif'));
D3 = imread(strcat(dirRoot, dirD, fileRootD, '_3.tif'));
D4 = imread(strcat(dirRoot, dirD, fileRootD, '_4.tif'));
D5 = imread(strcat(dirRoot, dirD, fileRootD, '_5.tif'));

E1 = imread(strcat(dirRoot, dirE, fileRootE, '_1.tif'));
E2 = imread(strcat(dirRoot, dirE, fileRootE, '_2.tif'));
E3 = imread(strcat(dirRoot, dirE, fileRootE, '_3.tif'));
E4 = imread(strcat(dirRoot, dirE, fileRootE, '_4.tif'));
E5 = imread(strcat(dirRoot, dirE, fileRootE, '_5.tif'));

end

function [GridA,GridB,GridC,GridD,GridE] = GetReferenceGrids()

% constant factors
num_x = 1024;
num_y = 1024;
num_x_sub = 512;
num_y_sub = 512;

% Grid A
GridA = zeros(num_x, num_y);
Sub1 = 2*ones(num_y_sub/2,num_x_sub);
Sub2 = 1*ones(num_y_sub/2,num_x_sub);
Startx = 0.5*(num_x - num_x_sub);
Starty = 0.5*(num_y - num_y_sub);
SubA = [Sub1;Sub2];
GridA(Startx:Startx+num_x_sub-1, Starty:Starty+num_y_sub-1) = SubA;

% Grid B
GridB = zeros(num_x, num_y);
SubB = SubA;
SubB(:,(0.5*num_x_sub)+1:end) = flipud(SubA(:,(0.5*num_x_sub)+1:end));
GridB(Startx:Startx+num_x_sub-1, Starty:Starty+num_y_sub-1) = SubB;

% Grid C
GridC = zeros(num_x, num_y);
Sub1 = 2*ones(num_y_sub/4,num_x_sub);
Sub2 = 1*ones(num_y_sub/4,num_x_sub);
SubC = [Sub1;Sub2;Sub1;Sub2];
GridC(Startx:Startx+num_x_sub-1, Starty:Starty+num_y_sub-1) = SubC;

% Grid D
GridD = zeros(num_x, num_y);
SubD = checkerboard(num_y_sub/4, 2, 2);
GridD(Startx:Startx+num_x_sub-1, Starty:Starty+num_y_sub-1) = flipud(SubD+1);
GridD = round(GridD);

% Grid E
GridE = zeros(num_x, num_y);
SubE = checkerboard(16, 16, 16);
GridE(Startx:Startx+num_x_sub-1, Starty:Starty+num_y_sub-1) = flipud(SubE+1);
GridE = round(GridE);

end
