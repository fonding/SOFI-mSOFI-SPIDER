% *************************************************************************
% Automated Interface to apply SParse deconvolution of hIgh-Density 
% supER-resolution images 
% *************************************************************************
%
% Reference to the publication:
%   Hugelier, S., de Rooi, J.J., Bernex, R., Duw?, S., Devos, O., Sliwa, M.,
%   Dedecker, P., Eilers, P.H.C. & Ruckebusch, C. (2016). Sparse 
%   deconvolution of high-density super-resolution images, Sci. Rep. 6, 
%   21413; doi: 10.1038/srep21413 (2016)
%
%
% Authors:  
%       - HUGELIER, S.      (1)
%       - DEVOS, O.         (1)
%       - RUCKEBUSCH, C.    (1)
%           
% (1):  LAboratoire de Spectrochimie Infrarouge et Raman (LASIR)
%       Universit? de Lille 1, UMR CNRS 8516
%       Bât C5, Cit? Scientifique
%       59655 Villeneuve d'Ascq - France
%
% *************************************************************************

function [SparseImage, Positions, subvector] = SpiderInterface(image, sigma, zoom, overlap, PatchSize, BaselineRemoval, kappa, subvector)


% *************************************************************************
% *************************************************************************
% 
% Automated process of performing Spider to localise emitters in blinking 
% data. Divides the image into different sub images (smaller size) to
% facilitate the calculations. The different images have a small overlap to
% compensate for the border effect. 
% When working with experimental data, certain 'patches' should not be 
% analyzed due to the absence of emitters. These patches can be selected to
% skip at the start of the calculation.
% Certain experimental data contain a background, which can be removed up
% selection, by using a method of assymetric penalised B-splines (de Rooi 
% et al.; 2013)
%
% Input: 
% 
%   Image:              The original image, in the form of a matrix cube, 
%                       to be analysed.
%   Sigma:              Sigma of spread function (expressed in pixels of 
%                       the original image).
%   Zoom:               Zoom factor for superresolution.
%   Overlap:            pixels of overlap between the different patches.
%   PatchSize:          The size of the (square) patches to split the
%                       original data into.
%                       Remark 1:   The original size of the matrix should 
%                                   be divisible by this size (ex.:
%                                   original size is 200x160, take patches 
%                                   of 20 and not 18).
%                       Remark 2:   The value for kappa changes according
%                                   to the size of the patch (as it varies 
%                                   along the total intensity of the
%                                   patch).
%                       Remark 3:   The memory used to do the calculations
%                                   depends heavily on the size of the 
%                                   patches (overlap included) and the zoom 
%                                   factor, so choose accordingly.
%   BaselineRemoval:    Selects whether a background needs to be removed or
%                       not (1: baseline removal; 0: no baseline removal).
%   Kappa:              Penalty factor on the least squares calculation
%                       that decides whether a localised map should be 
%                       sparse or not.
%                       Remark 1:   Kappa depends on the total intensity of
%                                   the image frame.
%                       Remark 2:   High kappa (e.g.: 2000) means a lot of
%                                   penalty and the result will be sparse. 
%                                   Low kappa (e.g.: 100) means less 
%                                   penalty, and thus the result will be 
%                                   less sparse.
%   Subvector:          A vector containing a 1 for the patches that 
%                       contain signal and thus need to be analyzed and a 0 
%                       for the other ones. If not provided, the interface 
%                       will show a screen to select the different 
%                       subimages to be analyzed. 
%                       Remark:     This vector should be provided as a 
%                                   columnvector
%  
% Output:
% 
%   SparseImage:        The calculated sparse image with the localised 
%                       emitters.
%   Positions:          (x,y) positions of the emitters and the intensity
%                       of these localised emitters.
%   Subvector:          A column vector containing 1's and 0's, depending
%                       on the fact whether a patch should be analysed or
%                       not.
%                       This vector is saved so it can be used in future
%                       calculations of the same data.
% 
% *************************************************************************
% *************************************************************************

% Parameters for the background removal. Change when needed!
PParam = 0.001;
SplineParam = 10;
LambdaParam = 0.01;

% Detecting if the graphics card is ATI brand (gives some visual errors
% when the patches are selected). This code compensates for this visual
% error.
gl = opengl('data');
if strfind(gl.Vendor,'ATI')
    
    opengl('software')
    
end

% Turn off all warnings (for readability)
warning('off','all')

% Getting the size of the original image.
m = size(image,1);
l = size(image,2);
n = size(image,3);

% Error messages when providing wrong input.
NumberOfImages = floor(m/PatchSize);
NumberOfImages1 = floor(l/PatchSize);

if PatchSize * NumberOfImages ~= m || PatchSize * NumberOfImages1 ~= l
    
    disp(' ');
    error('The given number of pixels for the patches cannot be used in combination with the size of the original data');
    
end

if nargin == 8
    
    if size(subvector,1) ~= NumberOfImages*NumberOfImages1
        disp(' ');
        error('The given subvector is not correct for this analysis');
        
    end
    
end

% Pre-allocate matrix to improve speed of the calculations.
SparseImage = zeros(m*zoom,l*zoom);
SubVectorRectangle = zeros(NumberOfImages*NumberOfImages1,1);

% Removes the baseline of the global image if necessary, by using a 2D 
% P-splines method with the parameters provided above).
if BaselineRemoval == 1
    
    for b=1:n
        
        if b == 1
            
            [Z,WeightBackground] = baselineSpider(image(:,:,b), PParam, SplineParam, LambdaParam);
            image(:,:,b) = image(:,:,b) - Z;
            
        else
            
            [Z,WeightBackground] = baselineSpider(image(:,:,b), PParam, SplineParam, LambdaParam, WeightBackground);
            image(:,:,b) = image(:,:,b) - Z;
            
        end
        
    end
    
end

% Removes an offset to the data.
for b=1:n
        
    image(:,:,b) = image(:,:,b) - min(min(image(:,:,b)));
    
end

% Determines the size of the patch to be analyzed and sets the counter to
% a value of 0 to follow the progress of the calculations.
o = PatchSize+2*overlap;
Counter = 0;

% Shows the image to decide which patches to calculate. Divides the image
% into patches of the given size and numbers them. This part is only
% launched when no subvector is given.
if nargin < 8
    
    % Create subvector to improve speed of calculations
    subvector = ones(NumberOfImages*NumberOfImages1,1);
    
    % Draw the mean image and show the numbered patches (divided by a grid)
    figure;imagesc(mean(image,3));hold on;

    for a = 0:PatchSize:m-PatchSize
        
        x = [0 l];
        y = [a a];
        plot(x,y,'Color','k');
        
    end
    
    for a = 0:PatchSize:l-PatchSize
        
        x = [a a];
        y = [0 m];
        plot(x,y,'Color','k');
        
    end
    
    xg = 0:PatchSize:m-PatchSize;
    yg = 0:PatchSize:l-PatchSize;
    [xlbl, ylbl] = meshgrid(xg+PatchSize/2, yg+PatchSize/2);
    lbl = strtrim(cellstr(num2str((1:numel(xlbl))')));
    text(ylbl(:), xlbl(:), lbl(:),'color','k','HorizontalAlignment','center','VerticalAlignment','middle');

    drawnow

    % Select the frames to be analyzed and create the vector that can be
    % used in future calculations.
    % Clicking a patch will deselect it (i.e. no calculation needed),
    % clicking it again will reselect it (i.e. calculation needed).
    % Remark:   Clicking outside the image will cause an error and the
    %           interface will stop (and it will have to be relaunched).
    disp('Select the patches that do not have to be analyzed and press any key when done.');
    
    while waitforbuttonpress ~= 1
        
        Position = get(gca,'CurrentPoint');
        YPos = floor(Position(1,2) / PatchSize);
        XPos = ceil(Position(1,1) / PatchSize);
        Sub = YPos*NumberOfImages1 + XPos;
        
        if subvector(Sub) == 1
            
            subvector(Sub) = 0;
            SubVectorRectangle(Sub) = patch([(XPos-1)*PatchSize XPos*PatchSize XPos*PatchSize (XPos-1)*PatchSize],[YPos*PatchSize YPos*PatchSize (YPos+1)*PatchSize (YPos+1)*PatchSize],'k');
            set(SubVectorRectangle(Sub),'FaceAlpha',0.25);
            
        else
            
            subvector(Sub) = 1;
            delete(SubVectorRectangle(Sub));
            
        end
        
    end
    
    % Show how many patches are skipped for the analysis.
    NumberSkipped = size(find(~subvector),1);
    disp(['Skipping ',num2str(NumberSkipped),' / ',num2str(NumberOfImages*NumberOfImages1),' subimages.']);
    disp(' ');
    
end

% Display starting time of calculations
% disp(['Calculations started: ',datestr(now),'.']);
% disp(' ');
% tic;

% Start of the calculations with the different patches. Calculates the
% localised positions of the emitters for a given patch in every frame.
% The procedure analyzes the patches in horizontal direction and then
% follows with the second 'set' of horizontal patches, and so on.
% 
% The method is designed so that the patches have a user-specified overlap
% to avoid missing information. 
% This means that the patches are a few pixels bigger in every direction, 
% but the final reconstructed image has its original size (multiplied by
% the zoom factor)

% Preparation of the calculations. This will calculate the necessary 
% matrices once for the patches and use them all the time (these matrices 
% don't change from patch to patch).
B = make_basis(o, zoom, sigma);
n2 = size(B, 2) ^2;
S = fast_inprod(B);
B0 = make_basis(o, 1, sigma);
n0 = size(B0, 2);
n02 = n0 * n0;
S0 = fast_inprod(B0);
R0 = kappa * speye(n02);
sel = 1:n02;
h = repmat(1:n0, zoom, 1);
h = h(:);
R = kappa * speye(n2);

for i=1:NumberOfImages
    
    % Calculates the vertical positions of the patches to be analyzed and
    % the vertical positions of the patches on the final image.
    Y1 = ((i-1)*PatchSize)+1;
    Y2 = (i*PatchSize)+2*overlap;
    YY1 = ((i-1)*PatchSize*zoom)+1;
    YY2 = i*PatchSize*zoom;
    
    for j=1:NumberOfImages1
        
        % Adds 1 to the counter to follow the progress (Counter works on a 
        % given patch, independent from the number of frames present in 
        % the matrix cube). 
        Counter = Counter + 1;
        
        % Calculated the horizontal positions of the patches to be
        % analyzed and the horizontal positions of the patches on the
        % final image.
        X1 = ((j-1)*PatchSize)+1;
        X2 = (j*PatchSize)+2*overlap;
        XX1 = ((j-1)*PatchSize*zoom)+1;
        XX2 = (j*PatchSize*zoom);
        
        % Pre-allocates the temporary sparse vector to improve speed of the
        % calculation.
        TempSparse = zeros((o*zoom)^2,n);
        
        reverseString = '';
        for t=1:n
      
            % Adds rows and columns of zeros (2x overlap) to the matrix.
            % Selects the patches of a given frame k.
            Frame = zeros(m+2*overlap,l+2*overlap);
            Frame(overlap+1:m+overlap,overlap+1:l+overlap) = image(:,:,t);
            
            SubImage = Frame(Y1:Y2,X1:X2);
            
            % Selection of the patches.
            if subvector(Counter) == 0
                break
            end
            
            % Localization calculations for the given patch of the given
            % frame k.
            [TempSparse(:,t)] = SpiderCalculations(SubImage, kappa, B, n2, S, B0, n0, n02, S0, sel, h, R0, R);
            
            % Printing the number of frame that was finished calculating
% % % % % % % % % % % % % % % % % % %             if t ~= 1
% % % % % % % % % % % % % % % % % % %                 msg = sprintf(['Calculating frame: ',num2str(t)]);
% % % % % % % % % % % % % % % % % % %             else
% % % % % % % % % % % % % % % % % % %                 msg = sprintf(['\nCalculating frame: ',num2str(t)]);
% % % % % % % % % % % % % % % % % % %             end
% % % % % % % % % % % % % % % % % % %             fprintf([reverseString, msg]);
% % % % % % % % % % % % % % % % % % %             reverseString = repmat(sprintf('\b'), 1, length(msg));
        end
               
        % Reshape the temporary sparse vector into a matrix.
        % Select the useful rows and columns of this matrix and constructs
        % the final sparse image.
        SparseMatrix = reshape(TempSparse, o*zoom, o*zoom, n);
        
        for t = 1:n

            SparseImage(YY1:YY2,XX1:XX2,t) = SparseMatrix((overlap*zoom)+1:(o-overlap)*zoom,(overlap*zoom)+1:(o-overlap)*zoom, t);
        
        end
        
        % Shows the progress.
% % % % %         if subvector(Counter) == 0
% % % % %             
% % % % %             fprintf([reverseString,'Skipping subimage: ',num2str(Counter),'\n']); 
% % % % %             
% % % % %         else
% % % % %             
% % % % %             fprintf([reverseString,'Finished with ',num2str(Counter),' / ',num2str(NumberOfImages*NumberOfImages1),' subimages.','\n']); 
% % % % %             
% % % % %         end
        
    end
    
end

% Localizes all pixels with an emitter.
% Add the values to the positions to correct between center of the pixel is
% integer for imagesc and corner of the grid is integer for plot
for i = 1:n
    
    [row,col,v] = find(SparseImage(:,:,i));
    
    if n ~= 1    
        Positions{i} = [row+zoom/2-0.5,col+zoom/2-0.5,v];
    else
        Positions = [row+zoom/2-0.5,col+zoom/2-0.5,v];
    end
    
end

% Display ending time of calculations .
% % % % % disp(['Calculations stopped: ',datestr(now),'.']);
% toc;

% Clear the unnecessary matrices from the memory.
clear B; clear nn; clear n2; clear S; clear B0; clear n0; clear n02; clear S0; clear R0; clear sel; clear h; clear R;
warning('on','all')

end