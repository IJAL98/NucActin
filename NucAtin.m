%% Header information
%%%%%%%%%%%%
% PRogram to quantify the structue of nuclea f-actin form fluorescence
% microscopy images. Supporting software for the publication "Nuclear
% fascin regulates cancer cell survival."
% Software authors Richard Marsh & Susan cox
% Takes uncrompressed single frame image of several cells as input
% Produces three metrics of filliment structure as output for each cell
%%%%%%%%%%%%%
% Instructions for use.
% Runing programe opens a file dialoge box to select input image
% select cells of interest by drawing a rectanglar ROI around each by left
% clicking top left & bottom right of each cell
% when all cell of interest selected right click to begin analysis.
% metric are output as a matric of numbers, each row represents a cell in
% the order they wher selected.
%%%%%%%%%%%%%


%% Get image file and info
progname = "NucleaActin"
uiwait(msgbox('Select a microscopy image file for analysis.',progname,'modal'));

% open a file  dialogue box
[filename pathname]=uigetfile('*.tif', 'Image sequence');
if (filename==0)
    'Aboting execution!'
    return
end

%  get  info  on selected file
try
    info=imfinfo(strcat(pathname,filename));
    numframes=numel(info);
    clear frame;
    frame=zeros(info(1).Height,info(1).Width,numframes);
    if (numframes>1)
        msgbox('please separate image frames - single image required',progname)
        return
    end
    frame(:,:) = imread(strcat(pathname,filename),1,'Info',info);
catch
    msgbox('An unexpected error occured opening the file. Single frame and channel image file expected.',progname)
end
    
% Section   - get list of cells by drawing box ROI aroud each one
% display image and initialise variables
close all
imgfig = figure;
imagesc(frame);
hold on
Xpoints=[];
Ypoints=[];
button=1;
numpoints = 1;
fig=[];

while (button==1) % loop each time left mouse button pressed
    [Xmouse,Ymouse,button]=ginput(1);
    if (button>1)
        continue %This brakes the loop if right mouse button pressed
    end
    
    %gets nouse coordinates when button in pressed
    Xpoints(numpoints,1)=Xmouse;
    Ypoints(numpoints,1)=Ymouse;
    %plot(Xpoints,Ypoints,'m')
    [Xmouse,Ymouse,button]=ginput(1);
    Xpoints(numpoints,2)=Xmouse;
    Ypoints(numpoints,2)=Ymouse;
    
    %Plot a square ROI using top left and botom right from points list
    plot([Xpoints(numpoints,1);Xpoints(numpoints,2)],[Ypoints(numpoints,1),Ypoints(numpoints,1)],'m.-');
    plot([Xpoints(numpoints,1);Xpoints(numpoints,2)],[Ypoints(numpoints,2),Ypoints(numpoints,2)],'m.-');
    plot([Xpoints(numpoints,1);Xpoints(numpoints,1)],[Ypoints(numpoints,1),Ypoints(numpoints,2)],'m.-');
    plot([Xpoints(numpoints,2);Xpoints(numpoints,2)],[Ypoints(numpoints,1),Ypoints(numpoints,2)],'m.-');
    text(Xpoints(numpoints,1),Ypoints(numpoints,1),int2str(numpoints),'Color','green');
    numpoints=numpoints+1;
end

% Section - get croped image from each ROI store in dataMat matrix
clear dataMat
for dcount = 1:numpoints-1
    %round mouse point to pixels in image
    xstart = min(round(Xpoints(dcount,1)),round(Xpoints(dcount,2)));
    xend = max(round(Xpoints(dcount,1)),round(Xpoints(dcount,2)));
    ystart = min(round(Ypoints(dcount,1)),round(Ypoints(dcount,2)));
    yend = max(round(Ypoints(dcount,1)),round(Ypoints(dcount,2)));
    %get region of image and display
    cropimg = frame(ystart:yend,xstart:xend);
    fig=figure;
    imagesc(cropimg)
    title(strcat("Cell No ",int2str(dcount)));
    fig.Name=strcat("Cell No ",int2str(dcount));
    dataMat{dcount}=cropimg;
end

% Section - Main loop through list of croped cell images
minint=0.10;                                     %provides floor for background estimate if to larger square selected of 
maxint=0.1;
bins=50;                                         % number of bin for histogram of pixel intensities
scoreM=[];

for dcount = 1:size(dataMat,2)-0
    data=double(dataMat{dcount});
    % get estimate of background level from histogram of pixel brightnesses
    datahist=hist(data(:),bins);
    maxval=max(data(:));
    minI=find((datahist(1:bins-1)-datahist(2:bins))<0);
    backest=minI(1);
    % set floor to background level incase ROI to large
    if (backest < 2)
        backest=minI(2);
    end
    if (backest < (bins*minint))
        backest = round(bins*minint);
    end
    %Sunbtact background and remove negative values
    backest=backest*maxval/bins;
    backdata=data-backest;
    maskeddata=(sign(backdata-0.01)+1).*backdata;
    maskeddata=maskeddata*0.5/(maxval-backest);
    
    %perform constrained grey scale errotion 
    filt=maskeddata;
    rlimits = [1.5 2];          %defines constraints
    for jcount = 1:50
        filtfun=greyerode(filt,1,rlimits);
        filt = filtfun;
    end
    
    
    % dialate and errode on binarised filter matrix - helps conect adjacent
    % end points
    bindata = ceil((filt./(max(filt(:))))-0.2);
    dilerdata = imerode(bindata,strel('square',2));
    dilerdata = imdilate(dilerdata,strel('square',2));
    score = sum(sum(dilerdata(:)))/sum(sum(bindata(:)));
    %scoreM(dcount,2) = score;
    
        
    % calc using skeliton - proportion of forground
    skeldata = bwmorph(bindata,'skeleton');
    score = sum(skeldata(:))/sum(bindata(:));
    score = (sum(bindata(:))-sum(skeldata(:)))/sum(bindata(:));
    scoreM(dcount,1) = score;
    
    
    % count number of statish line components
    stringdata=bwmorph(skeldata,'clean');
    stringdata=bwmorph(stringdata,'endpoints');
    score = sum(bindata(:))-sum(stringdata(:));
    score = score/sum(bindata(:));
    scoreM(dcount,2) = score;
    
    
    % count proportion of filterd points in bin data threshoild
    temp = ceil((maskeddata/max(maskeddata(:)))-0.01);
    score = sum(bindata(:))/sum(temp(:));
    scoreM(dcount,3) = score;
end

%Constuct table of results
scoreM=scoreM*1000;
Param1=scoreM(:,1);
Param2=scoreM(:,2);
Param3=scoreM(:,3);
cellNo=rot90(1:size(scoreM,1),3);
Results=table(cellNo,Param1,Param2,Param3)

%Format and copy results to clipboard
cliptab=convertCharsToStrings(char(9));
clipCR=convertCharsToStrings(char(13));
clipmsg=[];
for icount = 1:size(scoreM,1)
    cliptext=strcat(num2str(scoreM(icount,1)),cliptab,num2str(scoreM(icount,2)),cliptab,num2str(scoreM(icount,3)),clipCR);
    clipmsg=strcat(clipmsg,cliptext);
end
clipboard('copy',clipmsg)
%msgbox(scoreM,progname)
%clipmsg
msgbox('Analysis complete - Results in command windo & copied to clipboard, just ctrl-V in excel!',progname)



%Section - functions to perform repeated calculations
%function to erode grey levels constrained to lines and bright points
function [retmat] = greyerode(inmat,rad,range)                  %inmat=input image, rad=radius range=constraint
    xrad = round(rad);
    yrad = round(rad);
    %generte matrix of complex values of unit circle
    circmat=[exp(i*pi*3/4),exp(i*pi*2/4),exp(i*pi*1/4);exp(i*pi*4/4),0,exp(i*pi*0/4);exp(i*pi*5/4),exp(i*pi*6/4),exp(i*pi*7/4)];
    
    
    retmat = zeros(size(inmat,1),size(inmat,2));
    retmat = inmat;
    
    %loop over x,y pixels of patch
    for icount = (1+yrad):(size(retmat,1)-((1*yrad)+1))
        for jcount = (1+xrad):(size(retmat,2)-((1*xrad)+1))
            %get image patch and rank pixels by intensity
            patch = inmat(icount-yrad:icount+yrad,jcount-xrad:jcount+xrad);
            minval = max(patch(:));
            maxvals = sort(patch(:));
            nonzeros = 0;
            %get min non zero value
            for xcount = 1:numel(patch)
                if (patch(xcount) > 0)
                    nonzeros=nonzeros+1;
                    if (patch(xcount) < minval)
                        minval = patch(xcount);
                    end
                end
            end
            
            %if ((nonzeros>3) && (inmat(icount,jcount) == minval))
            if ((nonzeros>3) && (inmat(icount,jcount) < maxvals(6)))
                
                retmat(icount,jcount) = 0;
            end
            
            %Calculate shape 0f 3 point structure by dot product with unit
            %circle matrix, if ouside specified range remove point provided
            %its the lowest intensity one
            if ((nonzeros==3) && (inmat(icount,jcount) == minval))
                shape = ceil(patch/max(patch(:)));
                product=abs(sum(sum(circmat.*shape)));
                if ((product>range(1)) && (product<range(2)))
                    'trim';
                    retmat(icount,jcount)=0;
                end
            end
        end
    end
end
