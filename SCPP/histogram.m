
%Author- ADITYA KUMAR

%--------------------------------------------------------------------------
%DIGITAL SIGNAL PROCESSING 
% To create histogram & perform brightness and contrast
% operations on an image
%--------------------------------------------------------------------------


%-------------------------------------------------------------------------
%   LIST OF VARIABLES AND ASSOCIATED FUNCTIONS
%-------------------------------------------------------------------------
%filename                       Store path for fetching image
%img                            Input image data
%bright_scale                   Any value b/w 0 to 85
%bright_choice                  Any one of [-3 -2 -1 0 1 2 3]
%bright_level                   bright_choice*bright*factor
%contrst_scale                  Any value greater than 1
%contrst_choice                 Any one of [-3 -2 -1 0 1 2 3]
%contrst_level                  1+((contrst_choice+.1)/contrst_scale)
%rhist, ghist, bhist            Bin to store correspondin levels of R G B
%x,y                            image size
%show_img                       create original image window
%show_hist                      create histogram window
%show_processd_img              create processed image window
%choice                         choose processing options
%-------------------------------------------------------------------------

%       PROGRAM

%-------------------------------------------------------------------------


clc;
close all;
bright_scale=50; %change it to change brightness factor 
                 %(should be less than 85)
                 
contrst_scale=6; %change to change contrast factor

filename = input('Enter link of image : ','s'); %e.g - rose2.jpg
img = imread(filename);
[x,y,z]=size(img);
show_img = figure('Name','original Image','NumberTitle','off');
figure(show_img);
imshow(img);

while 1

%clc; %clear command window at the start of each loop

%-----Initiallising RGB bins--------------------------------------------
rhist = ones(1,256);
if(z>1)
ghist = ones(1,256);
bhist = ones(1,256);
end
%rhist,ghist,bhist has been initiallised wid 1's instead of 0's so that
%taking log of unfilled levels(after scanning image) would not yield
%minus infinity (log of zero is minus infinity) while plotting them

%if(z>1) condition used at many places bcoz greyscale image will not have
%RGB component. It'll have greyscale values of range 0 to 255
%see line 44. z returns 1 for black n white , 3 for coloured images.


%-----Scanning image pixels----------------------------------------------
    for i = 1:x
        for j = 1:y
            rhist(img(i,j,1)+1)=rhist(img(i,j,1)+1)+1;
            if(z>1)
            ghist(img(i,j,2)+1)=ghist(img(i,j,2)+1)+1;
            bhist(img(i,j,3)+1)=bhist(img(i,j,3)+1)+1;
            end
        end
    end

%------------------------------------------------------------------------
   
   rhist = double(log(rhist));
   if(z>1)
   ghist = double(log(ghist));
   bhist = double(log(bhist));
   end
%-------Create Histogram window---------------------------------------
   show_hist = figure('Name',' Histogram','NumberTitle','off');
   figure(show_hist);
   
%-------Plotting Histogram--------------------------------------------
    
   if(z>1)
    level=0:1:255;
    bar(level,rhist,'Barwidth',1,'Facecolor',[1 0 0],'Edgecolor',[1 0 0]);
   
    hold on;
    bar(level,ghist,'Barwidth',1,'Facecolor',[0 1 0],'Edgecolor',[0 1 0]);
    hold on;
    
    bar(level,bhist,'Barwidth',1,'Facecolor',[0 0 1],'Edgecolor',[0 0 1]);
    hold off;
   else
    level=0:1:255;
    bar(level,rhist,'Barwidth',1,'Facecolor',[0 0 0],'Edgecolor',[1 1 1]);
   
   end
    axis tight;
%-------END OF HISTOGRAM FUNCTION-------%
   
    
    
%--------Start Various image processing operations-------------------
    display('****Image Processing Functions****');
    display('Press.........');
    display('1 To set Brightness');
    display('2 To set Contrast');
    display('0 To Exit');
    display('     ');
    choice = input('Enter Your Choice : ','s');
    switch choice
        case '1'
            bright_choice=input('Enter Brightness level( Any value b/w -3 to 3): ');
            bright_level=bright_choice*bright_scale;
            img=img+bright_level;
            
        case '2'
            contrst_choice=input('Enter Contrast level( Any value b/w -3 to 3): ');
            contrst_level= 1+((contrst_choice+0.1)/contrst_scale);
            img = img*contrst_level;
    
        
        case '0'
            clc;
            close all;
            return;

        otherwise
            continue;
    end
%----------end of switch----------------------------------------------
 

    close all;
    
    show_processd_img= figure('Name','Processed Image','NumberTitle','off');
    figure(show_processd_img);
    image(img);
end
