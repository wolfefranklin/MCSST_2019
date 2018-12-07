
%Franklin Wolfe
%Harvard University
%Department of Earth and Planetary Sciences
%Doctoral Student
%wolfe_franklin@g.harvard.edu (primary contact)

%Modified from Tim Stahl
%University of Canterbury
%Department of Geological Sciences
%Lecturer in Tectonics
%timothy.stahl@canterbury.ac.nz

%-------------------------------------------------------------------------
%Function to allow calculation of a distributiom of net slips and slip 
%statistics (ex: rates) from input matrix X, Age vector, and Fault geometry 
%matrix.
%-------------------------------------------------------------------------

%your csv files must be in the following format:
%column 1: x values for hanging wall
%column 2: z values for hanging wall
%column 3: x values for foot wall
%column 4: z values for foot wall
%column 5: x values for scarp
%column 6: z values for scarp
%column 7: age min
%column 8: age max

%your csv files must be in a subfolder, or multiple subfolders within the
%parent directory

%-------------------------------------------------------------------------

%X is arranged as follows:
%Row 1: Hwall intercept,Hwall dev., Hangingwall slope,Slope deviation 
%Row 2: Fwall intercept,Fwall dev, Footwall slope,Fwall slope dev. 
%Row 3: Scarp intercept, Scarp dev., Scarp slope, Scarp slop dev.

%Age is a max and min value

%Fault is a matrix arranged as follows:
%Row 1: Dip mean, Dip deviation, Position mean, Position deviation
%Row 2: Dip min, max1, max2, end (trapezoidal distribution of dips)
%Row 3: Position min, max1, max2, end

%Section 0
%--------------------------------------------------------------------------
%X - Create Table
%Fault - Create Fault
%Fault has a trapezoidal distribution and the fault tip is located within
%the lower 10% of the scarp.
%nsamples - Number of Samples to run
%netslipmatrix - Create matrix to accept slip for each transect
%slipstatistics - create matrix to accept slip statistics for each transect
%filenames - create a dictionary to store file names
%myDir - set the working directory
%d,isub,foldernames - find subfolders within the working directory
X = [0,0,0,0;0,0,0,0;0,0,0,0];
Fault=[0,0,0,0;60,75,85,90;0,0,0.50,0.10];
nsamples = 10000;
netslipmatrix=[];
slipstatistics=[];
filenames={};
myDir = '/Users/franklinwolfe/Google Drive/Drum Mts Fault Zone/GIS and other data/Profiles/Extracted_Elevations_Ages/processed/';
addpath(genpath(myDir))
outputdirectory = '/Users/franklinwolfe/Google Drive/Drum Mts Fault Zone/GIS and other data/Profiles/Extracted_Elevations_Ages/Output_Tables/';
%get sub folder names in myDir
d = dir(myDir);
isub = [d(:).isdir];
foldernames = {d(isub).name}';
foldernames(ismember(foldernames,{'.','..'})) = [];


%Section 1
%--------------------------------------------------------------------------
%Code loops through subfolders of parent directory and outputs a separate
%csv file for each subfolder. Each subfolder is made up of different sets
%of profiles that you want to be analyzed and displayed separately (ex:
%profiles collected accross gulleys vs. ridges)

for p = 1:length(foldernames) %loop through each subfolder
    myDir1 = '/Users/franklinwolfe/Google Drive/Drum Mts Fault Zone/GIS and other data/Profiles/Extracted_Elevations_Ages/processed/'; %reset for each subfolder
    myDir = [myDir1 char(foldernames(p))]; %add subfolder name to path
    myFiles = dir(fullfile(myDir,'*.csv')); %find csv files within subfolder
    foldername = [foldernames(p)]; %output folder name so it can be used as column header in later table
for k = 1:length(myFiles) %loop through each csv file
    filename = myFiles(k).name; %a csv file
    %tbl = csvread(filename,1,1);
    tbl = readtable(filename); %read the csv file
    
    if isa(tbl.Var1,'double') == 0
    tbl.Var1(strcmp(tbl.Var1,'﻿')) = {''};
    tbl.Var1 = str2double(tbl.Var1);
    %tbl.Var1(1) = tbl.Var1(2) - 5;
    end
      
%clean input - rename headers
    tbl.Properties.VariableNames = {'HWx','HWy','FWx','FWy','Scx','Scy','age_min','age_max'};
    %filename(regexp(filename,'[.csv]'))=[ ];
    filename(regexp(filename,'[".csv"]'))=[ ];
    filenames = horzcat(filenames, filename);

%Determine position of hanging wall and footwall
%Determine direction of motion (down to east or down to west
%If the hanging wall X values as defined in excel are less than the
%footwall then they are to the west... ``````\,,,,, vs ,,,,,/'''''

% keep in mind here that the rest of the code defines the 'HW' as the
% actual real world footwall... not sure why this is...? But that is why
% this is reversed here
FWx = tbl.HWx; FWx(isnan(FWx))=[]; FWmeanX = mean(FWx);
HWx = tbl.FWx; HWx(isnan(HWx))=[]; HWmeanX = mean(HWx);
DownToTheEast = HWmeanX > FWmeanX;

%Run linear regression on hanging wall, footwall, and scarp values that
%were inputted in the csv files (see above on how to format csv files)
lmHW = fitlm(tbl,'HWy~HWx');
lmFW = fitlm(tbl,'FWy~FWx');
lmSc = fitlm(tbl,'Scy~Scx');

%Retrieve important statistics (intercept, standard error, slope, standard error)
%Place them in a 3X4 Table (X)

% Grab statistics
%%X = [508.0957,2.365813,0.097168,0.013109;478.0876,0.480601,-0.00648,0.006898;412.3754,0.0941179,0.710358,0.007389];
X(1,1) = lmHW.Coefficients.Estimate(1);
X(1,2) = lmHW.Coefficients.SE(1);
X(1,3) = lmHW.Coefficients.Estimate(2);
X(1,4) = lmHW.Coefficients.SE(2);
X(2,1) = lmFW.Coefficients.Estimate(1);
X(2,2) = lmFW.Coefficients.SE(1);
X(2,3) = lmFW.Coefficients.Estimate(2);
X(2,4) = lmFW.Coefficients.SE(2);
X(3,1) = lmSc.Coefficients.Estimate(1);
X(3,2) = lmSc.Coefficients.SE(1);
X(3,3) = lmSc.Coefficients.Estimate(2);
X(3,4) = lmSc.Coefficients.SE(2);

%Generate all 'distributions'--> vectors to be sampled from directly

%Create normal distribution for slope (s) and intercept (i) terms for each
%Hangingwall(H), Footwall (F), and Scarp(S)

   Hs = X(1,3) + X(1,4).*randn(nsamples,1);
   Hi = X(1,1) + X(1,2).*randn(nsamples,1);
   Fs = X(2,3) + X(2,4).*randn(nsamples,1);
   Fi = X(2,1) + X(2,2).*randn(nsamples,1);
   Ss = X(3,3) + X(3,4).*randn(nsamples,1);
   Si = X(3,1) + X(3,2).*randn(nsamples,1);
   
%Create fault geometry distributions
if Fault(1,1)~= 0 
    Dip = Fault(1,1) + Fault(1,2).*randn(nsamples,1);
else
    
%Trapezoidal distribution
%a = min;
%b1 = startofmax;
%b2 = endofmax;
%c = max;

a=Fault(2,1);
b1=Fault(2,2);
b2=Fault(2,3);
c=Fault(2,4);

for i=1:nsamples
    
h = 1/(b2-b1 + 0.5*(b1-a+c-b2));
u1 = rand(1);
if u1 < (b2-b1)*h
   x = b1 + u1/h;
else
   u1 = (u1-(b2-b1)*h)/(1-(b2-b1)*h);
   if u1 < (b1-a)/(c-b2+b1-a)
       x = a + sqrt(u1*(b1-a)*(c-b2+b1-a));
   else
       x = c - sqrt((1-u1)*(c-b2)*(c-b2+b1-a));
   end
end
 Dip(i,1) = x;
end
end   

%Fault position on scarp

%Normal
if Fault(1,3)~= 0
    Pos = Fault(1,3)+Fault(1,4).*randn(nsamples,1);
else
%Trapezoidal or Uniform
d=Fault(3,1);
e1=Fault(3,2);
e2=Fault(3,3);
f=Fault(3,4);

for j=1:nsamples
    
h = 1/(e2-e1 + 0.5*(e1-d+f-e2));
u1 = rand(1);
if u1 < (e2-e1)*h
   y = e1 + u1/h;
else
   u1 = (u1-(e2-e1)*h)/(1-(e2-e1)*h);
   if u1 < (e1-d)/(f-e2+e1-d)
       y = d + sqrt(u1*(e1-d)*(f-e2+e1-d));
   else
       y = f - sqrt((1-u1)*(f-e2)*(f-e2+e1-d));
   end
end
 Pos(j) = y;
end
end
%Find position along fault (numerical value) scarp for given propotion
%along scarp and...

%Calculate net slip using sampling with replacement from distributions
%generated
   
for k=1:nsamples
    
      
  mf(k,1)=datasample(Fs,1);
  bf(k,1)=datasample(Fi,1);
  mh(k,1)=datasample(Hs,1);
  bh(k,1)=datasample(Hi,1);
  dip(k,1)=datasample(Dip,1);
  ms(k,1)=datasample(Ss,1);
  bs(k,1)=datasample(Si,1);
  
  %Check that scarp is in correct position (i.e. that hangingwall is above
  %footwall for given scarp location)
  
  %SCENARIO 1: Downhill facing scarp, Hangingwall slope greater than
  %footwall
  
  if (abs(mh)>abs(mf))&(mh<=0)
      
  checkx1(k,1)=((bh(k,1)-bf(k,1))/(mf(k,1)-mh(k,1))); %x coordinate of hangingwall/footwall int.
  checkx2(k,1)=((bh(k,1)-bs(k,1))/(ms(k,1)-mh(k,1))); %x coordinate of hangingwall/scarp int.
  
  if checkx2(k,1)>checkx1(k,1)
     netslip(k,1)=NaN;
  else
      %Position Along scarp
  x1(k,1)=((bh(k,1)-bs(k,1))/(ms(k,1)-mh(k,1)));
  y1(k,1)=(ms(k,1)*x1(k,1))+bs(k,1);
  x2(k,1)=((bf(k,1)-bs(k,1))/(ms(k,1)-mf(k,1)));
  y2(k,1)=(ms(k,1)*x2(k,1))+bs(k,1);
  slopedistance(k,1)=sqrt(((x2(k,1)-x1(k,1))^2)+((y2(k,1)-y1(k,1))^2));
  distance(k,1)=slopedistance(k,1)*datasample(Pos,1);
        if x2>=x1
      angle(k,1)=asind(((abs(x2(k,1)-x1(k,1))))/(slopedistance(k,1)));
      x3(k,1)=x1(k,1)+distance(k,1)*sind(angle(k,1));
      pos(k,1)= x3(k,1); 
         else
      angle(k,1)=asind(((abs(x1(k,1)-x2(k,1))))/(slopedistance(k,1)));
      x3(k,1)=x1(k,1)-distance(k,1)*sind(angle(k,1));
      pos(k,1)= x3(k,1); 
        end
      
  %Net slip  
  fslip(k,1)= (pos(k,1)*(ms(k,1)-mf(k,1))+bs(k,1)-bf(k,1))/(sind(dip(k,1))+mf(k,1)*(cosd(dip(k,1))));
  
  hslip(k,1)= (pos(k,1)*(mh(k,1)-ms(k,1))+bh(k,1)-bs(k,1))/(sind(dip(k,1))+mh(k,1)*(cosd(dip(k,1))));
  
  netslip(k,1)= fslip(k,1) + hslip(k,1);
  end
   
  %SCENARIO 2: Downhill-facing scarp, Left-to-right, Footwall slope greater than
  %hangingwall
  elseif abs(mf)>abs(mh)&(mh<0)
  checkx1(k,1)=((bh(k,1)-bf(k,1))/(mf(k,1)-mh(k,1))); %x coordinate of hangingwall/footwall int.
  checkx2(k,1)=((bh(k,1)-bs(k,1))/(ms(k,1)-mh(k,1))); %x coordinate of hangingwall/scarp int.
  
        if checkx1(k,1)>checkx2(k,1)
        netslip(k,1)=NaN;
        else
  
  %Position Along scarp
  x1(k,1)=((bh(k,1)-bs(k,1))/(ms(k,1)-mh(k,1)));
  y1(k,1)=(ms(k,1)*x1(k,1))+bs(k,1);
  x2(k,1)=((bf(k,1)-bs(k,1))/(ms(k,1)-mf(k,1)));
  y2(k,1)=(ms(k,1)*x2(k,1))+bs(k,1);
  slopedistance(k,1)=sqrt(((x2(k,1)-x1(k,1))^2)+((y2(k,1)-y1(k,1))^2));
  distance(k,1)=slopedistance(k,1)*datasample(Pos,1);
        if x2>=x1
      angle(k,1)=asind(((abs(x2(k,1)-x1(k,1))))/(slopedistance(k,1)));
      x3(k,1)=x1(k,1)+distance(k,1)*sind(angle(k,1));
      pos(k,1)= x3(k,1); 
        else
      angle(k,1)=asind(((abs(x1(k,1)-x2(k,1))))/(slopedistance(k,1)));
      x3(k,1)=x1(k,1)-distance(k,1)*sind(angle(k,1));
      pos(k,1)= x3(k,1); 
        end
      
  %Net slip  
  fslip(k,1)= (pos(k,1)*(ms(k,1)-mf(k,1))+bs(k,1)-bf(k,1))/(sind(dip(k,1))+mf(k,1)*(cosd(dip(k,1))));
  
  hslip(k,1)= (pos(k,1)*(mh(k,1)-ms(k,1))+bh(k,1)-bs(k,1))/(sind(dip(k,1))+mh(k,1)*(cosd(dip(k,1))));
  
  netslip(k,1)= fslip(k,1) + hslip(k,1);
        end
        
  %SCENARIO 3: Downhill facing scarp, Right-to-left, Hangingwall slope
  %greater than footwall
   elseif abs(mh)>abs(mf)&(mh>0)
  checkx1(k,1)=((bh(k,1)-bf(k,1))/(mf(k,1)-mh(k,1))); %x coordinate of hangingwall/footwall int.
  checkx2(k,1)=((bh(k,1)-bs(k,1))/(ms(k,1)-mh(k,1))); %x coordinate of hangingwall/scarp int.
  
        if checkx1(k,1)>checkx2(k,1)
        netslip(k,1)=NaN;
        else
  
  %Position Along scarp
  x1(k,1)=((bh(k,1)-bs(k,1))/(ms(k,1)-mh(k,1)));
  y1(k,1)=(ms(k,1)*x1(k,1))+bs(k,1);
  x2(k,1)=((bf(k,1)-bs(k,1))/(ms(k,1)-mf(k,1)));
  y2(k,1)=(ms(k,1)*x2(k,1))+bs(k,1);
  slopedistance(k,1)=sqrt(((x2(k,1)-x1(k,1))^2)+((y2(k,1)-y1(k,1))^2));
  distance(k,1)=slopedistance(k,1)*datasample(Pos,1);
        if x2>=x1
      angle(k,1)=asind(((abs(x2(k,1)-x1(k,1))))/(slopedistance(k,1)));
      x3(k,1)=x1(k,1)+distance(k,1)*sind(angle(k,1));
      pos(k,1)= x3(k,1); 
        else
      angle(k,1)=asind(((abs(x1(k,1)-x2(k,1))))/(slopedistance(k,1)));
      x3(k,1)=x1(k,1)-distance(k,1)*sind(angle(k,1));
      pos(k,1)= x3(k,1); 
        end
      
  %Net slip  
  fslip(k,1)= (pos(k,1)*(ms(k,1)-mf(k,1))+bs(k,1)-bf(k,1))/(sind(dip(k,1))+mf(k,1)*(cosd(dip(k,1))));
  
  hslip(k,1)= (pos(k,1)*(mh(k,1)-ms(k,1))+bh(k,1)-bs(k,1))/(sind(dip(k,1))+mh(k,1)*(cosd(dip(k,1))));
  
  netslip(k,1)= fslip(k,1) + hslip(k,1);
        end
    
  %SCENARIO 4: Downhill-facing scarp, Right-to-left, Footwall slope greater
  %than hangingwall
  
     elseif abs(mf)>abs(mh)&(mh>0)
  checkx1(k,1)=((bh(k,1)-bf(k,1))/(mf(k,1)-mh(k,1))); %x coordinate of hangingwall/footwall int.
  checkx2(k,1)=((bh(k,1)-bs(k,1))/(ms(k,1)-mh(k,1))); %x coordinate of hangingwall/scarp int.
  
        if checkx2(k,1)>checkx1(k,1)
        netslip(k,1)=NaN;
        else
  
  %Position Along scarp
  x1(k,1)=((bh(k,1)-bs(k,1))/(ms(k,1)-mh(k,1)));
  y1(k,1)=(ms(k,1)*x1(k,1))+bs(k,1);
  x2(k,1)=((bf(k,1)-bs(k,1))/(ms(k,1)-mf(k,1)));
  y2(k,1)=(ms(k,1)*x2(k,1))+bs(k,1);
  slopedistance(k,1)=sqrt(((x2(k,1)-x1(k,1))^2)+((y2(k,1)-y1(k,1))^2));
  distance(k,1)=slopedistance(k,1)*datasample(Pos,1);
        if x2>=x1
      angle(k,1)=asind(((abs(x2(k,1)-x1(k,1))))/(slopedistance(k,1)));
      x3(k,1)=x1(k,1)+distance(k,1)*sind(angle(k,1));
      pos(k,1)= x3(k,1); 
        else
      angle(k,1)=asind(((abs(x1(k,1)-x2(k,1))))/(slopedistance(k,1)));
      x3(k,1)=x1(k,1)-distance(k,1)*sind(angle(k,1));
      pos(k,1)= x3(k,1); 
        end
      
  %Net slip  
  fslip(k,1)= (pos(k,1)*(ms(k,1)-mf(k,1))+bs(k,1)-bf(k,1))/(sind(dip(k,1))+mf(k,1)*(cosd(dip(k,1))));
  
  hslip(k,1)= (pos(k,1)*(mh(k,1)-ms(k,1))+bh(k,1)-bs(k,1))/(sind(dip(k,1))+mh(k,1)*(cosd(dip(k,1))));
  
  netslip(k,1)= fslip(k,1) + hslip(k,1);
        end   
  
  else
      
   %Scenario 5: All other configurations:
   
   %Position Along scarp
  x1(k,1)=((bh(k,1)-bs(k,1))/(ms(k,1)-mh(k,1)));
  y1(k,1)=(ms(k,1)*x1(k,1))+bs(k,1);
  x2(k,1)=((bf(k,1)-bs(k,1))/(ms(k,1)-mf(k,1)));
  y2(k,1)=(ms(k,1)*x2(k,1))+bs(k,1);
  slopedistance(k,1)=sqrt(((x2(k,1)-x1(k,1))^2)+((y2(k,1)-y1(k,1))^2));
  distance(k,1)=slopedistance(k,1)*datasample(Pos,1);
        if x2>=x1
      angle(k,1)=asind(((abs(x2(k,1)-x1(k,1))))/(slopedistance(k,1)));
      x3(k,1)=x1(k,1)+distance(k,1)*sind(angle(k,1));
      pos(k,1)= x3(k,1); 
        else
      angle(k,1)=asind(((abs(x1(k,1)-x2(k,1))))/(slopedistance(k,1)));
      x3(k,1)=x1(k,1)-distance(k,1)*sind(angle(k,1));
      pos(k,1)= x3(k,1); 
        end
      
  %Net slip  
  fslip(k,1)= (pos(k,1)*(ms(k,1)-mf(k,1))+bs(k,1)-bf(k,1))/(sind(dip(k,1))+mf(k,1)*(cosd(dip(k,1))));
  
  hslip(k,1)= (pos(k,1)*(mh(k,1)-ms(k,1))+bh(k,1)-bs(k,1))/(sind(dip(k,1))+mh(k,1)*(cosd(dip(k,1))));
  
  netslip(k,1)= fslip(k,1) + hslip(k,1);
  end
      
end
%remove NaNs and replace with a random pick from the slip vector.
netslip(isnan(netslip))=datasample(netslip,1);
%netslip(isnan(netslip))=[];  %Remove NaNs from dataset to perform statistics

%Append net slip distribution to a matrix that will contain the slip values for each transect
netslipmatrix = [netslipmatrix netslip];

%--------------------------------------------------------------------------
%Slip Statistics
slipmean = mean(netslip);
slipmedian = median(netslip);
slipmode = mode(netslip);
SEM = std(netslip)/sqrt(length(netslip)); %Standard Error
ts = tinv(0.975,length(netslip)-1); %T Statistic
CIhigh = mean(netslip) + ts*SEM; %Confidence Intervals
CIlow = mean(netslip) - ts*SEM;  

%Create age vector
age_max = tbl.age_max;
age_min = tbl.age_min;
    
age_max(isnan(age_max))=[];
age_min(isnan(age_min))=[];
age = (age_min+rand(1,nsamples)*(age_max-age_min)).';

%Create fault angle vector
fault_angle = (Fault(2,1)+rand(1,nsamples)*(Fault(2,4)-Fault(2,1))).';

%Slip rate statistics
slip_rate = netslip ./ (age./1000);
slip_rate_mean = mean(slip_rate);
slip_rate_SEM = std(slip_rate)/sqrt(length(slip_rate)); %Standard Error
slip_rate_ts = tinv(0.975,length(slip_rate)-1); %T Statistic
slip_rate_CIhigh = mean(slip_rate) + slip_rate_ts*slip_rate_SEM; %Confidence Intervals
slip_rate_CIlow = mean(slip_rate) - slip_rate_ts*slip_rate_SEM;  

%Slip rate - horizontal component 
slip_rate_horiz = slip_rate .* cos(deg2rad(fault_angle));
slip_rate_horiz_mean = mean(slip_rate_horiz);
slip_rate_horiz_SEM = std(slip_rate_horiz)/sqrt(length(slip_rate_horiz)); %Standard Error
slip_rate_horiz_ts = tinv(0.975,length(slip_rate_horiz)-1); %T Statistic
slip_rate_horiz_CIhigh = mean(slip_rate_horiz) + slip_rate_horiz_ts * slip_rate_horiz_SEM; %Confidence Intervals
slip_rate_horiz_CIlow = mean(slip_rate_horiz) - slip_rate_horiz_ts * slip_rate_horiz_SEM; 

%Slip rate - vertical component 
slip_rate_vert = slip_rate .* sin(deg2rad(fault_angle));
slip_rate_vert_mean = mean(slip_rate_vert);
slip_rate_vert_SEM = std(slip_rate_vert)/sqrt(length(slip_rate_vert)); %Standard Error
slip_rate_vert_ts = tinv(0.975,length(slip_rate_vert)-1); %T Statistic
slip_rate_vert_CIhigh = mean(slip_rate_vert) + slip_rate_vert_ts * slip_rate_vert_SEM; %Confidence Intervals
slip_rate_vert_CIlow = mean(slip_rate_vert) - slip_rate_vert_ts * slip_rate_vert_SEM; 

%slip mode
slipmode = mode(round(netslip,2));

%create a vector of statistics
statistics = [slipmean; slipmedian; slipmode; CIlow; CIhigh; slip_rate_mean; slip_rate_CIlow; slip_rate_CIhigh; slip_rate_horiz_mean; slip_rate_horiz_CIlow; slip_rate_horiz_CIhigh; slip_rate_vert_mean; slip_rate_vert_CIlow; slip_rate_vert_CIhigh; DownToTheEast; age_min; age_max];

%Append statistics for one transectto a matrix that will contain the
%statistics for all transects
slipstatistics = [slipstatistics statistics];


% % % %--------------------------------------------------------------------------
% %Histograms of slip distributions for each section
% %
% %Beautify
figurename = ['Slip Distributions for ', char(filename), ' (file name)' ,' in ' , char(foldername), ' (folder name)'];
figure
set( gcf, 'Color', 'White', 'Unit', 'Normalized', 'Position', [0.1,0.1,0.6,0.6] ) ;
axes( 'Position', [0, 0.95, 1, 0.05] ) ;
set( gca, 'Color', 'None', 'XColor', 'White', 'YColor', 'White' ) ;
text(0.5,0,figurename,'FontSize', 20, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' )

%Plot
subplot(3,2,[1,2])
histogram(netslip)
title('Slip')
subplot(3,2,3)
histogram(age)
title('Age')
subplot(3,2,4)
histogram(slip_rate)
title('Slip Rate Total')
subplot(3,2,5)
histogram(slip_rate_horiz)
title('Slip - Horizontal Component')
subplot(3,2,6)
histogram(slip_rate_vert)
title('Slip - Vertical Component')

end
%--------------------------------------------------------------------------


%Section 2
%--------------------------------------------------------------------------
%post processing for table of slip values to beautify output
netslipmatrix_table = array2table(netslipmatrix,'VariableNames',filenames);

%post processing for table of stats
slipstatistics = slipstatistics;
slipstatistics_table = [];

%assign headers
header = {'SlipMean','SlipMedian','SlipMode','CILow','CIHigh','SlipRate_mean','SlipRate_CIlow', 'SlipRate_CIhigh', 'SlipRate_avg_horiz', 'SlipRate_CILow_horiz', 'SlipRate_CIhigh_horiz', 'SlipRate_avg_vert', 'SlipRate_CIlow_vert', 'SlipRate_CIhigh_vert','Direction (E = 1)','AgeMin','AgeMax'};
slipstatistics_table = array2table(slipstatistics,'RowNames',header,'VariableNames',filenames)%,'VariableNames',header)%,'RowNames',filenames);

%Creates a new folder to store the output tables
filenameoutput = [outputdirectory, 'Slip_Statistics_',char(foldername),'.csv'];

%write output csv files for each subfolder within Profiles and stores them
%in new folder
writetable(slipstatistics_table,filenameoutput)
    
%reset tables and dictionary
netslipmatrix=[];
slipstatistics=[];
filenames={};
ages=[];
end
          