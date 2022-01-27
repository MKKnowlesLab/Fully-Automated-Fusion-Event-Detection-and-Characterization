
function  pro_cess_fusion_data(imdata,imdata2)
% pro_cess(imdata)
% imdata is the imaging data in tiff file format
%
% PURPOSE: 1) to process imaging data and detect fusion event. 2)Analyze
% detected fusion events to obtain various kinetic parameters.
% movie
% INPUT:
% imdata: movie stack, in tiff format. 
% OUTPUT:  an excel containing: 
% 1)Expression level of cell (background corrected)  
% 2)fusion event traces with fusion event onset aligned to t=0,
% 3)Normalized event traces
% 4)Slope of each fusion event
% 5)% intensity loss at different time points
% 6)Decay curve fit and rate of decay
% 7)Radial plots of each curve starting t=0
% 8)Full width half max at different time points starting with t=0
% 9)Area under curve 
%
% PROCESSING OVERVIEW(please look at each program below for details):
% I.   Read tiff movie
% II.  Generate Cell Mask, Background Mask and Relative Expression Level
% III. Generate Max projection of difference movie and find fusion event
% IV. Seperate each fusion event and obtain fluorescence intensity of the
% fusion event (circle) and of the background (annulus)
% V. Obtain decay rates, radial plots, HWHM, Area under RP.
% VI. Generate time aligned movies and montages for each fusion event

%
%
% I. Read tiff movie
    info = imfinfo(imdata);
    numframe = length(info);
    for K = 1 : numframe
        rawframes(:,:,:,K) = imread(imdata, K);  
    end
    %convert data to 3D
    rawframes = rawframes(:,:,:);

    % for dual channel input, nargin>1, process second input movie
    if nargin >1
    info2 = imfinfo(imdata2);
    numframe2 = length(info2);
    for K = 1 : numframe2
        rawframes2(:,:,:,K) = imread(imdata2, K);  
    end
    %convert data to 3D
    rawframes2 = rawframes2(:,:,:);
    end

    
% II. Generate Cell Mask, Background Mask and Relative Expression Level 
 % Use first frame to obtain masks and write out image
    CFimage = rawframes(:,:,1);
    imwrite(CFimage,'CFimage.tif','compression','none'); 
    %Obtain cell and background masks
    [~,cellmask,~,backgroundmask]=obtaincellmasks('CFimage.tif');
    %calculate relative expression level(RExLevel) =cell intensity in the first
    %image minus background
    CFimage=double(CFimage);
    % multiply by cell mask so intensity value outside cell will become
    % zero
    CF2=CFimage.*cellmask;
    % Calculate Avg cell intensity in first image of the cell
    CFavg = mean(nonzeros(CF2));
    % Calculate Avg background in first image of the cell using background 
    %mask
    CFb=CFimage.*backgroundmask;
    CFbAvg = mean(nonzeros(CFb));
    % Relative Expression level = Average Cell intensity minus average
    % background intensity
    RExLevel=CFavg-CFbAvg;
    % Create excel to output data
  	curr_directory =  pwd;current_folder_name = curr_directory(end-3:end);
    % current_folder_name = name of folder 
    output_excel_file = sprintf('Analyzed Data_%s.xls',current_folder_name);
    % Data is output in an excel file, and the data is converted to excel
    % compatible format
    RExLevel=array2table(RExLevel);
    writetable(RExLevel, output_excel_file,'Sheet', 'RExLevel');
    
% III. Generate Max projection of difference movie and find fusion event
    % Generate Max projection image
    [~,diff_max_proj] = CalcDifferenceMovie(rawframes,0,25);
    % Multiply Max projection and cell mask to eliminate background
    % intensity, this converts all background values to 0
    diff_max_proj=diff_max_proj.*cellmask;
    %Intensity values within cell perimeter are >0
    cell_Intensity_values = (diff_max_proj)>0;
    %Increase brightness of the bright spots for enhanced contrast and
    %easier peak selection
    max_proj_2=diff_max_proj.*cell_Intensity_values;
    % enanced intensity values located within cell perimeter
    cell_Intensity_enhanced = (nonzeros(max_proj_2));
    %Find threshold to detect fusion event locations using pkfnd 
    %Set threshold using average brigntness of the cell
    pk_thresh=(mean(cell_Intensity_enhanced));
    %find location of brigt peaks
    pk_loc = pkfnd (diff_max_proj,pk_thresh,30);
  
% IV. Seperate each fusion event and obtain fluorescence intensity of the
    % fusion event (circle) and of the background (annulus)

    %Crop out 25 x 25 pixel regions witht fusion events centered for indiviual
    %analysis of each event 
    %Create 7x7 pixel circular mask 
    Indiv_events = ministk_movie(rawframes,pk_loc,25,0);
    cir_mask = create_mask(25,7);
    outcir_mask = ~(cir_mask);
    outcir_mask(1:5,:)=0;outcir_mask(:,1:5)=0;
	outcir_mask(:,20:25)=0;outcir_mask(20:25,:)=0;

if nargin>1
 %Obtain array with timecourse average intensity
    %Array containing fusion event intensity and protein event intensity
    Indiv_eventsP = ministk_movie(rawframes2,pk_loc,25,0);
     Int_ArrayP=ministk_intensityavg(Indiv_eventsP,cir_mask);
    Int_ArraPO=ministk_intensityavg(Indiv_eventsP,outcir_mask);
    %Array containing fusion event intensity     
    Int_DiffP=Int_ArrayP-Int_ArraPO;
end

    %Obtain array with timecourse average intensity
    %Array containing fusion event intensity and protein event intensity
    Int_Array=ministk_intensityavg(Indiv_events,cir_mask);
    Int_ArraO=ministk_intensityavg(Indiv_events,outcir_mask);
  
    %Array containing fusion event intensity     
    Int_Diff=Int_Array-Int_ArraO;   
%Filter unreal peaks and 
%Align Fusion events by setting onset of fusion as T=0 
  if nargin <2
      [Aligned_Pks] = pk_align(Int_Diff,50,CFbAvg);
      [Aligned_PksC] = compile_data(Aligned_Pks);
      [ND1] = normalizedata2(Aligned_PksC);
  end
    %Filter unreal peaks and 
    %Align Fusion events by setting onset of fusion as T=0 
if nargin >1 
    %Compiled Aligned pks-filter out 0s
    %Normalize Pks 
    %Output intensity of fusion events in excel 
    [Aligned_Pks, Aligned_Pks_P,FoverS,Aligned_Pks_I] = pk_align(Int_Diff,50,Int_DiffP,Int_Array,CFbAvg);
    [Aligned_PksC,Aligned_PksC_P] = compile_data(Aligned_Pks,Aligned_Pks_P,FoverS,Aligned_Pks_I);
    [ND1,ND_P] = normalizedata2(Aligned_PksC,Aligned_PksC_P);
end
    %Compiled Aligned pks-filter out 0s
    %Normalize Pks 
    %Output intensity of fusion events in excel 
    
    
% V. Obtain decay rates, radial plots, HWHM, Area under RP
    %Fit Decay Curves to a double exp function and
    %Obtain fit parameters
        decayfit2(ND1);
    %MAKE Radial Plots, and calculate HWHM,and Area under radial plot
        %radial_plot(Aligned_Pks,Int_Diff,Indiv_events,Aligned_PksC);
        
% VI. Generate time aligned movies and montages for each fusion event
        
%Create a 25 x25 pixel movies for each fusion event with fusion event of
%middle of the movie
if nargin ==1
     align_movies(Aligned_Pks,Indiv_events);
     create_montage(Aligned_Pks,Indiv_events);
elseif nargin > 1
 align_movies(Aligned_Pks,Indiv_events,Indiv_eventsP);
 create_montage(Aligned_Pks,Indiv_events,Indiv_eventsP);
end
%Create 1s interval montages
     create_montage(Aligned_Pks,Indiv_events,Indiv_eventsP);
% Save Matlab Workspace 
    curr_directory =  pwd;current_folder_name = curr_directory(end-3:end);
    save_workspace= sprintf('%s Analysis.mat',current_folder_name);
    save(save_workspace);
end

%DIFFERENCE MOVIE
function [mov_out,maxproj,max_sub] = CalcDifferenceMovie( mov, filter, dt )
%PURPOSE: to calculate the difference in intensity from one frame to the
%next, output is a movie of the difference images.
%   Call:mot=CalcDifferenceMovie(a, 9, 1);
%   Variables: a is a movie.  The format is usually (256,256, nframes)
%   filter is the size of the bandpass filter, and dt is the time lag in
%   frames over which the difference movie is calculated 
%   Output: this calculates the percent of pixels that changed from one frame to the next
%   Steps - calculate difference images (frame t+1) - (frame t).
%  MKK Dec 2 2019
%
% max_sub=zeros(256,512);
nframes = max(size(mov(1,1,:)));
x = max(size(mov(:,1,1)));
y = max(size(mov(1,:,1)));


if nargin <2
    filter=7;
end    
if nargin <3
     dt=25;
%dt=floor(nframes*(0.025));
end   
 
mov_out = zeros(x, y, nframes-dt);
for i= 1:nframes-dt
    %filter images at frame t and t+dt
    im = bpass(mov(:,:,i),1,filter);
    im2 = bpass(mov(:,:,i+dt),1,filter);

    % take the difference image from frame t+1 and frame t to determine new
    % pixels with values. This is where signal moved TO. 
    delta = im2 +1000 - im;
    %store the image in a movie file
    mov_out(:,:,i) = delta;

end
%should backgrd always be 35 frames long? 
%Try backgrd with fewer frames to check how to
%automate length of bck
max_bckgrd = max(mov_out(:,:,1:35), [], 3);
max_rem= max(mov_out(:,:,36:end), [], 3);
max_sub=max_rem-max_bckgrd;

maxproj = max(mov_out, [], 3);
%figure;imshow(maxproj, [1000,1500]); 
end


%CIRCLE MASK
function [cir_mask] = create_mask(sz,d)
 %[cir_mask] = create_mask(25,7)
% PURPOSE:  to make a mask of size sz and diameter d
% INPUT:
% sz: in units of pixels,size of mask
% d:the diameter of the circular mask
% OUTPUT:  an array that is very long (number of spots found in each frame)
% that has X, Y, time 

mask = zeros(sz,sz);
cntr =ceil(sz/2);
% total rows/colums -1
trc=d-1; 
% rows/colums above/below center-halof of totakl rows/col=3
h=trc/2;
%center=13
for n=cntr-3:cntr+3
mask(n,cntr-1:cntr+1)=1;
end
for n=cntr-2:cntr+2
mask(n,cntr-2:cntr+2)=1;
end
for n=cntr-1:cntr+1
mask(n,cntr-3:cntr+3)=1;
end
[cir_mask] =mask;
end

%AVERAGE INTENSITY
function out=ministk_intensityavg(mov,msk)
% out=ministk_intensity(images, mask)
% 
% PURPOSE:  to measure the intensity in the center of a 25 x 25 image
% INPUT:
% mov: movie stack to measure from
% msk: an image that is 1 in the center and 0 outside the center
%
% OUTPUT:  an array that is N_frames long by the number of spots. Each
% column is a new location.
%           
% CREATED: Michelle Knowles Jan 2018

nrgn = length(mov(1,1,1,:));
nframes = length(mov(1,1,:,1));
intensity = zeros(nframes,nrgn);
size(intensity)
%loop through all spots
for i =1:nrgn;
    %loop through all frames - there's got to be a faster way!
    for j =1:nframes;
        a=double(mov(:,:, j, i)).*msk;
        intensity(j,i)=sum(sum(a))/sum(msk(:));
    end
end

out=intensity;
end


%CROP 25 x 25 pixel regions
function out=ministk_movie(im,rgn,sz,sepdist)
% out=ministk(im,rgn,sz)
% 
% PURPOSE:  to cut out small regions of an image based on spots found in
% the image or a corresponding image of a different color. This is used to
% measure colocalization based on the work of Knowles and Barg in two PNAS
% 2011 papers. Spots are found using the work of tracking routines
% available on Eric Weeks' website (Emory University) and made into Matlab
% by Eric Dufrense.
% INPUT:
% im: movie stack to cut from
% rgn: spots (x,y) about which regions should be cut
% sz: size of cut out (a square of sz by sz pixels)
% sepdist: is the minimum separation distance between two spots. If two
% spots are within this distance of one another, neither are counted.
%
% OUTPUT:  a sz x sz x N_frames array for each spot
%           
% CREATED: Michelle Knowles May 2012
% EDITED: to crop movies rather than single images by Michelle Knowles and
% Aubrie Blevins December 2017, note that X and Y from the peak finder seem
% to be switched when cropping in this way.

%if sz/2 == floor(sz/2)
%warning('sz must be even so that the spots can be centered on a pixel: 1 pixel added');
%sz = sz+1
%end

%scott's code for a mask
pix=(sz+1)/2;
dimy = length(im(1,:,1));
dimx = length(im(:,1,1));
nrgn = length(rgn(:,1));
nframes = length(im(1,1,:));
sepdist =sepdist;
%create a blank image array that you can fill
msk=zeros([sz,sz,nframes]);
%loop through all regions that locate spots in an image
for i =1:nrgn;
    x = rgn(i,1);
    y = rgn(i,2);
    %don't include regions within pix distance from the edge of the image.
    if ((x>pix) && ((x+pix)<dimy) && (y>pix) && ((y+pix)<dimx))
        % don't include regions within "sepdist" pixels of another region
        % calculate an array that contains the difference between the
        % current particle's postion and all others.
        if sepdist < 0
        diffy=rgn(:,2)-rgn(i,2);
        diffx=rgn(:,1)-rgn(i,1);
        mag=((diffx.*diffx)+(diffy.*diffy).^0.5);
        % find all the locations in the magnitude array that are non-zero.
        % This should remove the comparison between particle i and itself,
        % which will always be zero. 
        w = find(mag);
        mag = mag(w);
        if (min(mag) > sepdist)  
            %how can we crop a MOVIE, do we need to loop over all images,
            %concatenate to make "cutout" then concatenate cutout to msk?
            cutout=im(y-pix+1:y+pix-1, x-pix+1:x+pix-1,:);
          msk = cat(4, msk, cutout);
        end
        
        elseif sepdist==0
        cutout=im(y-pix+1:y+pix-1, x-pix+1:x+pix-1,:);
          msk = cat(4, msk, cutout);
        end
    out=msk;
    end



end
end


%Align FusionEvent Traces

function [out_Aligned_PksC,out_Aligned_PksC_P,outFoverS,out_Aligned_PksC_I]= pk_align(Int_Array,exposure_n,Int_DiffP,Int_2Array,CFbAvg)

%PURPOSE: to determine the location of fusion event onset and artifically set the determined
%fusion onset point to Time = 0s

% When calling function pk_align(Int_Array,exposure_n,Int_DiffP,CFbAvg), each INPUT
% file is as follows:
%       1) Int_Array: An array containing Average Intensity values, each
%          fusion event as one column 
%       2) exposure_n = exposure value (in milliseconds) used when collecting data, this
%          represents te time interval used to collect each data point and this value is
%          used to calculate the time range and interval for the
%          fusion event
%       3) Int_DiffP:  An array containing Average Intensity values from the second
%          color channel, if any.
%       4) CFbAvg= Cell background value calculated before function is
%          called and used for FoverS calculation.
% OUTPUT includes:
%    1) out_Aligned_PksC = fusion event intensity data, with onset of fusion aligned to T=0s 
%    2) out_Aligned_PksC_P= for two channel data, this output aligns the
%    intensity from te second channel to the onset of fusion aligned to T=0s   
%    3)  outFoverS = F over S is divides Fluorescence intesnity values by the Cell background 
%    value (CFbAvg-if provided when calling the funcion)  
%       
%   
%Determine number of columns (ncol) and number of rows (nr) of the input
%file
ncol = length(Int_Array(1,:));
[nr,~] = size(Int_Array);
%
% Create empty arrays to fill with output data
out_Aligned_PksC = zeros(nr,ncol);
out_Aligned_PksC_P = zeros(nr,ncol);
out_Aligned_PksC_I = zeros(nr,ncol);
outFoverS=zeros(nr,ncol);

%Loop through all fusion events, first column has zeros, so loop begins at
%column (c) =2

for c = 2:ncol
   
%Isolate the fusion event being aligned &
%Find the max value (Max_value), and location of the max value (Loc_max) for the isolated
%event
Fusion_event = Int_Array (:,c);  
[Max_value,Loc_max] = max (Fusion_event);
 
%keep peaks with at least a few points (0.99* nr (total number of rows))
%if not enough points before fusion event, do not include event in analysis
%and continue analysis of the next column

if Loc_max >0.99*nr
continue 
end

% Find the range of data points to calculate background average for the
% fusion event
% Loc_max = location of max, st= starting point for background avg. calculation
%en = end point for background avg. calculation
if Loc_max<20 
start_time =1;
en=5;
else 
start_time=1;
en=20;
end
% Calculate background average for fusion event (Bck_avg_FE)
%and only keep fusion events with max value > 1.4x Bck_avg_FE 

Bck_avg_FE = mean(Fusion_event(start_time:en,1));

if Bck_avg_FE>1
Bck_avg_FE=Bck_avg_FE;
else
Bck_avg_FE=1;
end

Max_Bckgrd_ratio= Max_value/Bck_avg_FE;

%Additonally, only keep peaks that maintain more that 60% of the Intensity at Max+1
% data point. This is to exclude pixel noise which generally only lasts one data
%point

if  Max_Bckgrd_ratio > 1.4 && Fusion_event(Loc_max+1,1)/Max_value > 0.6

%For peaks that meet conditions, begin alinging of peak onset at T=0 
%Divide all Intensity values for fusion event (Fusion_event(1:end)) by the background average
%(Bck_avg_FE)
%Search for point where 

Fusion_Int_Bckgrd_ratio= Fusion_event(1:Loc_max)./Bck_avg_FE;
Max_Int_ratio=Fusion_event(1:Loc_max)./Max_value;


%Find te point were fusion begins(Fusion_onset_Point)

Fusion_onset_Point = find(Fusion_Int_Bckgrd_ratio >= (Max_Bckgrd_ratio-0.5) & Max_Int_ratio >= 0.5 ,1,'first');

%Align the Intensity data with fusion onset point at T=0, T=0 is located at
%row 200
%If fusion onset point is located in a row < 200, data will be pushed
%forward so onset point is at T=0 and empty cells will be filled with NaNs

    
if Fusion_onset_Point < 200
        %Figure out how many data points data needs to be shifted by so onset of
        %fusion is at T=0
        % Solve for x, were x = the number of data points that need the data needs
        % to be pushed forward so Fusion onset  is at row 200 (T=0s).
        %Find x
        syms x
        x_shift =  solve(x + Fusion_onset_Point == 200, x);
        %x_shift = x, the number of points data needs to be pushed forward
        x_shift = double(x_shift);
        %convert x_shift to double
        %create an array of dimensions [S 1], filled with NaNs 
        array_nans = NaN([x_shift 1],'double');
        %Create an array containing fusion event intensity data
        %Concatenate the array containing NaNs to the beginning of the array containing
        %fusion event intensity data
        array_Fusionevent = Fusion_event;
        Concatenate_arrays= [array_nans;array_Fusionevent];
        %add concatenated array into output file
        out_Aligned_PksC(:,c) = Concatenate_arrays(1:nr);

        %Convert row numbers to Time points, using exposure value provided by user
        % convert exposure value from milliseconds to seconds, exposure value =
        % time interval betwen each frame recorded
        expo = exposure_n/1000;
        %total duration of data collection 
        total_time = nr*expo; 
        start_time=-9.95;
        end_time=total_time+(start_time-expo);
        %create array with time interval from expo
        time_array = (start_time:expo:end_time); time_array_h = time_array';
        out_Aligned_PksC(:,1) = time_array_h(:,1);

        % For data with two channels i.e simultaneous analysis of two diff. fluo. 
        % molecules,the second fluorescent molecule is referred to as protein in this data
        % analysis
        %And we need to shift the protein data from the second channel so it
        % matches with our aligned data
        if nargin > 3

        Protein_event = Int_DiffP(:,c);
        %FoverSP= Protein event intensity data divided by Average Background
        %Value(CFbAvg) for the data 
        FoverSP=Protein_event/CFbAvg;
        %Create an array containing protein event intensity data
        array_proteinevent = Protein_event;
        %Concatenate the array containing NaNs to the beginning of the array containing
        %protein event intensity data
        Concatenate_p_array= [array_nans;array_proteinevent];
        out_Aligned_PksC_P(:,c) = Concatenate_p_array(1:nr);
        Concatenate_FoverS_array= [array_nans;FoverSP];
        outFoverS(:,c) = Concatenate_FoverS_array(1:nr);

        %Convert row numbers to Time points, using exposure value provided by user
        % convert exposure value from milliseconds to seconds, exposure value =
        % time interval betwen each frame recorded
        expo = exposure_n/1000;
        %total duration of data collection 
        total_time = nr*expo; 
        start_time=-9.95;
        end_time=total_time+(start_time-expo);
        %create array with time interval from expo
        time_array = (start_time:expo:end_time); time_array_h = time_array';

        out_Aligned_PksC_P(:,1) = time_array_h(:,1); 
        outFoverS(:,1) = time_array_h(:,1);
        end
        
        if nargin > 4

        Int_array_event = Int_2Array(:,c);
        Concatenate_I_array= [array_nans;Int_array_event];
        out_Aligned_PksC_I(:,c) = Concatenate_I_array(1:nr);
    
       %Convert row numbers to Time points, using exposure value provided by user
        % convert exposure value from milliseconds to seconds, exposure value =
        % time interval betwen each frame recorded
        expo = exposure_n/1000;
        %total duration of data collection 
        total_time = nr*expo; 
        start_time=-9.95;
        end_time=total_time+(start_time-expo);
        %create array with time interval from expo
        time_array = (start_time:expo:end_time); time_array_h = time_array';
        out_Aligned_PksC_I(:,1) = time_array_h(:,1); 

        end

elseif Fusion_onset_Point > 200
    
        %Figure out how many data points data needs to be shifted by so onset of
        %fusion is at T=0
        % Solve for x_shift, were x = the number of data points that need the data needs
        % to be pushed backward so Fusion onset  is at row 200 (T=0s).
        %Find x
        x_shift = Fusion_onset_Point - 200 ;
        %x_shift = x, the number of points data needs to be pulled backward
        %create an array of dimensions [x_shift 1], filled with NaNs 
        array_nans = NaN([x_shift 1]);
        shifted_fusionevent = Fusion_event(x_shift+1:end);
        Concatenate_array = [shifted_fusionevent;array_nans] ; 
        out_Aligned_PksC(:,c) = Concatenate_array(1:nr);

        %Convert row numbers to Time points, using exposure value provided by user
        % convert exposure value from milliseconds to seconds, exposure value =
        % time interval betwen each frame recorded
        expo = exposure_n/1000;
        %total duration of data collection 
        total_time = nr*expo; 
        start_time=-9.95;
        end_time=total_time+(start_time-expo);
        %create array with time interval from expo
        time_array = (start_time:expo:end_time); time_array_h = time_array';

        out_Aligned_PksC(:,1) = time_array_h(:,1);
        % For data with two channels i.e simultaneous analysis of two diff. fluo. 
        % molecules,the second fluorescent molecule is referred to as protein in this data
        % analysis
        %And we need to shift the protein data from the second channel so it
        % matches with our aligned data

        if nargin > 3

        Protein_event = Int_DiffP(:,c);
        sifted_proteinevent = Protein_event(x_shift+1:end);
        Concatenate_p_array= [sifted_proteinevent;array_nans];
        out_Aligned_PksC_P(:,c) = Concatenate_p_array(1:nr);
        FoverSP=Protein_event/CFbAvg;
        shifted_FoverSP = FoverSP(x_shift+1:end);
        Concatenate_FoverSP = [shifted_FoverSP;array_nans];
        outFoverS(:,c) = Concatenate_FoverSP(1:nr); 

       %Convert row numbers to Time points, using exposure value provided by user
        % convert exposure value from milliseconds to seconds, exposure value =
        % time interval betwen each frame recorded
        expo = exposure_n/1000;
        %total duration of data collection 
        total_time = nr*expo; 
        start_time=-9.95;
        end_time=total_time+(start_time-expo);
        %create array with time interval from expo
        time_array = (start_time:expo:end_time); time_array_h = time_array';
        out_Aligned_PksC_P(:,1) = time_array_h(:,1); 
        outFoverS(:,1) = time_array_h(:,1);

        end
%Output of ALigned Int_Array to check if any peaks were misses
        if nargin > 4

        Int_array_event = Int_2Array(:,c);
        sifted_intensityevent = Int_array_event(x_shift+1:end);
        Concatenate_I_array= [sifted_intensityevent;array_nans];
        out_Aligned_PksC_I(:,c) = Concatenate_I_array(1:nr);
    
       %Convert row numbers to Time points, using exposure value provided by user
        % convert exposure value from milliseconds to seconds, exposure value =
        % time interval betwen each frame recorded
        expo = exposure_n/1000;
        %total duration of data collection 
        total_time = nr*expo; 
        start_time=-9.95;
        end_time=total_time+(start_time-expo);
        %create array with time interval from expo
        time_array = (start_time:expo:end_time); time_array_h = time_array';
        out_Aligned_PksC_I(:,1) = time_array_h(:,1); 

        end
end
end
end
end


%CURVE FITTING
function [fitparam] = decayfit2(file)
% 
% PURPOSE:  to obtain decay fit parameters for fusion event peak onwards
% INPUT:
% file: fusion event intensity data peak onwards
% OUTPUT:fit parameters

    ncol = length(file(1,:));
    fitparam = zeros(ncol,5);
    time_x = file(:,1);

    ft = fittype( 'p+a*exp(-b*x)+c*exp(-d*x)', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.StartPoint = [1 1 1 0 1];
    opts.Lower=[0 0 0 0 0]; opts.Upper=[100 100 100 100 100];
%loop through all decay curves
    for c = 2 : ncol
    decay_curve = file(:,c);  
%isolate fusion event decay curve to onwards
%find location of peak (I)
    [~,I] = max(decay_curve(:,1));
% Time-aligned fusion event traces have Nans either before or after the
% onset of fuiosn event. Nans help adjust the onset of fuion event to t=0
%Check if Nan's were added before or after fusion event in order to exclude
%Nans from decay curve otherwise fit does not work
        y1 = isnan(decay_curve(1,1));
%if the first value of fusion event is a Nan we can fit the entire decay
%curve from peak to end
        if y1 == 1
        [~,I] = max(decay_curve(:,1));
        decaycurve = decay_curve(I:end,:);
% length of x (time) and decay curve need to be equal
        length_y= length(decaycurve(:,1));
        time_x2 = time_x(200:length_y+199,1);
        xandy=horzcat(time_x2,decaycurve);
        x2=xandy(:,1);y2=xandy(:,2);
        fp = fit( x2, y2, ft, opts );
        fitparam(c,:)=coeffvalues(fp);
%if the first value of fusion event is not Nan we need to find the fist Nan
%value after fusion onset (peak) so Nans can be excluded from the decay
%curve fit
        elseif y1==0
%find first Nan value post peak instensity
        fnd_nan_val=find(isnan(decay_curve),1,'first');
%lastpoint
        last_point_dcurve=fnd_nan_val-1;
        decaycurve = decay_curve(I:last_point_dcurve,:);
% length of x (time) and decay curve need to be equal
        length_y= length(decaycurve(:,1));
        time_x2 = time_x(1:length_y,1);
        xandy=horzcat(time_x2,decaycurve);
        x2=xandy(:,1);y2=xandy(:,2);
% f = fit(x2,y2,'Exp2'); 
        fp = fit( x2, y2, ft, opts );
        fitparam(c,:)=coeffvalues(fp);
        end
        
    end
% name each column and write out fit parameters for each decay curve    
fitparam=fitparam(2:end,:);
curr_directory =  pwd;current_folder_name = curr_directory(end-3:end);
C = cell(ncol,1);
for n2 = 2: ncol
filnam = sprintf('%s_%d',current_folder_name,n2);
C(n2,:) = {filnam};
end
rowNames = C(2:end,:);
fitparam = array2table(fitparam,'RowNames',rowNames);
        curr_directory =  pwd;current_folder_name = curr_directory(end-3:end);
       
output_excel_file = sprintf('Analyzed Data_%s.xls',current_folder_name);
writetable(fitparam, output_excel_file,'Sheet','Decayfit Data','WriteRowNames',true);
end
    

%COMPILE DATA
function [out_Aligned_PksC2, out2_Aligned_PksC_P2,out_slope,out_percentloss,out2Aligned_PksC_I]= compile_data(Aligned_Pks,Aligned_Pks_P,FoverS,Aligned_PksC_I)
% 
% PURPOSE: 
% 1) To refine data by excluding non-real curves. Non real curves
%    have all rows ==0
% 2) Calculate Slope, and % Loss in intensity at Time points = 1s and 5s
%     for each decay curve
% When calling function pk_align(Int_Array,exposure_n,Int_DiffP,CFbAvg), each INPUT
% file is as follows:
%       1) Aligned_Pks: An array containing Average Intensity values with onset of fusion event
%          aligned to T=0, each  fusion event as one column 
%       2) Aligned_Pks_P = An array containing Average Intensity values for second color channel (if any) with onset of fusion event
%          aligned to T=0, each fusion event as one column 
%       3) FoverS:  An array containing BAckground corrected Average Intensity values for second color channel (if any) with onset of fusion event
%          aligned to T=0, each fusion event as one column 
% OUTPUT includes:
%    1) out_Aligned_PksC = fusion event intensity data, with non real
%    curves (=0) excluded.
%    2) out2_Aligned_PksC_P= for two channel data, this output as intensity data from the second cannel, with non real
%    curves (=0) excluded.
%    3)  outFoverS = F over S is divides Fluorescence intesnity values by the Cell background 
%    value (CFbAvg-if provided when calling the funcion)  
% file: Aligned peaks intensity data
% OUTPUT:  Aligned peaks intensity data excluding non real peaks


%   
%Determine number of columns (ncol) and number of rows (nr) of the input
%file
ncol = length(Aligned_Pks(1,:));
nr = length(Aligned_Pks(:,1));
%
% Create empty arrays to fill with output data
out_Aligned_PksC = Aligned_Pks(:,1);
out_slope = Aligned_Pks(1,1);
out_percentloss = Aligned_Pks(1:2,1);
outAligned_PksC_I=Aligned_PksC_I(:,1);
%for two cannel data, nargin >1
%if nargin>1 create empty arrays to fill second channel data

if nargin > 1
out2_Aligned_PksC_P = Aligned_Pks_P(:,1);
outFoverS = FoverS(:,1);
end


%Loop through all fusion events, first column has time, so loop begins at
%column (c) =2
for c = 2 : ncol
% Sum of column with fusion event data >0 if column contains data
Sum_Fusionevent(:,c) = sum(Aligned_Pks(:,c),'omitnan');
%if Sum_Fusionevent >0, fusion event is output and further analyzed via calculation of
%slope and % loss
    if Sum_Fusionevent (:,c) > 0 
    out_Aligned_PksC = [out_Aligned_PksC Aligned_Pks(:,c)];
    outAligned_PksC_I=[outAligned_PksC_I Aligned_PksC_I(:,c)];
    %Find Max and location of max to calculate slope and % loss
    [Max_value,Max_location] = max(Aligned_Pks(:,c));

        %Calculate slope 

        if Max_location+20 <nr
            avg_time= mean(Aligned_Pks(Max_location:Max_location+20,1));
            avg_decaycurve = mean(Aligned_Pks(Max_location:Max_location+20,c));
            for i =  Max_location:Max_location+20 
                slope_decaycurve(i,1) = ((Aligned_Pks(i,1) - avg_time)*(Aligned_Pks(i,c) - avg_decaycurve));
                slope_decaycurve(i,2) = ((Aligned_Pks(i,1) - avg_time)^2);
            end 
            slope_decaycurve = slope_decaycurve(Max_location:Max_location+20,:);
            slope_decaycurve_y = sum(slope_decaycurve(:,1));
            slope_decaycurve_x=sum(slope_decaycurve(:,2));slope_rise_overrun = slope_decaycurve_y/slope_decaycurve_x;
            slope_final = (slope_rise_overrun/Max_value);
            %Calculate % intensity lost after 1s data points, baseline (first 150 points) substracted
            intensity_initial= Max_value;
            intensity_final= (Aligned_Pks(Max_location+20,c));
            baseline= mean(Aligned_Pks(1:150,c),'omitnan');
            percent_loss_1s= ((intensity_initial-intensity_final)/(intensity_initial-baseline))*100;
        else
            %if the the maximum value is not followed by 20 data points, we cannot
            %calculate the slope and an output value of 0 is generated for the slope
            %and for % loss
            slope_final=0;
            percent_loss_1s=0;
        end

        %if the the maximum value is followed by 100 data points, we can
        %calculate the % intensity loss after 5s of peak value
        if Max_location+100 < nr
            % %Calculate % intensity lost after 5s data points, baseline substracted
            intensity_initial= Max_value;
            intensity_final= (Aligned_Pks(Max_location+100,c));
            baseline= mean(Aligned_Pks(1:150,c),'omitnan');
            percent_loss_5s= ((intensity_initial-intensity_final)/(intensity_initial-baseline))*100;
        else 
            percent_loss_5s=0; 
        end
compiled_percent_loss=[percent_loss_1s;percent_loss_5s];
out_percentloss =[out_percentloss compiled_percent_loss];
out_slope=[out_slope slope_final];

%for two channel data, compile second channel data
if nargin > 1
out2_Aligned_PksC_P = [out2_Aligned_PksC_P Aligned_Pks_P(:,c)];
outFoverS = [outFoverS FoverS(:,c)];
end
end
end
%Folder name is the same as data name (#), so file names are assigned using the folder name 
curr_directory =  pwd; current_folder_name = curr_directory(end-3:end);
C = cell(1, ncol);
for n = 1: ncol
filnam = sprintf('%s_%d',current_folder_name,n);
C(:,n) = {filnam};
end

% Data is output in an excel file, and the data is converted to excel
% compatible format
    out_Aligned_PksC2=out_Aligned_PksC;
    
    out_Aligned_PksC = mat2dataset(out_Aligned_PksC);
    out_Aligned_PksC = dataset2table(out_Aligned_PksC); 
    out_slope = mat2dataset(out_slope);
    out_slope = dataset2table(out_slope); 
    out_percentloss = mat2dataset(out_percentloss);
    out_percentloss = dataset2table(out_percentloss); 
    out2Aligned_PksC_I=out_Aligned_PksC_I;
    out_Aligned_PksC_I = mat2dataset(out_Aligned_PksC_I);
    out_Aligned_PksC_I = dataset2table(out_Aligned_PksC_I); 
% Data is output in an excel file, and the data is converted to excel
% compatible format for second channel
if nargin > 1
    out2_Aligned_PksC_P2=out2_Aligned_PksC_P;
    out2_Aligned_PksC_P = mat2dataset(out2_Aligned_PksC_P);
    out2_Aligned_PksC_P = dataset2table(out2_Aligned_PksC_P);
    outFoverS = mat2dataset(outFoverS);
    outFoverS = dataset2table(outFoverS);
end    
 
%Folder name is the same as data name (#), so column names are assigned to both 
% te fusion event and the corresponding slope/percentloss data using the folder name 
ll = size(out_Aligned_PksC(1,:));
ll = ll(:,2);
for l = 1 : ll
 out_Aligned_PksC.Properties.VariableNames{l} = C{l};
 out2_Aligned_PksC_P.Properties.VariableNames{l} = C{l}; 
 outFoverS.Properties.VariableNames{l} = C{l}; 
 out_slope.Properties.VariableNames{l} = C{l};
 out_percentloss.Properties.VariableNames{l} = C{l};
end 

%All data is written out in excel that is named according to the folder
%name, Fusion, Slope, Protein loss etc are all outout to different tabs in
%excel with labels

output_excel_file = sprintf('Analyzed Data_%s.xls',current_folder_name);
writetable(out_Aligned_PksC, output_excel_file,'Sheet','Fusion');
writetable(out_slope, output_excel_file,'Sheet','Slope');
writetable(out_percentloss, output_excel_file,'Sheet','% Loss');
writetable(out_Aligned_PksC_I, output_excel_file,'Sheet','Fusion_RawData');
if nargin > 1
writetable(out2_Aligned_PksC_P, output_excel_file,'Sheet','Protein');
writetable(outFoverS, output_excel_file,'Sheet','Protein_FoverS');
end

end

%NORMALIZE DATA
function[out_Norm_Data_Fusion,out_Norm_Data_Protein] = normalizedata2(Aligned_PksC,Aligned_PksC_P)

%PURPOSE: To normalize fusion and protein data using minimum value from
%data points before te onset of fusion event. This is to avoid changing the
%shape and kinetics of the decay curve

% When calling function normalizedata2(Aligned_PksC,fiAligned_PksC_Ple2), each INPUT
% file is as follows:
%       1) Aligned_Pks: An array containing Average Intensity values with onset of fusion event
%          aligned to T=0, each  fusion event as one column 
%       2) Aligned_Pks_P = An array containing Average Intensity values for second color channel (if any) with onset of fusion event
%          aligned to T=0, each fusion event as one column 
% OUTPUT includes:
%    1) Norm_Data_Fusion = Normalized fusion event intensity data, with onset of fusion aligned to T=0s 
%    2) Norm_Data_Protein= for two channel data, this output normalizes the
%    intensity from the second channel 
     
%   
%Determine number of columns (ncol) and number of rows (nr) of the input
%file

nframes = length(Aligned_PksC(:,1));
ncols = length(Aligned_PksC(1,:));
%Determine number of columns (ncol) and number of rows (nr) of the input
% in multi cannel data
if nargin>1 
ncols2 = length(Aligned_PksC_P(1,:));
end
% Create empty arrays to fill with output data
Norm_Data_Fusion = zeros(nframes,ncols);
Norm_Data_Protein = zeros(nframes,ncols);

%Loop through all fusion events, first column has time, so loop begins at
%column (c) =2
for c = 2:ncols

%Find the Max value of the fusion event i.e. peak
%Find a minimum value from frames before fusion event
max_data_fusion = max(Aligned_PksC(:,c));
min_data_fusion = min(Aligned_PksC(1:199,c));

%Loop through each data point and normalize using the min and max values
%obtained above nframes=number of rows
    for i = 1:nframes
    i_intensity_value = Aligned_PksC(i,c);
    Norm_Data_Fusion(i,c) = ((i_intensity_value - min_data_fusion)/(max_data_fusion-min_data_fusion));
%Copy the Time column from input file
    Norm_Data_Fusion(:,1) = Aligned_PksC(:,1);

%for two cannel data, nargin >1
%if nargin>1 cnormalize data from second channel 
if nargin > 1 
%Loop through all columns of second data set
for k = 1:ncols2   

%Find the min and max values for each column of second data set

max_data_protein = max(Aligned_PksC_P(:,k));
min_data_protein = min(Aligned_PksC_P(:,k));
    
k_intensity_value = Aligned_PksC_P(i,k);
Norm_Data_Protein(i,k) = ((k_intensity_value - min_data_protein)/(max_data_protein-min_data_protein));

Norm_Data_Protein(:,1) = Aligned_PksC_P(:,1);

end
end
    end

end

% Data is output in an excel file, and the data is converted to excel
% compatible format
 out_Norm_Data_Fusion=Norm_Data_Fusion;
 Norm_Data_Fusion = mat2dataset(Norm_Data_Fusion);
 Norm_Data_Fusion = dataset2table(Norm_Data_Fusion); 
% Data is output in an excel file, and the data is converted to excel
% compatible format for second channel
 if nargin > 1
out_Norm_Data_Protein = Norm_Data_Protein;
 Norm_Data_Protein = mat2dataset(Norm_Data_Protein);
 Norm_Data_Protein = dataset2table(Norm_Data_Protein); 
 end
%Folder name is the same as data name (#), so file (column) names are assigned using the folder name 
curr_directory =  pwd;current_folder_name = curr_directory(end-3:end);
C = cell(1, ncols);
for n2 = 1: ncols
filnam = sprintf('%s_%d',current_folder_name,n2);
C(:,n2) = {filnam};
end

%Column names for second data set 
ll = size(Norm_Data_Fusion(1,:));
ll = ll(:,2);
for l = 1 : ll
 Norm_Data_Fusion.Properties.VariableNames{l} = C{l};
 Norm_Data_Protein.Properties.VariableNames{l} = C{l};
end 

%All data is written out in excel that is named according to the folder
%name, and data sets are put in tabs wit their correspomnding names

output_excel_file = sprintf('Analyzed Data_%s.xls',current_folder_name);
writetable(Norm_Data_Fusion,output_excel_file,'Sheet','Normalized Fusion');

if nargin > 1
writetable(Norm_Data_Protein,output_excel_file,'Sheet','Normalized Protein');  

end
end


%ALIGN MINI MOVIES
function [sh_movie1,sh_movie2]= align_movies(Aligned_Pks,movie1,movie2)
% PURPOSE: 
% 1) To create mini movies with fusion events in the middle of the mini
% movie e.g for a 200 frame lengt mini move the onset of fusion event is at
% approximately frame 100
% 2) Create montages of input movies
% When calling functionalign_movies(Aligned_Pks,movie1,movie2), each INPUT
% file is as follows:
%       1) Aligned_Pks: An array containing Average Intensity values with onset of fusion event
%          aligned to T=0, each  fusion event as one column. This file is
%          used to create minis with the onset of fusion event centered. 
%       2) movie1 = A 25x25 pixels 4D array containing all fusion events.
%       3) movie2 = A 25x25 pixels 4D array containing all second channel frames (if any) accompanying the fusion event.         
% OUTPUT includes:
%    1) sh_movie1 =  A 25x25 pixels 4D array containing all fusion event with fusion
%    event aligned to the center of each mini for the respective fusion
%    event.
%    2) sh_movie2=  A 25x25 pixels 4D array containing all second channel frames (if any) accompanying the fusion event. 
%   Each event is aligned to the center of each mini for the respective fusion
%   event.
% 
%
%
%Determine number of columns (ncol) and number of rows (nr) of the input
%file 
ncol = length(Aligned_Pks(1,:));
[nr,~] = size(Aligned_Pks);
%
%figure out how big the output arrays will be
n_movies = ncol-(length(find(Aligned_Pks(1,:) == 0))+1);

%find all the time shifts, note that the first column
%contains time info (just 0s).
idx = 1;
%find the peak onset time (T=0s) that the peaks are aligned to in the array
%if data includes no fusion event Aligned_Pks(:,1)==0
if Aligned_Pks(:,1)==0
    sh_movie1=[];
    disp('no events to make minis');
    return;
end

%Loop through all fusion events, first column has time, so loop begins at
%column (c) =2
for c = 2 : ncol
%isolate each fusion event and check if the column contains a real event
%i.e sum of column >0
    Fusion_Event_x(:,c)= Aligned_Pks(:,c);
    sum_Fusion_Event_x(:,c) = sum(Aligned_Pks(:,c),'omitnan');
    if sum_Fusion_Event_x(:,c) > 0 
%Determine length of output mini for the fusion event
%Determine if fusion event has 100 frames before the onset of the event,
%for events with less than 200 frames before onset, the first row is nan
        real_fusion_event = Fusion_Event_x(1,c);
if isnan(real_fusion_event) 
%Determine the last nan location to determine how many frames of data are
%present before the fusion onset 
            last_nan = find(isnan(Fusion_Event_x(:,c)),1,'last');
            if last_nan < 100 
%Determine location of the first frame for output mini for the fusion event
            start_frame = 200-last_nan;
%Create empty arrays to store the mini, and to create the mini, these files
%are output 
            sh_movie1 = uint16(zeros(25,25,200,idx));
            shmovie1 = uint16(zeros(25,25,200));
%write mini into the output file
            sh_movie1(:,:,:,idx) = movie1(:,:,start_frame-100:start_frame+99,c);
            shmovie1(:,:,:) = sh_movie1(:,:,:,idx);
if nargin >2
%Create empty arrays to store the second channel mini, and to create the mini, these files
%are output            
            sh_movie2 = uint16(zeros(25,25,200,idx));
            shmovie2 = uint16(zeros(25,25,200));
             sh_movie2(:,:,:,idx) = movie2(:,:,start_frame-100:start_frame+99,c);
            shmovie2(:,:,:) = sh_movie2(:,:,:,idx);
end
%loop through each frame to write out mini as a tif file
%each mini is named according to the column number of the fusion event
%for 2 channel data the second mini is also created

                    for frame_num = 1:200
                    shmovie1(:,:,frame_num) = shmovie1(:,:,frame_num);
                    mini_fusion= sprintf('Spot%dmini.tif',c);
                    imwrite(shmovie1(:,:,frame_num),mini_fusion,'WriteMode','append','compression','none');
                    if nargin>2
                    shmovie2(:,:,frame_num) = shmovie2(:,:,frame_num);
                    mini_protein = sprintf('Spot%dminip.tif',c);
                    imwrite(shmovie2(:,:,frame_num),mini_protein,'WriteMode','append','compression','none');
                    end
                    end
%Create montages for the fusion event and the accompanying second channel
%data, if any
%Write montages out as tif 
               
                    shmovie1_montage=(shmovie1(:,:,1:20:end));
                    if nargin>2
                    shmovie2_montage=(shmovie2(:,:,1:20:end));
                    end
                     mont_size=size(shmovie1_montage,3);                 
                     for num_frame = 1:mont_size
                     shmovie1_montage(:,:,num_frame) =  shmovie1_montage(:,:,num_frame);
                     mininmontage = sprintf('Spot%dminimontage.tif',c);
                     imwrite(shmovie1_montage(:,:,num_frame),mininmontage,'WriteMode','append','compression','none');
                     if nargin>2
                     shmovie2_montage(:,:,num_frame) =  shmovie2_montage(:,:,num_frame);
                     mininpmontage = sprintf('Spot%dminipmontage.tif',c);
                     imwrite (shmovie2_montage(:,:,num_frame),mininpmontage,'WriteMode','append','compression','none');
                     end
                     end
            end
%Determine length of minis when fusion event has less than 100 frame length
%data before onset of fusion (T=0s,row=200)
            if last_nan > 100 && last_nan < 200
            frames_before_fusion = 200-last_nan;
            strt = frames_before_fusion -1;
            pkstrt = frames_before_fusion-strt;
            mini_length = frames_before_fusion*2;
%Create empty arrays to store the mini, and to create the mini, these files
%are output 
            sh_movie1 = uint16(zeros(25,25,mini_length,idx));    
            shmovie1 = uint16(zeros(25,25,mini_length));      
            sh_movie1(:,:,:,idx) = movie1(:,:,pkstrt:mini_length,c);
            shmovie1(:,:,:) = sh_movie1(:,:,:,idx);  
            if nargin >2
            sh_movie2 = uint16(zeros(25,25,mini_length,idx));    
            shmovie2 = uint16(zeros(25,25,mini_length));      
            sh_movie2(:,:,:,idx) = movie2(:,:,pkstrt:mini_length,c);
            shmovie2(:,:,:) = sh_movie2(:,:,:,idx);
            end
%loop through each frame to write out mini as a tif file
%each mini is named according to the column number of the fusion event
%for 2 channel data the second mini is also created

                    for frame_num = 1:mini_length
                    shmovie1(:,:,frame_num) = shmovie1(:,:,frame_num);
                    mini_fusion= sprintf('Spot%dmini.tif',c);
                    imwrite(shmovie1(:,:,frame_num),mini_fusion,'WriteMode','append','compression','none');
                    if nargin>2
                    shmovie2(:,:,frame_num) = shmovie2(:,:,frame_num);
                    mini_protein = sprintf('Spot%dminip.tif',c);
                    imwrite(shmovie2(:,:,frame_num),mini_protein,'WriteMode','append','compression','none');
                    end
                    end

%Create montages for the fusion event and the accompanying second channel
%data, if any
%Write montages out as tif 
                    
                     shmovie1_montage=(shmovie1(:,:,1:20:end));
                     if nargin>2
                     shmovie2_montage=(shmovie2(:,:,1:20:end));
                     end
                     mont_size=size(shmovie1_montage,3);
                     for num_frame = 1:mont_size
                     shmovie1_montage(:,:,num_frame) =  shmovie1_montage(:,:,num_frame);
                     mininmontage = sprintf('Spot%dminimontage.tif',c);                    
                     imwrite(shmovie1_montage(:,:,num_frame),mininmontage,'WriteMode','append','compression','none');
                     
                     if nargin>2
                     shmovie2_montage(:,:,num_frame) =  shmovie2_montage(:,:,num_frame);
                     mininpmontage = sprintf('Spot%dminipmontage.tif',c);
                     imwrite (shmovie2_montage(:,:,num_frame),mininpmontage,'WriteMode','append','compression','none');
                     end
                     end
            end
%Determine length of output mini for the fusion event
%Determine how many frames of data are present post onset of the event,
%for events with greater than 200 frames before onset, the first row is not nan
 elseif ~isnan(real_fusion_event)  
             first_nan = find(isnan(Fusion_Event_x(:,c)),1,'first');
             total_nan = nr-first_nan+201 ;
             if total_nan < nr-99
             lengthvid = 200;
%Create empty arrays to store the mini, and to create the mini, these files
%are output 
             sh_movie1 = uint16(zeros(25,25,200,idx));
             shmovie1 = uint16(zeros(25,25,200));
             sh_movie1(:,:,:,idx) = movie1(:,:,total_nan-100:total_nan+99,c);
             shmovie1(:,:,:) = sh_movie1(:,:,:,idx);
              if nargin>2
             sh_movie2 = uint16(zeros(25,25,200,idx));
             shmovie2 = uint16(zeros(25,25,200));
             sh_movie2(:,:,:,idx) = movie2(:,:,total_nan-100:total_nan+99,c);
             shmovie2(:,:,:) = sh_movie2(:,:,:,idx);
              end
%loop through each frame to write out mini as a tif file
%each mini is named according to the column number of the fusion event
%for 2 channel data the second mini is also created                  
                for frame_num = 1:lengthvid
                shmovie1(:,:,frame_num) = shmovie1(:,:,frame_num);
                mini_fusion= sprintf('Spot%dmini.tif',c);
                imwrite(shmovie1(:,:,frame_num),mini_fusion,'WriteMode','append','compression','none');  
                if nargin>2
                shmovie2(:,:,frame_num) = shmovie2(:,:,frame_num);   
                mini_protein = sprintf('Spot%dminip.tif',c);
                imwrite(shmovie2(:,:,frame_num),mini_protein,'WriteMode','append','compression','none');
                end
                end
 %Create montages for the fusion event and the accompanying second channel
%data, if any
%Write montages out as tif  
             shmovie1_montage=(shmovie1(:,:,1:20:end));
                    shmovie2_montage=(shmovie2(:,:,1:20:end));
                     mont_size=size(shmovie1_montage,3);
                     for num_frame = 1:mont_size
                     shmovie1_montage(:,:,num_frame) =  shmovie1_montage(:,:,num_frame);
                     shmovie2_montage(:,:,num_frame) =  shmovie2_montage(:,:,num_frame);
                     mininmontage = sprintf('Spot%dminimontage.tif',c);
                     mininpmontage = sprintf('Spot%dminipmontage.tif',c);
             
                    imwrite(shmovie1_montage(:,:,num_frame),mininmontage,'WriteMode','append','compression','none');
                     imwrite (shmovie2_montage(:,:,num_frame),mininpmontage,'WriteMode','append','compression','none');
                    end
             end
%Determine length of minis for fusion events with less than 100 frames post
%fusion
                if total_nan> nr-99        
                mini_length = nr-total_nan;
                start_frame = total_nan-mini_length;
                lengthvid= nr-start_frame;
%Create empty arrays to store the mini, and to create the mini, these files
%are output 
            sh_movie1 = uint16(zeros(25,25,lengthvid,idx));
            shmovie1 = uint16(zeros(25,25,lengthvid));
            sh_movie1(:,:,:,idx) = movie1(:,:,start_frame:total_nan+mini_length-1,c);
            shmovie1(:,:,:) = sh_movie1(:,:,:,idx);
            if nargin>2
            sh_movie2 = uint16(zeros(25,25,lengthvid,idx));
            shmovie2 = uint16(zeros(25,25,lengthvid));
            sh_movie2(:,:,:,idx) = movie2(:,:,start_frame:total_nan+mini_length-1,c);
            shmovie2(:,:,:) = sh_movie2(:,:,:,idx);
            end
%loop through each frame to write out mini as a tif file
%each mini is named according to the column number of the fusion event
%for 2 channel data the second mini is also created                
            for frame_num = 1:lengthvid
            shmovie1(:,:,frame_num) = shmovie1(:,:,frame_num);
            mini_fusion= sprintf('Spot%dmini.tif',c);
            imwrite(shmovie1(:,:,frame_num),mini_fusion,'WriteMode','append','compression','none'); 
            if nargin>2
            shmovie2(:,:,frame_num) = shmovie2(:,:,frame_num);
            mini_protein = sprintf('Spot%dminip.tif',c);
            imwrite(shmovie2(:,:,frame_num),mini_protein,'WriteMode','append','compression','none');
            end
            end   
%Create montages for the fusion event and the accompanying second channel
%data, if any
%Write montages out as tif 
                    shmovie1_montage=(shmovie1(:,:,1:20:end));
                    shmovie2_montage=(shmovie2(:,:,1:20:end));
                    mont_size=size(shmovie1_montage,3);
                    for num_frame = 1:mont_size
                    shmovie1_montage(:,:,num_frame) =  shmovie1_montage(:,:,num_frame);
                    mininmontage = sprintf('Spot%dminimontage.tif',c);
                    imwrite(shmovie1_montage(:,:,num_frame),mininmontage,'WriteMode','append','compression','none');
                    if nargin>2
                    shmovie2_montage(:,:,num_frame) =  shmovie2_montage(:,:,num_frame);
                    mininpmontage = sprintf('Spot%dminipmontage.tif',c);
                    imwrite (shmovie2_montage(:,:,num_frame),mininpmontage,'WriteMode','append','compression','none');
                    end
                    end
                    end
end    
    end
   end
end

%CREATE_MONTAGE
function [sh_movie1,sh_movie2]= create_montage(Aligned_Pks,movie1,movie2)
% PURPOSE: 
% 1) To create mini movies with fusion events in the middle of the mini
% movie e.g for a 200 frame lengt mini move the onset of fusion event is at
% approximately frame 100
% 2) Create montages of input movies
% When calling functionalign_movies(Aligned_Pks,movie1,movie2), each INPUT
% file is as follows:
%       1) Aligned_Pks: An array containing Average Intensity values with onset of fusion event
%          aligned to T=0, each  fusion event as one column. This file is
%          used to create minis with the onset of fusion event centered. 
%       2) movie1 = A 25x25 pixels 4D array containing all fusion events.
%       3) movie2 = A 25x25 pixels 4D array containing all second channel frames (if any) accompanying the fusion event.         
% OUTPUT includes:
%    1) sh_movie1 =  A 25x25 pixels 4D array containing all fusion event with fusion
%    event aligned to the center of each mini for the respective fusion
%    event.
%    2) sh_movie2=  A 25x25 pixels 4D array containing all second channel frames (if any) accompanying the fusion event. 
%   Each event is aligned to the center of each mini for the respective fusion
%   event.
% 
%
%
%Determine number of columns (ncol) and number of rows (nr) of the input
%file 
ncol = length(Aligned_Pks(1,:));
[nr,~] = size(Aligned_Pks);
%
%figure out how big the output arrays will be
n_movies = ncol-(length(find(Aligned_Pks(1,:) == 0))+1);

%find all the time shifts, note that the first column
%contains time info (just 0s).
idx = 1;
%find the peak onset time (T=0s) that the peaks are aligned to in the array
%if data includes no fusion event Aligned_Pks(:,1)==0
if Aligned_Pks(:,1)==0
    sh_movie1=[];
    disp('no events to make minis');
    return;
end

%Loop through all fusion events, first column has time, so loop begins at
%column (c) =2
for c = 2 : ncol
%isolate each fusion event and check if the column contains a real event
%i.e sum of column >0
    Fusion_Event_x(:,c)= Aligned_Pks(:,c);
    sum_Fusion_Event_x(:,c) = sum(Aligned_Pks(:,c),'omitnan');
    if sum_Fusion_Event_x(:,c) > 0 
%Determine length of output mini for the fusion event at 20 frame interval

        length_montage= floor(nr/20);
%Convert the movie at c to 3D
mont_at_c = uint16(zeros(25,25,nr,idx));
mont_at_c(:,:,:,idx) = movie1(:,:,:,c);
mont_at_c=mont_at_c(:,:,:);
if nargin >2
mont_at_c_p = uint16(zeros(25,25,nr,idx));
mont_at_c_p(:,:,:,idx) = movie2(:,:,:,c);
mont_at_c_p=mont_at_c_p(:,:,:);
end
%Create empty arrays to store the mini, and to create the mini, these files
%are output 
            sh_movie1 = uint16(zeros(25,25,length_montage));
            shmovie1 = uint16(zeros(25,25,length_montage));
if nargin >2
            sh_movie2 = uint16(zeros(25,25,length_montage));
            shmovie2 = uint16(zeros(25,25,length_montage));
end
%write mini into the output file
        montage_interval= 1:20:nr;
     for i = 1: length_montage
       first_frame=montage_interval(1,i);
       montage_frame_mean= mean(mont_at_c(:,:,first_frame:first_frame+19),3);
       montage_frame_mean = uint16(montage_frame_mean);
       shmovie1(:,:,i) = montage_frame_mean;

if nargin >2
%Create empty arrays to store the second channel mini, and to create the mini, these files
%are output            

            montage_frame_meanp= mean(mont_at_c_p(:,:,first_frame:first_frame+19),3);
            montage_frame_meanp = uint16(montage_frame_meanp);
            shmovie2(:,:,i) = montage_frame_meanp;
end
    end
%loop through each frame to write out mini as a tif file
%each mini is named according to the column number of the fusion event
%for 2 channel data the second mini is also created

                    for frame_num = 1:length_montage
                    shmovie1(:,:,frame_num) = shmovie1(:,:,frame_num);
                    mini_fusion= sprintf('Spot%dmini_fullmontage.tif',c);
                    imwrite(shmovie1(:,:,frame_num),mini_fusion,'WriteMode','append','compression','none');
                    if nargin>2
                    shmovie2(:,:,frame_num) = shmovie2(:,:,frame_num);
                    mini_protein = sprintf('Spot%dminip_fullmontage.tif',c);
                    imwrite(shmovie2(:,:,frame_num),mini_protein,'WriteMode','append','compression','none');
                    end
                    end

   end
end
end
%RADIAL PLOTS
function [r_avg3,area_rp2,HWHMO]=radial_plot(Aligned_Pks,Int_Diff,Indiv_events,totalpks)
% 
% PURPOSE: loop though entire input movie to make calculate Half width full max from
% onset of fusion event to the end of movie using better_colo_v2(a) program included below
% 
% INPUT:
% mov: movie stack, in the format of (512,512,300). It can be smaller in X
% and Y dimensions. 
% time: 
% OUTPUT:  Half width half max in nanometers)and area under radial plot
% over time from fusion onset onwards. 
%
ncol = length(Aligned_Pks(1,:));
nr = length(Aligned_Pks(:,1));
totalpks1=length(totalpks(1,:));
tp2=(totalpks1-1);
totalpks=(totalpks1-1)*7; 
area_rp2=zeros(ncol,1000);
radial_avg=zeros(83,totalpks);
HWHMO=zeros(1,1001);
AREAO=zeros(1,1000);

%convert from pixels to nm
conversion = [0,1,1.414213562,2,2.236067977,2.828427125,3,3.16227766,3.605551275,4,4.123105626,4.242640687,4.472135955,5,5.099019514,5.385164807,5.6568542490,5.83095189500000,6,6.08276253000000,6.32455532000000,6.40312423700000,6.70820393200000,7,7.07106781200000,7.21110255100000,7.28010988900000,7.61577310600000,7.81024967600000,8,8.06225774800000,8.24621125100000,8.48528137400000,8.54400374500000,8.60232526700000,8.94427191000000,9,9.05538513800000,9.21954445700000,9.43398113200000,9.48683298100000,9.84885780200000,9.89949493700000,10,10.0498756200000,10.1980390300000,10.2956301400000,10.4403065100000,10.6301458100000,10.7703296100000,10.8166538300000,11,11.0453610200000,11.1803398900000,11.3137085000000,11.4017542500000,11.6619037900000,11.7046999100000,12,12.0415945800000,12.0830459700000,12.1655250600000,12.2065556200000,12.3693168800000,12.5299640900000,12.6491106400000,12.7279220600000,12.8062484700000,13,13.0384048100000,13.4164078600000,13.4536240500000,13.6014705100000,13.8924439900000,14.1421356200000,14.2126704000000,14.4222051000000,14.8660687500000,15,15.5563491900000,15.6204993500000,16.2788206000000,16.97056275];
x_rad=conversion';
x_rad=x_rad*102.7;

r_avg2=x_rad;
r_avg3=x_rad;
% loop through all frames of fusion event mini movie if the movie
% corrensponding to c represents a real fusion event
for c =2:ncol
HWHM=zeros(1,1001);
radial_avg=zeros(83,totalpks);
AS = sum(Aligned_Pks(:,c),'omitnan');
AAS2 = Int_Diff(:,c);
    [~,I] = max(AAS2);
if AS  < 0 
continue
end
if AS>0
   
    r_avg2=x_rad;
    AAS2 = Int_Diff(:,c);
    [~,I] = max(AAS2);
    mini = Indiv_events(:,:,:,c);
    mini=mini(:,:,:);
    mini = mini(:,:,I-1:end);
    lenmini=size(mini,3);
   
%writeout mini for radial plot
    for i=1:lenmini 
    radial_out(:,i) = better_colo_V2(mini(:,:,i));   
    end
if lenmini<20
continue
end
    for n=1:20:lenmini-20  
    radial_avg(:,n)=mean(radial_out(:,n:n+2),2);    
    r_avg(:,n)=radial_avg(:,n);
    maxn=max(r_avg(:,n));
    bckn=r_avg(83,1);
    bck=(maxn-bckn);
    hmaxya1=bck/2;
    hmaxya=maxn-hmaxya1;
    uAvgn=unique(r_avg(:,n), 'stable');
    lenuAvg=length(uAvgn);
    xn= x_rad (1:lenuAvg,1);
    HWHM(1,n)= interp1(uAvgn,xn,hmaxya, 'makima');
    r_avgx= [r_avg r_avg(:,n)];
    r_avgx= r_avgx(:,1:end);
    end
    r_avgx=r_avgx(:,1:end-1);
    HWHMO=[HWHMO;HWHM];
r_avgx2= r_avgx(:,1:20:end);
r_avg3=[r_avg3 r_avgx2];


clear mini;


 area_rp=zeros(1,1000);
 for j = 1:20:size(r_avgx,2)
s=r_avgx(:,j);
%This half of the program finds the area of the central peak. The first few
%lines are not necessary, but
ydatatemp=s(7:38,:);
ypeak=s(1:7,:);
%limits the yvalues to just the values past 3 pixels from the center of the
%image. The nanoparticle seems to only be 3 pixels in radius.
base=mean(ydatatemp);
%finds the average of all of the data past the 3rd pixel for the height of
%the rectangle
rect=base*7;
%7 is the width of the rectangle because thats how many data points there
%are under the peak
clear sum
x=sum(ypeak);
area_rp2(c,j)=x-rect;
% area_rp = area_rp2;
%this gives a number for the area of the peak from the colocalization
%as the answer
% area_rp2=[area_rp2 area_rp(:,j)]
area_rp(:,:)=area_rp2(c,:);

end

AREAO=[AREAO;area_rp];

clear radial_avg
clear radial_out
clear r_avg
clear HWHM
clear r_avgx
clear area_rp
end

end

% Compile calculated HWHM and Area data and output both as matfile and in excel 
HWHMO= HWHMO(:,1:20:end);  
AREAO= AREAO(:,1:20:end); 
sv=size(HWHMO,1);
curr_directory =  pwd;current_folder_name = curr_directory(end-3:end);
C3 = cell(sv,1);
for n2 = 1: sv
filnam = sprintf('%s_%d',current_folder_name,n2);
C3(n2,:) = {filnam};
end
rowNames3 = C3(1:end,:);
% rowNames4 = C(2:end,:);
    HWHMO = array2table(HWHMO,'RowNames',rowNames3);
sva=size(AREAO,1);
curr_directory =  pwd;current_folder_name = curr_directory(end-3:end);
C3a = cell(sva,1);
for n2 = 1: sva
filnam = sprintf('%s_%d',current_folder_name,n2);
C3a(n2,:) = {filnam};
end
rowNames3a = C3a(1:end,:);

    AREAO = array2table(AREAO,'RowNames',rowNames3a);
    output_excel_file = sprintf('Analyzed Data_%s.xls',current_folder_name);
    writetable(AREAO, output_excel_file,'Sheet','area_rp','WriteRowNames',true);
    writetable(HWHMO, output_excel_file,'Sheet','HWHM,','WriteRowNames',true);
end

function  [ radial_out ] = better_colo_V2(a)
%better_colo plots intensities of a location guided average image as a 
%function of distance from the center of the image. 
% The second half of the
%program finds the area under the peak of the radial intensity plot. It
%does this by subtracting the area of the rectangle under the peak (the 
%peak ends 3 pixels away from the center) using the mean of the data after 
%the third pixel as the height of the rectangle. 
%Written by Mitch Alton on 7/1/15
%adapted for no user input and mass throughput of images by MKK 9/27/2016

%finds size of the image
[m,n]=size(a);

%creates matrix of repeating columns of 1 through how ever many columns
%there are in the image
x=[1:n];
l=repmat(x,n,1);

%creates matrix of repeating rows
y=[1:m];
y';
p=repmat(y,m,1);
w=p';
%makes matrixes -# to +#. Not sure how this will work with even numbers
%though
x=l-(((n-1)*0.5)+1);
y=w-(((n-1)*0.5)+1);
%creates a matrix of intensities from image
z=a(1:m,1:n);

%converts three matrixes to a three dimentional polar coordinate matrix
[THETA,RHO,Z]=cart2pol(x,y,z);

%turns matrix into a list of values
P=reshape(Z,1,[]);
H=reshape(RHO,1,[]);

%this averages the intensities (Z) for single distance (RHO) values
[plt,I,J] = unique(H);
s = zeros(size(plt));
frequencies = zeros(size(plt));
for i = 1:max(J)
    I = find(J==i);
    s(i) = mean(P(I));
    frequencies(i) = length(I);
end
radial_out = s.';
end

%CELL MASK
function [cellperimeter,cellmask,bckcormsk,backgroundmask]=obtaincellmasks(I);
% 
% PURPOSE:  to create cell mask
% INPUT:
% OUTPUT:  Cellperimeter,cellmask,bckcormsk, and backgroundmask. 
% OUTPUT:  an array that is very long (number of spots found in each frame)
% that has X, Y, time 

   
    I = imread(I);
    fudgeFactor = 0.4;
    [~, threshold] = edge(I, 'sobel');
    BWs = edge(I,'sobel',threshold * fudgeFactor);
    se90 = strel('line',3,90);
    se0 = strel('line',3,0);
    BWsdil = imdilate(BWs,[se90 se0]);
    BWdfill = imfill(BWsdil,8,'holes');
    BWnobord = imclearborder(BWdfill);
    bckcormsk=~BWnobord(:,:);
    cellmask = BWdfill-BWnobord;
    backgroundmask(:,:) = ~cellmask(:,:);
    seD = strel('diamond',1);
    BWfinal1 = imerode(cellmask,seD);
    cellperimeter = bwperim(BWfinal1);
    Segout = I;
    Segout(cellperimeter) = 255;
% 
end

%REMOVE VESICLES THAT MOVE
function newpkfile = n_pkfile(pk,msd_indiv)
nr = length(pk(:,1));
ncol = length(msd_indiv(1,:));
newpkfile = zeros(nr,2);
for n=2:ncol
MSDs= sum(msd_indiv(:,n));
if MSDs == 0
newpkfile(n-1,:)=pk(n-1,:);
else
newpkfile(n-1,:)=[];
end
end
end

