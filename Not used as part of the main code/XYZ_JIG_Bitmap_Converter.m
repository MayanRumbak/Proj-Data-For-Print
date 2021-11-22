function[]=XYZ_JIG_Bitmap_Converter(PathName,BitMapStruct)
FileType.png=1;
FileType.bmp=0;

FileType.type=FileType.bmp;
h= waitbar(0,'Pleas wait');
% Z shift in mictometer per interlace
Zshift=[100 13.7 100 27];
fnames={};
FileListArray={};
HeadIndex=0; %str2num(BitMapStruct.BitMapHeadIndex);
HeadShift=0;%str2num(BitMapStruct.BitMapHeadShift);
StartRaw=1;%str2num(BitMapStruct.BitMapStartRaw);
EndRaw=192;%str2num(BitMapStruct.BitMapEndRaw);
StartFireIndex=1%str2num(BitMapStruct.BitMapStartFireIndex);
EndFireIndex=28000;%str2num(BitMapStruct.BitMapEndFireIndex);
DataMultiplier=4;% data multiplier is the numbers of fires per pixel
%fnames = dir([PathName '\*.bmp']);
numfids = length(fnames);
vals = cell(1,numfids);
% interlace parameters
a=5;
b=10;

NOfTravels=4;  %number of travels


TotalNoOfNozzles=192; % number of m\nozzles in print head
NumberOfPrintHeads= 4; % number of print heads
NumberOfNozzleLines= 2; % numbers of nozzles line per print head
F=fopen([PathName '\BitmapList.gis'],'w'); % this file contains a list of the bitmaps full path

%%%%
%  a=3;
%  b=5;
%  TotalNoOfNozzles=100;
%  NOfTravels=4;
% Scattering=100-66;
% NumberOfPrintHeads=1
% NumberOfNozzleLines=1
%%%%

travelShift=[0 -b b-a -a]; % shift in pixels to the bitmaps - interlace parameters



fnames = dir([PathName '\*.png']);
numfids = length(fnames);


p=NumberOfPrintHeads*NumberOfNozzleLines+1;


for K = 1:numfids
    %      profile on
    Scattering=round(100-50*rand); %sactering parameters 100 and 50 - there is a random shift of 50-100 pixels
    
    InIm=imread([PathName '\' sprintf('finished_Slice%d.png',K)]);
    
    RGBImage=d2b(InIm,K==1);
    T1=K-1;%sscanf(fnames(K).name,'finished_Slice%d.png');
    
    if mod(T1,2)
        RGBImage(:,:,:)=RGBImage(end:1,:,:);
    end
    p=8;
    %[path name ext]=fileparts([fnames(K).name]);
    
    for i=0:NumberOfPrintHeads-1
        for w=1:NumberOfNozzleLines
            [rows columns numberOfColorChannels] = size(RGBImage{i*NumberOfNozzleLines+w});
            newWidth = [DataMultiplier * columns];
            
            RGBImage{i*NumberOfNozzleLines+w} = imresize(RGBImage{i*NumberOfNozzleLines+w}, [rows newWidth]);
            %%%
            [Xr Yr]=size(RGBImage{i*NumberOfNozzleLines+w});
            if Scattering~=0
                RGBImage{i*NumberOfNozzleLines+w}=[zeros(Scattering,Yr);RGBImage{i*NumberOfNozzleLines+w}(:,:)];
            end
            
            %RGBImage=RGBImage(:150,50:150);
            %%%
            [Xr Yr]=size(RGBImage{i*NumberOfNozzleLines+w});
            RGBImage{i*NumberOfNozzleLines+w}(:,end)=zeros(Xr,1);
            %imshow(RGBImage);
            %EndFireIndex=28000;%str2num(BitMapStruct.BitMapEndFireIndex);
            
            if Yr<EndFireIndex
                EndFireIndex=Yr;
            end
            %
            NOPass=ceil(Xr/(TotalNoOfNozzles*NOfTravels));
            
            for PassIndex=1:NOPass
                PassShift=NOfTravels*TotalNoOfNozzles*(PassIndex-1);
                for TravelIndex=1:NOfTravels
                    
                    RawShift=PassShift+travelShift(TravelIndex);
                    StratRaw=RawShift+1;
                    EndRaw=RawShift+NOfTravels*TotalNoOfNozzles;
                    RawFilter=StratRaw:NOfTravels:EndRaw;
                    PreZero=find(RawFilter<1);
                    DeleteFilt=find(RawFilter>Xr | RawFilter<1);
                    RawFilter(DeleteFilt)=[];
                    FireFilter=StartFireIndex:EndFireIndex;
                    LenRaw=length(RawFilter);
                    LenFire=length(FireFilter);
                    pattern=zeros(TotalNoOfNozzles,LenFire);
                    LenPreZero=length(PreZero);
                    pattern((TotalNoOfNozzles*HeadIndex+1+LenPreZero):(TotalNoOfNozzles*HeadIndex+LenRaw+LenPreZero),StartFireIndex:EndFireIndex)=RGBImage{i*NumberOfNozzleLines+w}(RawFilter,FireFilter);
                    
                    % imshow(pattern);
                    pattern=logical(pattern);
                    %T1= sscanf(fnames(K).name,'finished_Slice%d.png');
                    
                    if i==p
                        FileListArray{T1*NOPass*NOfTravels+(PassIndex-1)*NOfTravels+TravelIndex}=sprintf('%s',[PathName '\' sprintf('Slice_%06d',T1*NOPass*NOfTravels+(PassIndex-1)*NOfTravels+TravelIndex-1) sprintf('_Z_%f',(T1)*Zshift(4)+Zshift((TravelIndex))) '_YShiftInPixels_' num2str(RawShift-Scattering) '_Travel_' num2str(TravelIndex) '_Pass_'  num2str(PassIndex) '_Scattering_' num2str(Scattering) '_PH.bmp']);
                        
                    end
                    if w==1
                        patternAll{T1*NOPass*NOfTravels*NumberOfPrintHeads+(PassIndex-1)*NOfTravels*NumberOfPrintHeads+(TravelIndex-1)*NumberOfPrintHeads+i+1}=ones(TotalNoOfNozzles*NumberOfNozzleLines*2,LenFire);
                    end
                    patternAll{T1*NOPass*NOfTravels*NumberOfPrintHeads+(PassIndex-1)*NOfTravels*NumberOfPrintHeads+(TravelIndex-1)*NumberOfPrintHeads+i+1}(w*2:4:end,1:end)=~pattern;
                    
                    if w==(NumberOfNozzleLines)
                        ImWt=logical(patternAll{T1*NOPass*NOfTravels*NumberOfPrintHeads+(PassIndex-1)*NOfTravels*NumberOfPrintHeads+(TravelIndex-1)*NumberOfPrintHeads+i+1});
                        DropsPerFireOdd=mean(sum(~ImWt([1:4:end 2:4:end],:)));
                        DropsPerFireEven=mean(sum(~ImWt([3:4:end 4:4:end],:)));
                        A=cumsum( sum(~ImWt([1:4:end 2:4:end],:)));
                        MeanDropCountOdd=mean(A)./A(end);
                        if(isnan(   MeanDropCountOdd))
                            MeanDropCountOdd=0.5;
                        end
                        A=cumsum( sum(~ImWt([3:4:end 4:4:end],:)));
                        MeanDropCountEven=mean(A)./A(end);
                        if(isnan(   MeanDropCountEven))
                            MeanDropCountEven=0.5;
                        end
                        %'_OFlow_' num2str(DropsPerFireOdd) '_EFlow_' num2str(DropsPerFireEven) '_OMCount_' num2str(MeanDropCountOdd) '_EMCount_' num2str(MeanDropCountEven);
                        
                        info{T1*NOPass*NOfTravels+(PassIndex-1)*NOfTravels+TravelIndex}.val{1,i+1}=num2str(DropsPerFireOdd);
                        info{T1*NOPass*NOfTravels+(PassIndex-1)*NOfTravels+TravelIndex}.name{1,i+1}='OFlow';
                        
                        info{T1*NOPass*NOfTravels+(PassIndex-1)*NOfTravels+TravelIndex}.val{2,i+1}=num2str(DropsPerFireEven);
                        info{T1*NOPass*NOfTravels+(PassIndex-1)*NOfTravels+TravelIndex}.name{2,i+1}='EFlow';
                        
                        info{T1*NOPass*NOfTravels+(PassIndex-1)*NOfTravels+TravelIndex}.val{3,i+1}=num2str(MeanDropCountOdd);
                        info{T1*NOPass*NOfTravels+(PassIndex-1)*NOfTravels+TravelIndex}.name{3,i+1}='OMCount';
                        
                        info{T1*NOPass*NOfTravels+(PassIndex-1)*NOfTravels+TravelIndex}.val{4,i+1}=num2str(MeanDropCountEven);
                        info{T1*NOPass*NOfTravels+(PassIndex-1)*NOfTravels+TravelIndex}.name{4,i+1}='EMCount';
                        info{T1*NOPass*NOfTravels+(PassIndex-1)*NOfTravels+TravelIndex}.path{i+1}=[PathName '\' sprintf('Slice_%06d',T1*NOPass*NOfTravels+(PassIndex-1)*NOfTravels+TravelIndex-1) sprintf('_Z_%f',(T1)*Zshift(4)+Zshift((TravelIndex))) '_YShiftInPixels_' num2str(RawShift-Scattering) '_Travel_' num2str(TravelIndex) '_Pass_'  num2str(PassIndex) '_Scattering_' num2str(Scattering) '_PH' '.bmp'];
                        
                        imwrite(ImWt,[PathName '\' sprintf('Slice_%06d',T1*NOPass*NOfTravels+(PassIndex-1)*NOfTravels+TravelIndex-1) sprintf('_Z_%f',(T1)*Zshift(4)+Zshift((TravelIndex))) '_YShiftInPixels_' num2str(RawShift-Scattering) '_Travel_' num2str(TravelIndex) '_Pass_'  num2str(PassIndex) '_Scattering_' num2str(Scattering) '_PH_' num2str( i) '.bmp']);
                        patternAll{T1*NOPass*NOfTravels*NumberOfPrintHeads+(PassIndex-1)*NOfTravels*NumberOfPrintHeads+(TravelIndex-1)*NumberOfPrintHeads+i+1}=[];
                    end
                    
                end
            end
        end
        
        
        %imwrite(Im1,[PathName '\' name '_' num2str(i-1) ext]);
        
        
        
        
    end
    set(0, 'currentfigure', h);
    waitbar(K/numfids,h) ;%%0.074
    drawnow;
    %profile off
    %  profile viewer
end




for j=1:length(info)
    Sinfo=[];
    for p=1:4
        for k=1:4
            Sinfo=[Sinfo info{j}.name{p,k} num2str(k) '_' num2str(info{j}.val{p,k}) '_'];
        end
    end
    
    fprintf(F,'%s|%s\n',info{j}.path{1},Sinfo);
    
    %PathPart=split( FileListArray{j},'.bmp');
    %imwrite(logical(patternAll{j}),[PathPart{1} '_0.bmp']);
    
end
fclose(F);
close(h);
%create([PathName '\' i], 'dir')
end


%function y = b2d(x)
%% Preface
% Convert a binary array to a decimal number
% Similar to bin2dec but works with arrays instead of strings and is found
% to be rather faster
% Source: https://www.mathworks.com/matlabcentral/fileexchange/ ...
% 26447-efficient-convertors-between-binary-and-decimal-numbers

%% Compute y
% z = 2.^(length(x)-1:-1:0);
% y = sum(x.*z);

%% Better Performance using polyval
%y=polyval(x,2);

function y = d2b(x,DefineBorder)
y={};
CollorMap=[127 % no material
    166
    255 % support
    240
    192
    200
    0
    26];
y{1}=x(:,:,1)==CollorMap(1);
global BXP;
global BYP;
if(DefineBorder)
    
    BXP= find(~y{1}(:,round(end/2),1));
    
    BYP= find(~y{1}(round(end/2),:,1));
end
for n=1:length(CollorMap)
    %y{n}=x(:,:,1)==CollorMap(n);
    
    y{n}=x(BXP,BYP,1)==-1;
    
    %Model Channel
end
y{3}=x(BXP,BYP,1)==255%Support Cahnell
y{4}=((x(BXP,BYP,1)~=255) & (x(BXP,BYP,1)~=127))%Support Cahnell
%y{1}(:,:)


end

function ColorMap=ColorMapSearch(PathName)
fnames = dir([PathName '\*.png']);
numfids = length(fnames);

%%check the ColorMap
ColorMap=[];
for j=1:numfids
    InIm=imread([PathName '\' fnames(j).name]);
    g=InIm(:,:,1);
    g=g(:);
    ColorMap=[ColorMap;g];
    ColorMap=unique(ColorMap);
end
end