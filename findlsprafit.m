%%
% # NanOpticS
% In-depth analysis of NANomaterials for OPTICal localized surface plasmon resonance Sensing
% 
% ## Problems Solved
% This software solves the problem of analysing the LSPR bands of plasmonic nanomaterials in monitoring processes, with hundreds of transmittance spectra, in a few seconds. It was conceived mainly for high-resolution LSPR spectroscopy systems. A statistical analysis of the LSPR bands with normalized spectral distributions is performed to compare different sensing platforms, and the signal-to-noise ratio (SNR) is calculated for each analysed parameter. Furthermore, all the spectra are fitted with a polynomial function, which enables a fast and direct analysis of the transmittance (LSPR band) minimum.
% 
% ## Files
% The instalation file - NANOPTICS24.exe
% A text file containning the raw data to try the software - raw_transmittance_data.zip
% This Readme file
% 
% ## Installation
% NANOPTICS must be installed under a 64-bit platform running Windows. 
% During installation the user will be prompted to download (free) and install MATLAB Runtime R2018a (v 9.4).
% 
% ## Running the software
% To run the software, a file containing the transmittance spectra acquired over time, must have the wavelength (in nm) in the first column and the transmittance (in %), of the consecutive spectra, ordered in the next columns.
% Some initial parameters must be input in the GUI:
% 1-The user starts by writing the "Monitoring name" that will be added to the produced figures and files to export;
% 2-The user must choose whether a "LSPR band" exists, or not;
% 3-If the LSPR band exists and is clearly seen, the user can choose to let the algorithm do an "Automatic wavelength range". If there is no clear observation of the LSPR band, or if the user prefers to choose another range, the wavelength range can be made manually;
% 4-The "Wavelength lower limit" should be set to include the relative maximum transmittance, while "Wavelength higher limit" should be set, when possible, as having the same y-axis coordinate of the wavelength lower limit;
% 5-The user may choose to "Save all spectra figures" for all the spectra analysed. It will show all transmittance spectra acquired during monitoring, as well as the analysed range with the correspondent fitting obtained using a polynomial function. If this option is not selected, only the first two spectra, along with last one, will be saved as figures. By default, this checkbox is not selected since this operation is time consuming;
% 6-The elapsed "Time between each spectrum" must be written, in milliseconds;
% 7-The user can choose to apply a linear "Drift correction" of the monitored parameters;
% 8- The total "Number of cycles" (each cycle with the reference and the test conditions) must be input. The last cycle of the sensing test must be done only with the reference condition, as it was programmed to be the “control” cycle. The time for each condition (half-cycle) should be the same, and the total cycle time should also be constant;
% 9-The user must enter the portion (percentage) of each "Half cycle analysis" of the monitoring to be considered for the average parameter shift calculation;
% 10-If the user is running a monitoring process to assess the refractive index sensitivity (RIS) of the sensor, the "Sensitivity" check box needs to be ticked, and then the refractive indices of the reference "P1" and test "P2" materials have to be input. 
% During the analysis, a folder named “figures” and another one named “results” are created. The figures containing the individual spectra are saved in the first folder, and the figures containing all the calculated parameters in the second folder. All the files with the analysis results are also saved in the “results” folder.
% 
% ## Version
% 2020a (2.4)
% 
% ## License
% This project is licensed under the MIT License
%%
function [FOM] = findlsprafit (path,file,savefig,scantime,cycles,lic,Peak,drift,range,minWL,maxWL,valey,system0name,cycleFrom,cycleTo,sensitiv,RI1,RI2)
%filename='CO2_variation\test300001.txt'%000-003,164,228

close all
% LSPRBandw=zeros(300,1);
% LSPRBandT=zeros(300,1);
h=msgbox('Running analysis');
%clear all

% valey = 1;
% savefig = 0;
% cycleFrom = 85;
% cycleTo = 95;
% sensitiv = 1;
% RI1 = 1;%1.000036;
% RI2 = 1.000277;%1.000281;
% drift=1;
% Peak=1;
% range=0;
% minWL=420;
% maxWL=750;
% scantime=500;
% cycles=5;
% system0name = 'Nano#000dep000-AuTiO2-400-Vac-O2-manual';
% lic=4349722;
% path='C:\MatLab\Transmittance\Marco\Gas\Artigo\testOx074_Artigo';
% file='ProcTransm_300_900nm_movmean21';%'Manuela_gases_TestOx107_ProcTransm_300_900nm_movmean21';


FOM = [0 0];
testmaterial='Test material';
tempfile=strcat(file,'.lvm');
scantime=scantime/1000;
smooth=21;

system0min = minWL;
system0max = maxWL;
snrMinNoise=1.8;

colorTime=zeros(cycles,2);

systemname=system0name;
systemname2=system0name;
wlminnm=system0min;
wlmaxnm=system0max;

if sensitiv == 1
    RIa=RI1;%config(1,1); %RIa=cell2mat(RIa);
    RIb=RI2;%config(1,2); %RIb=cell2mat(RIb);
    deltaRI=RIb-RIa;
end

CentMom=1;
smethod='movmean';
smooth10=61;

basefoldername=system0name;

%Create folders if inexistent
status = mkdir(path,'figures');
status3 = mkdir(path,'results');

%% Date usage Limit
Date = today;
formatIn = 'dd-mmm-yyyy';
if lic == 4349722
    DateString = '31-dec-2020';
end

DateLimit=datenum(DateString,formatIn);

%%
%Look for pressure and temperature file
existspressurefile=~isempty(dir(fullfile(path,'pressure.csv')));
existstemperaturefile=~isempty(dir(fullfile(path,tempfile)));
pp=existspressurefile*1;
tt=existstemperaturefile*2;
pptt=pp+tt;


%%

if status==1
    tic
    
    %% Loop every column (every spectra in one txt file)
    
    
    a='Calculating LSPR band minimum position'
    
    lsprmin=single(wlminnm);
    lsprmax=single(wlmaxnm);
    filenumber=2;
    
    
    
    %Import data
    delimiterIn= '\t';
    headerlinesIn= 1;
    
    %linefilename = fullfile(path,'results',strcat(basefoldername,'_ProcTransm_300_900nm_',smethod,num2str(smooth),'.txt'));
    linefilename = fullfile(path,strcat(file,'.txt'));
    
    A=importdata(linefilename,delimiterIn,headerlinesIn);
    
    %Wavelength vector creation
    wavelength=A.data(:,1);
    %wavelength2=wavelength(378:end,:);
    
    %Transmittance matrix creation
    T=A.data(:,2:end);
    [~,ScanN]=size(T);
      
    if pptt ==2 || pptt==3
        %Import Temperature measurements
        tempfilename=fullfile(path,tempfile);
        
        %substitute commas by dots
        Data = fileread(tempfilename);
        Data = strrep(Data, ',', '.');
        FID = fopen(tempfilename, 'w');
        fwrite(FID, Data, 'char');
        fclose(FID);
        
        G=importdata(tempfilename,'\t',23);
        tempTime=G.data(1:end-2,1);
        tempData=G.data(3:end,2);
    end
    
    if pptt ==1 || pptt==3
        %Import Pressure
        pressurefilename=fullfile(path,'pressure.csv');
        B=importdata(pressurefilename,'\t',3);
        timepressure=B.data(:,1);
        pressure=B.data(:,2);
        
    end
    
    %Shift vector creation n cycles, x min, from 50% to 80%n of each half
    %cycle
    cycleFrom0=1-(cycleFrom/100);
    cycleTo0=1-(cycleTo/100);
    Shiftarray=zeros(cycles*2,2);
    for k=1:cycles
        i=k*2-1;
        j=k*2;
        Shiftarray(i,:)=[round(ScanN*(1/(cycles*2)*(i-cycleFrom0))) round(ScanN*(1/(cycles*2)*(i-cycleTo0)))];
        Shiftarray(j,:)=[round(ScanN*(1/(cycles*2)*(j-cycleFrom0))) round(ScanN*(1/(cycles*2)*(j-cycleTo0)))];
        
    end
    
    %Discrete Wavelength transmittance analysis
    wl=linspace(450,850,9);
    
    %Discrete Wavelength transmittance analysis file creation
    wlfilename = fullfile(path,'results',strcat(basefoldername,'_9Wavelengths_',smethod,num2str(smooth),'.txt'));
    fileID = fopen(wlfilename,'at');
    fprintf(fileID,'%2s %4f %4f %4f %4f %4f %4f %4f %4f %4f %4f %4f %4f %4f %4f %4f %4f %4f %4f\n','WL', wl(1), wl(2), wl(3), wl(4), wl(5), wl(6), wl(7), wl(8), wl(9), wl(1), wl(2), wl(3), wl(4), wl(5), wl(6), wl(7), wl(8), wl(9));
    fclose(fileID);
    
    
    %Discrete Wavelength transmittance calculation in 9 different wavelengths wl(i)
    %choosen for each system
    
    Iwl=zeros(1,9);
    Twl=zeros(ScanN,9);
    deltaTwl=zeros(ScanN,9);
    for k=1:ScanN
        for i =1:9
            [a(i),Iwl(i)]=min(abs(wl(i)-wavelength));
            Twl(k,i)=T(Iwl(i),k);
            deltaTwl(k,i)=Twl(k,i)-Twl(1,i);
        end
        
        
        wlfilename = fullfile(path,'results',strcat(basefoldername,'_9Wavelengths_',smethod,num2str(smooth),'.txt'));
        fileID = fopen(wlfilename,'at');
        fprintf(fileID,'%1d %5f %5f %5f %5f %5f %5f %5f %5f %5f %5f %5f %5f %5f %5f %5f %5f %5f %5f\n', k, Twl(k,1), Twl(k,2), Twl(k,3), Twl(k,4), Twl(k,5), Twl(k,6), Twl(k,7), Twl(k,8), Twl(k,9), deltaTwl(k,1), deltaTwl(k,2), deltaTwl(k,3), deltaTwl(k,4), deltaTwl(k,5), deltaTwl(k,6), deltaTwl(k,7), deltaTwl(k,8), deltaTwl(k,9));
        fclose(fileID);
        
    end
    
    %Select spectrum portion to search for raw LSPR edge and peak
    
    [~,wcutminL]=min(abs(lsprmin-wavelength));
    [~,wcutmaxL]=min(abs(lsprmax-wavelength));
    
    Tcut=T(wcutminL:wcutmaxL,:);
    wcut=wavelength(wcutminL:wcutmaxL);
    
    if range == 1 & Peak == 1
        TminMatr = islocalmin(T(:,1),'MinProminence',valey);
        TminMatrWL = wavelength(TminMatr);
        TminVect=T(:,1);
        minT=TminVect(TminMatr);
        minT=minT(end,1);
        
        TmaxMatr = islocalmax(T(:,1),'MinProminence',valey);
        TmaxMatrWL = wavelength(TmaxMatr);
        
        WLmin=TminMatrWL(end,1);
        WLmax=TmaxMatrWL(end,1);
        
        if WLmin<WLmax
            Peak = 0;
            range= 0;
            CreateStruct.Interpreter = 'tex';
            CreateStruct.WindowStyle = 'modal';
            msgbox('\fontsize{14} 1LSPR band was not detected. Select the wavelength range manually or run app without LSPR band Box selected.','Finished',CreateStruct);
            return
        end
    else
        %Calculate minimum transmittance point
        [minT,I]=min(Tcut(:,1));
        
        WLmin=wcut(I);
        
        %Calculate maximum transmittance point
        [~,Imax]=max(Tcut(1:I,1));
        
        WLmax=wcut(Imax);
    end
    
    %Shift real spectra down to LSPR peak at 0 transmittance and cut left
    %LSPR peak
    [~,wcutminL2]=min(abs(WLmax-wavelength));
    [~,wcutmaxL2]=min(abs(WLmin-wavelength));
    
    Tcut2=T(wcutminL2:wcutmaxL2,:);
    wcut2=wavelength(wcutminL2:wcutmaxL2);
    integralLeft=trapz(wcut2,Tcut2);
    
    %Shift real spectra down to LSPR peak at 0 transmittance and cut right
    %LSPR peak
    WlR=WLmin-WLmax+WLmin;
    [~, wcutminL3]=min(abs(WLmin-wavelength));
    [~, wcutmaxL3]=min(abs(WlR-wavelength));
    Tcut3=T(wcutminL3:wcutmaxL3,:);
    wcut3=wavelength(wcutminL3:wcutmaxL3);
    integralRight=trapz(wcut3,Tcut3);
    
    
    
    %% Time vector
    xx=(linspace(1,ScanN,ScanN))';
    xx=xx*scantime;
    maxTime=max(xx);
    cycleTime=maxTime/cycles;
    
    
    for k=1:cycles
        i=k*2-1;
        colorTime(k,1)=cycleTime*i/2;
        colorTime(k,2)=cycleTime*i/2+cycleTime/2;
    end
    
    
    
    %% Build LSPR band before LSPR edge and after LSPR peak
    %
    if range == 1 & Peak == 1%automatic or not wavelength range detection
        WLmin1=WLmin;
        WLmax1=WLmax;
        deltaWL=WLmin-WLmax;
        if WLmax1<350
            WLminL=340;
        else
            WLminL=WLmax1-10;
        end
        
        if deltaWL<180
            WLminR=WLmin1+deltaWL/3*5;
        else
            WLminR=850;
        end
        
    else
        WLminR=maxWL;
        WLminL=minWL;
    end
    
    %% Check if there is a peak
    if WLmin==WLmax || WLmin==WLminR  %with this conditions there is no peak
        
        
        Peak=0;
        
        if range == 1 & Peak == 1
            
            WLminL=350;
            WLminR=850;
            WLmiddle=600;
        else
            
            WLmiddle=(WLminR-WLminL)/2;
        end
        
        [~, BandWLleft1]=min(abs(WLminL-wavelength));
        [~, BandWLright1]=min(abs(WLminR-wavelength));
        [~, BandWLmiddle]=min(abs(WLmiddle-wavelength));
        
        %Band without peak
        TBandcutPeak1=T(BandWLleft1:BandWLright1,:);
        wBandcutPeak1=wavelength(BandWLleft1:BandWLright1);
        
        %Real data Integration
        realDataIntegral=trapz(wBandcutPeak1,TBandcutPeak1);
        %Left
        Tcut2=T(BandWLleft1:BandWLmiddle,:);
        wcut2=wavelength(BandWLleft1:BandWLmiddle);
        integralLeft=trapz(wcut2,Tcut2);
        %right
        Tcut3=T(BandWLmiddle:BandWLright1,:);
        wcut3=wavelength(BandWLmiddle:BandWLright1);
        
        integralRight=trapz(wcut3,Tcut3);
        
    else
        
        
        
        [~, BandWLleft1]=min(abs(WLminL-wavelength));
        [~, BandWLright1]=min(abs(WLminR-wavelength));
        
        %Band around peak
        TBandcutPeak1=T(BandWLleft1:BandWLright1,:);
        wBandcutPeak1=wavelength(BandWLleft1:BandWLright1);
        
        %Real data Integration
        realDataIntegral=trapz(wBandcutPeak1,TBandcutPeak1);
        
                
    end
    
    if drift==1
        Twlcorrect=detrend(Twl);
        realDataIntegralcorrect=detrend(realDataIntegral);
        integralLeftcorrect=detrend(integralLeft);
        integralRightcorrect=detrend(integralRight);
    else
        realDataIntegralcorrect=realDataIntegral;
        Twlcorrect=Twl;
        integralLeftcorrect=integralLeft;
        integralRightcorrect=integralRight;
    end
    
    %% Central Moments
    if CentMom==1 & Peak ==1
        %Prepare data for Central Moments analysis
        wBandcutPeak1CMom=wBandcutPeak1/wBandcutPeak1(1,1); %Wavelength normalization
        TBandcutPeak1CMom=1-TBandcutPeak1/100; %from 100 % to 1 and then "1-T" normalization
        %M0 calculation - 1-T data trapezoidal method integral
        Mzero=trapz(wBandcutPeak1CMom,TBandcutPeak1CMom);
        %Transmittance data (1-T) normalization using M0
        pseusopdf=TBandcutPeak1CMom./Mzero; %defines "distribution function" "P(?)"
        Mzero=Mzero/Mzero(1); %Mzero normalization to first integral
        %M1 calculation - Expectation value
        M1Notcentred=wBandcutPeak1CMom.*pseusopdf;
        M1=trapz(wBandcutPeak1CMom,M1Notcentred);
        %M2 calculation - Variance or deviation
        M2Notcentred=(wBandcutPeak1CMom-M1).^2.*pseusopdf;
        M2=trapz(wBandcutPeak1CMom,M2Notcentred);
        %M3 calculation - Skewness - lack of simetry. If positive left
        %shifted, if negative, right shifted
        M3Notcentred=(wBandcutPeak1CMom-M1).^3.*pseusopdf;
        M3=trapz(wBandcutPeak1CMom,M3Notcentred)./M2.^(3/2);
        %M4 calculation - Kurtosis - tail shape. If positive, heavy tailed
        %distribution, more than an exponential, if negative, the opposite
        M4Notcentred=(wBandcutPeak1CMom-M1).^4.*pseusopdf;
        M4=trapz(wBandcutPeak1CMom,M4Notcentred)./M2.^(2);
    
        %Linear drift correction
        if drift==1
            Mzerocorrect=detrend(Mzero);
            M1correct=detrend(M1);
            M2correct=detrend(M2);
            M3correct=detrend(M3);
            M4correct=detrend(M4);
        else
            Mzerocorrect=Mzero;
            M1correct=M1;
            M2correct=M2;
            M3correct=M3;
            M4correct=M4;
        end
    else
        
    end
    %% LSPR band fitting and calculation of the maximum and minimum
    %transmittance wavelength
    
    coeffitPeakn=zeros(ScanN,10);
    coeffitPeakn4=zeros(ScanN,6);
    minwlfitn=zeros(ScanN,1);
    minTfitn=zeros(ScanN,1);
    maxwlfitn=zeros(ScanN,1);
    maxTfitn=zeros(ScanN,1);
    fullH=zeros(ScanN,1);
    halfH=zeros(ScanN,1);
    FWHHL=zeros(ScanN,1);
    FWHHR=zeros(ScanN,1);
    FWHH=zeros(ScanN,1);
    
    delete(h);
    
    for k=1:ScanN
        if Peak==0
            
        else
            
            BandfitPeak9=fit(wBandcutPeak1,TBandcutPeak1(:,k),'poly9'); %Poly4
            BandfitPeak4=fit(wBandcutPeak1,TBandcutPeak1(:,k),'poly5');
            coeffitPeakn(k,:)=coeffvalues(BandfitPeak9);
            coeffitPeakn4(k,:)=coeffvalues(BandfitPeak4);
            
            BandfitPeakfunctn=@(x)coeffitPeakn(k,1)*x.^9+coeffitPeakn(k,2)*x.^8+coeffitPeakn(k,3)*x.^7+...
                coeffitPeakn(k,4)*x.^6+coeffitPeakn(k,5)*x.^5+coeffitPeakn(k,6)*x.^4+...
                coeffitPeakn(k,7)*x.^3+coeffitPeakn(k,8)*x.^2+coeffitPeakn(k,9)*x+coeffitPeakn(k,10);
            
            BandfitPeakfunctn4=@(x)coeffitPeakn4(k,1)*x.^5+coeffitPeakn4(k,2)*x.^4+...
                coeffitPeakn4(k,3)*x.^3+coeffitPeakn4(k,4)*x.^2+coeffitPeakn4(k,5)*x+coeffitPeakn4(k,6);
            
            BandfitPeakfunctnminus=@(x)-BandfitPeakfunctn(x);
            [minwlfitn(k,1),minTfitn(k,1)] = fminbnd(BandfitPeakfunctn,WLmax,WLminR);
            [maxwlfitn(k,1),maxTfitn(k,1)] = fminbnd(BandfitPeakfunctnminus,WLminL,WLmin);
            
            maxTfitn(k,1)=maxTfitn(k,1).*-(1);
            fullH(k,1)=maxTfitn(k,1)-minTfitn(k,1);
            halfH(k,1)=fullH(k,1)/2;
            TShift=halfH(k,1)+minTfitn(k,1);
            
            BandfitPeakfunctnTShift=@(x) BandfitPeakfunctn(x)-TShift;
            x01 = [WLmax WLmin]; % initial point
            
           
            try
                TFW = fzero(BandfitPeakfunctnTShift,x01);
                FWHHL(k,1) = TFW;  % If FZERO errors, this assignment won't happen.
                
            catch
                fprintf('FWHHL No root found for parameter %f\n',k); %optional
                
                FWHHL(k,1) = 0;
               
            end
            
            x02try = [WLmin WLminR];
            
            try
                TFW = fzero(BandfitPeakfunctnTShift,x02try);
                FWHHR(k,1) = TFW;  % If FZERO errors, this assignment won't happen.
                
            catch
                fprintf('FWHHR No root found for parameter %f\n',k); %optional
                
                FWHHR(k,1) = 0;
                
            end
            
            % full width at half height calculation
            if FWHHL(k,1) == 0 || FWHHR(k,1) == 0;
                if k == 1
                    FWHH(k,1) = 0; %full width at half height calculation
                    LWHH(k,1) = 0; %left width at half height calculation
                    RWHH(k,1) = 0; %right width at half height calculation
                else
                    ki = k-1;
                    FWHH(k,1) = FWHH(ki,1); %full width at half height calculation
                    LWHH(k,1) = LWHH(ki,1); %left width at half height calculation
                    RWHH(k,1) = RWHH(ki,1); %right width at half height calculation
                end
                
            else
                FWHH(k,1) = FWHHR(k,1)-FWHHL(k,1); %full width at half height calculation
                LWHH(k,1) = minwlfitn(k,1)-FWHHL(k,1); %left width at half height calculation
                RWHH(k,1) = FWHHR(k,1)-minwlfitn(k,1); %right width at half height calculation
            end
            
            
            
        end
        
        %% Save Some Spectra Analysis Figures
        if savefig == 2 || k == 3 || k == ScanN-1
            
            
            xmarkers=WLmin;
            ymarkers=linspace(2,minT+2,2);
            figure(10)
            
            if Peak==0
            else
                subplot(1,2,1);
            end
            plot(wavelength,T(:,k),'r',xmarkers,ymarkers,'bv')
            set(gca,'fontsize',14);
            axis([350 850 0 100]);
            ylabel('Transmittance (%)', 'Fontsize',14);
            xlabel('Wavelength (nm)', 'Fontsize',14);
            
            
            hold on
            
            plot(wcut,Tcut(:,k),'g',xmarkers,ymarkers,'bv')
            text(xmarkers-60,minT+10,strcat('WLmin = ',num2str(WLmin,'%0.2f'),' nm'))
            
             
            
            if Peak==0
                title(systemname, 'Fontsize',16);
            else
                subplot(1,2,2);
                plot(wBandcutPeak1,TBandcutPeak1(:,k),'*')
                hold on
                fplot(BandfitPeakfunctn,'r')
                hold on
                fplot(BandfitPeakfunctn4,'y')
                hold off
                set(gca,'fontsize',14);
                axis([WLminL WLminR 0 100]);
                ylabel('Transmittance (%)', 'Fontsize',14);
                xlabel('Wavelength (nm)', 'Fontsize',14);
                suptitle(systemname);
%                 
            end
            set(gcf, 'Position', get(0, 'Screensize'));
            
            savefigname=fullfile(path,'figures',strcat(basefoldername,'transmittance_Fittin_LSPRband_',smethod,num2str(smooth),'_',sprintf('%05d',k),'.png'));
            saveas(gcf,savefigname)
            savefigname=fullfile(path,'figures',strcat(basefoldername,'transmittance_Fittin_LSPRband_',smethod,num2str(smooth),'_',sprintf('%05d',k),'.fig'));
            saveas(gcf,savefigname)
           
        end
    end
    
    if Peak == 1 && sensitiv == 1
        
        %FWHH without zeros...
        FWHHnonzero = find(FWHH);
        FWHHnonzerofirst = FWHHnonzero(1,1);
        FWHHFOM = FWHH(FWHHnonzerofirst,1);
    end
    
    
    
    %% Drift correction of LSPR Band fit
    %Linear drift correction
    if Peak==0
    else
        if drift==1
            minwlfitncorrect=detrend(minwlfitn);
            minTfitncorrect=detrend(minTfitn);
            FWHHcorrect=detrend(FWHH);
            LWHHcorrect=detrend(LWHH);
            RWHHcorrect=detrend(RWHH);
        else
            minwlfitncorrect=minwlfitn;
            minTfitncorrect=minTfitn;
            FWHHcorrect=FWHH;
            LWHHcorrect=LWHH;
            RWHHcorrect=RWHH;
        end
        
        
    end
    
    
    %% DELTA Transmitance - Optical Transmitance Change (OTC)
    if cycles==0
    else
        %
        
        %%Calculate full optical spectra diference at each 10% (5 cycles, 2 min in each condition)
        %%20 min after the drift correction
        legendOTC=cell(1,cycles);
        
        
        for k=1:cycles
            i=k*2-1;
            j=k*2;
            
            A=smoothdata(mean(T(:,Shiftarray(j,1):Shiftarray(j,2)),2),smethod,smooth10);%Partial OTC Data with part 2 of each cycle
            B=smoothdata(mean(T(:,Shiftarray(i,1):Shiftarray(i,2)),2),smethod,smooth10);%Partial OTC Data with part 1 of each cycle
            
            OTC(:,k)=A-B;
            
            C=strcat('cycle',string(k));
            legendOTC{k}=sprintf(C);
            
        end
    end
    
    %% Shift calculation of each cycle
    if cycles ==0
    else
        %% Shift at each cycle of RAW data and SNR
        
        % 450 nm
        % Calculate the Shift at each cycle
        % after the drift correction if "drift" is set to "1"
        
        for k=1:cycles
            i=k*2-1;
            j=k*2;
            A=Twlcorrect(Shiftarray(j,1):Shiftarray(j,2),1);%Partial Data with part 2 of each cycle
            B=Twlcorrect(Shiftarray(i,1):Shiftarray(i,2),1);%Partial Data with part 1 of each cycle
            
            Asize=length(A);
            Bsize=length(B);
            
            if Asize == Bsize
                shift=A-B;
            else
                if Asize < Bsize
                    B = B(1:Asize,1);
                else
                    A = A(1:Bsize,1);
                end
                shift=A-B;
            end
            
            

            SNRnoise=shift-mean(shift);
            SNRSigrms=rms(shift);
            SNRNoirms=rms(SNRnoise);
            SNRall(k)=SNRSigrms/SNRNoirms; %SNR of signal rms divided by noise rms
            
            
            
            shift450mean(k,1)=mean(shift); %Shifts of each cycle
            shift450mean(k,2)=std(shift)/sqrt(Shiftarray(j,2)-Shiftarray(j,1)); %ST error of each cycle
        end
        
        
        meanshift450all(:,1)=mean(shift450mean(1:end-1,1)); %Mean
        
        meanshift450all(:,2) = mad(shift450mean(1:end-1,1)); %Mean deviation
        
        if abs(meanshift450all(:,1)/meanshift450all(:,2))<snrMinNoise
            SNR450='Unstable Signal';
        else
            SNR450=mean(SNRall());
        end
        
        if sensitiv ==1
            RISX450=shift450mean/deltaRI;
            RISPeak450=meanshift450all/deltaRI; %Sensitivity if RIS file exists
        end
        
        % 500 nm
        % Calculate the Shift at each cycle
        % after the drift correction if "drift" is set to "1"
        for k=1:cycles
            i=k*2-1;
            j=k*2;
            A=Twlcorrect(Shiftarray(j,1):Shiftarray(j,2),2);%Partial Data with part 2 of each cycle
            B=Twlcorrect(Shiftarray(i,1):Shiftarray(i,2),2);%Partial Data with part 1 of each cycle
            
            Asize=length(A);
            Bsize=length(B);
            if Asize == Bsize
                shift=A-B;
            else
                if Asize < Bsize
                    B = B(1:Asize,1);
                else
                    A = A(1:Bsize,1);
                end
                shift=A-B;
            end
            
            SNRnoise=shift-mean(shift);
            SNRSigrms=rms(shift);
            SNRNoirms=rms(SNRnoise);
            SNRall(k)=SNRSigrms/SNRNoirms; %SNR of signal rms divided by noise rms
            
            shift500mean(k,1)=mean(shift); %Shifts of each cycle
            shift500mean(k,2)=std(shift)/sqrt(Shiftarray(j,2)-Shiftarray(j,1)); %ST error of each cycle
        end
        
        meanshift500all(:,1)=mean(shift500mean(1:end-1,1)); %Mean
        meanshift500all(:,2) = mad(shift500mean(1:end-1,1)); %Mean deviation
        if abs(meanshift500all(:,1)/meanshift500all(:,2))<snrMinNoise
            SNR500='Unstable Signal';
        else
            SNR500=mean(SNRall());
        end
        if sensitiv ==1
            RISX500=shift500mean/deltaRI;
            RISPeak500=meanshift500all/deltaRI; %Sensitivity if RIS file exists
        end
        
        % 550 nm
        % Calculate the Shift at each cycle
        % after the drift correction if "drift" is set to "1"
        for k=1:cycles
            i=k*2-1;
            j=k*2;
            A=Twlcorrect(Shiftarray(j,1):Shiftarray(j,2),3);%Partial Data with part 2 of each cycle
            B=Twlcorrect(Shiftarray(i,1):Shiftarray(i,2),3);%Partial Data with part 1 of each cycle
            
            Asize=length(A);
            Bsize=length(B);
            if Asize == Bsize
                shift=A-B;
            else
                if Asize < Bsize
                    B = B(1:Asize,1);
                else
                    A = A(1:Bsize,1);
                end
                shift=A-B;
            end
            
            SNRnoise=shift-mean(shift);
            SNRSigrms=rms(shift);
            SNRNoirms=rms(SNRnoise);
            SNRall(k)=SNRSigrms/SNRNoirms; %SNR of signal rms divided by noise rms
            
            shift550mean(k,1)=mean(shift); %Shifts of each cycle
            shift550mean(k,2)=std(shift)/sqrt(Shiftarray(j,2)-Shiftarray(j,1)); %ST error of each cycle
        end
        
        meanshift550all(:,1)=mean(shift550mean(1:end-1,1)); %Mean
        meanshift550all(:,2) = mad(shift550mean(1:end-1,1)); %Mean deviation
        if abs(meanshift550all(:,1)/meanshift550all(:,2))<snrMinNoise
            SNR550='Unstable Signal';
        else
            SNR550=mean(SNRall());
        end
        if sensitiv ==1
            RISX550=shift550mean/deltaRI;
            RISPeak550=meanshift550all/deltaRI; %Sensitivity if RIS file exists
        end
        
        % 600 nm
        % Calculate the Shift at each cycle
        % after the drift correction if "drift" is set to "1"
        for k=1:cycles
            i=k*2-1;
            j=k*2;
            A=Twlcorrect(Shiftarray(j,1):Shiftarray(j,2),4);%Partial Data with part 2 of each cycle
            B=Twlcorrect(Shiftarray(i,1):Shiftarray(i,2),4);%Partial Data with part 1 of each cycle
            
            Asize=length(A);
            Bsize=length(B);
            if Asize == Bsize
                shift=A-B;
            else
                if Asize < Bsize
                    B = B(1:Asize,1);
                else
                    A = A(1:Bsize,1);
                end
                shift=A-B;
            end
            
            SNRnoise=shift-mean(shift);
            SNRSigrms=rms(shift);
            SNRNoirms=rms(SNRnoise);
            SNRall(k)=SNRSigrms/SNRNoirms; %SNR of signal rms divided by noise rms
            
            shift600mean(k,1)=mean(shift); %Shifts of each cycle
            shift600mean(k,2)=std(shift)/sqrt(Shiftarray(j,2)-Shiftarray(j,1)); %ST error of each cycle
        end
        
        meanshift600all(:,1)=mean(shift600mean(1:end-1,1)); %Mean
        meanshift600all(:,2) = mad(shift600mean(1:end-1,1)); %Mean deviation
        if abs(meanshift600all(:,1)/meanshift600all(:,2)) < snrMinNoise
            SNR600='Unstable Signal';
        else
            SNR600=mean(SNRall());
        end
        if sensitiv ==1
            RISX600=shift600mean/deltaRI;
            RISPeak600=meanshift600all/deltaRI; %Sensitivity if RIS file exists
        end
        
        % 650 nm
        % Calculate the Shift at each cycle
        % after the drift correction if "drift" is set to "1"
        for k=1:cycles
            i=k*2-1;
            j=k*2;
            A=Twlcorrect(Shiftarray(j,1):Shiftarray(j,2),5);%Partial Data with part 2 of each cycle
            B=Twlcorrect(Shiftarray(i,1):Shiftarray(i,2),5);%Partial Data with part 1 of each cycle
            
            Asize=length(A);
            Bsize=length(B);
            if Asize == Bsize
                shift=A-B;
            else
                if Asize < Bsize
                    B = B(1:Asize,1);
                else
                    A = A(1:Bsize,1);
                end
                shift=A-B;
            end
            
            SNRnoise=shift-mean(shift);
            SNRSigrms=rms(shift);
            SNRNoirms=rms(SNRnoise);
            SNRall(k)=SNRSigrms/SNRNoirms; %SNR of signal rms divided by noise rms
            
            shift650mean(k,1)=mean(shift); %Shifts of each cycle
            shift650mean(k,2)=std(shift)/sqrt(Shiftarray(j,2)-Shiftarray(j,1)); %ST error of each cycle
        end
        
        meanshift650all(:,1)=mean(shift650mean(1:end-1,1)); %Mean
        meanshift650all(:,2) = mad(shift650mean(1:end-1,1)); %Mean deviation
        if abs(meanshift650all(:,1)/meanshift650all(:,2)) < snrMinNoise
            SNR650='Unstable Signal';
        else
            SNR650=mean(SNRall());
        end
        if sensitiv ==1
            RISX650=shift650mean/deltaRI;
            RISPeak650=meanshift650all/deltaRI; %Sensitivity if RIS file exists
        end
        % 700 nm
        % Calculate the Shift at each cycle
        % after the drift correction if "drift" is set to "1"
        for k=1:cycles
            i=k*2-1;
            j=k*2;
            A=Twlcorrect(Shiftarray(j,1):Shiftarray(j,2),6);%Partial Data with part 2 of each cycle
            B=Twlcorrect(Shiftarray(i,1):Shiftarray(i,2),6);%Partial Data with part 1 of each cycle
            
            Asize=length(A);
            Bsize=length(B);
            if Asize == Bsize
                shift=A-B;
            else
                if Asize < Bsize
                    B = B(1:Asize,1);
                else
                    A = A(1:Bsize,1);
                end
                shift=A-B;
            end
            
            SNRnoise=shift-mean(shift);
            SNRSigrms=rms(shift);
            SNRNoirms=rms(SNRnoise);
            SNRall(k)=SNRSigrms/SNRNoirms; %SNR of signal rms divided by noise rms
            
            shift700mean(k,1)=mean(shift); %Shifts of each cycle
            shift700mean(k,2)=std(shift)/sqrt(Shiftarray(j,2)-Shiftarray(j,1)); %ST error of each cycle
        end
        
        meanshift700all(:,1)=mean(shift700mean(1:end-1,1)); %Mean
        meanshift700all(:,2) = mad(shift700mean(1:end-1,1)); %Mean deviation
        if abs(meanshift700all(:,1)/meanshift700all(:,2)) < snrMinNoise
            SNR700='Unstable Signal';
        else
            SNR700=mean(SNRall());
        end
        if sensitiv ==1
            RISX700=shift450mean/deltaRI;
            RISPeak700=meanshift700all/deltaRI; %Sensitivity if RIS file exists
        end
        % 750 nm
        % Calculate the Shift at each cycle
        % after the drift correction if "drift" is set to "1"
        for k=1:cycles
            i=k*2-1;
            j=k*2;
            A=Twlcorrect(Shiftarray(j,1):Shiftarray(j,2),7);%Partial Data with part 2 of each cycle
            B=Twlcorrect(Shiftarray(i,1):Shiftarray(i,2),7);%Partial Data with part 1 of each cycle
            
            Asize=length(A);
            Bsize=length(B);
            if Asize == Bsize
                shift=A-B;
            else
                if Asize < Bsize
                    B = B(1:Asize,1);
                else
                    A = A(1:Bsize,1);
                end
                shift=A-B;
            end
            
            SNRnoise=shift-mean(shift);
            SNRSigrms=rms(shift);
            SNRNoirms=rms(SNRnoise);
            SNRall(k)=SNRSigrms/SNRNoirms; %SNR of signal rms divided by noise rms
            
            shift750mean(k,1)=mean(shift); %Shifts of each cycle
            shift750mean(k,2)=std(shift)/sqrt(Shiftarray(j,2)-Shiftarray(j,1)); %ST error of each cycle
        end
        
        meanshift750all(:,1)=mean(shift750mean(1:end-1,1)); %Mean
        meanshift750all(:,2) = mad(shift750mean(1:end-1,1)); %Mean deviation
        if abs(meanshift750all(:,1)/meanshift750all(:,2)) < snrMinNoise
            SNR750='Unstable Signal';
        else
            SNR750=mean(SNRall());
        end
        if sensitiv ==1
            RISX750=shift750mean/deltaRI;
            RISPeak750=meanshift750all/deltaRI; %Sensitivity if RIS file exists
        end
        % 800 nm
        % Calculate the Shift at each cycle
        % after the drift correction if "drift" is set to "1"
        for k=1:cycles
            i=k*2-1;
            j=k*2;
            A=Twlcorrect(Shiftarray(j,1):Shiftarray(j,2),8);%Partial Data with part 2 of each cycle
            B=Twlcorrect(Shiftarray(i,1):Shiftarray(i,2),8);%Partial Data with part 1 of each cycle
            
            Asize=length(A);
            Bsize=length(B);
            if Asize == Bsize
                shift=A-B;
            else
                if Asize < Bsize
                    B = B(1:Asize,1);
                else
                    A = A(1:Bsize,1);
                end
                shift=A-B;
            end
            
            SNRnoise=shift-mean(shift);
            SNRSigrms=rms(shift);
            SNRNoirms=rms(SNRnoise);
            SNRall(k)=SNRSigrms/SNRNoirms; %SNR of signal rms divided by noise rms
            
            shift800mean(k,1)=mean(shift); %Shifts of each cycle
            shift800mean(k,2)=std(shift)/sqrt(Shiftarray(j,2)-Shiftarray(j,1)); %ST error of each cycle
        end
        
        meanshift800all(:,1)=mean(shift800mean(1:end-1,1)); %Mean
        meanshift800all(:,2) = mad(shift800mean(1:end-1,1)); %Mean deviation
        if abs(meanshift800all(:,1)/meanshift800all(:,2)) < snrMinNoise
            SNR800='Unstable Signal';
        else
            SNR800=mean(SNRall());
        end
        if sensitiv ==1
            RISX800=shift800mean/deltaRI;
            RISPeak800=meanshift800all/deltaRI; %Sensitivity if RIS file exists
        end
        
        % 850 nm
        % Calculate the Shift at each cycle
        % after the drift correction if "drift" is set to "1"
        for k=1:cycles
            i=k*2-1;
            j=k*2;
            A=Twlcorrect(Shiftarray(j,1):Shiftarray(j,2),9);%Partial Data with part 2 of each cycle
            B=Twlcorrect(Shiftarray(i,1):Shiftarray(i,2),9);%Partial Data with part 1 of each cycle
            
            Asize=length(A);
            Bsize=length(B);
            if Asize == Bsize
                shift=A-B;
            else
                if Asize < Bsize
                    B = B(1:Asize,1);
                else
                    A = A(1:Bsize,1);
                end
                shift=A-B;
            end
            
            SNRnoise=shift-mean(shift);
            SNRSigrms=rms(shift);
            SNRNoirms=rms(SNRnoise);
            SNRall(k)=SNRSigrms/SNRNoirms; %SNR of signal rms divided by noise rms
            
            shift850mean(k,1)=mean(shift); %Shifts of each cycle
            shift850mean(k,2)=std(shift)/sqrt(Shiftarray(j,2)-Shiftarray(j,1)); %ST error of each cycle
        end
        
        meanshift850all(:,1)=mean(shift850mean(1:end-1,1)); %Mean
        meanshift850all(:,2) = mad(shift850mean(1:end-1,1)); %Mean deviation
        SNR850=snr(meanshift850all(:,1),meanshift850all(:,2)); %SNR
        if abs(meanshift850all(:,1)/meanshift850all(:,2)) < snrMinNoise
            SNR850='Unstable Signal';
        else
            SNR850=mean(SNRall());
        end
        if sensitiv ==1
            RISX850=shift850mean/deltaRI;
            RISPeak850=meanshift850all/deltaRI; %Sensitivity if RIS file exists
        end
        
        % realDataIntegral
        % Calculate the Shift at each cycle
        % after the drift correction if "drift" is set to "1"
        for k=1:cycles
            i=k*2-1;
            j=k*2;
            A=realDataIntegralcorrect(Shiftarray(j,1):Shiftarray(j,2));%Partial Data with part 2 of each cycle
            B=realDataIntegralcorrect(Shiftarray(i,1):Shiftarray(i,2));%Partial Data with part 1 of each cycle
            
            Asize=length(A);
            Bsize=length(B);
            if Asize == Bsize
                shift=A-B;
            else
                if Asize < Bsize
                    B = B(1,1:Asize);
                else
                    A = A(1,1:Bsize);
                end
                shift=A-B;
            end
            
            SNRnoise=shift-mean(shift);
            SNRSigrms=rms(shift);
            SNRNoirms=rms(SNRnoise);
            SNRall(k)=SNRSigrms/SNRNoirms; %SNR of signal rms divided by noise rms
            
            shiftrealDataIntegralmean(k,1)=mean(shift); %Shifts of each cycle
            shiftrealDataIntegralmean(k,2)=std(shift)/sqrt(Shiftarray(j,2)-Shiftarray(j,1)); %ST error of each cycle
        end
        
        meanshiftrealDataIntegralall(:,1)=mean(shiftrealDataIntegralmean(1:end-1,1)); %Mean
        meanshiftrealDataIntegralall(:,2) = mad(shiftrealDataIntegralmean(1:end-1,1)); %Mean deviation
        if abs(meanshiftrealDataIntegralall(:,1)/meanshiftrealDataIntegralall(:,2)) < snrMinNoise
            SNRrealDataIntegral='Unstable Signal';
        else
            SNRrealDataIntegral=mean(SNRall());
        end
        if sensitiv ==1
            RISXrealDataIntegral=shiftrealDataIntegralmean/deltaRI;
            RISPeakrealDataIntegral=meanshiftrealDataIntegralall/deltaRI; %Sensitivity if RIS file exists
        end
        
        % integralLeft
        % Calculate the Shift at each cycle
        % after the drift correction if "drift" is set to "1"
        for k=1:cycles
            i=k*2-1;
            j=k*2;
            A=integralLeftcorrect(Shiftarray(j,1):Shiftarray(j,2));%Partial Data with part 2 of each cycle
            B=integralLeftcorrect(Shiftarray(i,1):Shiftarray(i,2));%Partial Data with part 1 of each cycle
            
            Asize=length(A);
            Bsize=length(B);
            if Asize == Bsize
                shift=A-B;
            else
                if Asize < Bsize
                    B = B(1,1:Asize);
                else
                    A = A(1,1:Bsize);
                end
                shift=A-B;
            end
            
            SNRnoise=shift-mean(shift);
            SNRSigrms=rms(shift);
            SNRNoirms=rms(SNRnoise);
            SNRall(k)=SNRSigrms/SNRNoirms; %SNR of signal rms divided by noise rms
            
            shiftintegralLeftmean(k,1)=mean(shift); %Shifts of each cycle
            shiftintegralLeftmean(k,2)=std(shift)/sqrt(Shiftarray(j,2)-Shiftarray(j,1)); %ST error of each cycle
        end
        
        meanshiftintegralLeftall(:,1)=mean(shiftintegralLeftmean(1:end-1,1)); %Mean
        meanshiftintegralLeftall(:,2) = mad(shiftintegralLeftmean(1:end-1,1)); %Mean deviation
        if abs(meanshiftintegralLeftall(:,1)/meanshiftintegralLeftall(:,2)) < snrMinNoise
            SNRintegralLeft='Unstable Signal';
        else
            SNRintegralLeft=mean(SNRall());
        end
        if sensitiv ==1
            RISXintegralLeft=shiftintegralLeftmean/deltaRI;
            RISPeakintegralLeft=meanshiftintegralLeftall/deltaRI; %Sensitivity if RIS file exists
        end
        
        % integralRight
        % Calculate the Shift at each cycle
        % after the drift correction if "drift" is set to "1"
        for k=1:cycles
            i=k*2-1;
            j=k*2;
            A=integralRightcorrect(Shiftarray(j,1):Shiftarray(j,2));%Partial Data with part 2 of each cycle
            B=integralRightcorrect(Shiftarray(i,1):Shiftarray(i,2));%Partial Data with part 1 of each cycle
            
            Asize=length(A);
            Bsize=length(B);
            if Asize == Bsize
                shift=A-B;
            else
                if Asize < Bsize
                    B = B(1,1:Asize);
                else
                    A = A(1,1:Bsize);
                end
                shift=A-B;
            end
            
            SNRnoise=shift-mean(shift);
            SNRSigrms=rms(shift);
            SNRNoirms=rms(SNRnoise);
            SNRall(k)=SNRSigrms/SNRNoirms; %SNR of signal rms divided by noise rms
            
            shiftintegralRightmean(k,1)=mean(shift); %Shifts of each cycle
            shiftintegralRightmean(k,2)=std(shift)/sqrt(Shiftarray(j,2)-Shiftarray(j,1)); %ST error of each cycle
        end
        
        meanshiftintegralRightall(:,1)=mean(shiftintegralRightmean(1:end-1,1)); %Mean
        meanshiftintegralRightall(:,2) = mad(shiftintegralRightmean(1:end-1,1)); %Mean deviation
        if abs(meanshiftintegralRightall(:,1)/meanshiftintegralRightall(:,2)) < snrMinNoise
            SNRintegralRight='Unstable Signal';
        else
            SNRintegralRight=mean(SNRall());
        end
        if sensitiv ==1
            RISXintegralRight=shiftintegralRightmean/deltaRI;
            RISPeakintegralRight=meanshiftintegralRightall/deltaRI; %Sensitivity if RIS file exists
        end
        
        %% Shift at each cycle of Central Moments data
        
        %Mzero
        
        if CentMom==1 & Peak ==1
            %%Calculate the Shift at each cycle
            %%after the drift correction if "drift" is set to "1"
            for k=1:cycles
                i=k*2-1;
                j=k*2;
                A=Mzerocorrect(Shiftarray(j,1):Shiftarray(j,2));%Partial Data with part 2 of each cycle
                B=Mzerocorrect(Shiftarray(i,1):Shiftarray(i,2));%Partial Data with part 1 of each cycle
                
                Asize=length(A);
                Bsize=length(B);
                if Asize == Bsize
                    shift=A-B;
                else
                    if Asize < Bsize
                        B = B(1,1:Asize);
                    else
                        A = A(1,1:Bsize);
                    end
                    shift=A-B;
                end
                
                SNRnoise=shift-mean(shift);
            SNRSigrms=rms(shift);
            SNRNoirms=rms(SNRnoise);
            SNRall(k)=SNRSigrms/SNRNoirms; %SNR of signal rms divided by noise rms
               
                shiftMzeromean(k,1)=mean(shift); %Shifts of each cycle
                shiftMzeromean(k,2)=std(shift)/sqrt(Shiftarray(j,2)-Shiftarray(j,1)); %STD of each cycle
            end
            
            meanshiftMzeroall(1,1)=mean(shiftMzeromean(1:end-1,1)); %Mean
            meanshiftMzeroall(1,2) = mad(shiftMzeromean(1:end-1,1)); %Mean deviation
            if abs(meanshiftMzeroall(:,1)/meanshiftMzeroall(:,2)) < snrMinNoise
                SNRMzero='Unstable Signal';
            else
                SNRMzero=mean(SNRall());
            end
            if sensitiv ==1
                RISXMzero=shiftMzeromean/deltaRI;
                RISPeakMzero=meanshiftMzeroall/deltaRI; %Sensitivity if RIS file exists
            end
            
            %M1
            
            %%Calculate the Shift at each cycle
            %%after the drift correction if "drift" is set to "1"
            for k=1:cycles
                i=k*2-1;
                j=k*2;
                A=M1correct(Shiftarray(j,1):Shiftarray(j,2));%Partial Data with part 2 of each cycle
                B=M1correct(Shiftarray(i,1):Shiftarray(i,2));%Partial Data with part 1 of each cycle
                
                Asize=length(A);
                Bsize=length(B);
                if Asize == Bsize
                    shift=A-B;
                else
                    if Asize < Bsize
                        B = B(1,1:Asize);
                    else
                        A = A(1,1:Bsize);
                    end
                    shift=A-B;
                end
                
                SNRnoise=shift-mean(shift);
            SNRSigrms=rms(shift);
            SNRNoirms=rms(SNRnoise);
            SNRall(k)=SNRSigrms/SNRNoirms; %SNR of signal rms divided by noise rms
               
                shiftM1mean(k,1)=mean(shift); %Shifts of each cycle
                shiftM1mean(k,2)=std(shift)/sqrt(Shiftarray(j,2)-Shiftarray(j,1)); %STD of each cycle
            end
            
            meanshiftM1all(1,1)=mean(shiftM1mean(1:end-1,1)); %Mean
            meanshiftM1all(1,2) = mad(shiftM1mean(1:end-1,1)); %Mean deviation
            if abs(meanshiftM1all(:,1)/meanshiftM1all(:,2)) < snrMinNoise
                SNRM1='Unstable Signal';
            else
                SNRM1=mean(SNRall());
            end
            if sensitiv ==1
                RISXM1=shiftM1mean/deltaRI;
                RISPeakM1=meanshiftM1all/deltaRI; %Sensitivity if RIS file exists
            end
            
            %M2
            
            %%Calculate the Shift at each cycle
            %%after the drift correction if "drift" is set to "1"
            for k=1:cycles
                i=k*2-1;
                j=k*2;
                A=M2correct(Shiftarray(j,1):Shiftarray(j,2));%Partial Data with part 2 of each cycle
                B=M2correct(Shiftarray(i,1):Shiftarray(i,2));%Partial Data with part 1 of each cycle
                
                Asize=length(A);
                Bsize=length(B);
                if Asize == Bsize
                    shift=A-B;
                else
                    if Asize < Bsize
                        B = B(1,1:Asize);
                    else
                        A = A(1,1:Bsize);
                    end
                    shift=A-B;
                end
                
                SNRnoise=shift-mean(shift);
            SNRSigrms=rms(shift);
            SNRNoirms=rms(SNRnoise);
            SNRall(k)=SNRSigrms/SNRNoirms; %SNR of signal rms divided by noise rms
               
                shiftM2mean(k,1)=mean(shift); %Shifts of each cycle
                shiftM2mean(k,2)=std(shift)/sqrt(Shiftarray(j,2)-Shiftarray(j,1)); %STD of each cycle
            end
            
            meanshiftM2all(1,1)=mean(shiftM2mean(1:end-1,1)); %Mean
            meanshiftM2all(1,2) = mad(shiftM2mean(1:end-1,1)); %Mean deviation
            if abs(meanshiftM2all(:,1)/meanshiftM2all(:,2)) < snrMinNoise
                SNRM2='Unstable Signal';
            else
                SNRM2=mean(SNRall());
            end
            if sensitiv ==1
                RISXM2=shiftM2mean/deltaRI;
                RISPeakM2=meanshiftM2all/deltaRI; %Sensitivity if RIS file exists
            end
            
            %M3
            
            %%Calculate the Shift at each cycle
            %%after the drift correction if "drift" is set to "1"
            for k=1:cycles
                i=k*2-1;
                j=k*2;
                A=M3correct(Shiftarray(j,1):Shiftarray(j,2));%Partial Data with part 2 of each cycle
                B=M3correct(Shiftarray(i,1):Shiftarray(i,2));%Partial Data with part 1 of each cycle
                
                Asize=length(A);
                Bsize=length(B);
                if Asize == Bsize
                    shift=A-B;
                else
                    if Asize < Bsize
                        B = B(1,1:Asize);
                    else
                        A = A(1,1:Bsize);
                    end
                    shift=A-B;
                end
                
                SNRnoise=shift-mean(shift);
            SNRSigrms=rms(shift);
            SNRNoirms=rms(SNRnoise);
            SNRall(k)=SNRSigrms/SNRNoirms; %SNR of signal rms divided by noise rms
               
                shiftM3mean(k,1)=mean(shift); %Shifts of each cycle
                shiftM3mean(k,2)=std(shift)/sqrt(Shiftarray(j,2)-Shiftarray(j,1)); %STD of each cycle
            end
            
            meanshiftM3all(1,1)=mean(shiftM3mean(1:end-1,1)); %Mean
            meanshiftM3all(1,2) = mad(shiftM3mean(1:end-1,1)); %Mean deviation
            if abs(meanshiftM3all(:,1)/meanshiftM3all(:,2)) < snrMinNoise
                SNRM3='Unstable Signal';
            else
                SNRM3=mean(SNRall());
            end
            if sensitiv ==1
                RISXM3=shiftM3mean/deltaRI;
                RISPeakM3=meanshiftM3all/deltaRI; %Sensitivity if RIS file exists
            end
            
            %M4
            
            %%Calculate the Shift at each cycle
            %%after the drift correction if "drift" is set to "1"
            for k=1:cycles
                i=k*2-1;
                j=k*2;
                A=M4correct(Shiftarray(j,1):Shiftarray(j,2));%Partial Data with part 2 of each cycle
                B=M4correct(Shiftarray(i,1):Shiftarray(i,2));%Partial Data with part 1 of each cycle
                
                Asize=length(A);
                Bsize=length(B);
                if Asize == Bsize
                    shift=A-B;
                else
                    if Asize < Bsize
                        B = B(1,1:Asize);
                    else
                        A = A(1,1:Bsize);
                    end
                    shift=A-B;
                end
                
                SNRnoise=shift-mean(shift);
            SNRSigrms=rms(shift);
            SNRNoirms=rms(SNRnoise);
            SNRall(k)=SNRSigrms/SNRNoirms; %SNR of signal rms divided by noise rms
                
                shiftM4mean(k,1)=mean(shift); %Shifts of each cycle
                shiftM4mean(k,2)=std(shift)/sqrt(Shiftarray(j,2)-Shiftarray(j,1)); %STD of each cycle
            end
            
            meanshiftM4all(1,1)=mean(shiftM4mean(1:end-1,1)); %Mean
            meanshiftM4all(1,2) = mad(shiftM4mean(1:end-1,1)); %Mean deviation
            if abs(meanshiftM4all(:,1)/meanshiftM4all(:,2)) < snrMinNoise
                SNRM4='Unstable Signal';
            else
                SNRM4=mean(SNRall());
            end
            if sensitiv ==1
                RISXM4=shiftM4mean/deltaRI;
                RISPeakM4=meanshiftM4all/deltaRI; %Sensitivity if RIS file exists
            end
        end
    end
    %%  Shift at each cycle of Polynomial fit
    %LSPRPeakwl
    if Peak==0
    else
        if cycles==0
        else
            % Calculate the Shift at each cycle
            % after the drift correction if "drift" is set to "1"
            for k=1:cycles
                i=k*2-1;
                j=k*2;
                A=minwlfitncorrect(Shiftarray(j,1):Shiftarray(j,2)); %Partial Data with part 2 of each cycle
                B=minwlfitncorrect(Shiftarray(i,1):Shiftarray(i,2)); %Partial Data with part 1 of each cycle
                
                Asize=length(A);
                Bsize=length(B);
                if Asize == Bsize
                    shift=A-B;
                else
                    if Asize < Bsize
                        B = B(1:Asize,1);
                    else
                        A = A(1:Bsize,1);
                    end
                    shift=A-B;
                end
                
                SNRnoise=shift-mean(shift);
            SNRSigrms=rms(shift);
            SNRNoirms=rms(SNRnoise);
            SNRall(k)=SNRSigrms/SNRNoirms; %SNR of signal rms divided by noise rms
                
                shiftWLmean(k,1)=mean(shift); %Shifts of each cycle
                shiftWLmean(k,2)=std(shift)/sqrt(Shiftarray(j,2)-Shiftarray(j,1)); %STD error of each cycle
            end
            
            meanshiftWLall(1,1) = mean(shiftWLmean(1:end-1,1)); %Mean
            meanshiftWLall(1,2) = mad(shiftWLmean(1:end-1,1)); %Mean deviation
            if abs(meanshiftWLall(:,1)/meanshiftWLall(:,2)) < snrMinNoise
                SNRWL='Unstable Signal';
            else
                SNRWL=mean(SNRall());
            end
            if sensitiv ==1
                RISXwl = shiftWLmean/deltaRI;
                RISPeakwl = meanshiftWLall/deltaRI; %Sensitivity if RIS file exists
                FOM = RISPeakwl/FWHHFOM;
            end
            
            %LSPRPeakTr
            
            %%Calculate the Shift at each 10% (5 cycles, 2 min in each condition)
            %%20 min after the drift correction
            for k=1:cycles
                i=k*2-1;
                j=k*2;
                A=minTfitncorrect(Shiftarray(j,1):Shiftarray(j,2));%Partial Data with part 2 of each cycle
                B=minTfitncorrect(Shiftarray(i,1):Shiftarray(i,2));%Partial Data with part 1 of each cycle
                
                Asize=length(A);
                Bsize=length(B);
                if Asize == Bsize
                    shift=A-B;
                else
                    if Asize < Bsize
                        B = B(1:Asize,1);
                    else
                        A = A(1:Bsize,1);
                    end
                    shift=A-B;
                end
                
                SNRnoise=shift-mean(shift);
            SNRSigrms=rms(shift);
            SNRNoirms=rms(SNRnoise);
            SNRall(k)=SNRSigrms/SNRNoirms; %SNR of signal rms divided by noise rms
                
                shiftTmean(k,1)=mean(shift);
                
                shiftTmean(k,2)=std(shift)/sqrt(Shiftarray(j,2)-Shiftarray(j,1));
            end
            
            
            
            
            meanshiftTall(1,1)=mean(shiftTmean(1:end-1,1)); %Mean
            meanshiftTall(1,2) = mad(shiftTmean(1:end-1,1)); %Mean deviation
            if abs(meanshiftTall(:,1)/meanshiftTall(:,2)) < snrMinNoise
                SNRT='Unstable Signal';
            else
                SNRT=mean(SNRall());
            end
            if sensitiv ==1
                RISXT=shiftTmean/deltaRI;
                RISPeakT=meanshiftTall/deltaRI;
            end
            
            
            
            
            %% Calculate the Shift at each 10% (5 cycles, 2 min in each condition)
            %%20 min after the drift correction
            for k=1:cycles
                i=k*2-1;
                j=k*2;
                A=FWHHcorrect(Shiftarray(j,1):Shiftarray(j,2));%Partial Data with part 2 of each cycle
                B=FWHHcorrect(Shiftarray(i,1):Shiftarray(i,2));%Partial Data with part 1 of each cycle
                
                Asize=length(A);
                Bsize=length(B);
                if Asize == Bsize
                    shift=A-B;
                else
                    if Asize < Bsize
                        B = B(1:Asize,1);
                    else
                        A = A(1:Bsize,1);
                    end
                    shift=A-B;
                end
                
                SNRnoise=shift-mean(shift);
            SNRSigrms=rms(shift);
            SNRNoirms=rms(SNRnoise);
            SNRall(k)=SNRSigrms/SNRNoirms; %SNR of signal rms divided by noise rms
                
                
                
                shiftFWHHmean(k,1)=mean(shift);
                shiftFWHHmean(k,2)=std(shift)/sqrt(Shiftarray(j,2)-Shiftarray(j,1));
            end
            
            
            meanshiftFWHHall(1,1)=mean(shiftFWHHmean(1:end-1,1)); %Mean
            meanshiftFWHHall(1,2) = mad(shiftFWHHmean(1:end-1,1)); %Mean deviation
            if abs(meanshiftFWHHall(:,1)/meanshiftFWHHall(:,2)) < snrMinNoise
                SNRFWHH='Unstable Signal';
            else
                SNRFWHH=mean(SNRall());
            end
            
            if sensitiv ==1
                RISXFWHH=shiftFWHHmean/deltaRI;
                RISPeakFWHH=meanshiftFWHHall/deltaRI;
            end
            
            %LSPRLWHH
            
            %%Calculate the Shift at each 10% (5 cycles, 2 min in each condition)
            %%20 min after the drift correction
            for k=1:cycles
                i=k*2-1;
                j=k*2;
                A=LWHHcorrect(Shiftarray(j,1):Shiftarray(j,2));%Partial Data with part 2 of each cycle
                B=LWHHcorrect(Shiftarray(i,1):Shiftarray(i,2));%Partial Data with part 1 of each cycle
                
                Asize=length(A);
                Bsize=length(B);
                if Asize == Bsize
                    shift=A-B;
                else
                    if Asize < Bsize
                        B = B(1:Asize,1);
                    else
                        A = A(1:Bsize,1);
                    end
                    shift=A-B;
                end
                
                SNRnoise=shift-mean(shift);
            SNRSigrms=rms(shift);
            SNRNoirms=rms(SNRnoise);
            SNRall(k)=SNRSigrms/SNRNoirms; %SNR of signal rms divided by noise rms
                
                                
                shiftLWHHmean(k,1)=mean(shift);
                shiftLWHHmean(k,2)=std(shift)/sqrt(Shiftarray(j,2)-Shiftarray(j,1));
            end
            
            
            meanshiftLWHHall(1,1)=mean(shiftLWHHmean(1:end-1,1)); %Mean
            meanshiftLWHHall(1,2) = mad(shiftLWHHmean(1:end-1,1)); %Mean deviation
            if abs(meanshiftLWHHall(:,1)/meanshiftLWHHall(:,2)) < snrMinNoise
                SNRLWHH='Unstable Signal';
            else
                SNRLWHH=mean(SNRall());
            end
            
            if sensitiv ==1
                RISXLWHH=shiftLWHHmean/deltaRI;
                RISPeakLWHH=meanshiftLWHHall/deltaRI;
            end
            
            %LSPRRWHH
            
            %%Calculate the Shift at each 10% (5 cycles, 2 min in each condition)
            %%20 min after the drift correction
            for k=1:cycles
                i=k*2-1;
                j=k*2;
                A=RWHHcorrect(Shiftarray(j,1):Shiftarray(j,2));%Partial Data with part 2 of each cycle
                B=RWHHcorrect(Shiftarray(i,1):Shiftarray(i,2));%Partial Data with part 1 of each cycle
                
                Asize=length(A);
                Bsize=length(B);
                if Asize == Bsize
                    shift=A-B;
                else
                    if Asize < Bsize
                        B = B(1:Asize,1);
                    else
                        A = A(1:Bsize,1);
                    end
                    shift=A-B;
                end
                
                SNRnoise=shift-mean(shift);
            SNRSigrms=rms(shift);
            SNRNoirms=rms(SNRnoise);
            SNRall(k)=SNRSigrms/SNRNoirms; %SNR of signal rms divided by noise rms
                
                               
                shiftRWHHmean(k,1)=mean(shift);
                shiftRWHHmean(k,2)=std(shift)/sqrt(Shiftarray(j,2)-Shiftarray(j,1));
            end
            
            
            meanshiftRWHHall(1,1)=mean(shiftRWHHmean(1:end-1,1)); %Mean
            meanshiftRWHHall(1,2) = mad(shiftRWHHmean(1:end-1,1)); %Mean deviation
            if abs(meanshiftRWHHall(:,1)/meanshiftRWHHall(:,2)) < snrMinNoise
                SNRRWHH='Unstable Signal';
            else
                SNRRWHH=mean(SNRall());
            end
            
            if sensitiv ==1
                RISXRWHH=shiftRWHHmean/deltaRI;
                RISPeakRWHH=meanshiftRWHHall/deltaRI;
            end
        end
        %         end
        
        
        
        %% Over time monitoring
        %% Polinomial fit over time monitorization
        
        figure(1)
        
        subplot(5,3,[1 2]);
        
        plot(xx,minwlfitn)
        title('LSPR band minimum Wavelength', 'Fontsize',14);
        axis([0 maxTime -inf inf])
        grid on
        hold on
        
        colorMin=min(minwlfitn);
        colorMax=max(minwlfitn);
        
        if cycles==0
        else
            for n=1:cycles
                k=n*4-3;
                i=k+1;
                j=k+2;
                l=k+3;
                v2(k,1)=colorTime(n,1);
                v2(k,2)=colorMin;
                v2(i,1)=colorTime(n,2);
                v2(i,2)=colorMin;
                v2(j,1)=colorTime(n,2);
                v2(j,2)=colorMax;
                v2(l,1)=colorTime(n,1);
                v2(l,2)=colorMax;
            end
            
            for n=1:cycles
                k=n*4-3;
                i=k+1;
                j=k+2;
                l=k+3;
                f2(n,1)=k;
                f2(n,2)=i;
                f2(n,3)=j;
                f2(n,4)=l;
            end
            
            patch('Faces',f2,'Vertices',v2,'FaceColor',[0.9 1 0.9],'FaceAlpha',.5,'EdgeColor','none')
        end
        hold off
        
        subplot(5,3,[10 11]);
        
        yyaxis left
        plot(xx,LWHH,'b','LineWidth',1.5)
        hold on
        yyaxis right
        plot(xx,RWHH,'r','LineWidth',1.5)
        title('LSPR Band Width L&R', 'Fontsize',14);
        axis([0 maxTime -inf inf])
        xlabel('Time (s)', 'Fontsize',14);
        grid on
        hold on
        
        
        colorMin=min(RWHH);
        colorMax=max(RWHH);
        if cycles==0
        else
            for n=1:cycles
                k=n*4-3;
                i=k+1;
                j=k+2;
                l=k+3;
                v2(k,1)=colorTime(n,1);
                v2(k,2)=colorMin;
                v2(i,1)=colorTime(n,2);
                v2(i,2)=colorMin;
                v2(j,1)=colorTime(n,2);
                v2(j,2)=colorMax;
                v2(l,1)=colorTime(n,1);
                v2(l,2)=colorMax;
            end
            
            for n=1:cycles
                k=n*4-3;
                i=k+1;
                j=k+2;
                l=k+3;
                f2(n,1)=k;
                f2(n,2)=i;
                f2(n,3)=j;
                f2(n,4)=l;
            end
            
            patch('Faces',f2,'Vertices',v2,'FaceColor',[0.9 1 0.9],'FaceAlpha',.5,'EdgeColor','none')
            lgd=legend('Left', 'Right',testmaterial);
            set(lgd,'FontSize',6, 'color', 'none','Box','off');
        end
        hold off
        
        %
        subplot(5,3,[7 8]);
        plot(xx,FWHH,'k')
        title('LSPR Band FWHH', 'Fontsize',14);
        axis([0 maxTime -inf inf])
        grid on
        hold on
        
        colorMin=min(FWHH);
        colorMax=max(FWHH);
        if cycles==0
        else
            for n=1:cycles
                k=n*4-3;
                i=k+1;
                j=k+2;
                l=k+3;
                v2(k,1)=colorTime(n,1);
                v2(k,2)=colorMin;
                v2(i,1)=colorTime(n,2);
                v2(i,2)=colorMin;
                v2(j,1)=colorTime(n,2);
                v2(j,2)=colorMax;
                v2(l,1)=colorTime(n,1);
                v2(l,2)=colorMax;
            end
            
            for n=1:cycles
                k=n*4-3;
                i=k+1;
                j=k+2;
                l=k+3;
                f2(n,1)=k;
                f2(n,2)=i;
                f2(n,3)=j;
                f2(n,4)=l;
            end
            
            patch('Faces',f2,'Vertices',v2,'FaceColor',[0.9 1 0.9],'FaceAlpha',.5,'EdgeColor','none')
        end
        hold off
        
        %         end
        
        if pptt ==3
            %Temperature Variation
            pTemperature=subplot(5,3,[13 14]);
            yyaxis left
            plot(pTemperature,tempTime,tempData,'r')
            axis([0 maxTime -inf inf])
            ylabel('T (ºC)', 'Fontsize',12)
            grid on
            
            hold on
            %Pressure Variation
            pPressure=subplot(5,3,[13 14]);
            yyaxis right
            plot(pPressure,timepressure,pressure,'k')
            lgd=legend(testmaterial,'Temperature','Pressure');
            set(lgd,'FontSize',6, 'color', 'none','Box','off');
            axis([0 maxTime -inf inf])
            xlabel('Time (s)', 'Fontsize',14);
            ylabel('P (mbar)', 'Fontsize',12)
            grid on
            hold off
            
        elseif pptt ==1
            %Pressure Variation
            pPressure=subplot(5,3,[13 14]);
            colorMin=min(pressure);
            colorMax=max(pressure);
            
            if cycles==0
            else
                for n=1:cycles
                    k=n*4-3;
                    i=k+1;
                    j=k+2;
                    l=k+3;
                    v2(k,1)=colorTime(n,1);
                    v2(k,2)=colorMin;
                    v2(i,1)=colorTime(n,2);
                    v2(i,2)=colorMin;
                    v2(j,1)=colorTime(n,2);
                    v2(j,2)=colorMax;
                    v2(l,1)=colorTime(n,1);
                    v2(l,2)=colorMax;
                end
                
                for n=1:cycles
                    k=n*4-3;
                    i=k+1;
                    j=k+2;
                    l=k+3;
                    f2(n,1)=k;
                    f2(n,2)=i;
                    f2(n,3)=j;
                    f2(n,4)=l;
                end
                
                patch('Faces',f2,'Vertices',v2,'FaceColor',[0.9 0.9 1],'FaceAlpha',.5,'EdgeColor','none')
            end
            hold on
            plot(pPressure,timepressure,pressure,'k')
            if cycles==0
            else
                lgd=legend(testmaterial,'Pressure');
                set(lgd,'FontSize',6, 'color', 'none','Box','off');
            end
            
            
            axis([0 maxTime -inf inf])
            xlabel('Time (s)', 'Fontsize',14);
            ylabel('P (mbar)', 'Fontsize',12)
            grid on
            hold off
            
            
        elseif pptt ==2
            %Temperature Variation
            pTemperature=subplot(5,3,[13 14]);
            plot(pTemperature,tempTime,tempData,'r')
            lgd=legend(testmaterial,'Temperature');
            set(lgd,'FontSize',6, 'color', 'none','Box','off');
            axis([0 maxTime -inf inf])
            ylabel('T (ºC)', 'Fontsize',12)
            grid on
        else
        end
        
        subplot(5,3,[4 5]);
        plot(xx,minTfitn)
        title('LSPR band minimum Transmittance', 'Fontsize',14);
        axis([0 maxTime -inf inf])
        grid on
        hold on
        
        colorMin=min(minTfitn);
        colorMax=max(minTfitn);
        
        if cycles==0
        else
            for n=1:cycles
                k=n*4-3;
                i=k+1;
                j=k+2;
                l=k+3;
                v2(k,1)=colorTime(n,1);
                v2(k,2)=colorMin;
                v2(i,1)=colorTime(n,2);
                v2(i,2)=colorMin;
                v2(j,1)=colorTime(n,2);
                v2(j,2)=colorMax;
                v2(l,1)=colorTime(n,1);
                v2(l,2)=colorMax;
            end
            
            for n=1:cycles
                k=n*4-3;
                i=k+1;
                j=k+2;
                l=k+3;
                f2(n,1)=k;
                f2(n,2)=i;
                f2(n,3)=j;
                f2(n,4)=l;
            end
            
            patch('Faces',f2,'Vertices',v2,'FaceColor',[0.9 1 0.9],'FaceAlpha',.5,'EdgeColor','none')
        end
        hold off
        
        
        
        subplot(5,3,[3 6 9 12 15])
        plot(wBandcutPeak1,TBandcutPeak1(:,end),'*')
        hold on
        fplot(BandfitPeakfunctn,'r')
        hold off
        set(gca,'fontsize',14);
        axis([WLminL WLminR 0 100]);
        ylabel('Transmittance (%)', 'Fontsize',14);
        xlabel('Wavelength (nm)', 'Fontsize',14);
        
        
        grid on
        set(gcf, 'Position', get(0, 'Screensize'));
        
        
        savefigname=fullfile(path,'results',strcat(basefoldername,'_Poly9LSPRPeak_',smethod,num2str(smooth),'.png'));
        saveas(gcf,savefigname)
        savefigname=fullfile(path,'results',strcat(basefoldername,'_Poly9LSPRPeak_',smethod,num2str(smooth),'.svg'));
        print('-f1',savefigname,'-dsvg')
        savefigname=fullfile(path,'results',strcat(basefoldername,'_Poly9LSPRPeak_',smethod,num2str(smooth),'.fig'));
        saveas(gcf,savefigname)
    end
    
    %% Central moments over time monitorization
    if CentMom==1 & Peak ==1
        figure(4)
        
        subplot(24,3,[1 2 4 5 7 8]);
        plot(xx,Mzero)
        title('M0 - spectra integral', 'Fontsize',14);
        axis([0 maxTime -inf inf])
        grid on
        hold on
        
        colorMin=min(Mzero);
        colorMax=max(Mzero);
        
        if cycles==0
        else
            for n=1:cycles
                k=n*4-3;
                i=k+1;
                j=k+2;
                l=k+3;
                v2(k,1)=colorTime(n,1);
                v2(k,2)=colorMin;
                v2(i,1)=colorTime(n,2);
                v2(i,2)=colorMin;
                v2(j,1)=colorTime(n,2);
                v2(j,2)=colorMax;
                v2(l,1)=colorTime(n,1);
                v2(l,2)=colorMax;
            end
            
            for n=1:cycles
                k=n*4-3;
                i=k+1;
                j=k+2;
                l=k+3;
                f2(n,1)=k;
                f2(n,2)=i;
                f2(n,3)=j;
                f2(n,4)=l;
            end
            
            patch('Faces',f2,'Vertices',v2,'FaceColor',[0.9 1 0.9],'FaceAlpha',.5,'EdgeColor','none')
            
        end
        hold off
        
        subplot(24,3,[13 14 16 17 19 20]);
        plot(xx,M1)
        title('M1 - Expectation value', 'Fontsize',14);
        axis([0 maxTime -inf inf])
        grid on
        
        hold on
        
        colorMin=min(M1);
        colorMax=max(M1);
        if cycles==0
        else
            for n=1:cycles
                k=n*4-3;
                i=k+1;
                j=k+2;
                l=k+3;
                v2(k,1)=colorTime(n,1);
                v2(k,2)=colorMin;
                v2(i,1)=colorTime(n,2);
                v2(i,2)=colorMin;
                v2(j,1)=colorTime(n,2);
                v2(j,2)=colorMax;
                v2(l,1)=colorTime(n,1);
                v2(l,2)=colorMax;
            end
            
            for n=1:cycles
                k=n*4-3;
                i=k+1;
                j=k+2;
                l=k+3;
                f2(n,1)=k;
                f2(n,2)=i;
                f2(n,3)=j;
                f2(n,4)=l;
            end
            
            patch('Faces',f2,'Vertices',v2,'FaceColor',[0.9 1 0.9],'FaceAlpha',.5,'EdgeColor','none')
            
        end
        hold off
        
        subplot(24,3,[25 26 28 29 31 32]);
        plot(xx,M2)
        title('M2 - Variance', 'Fontsize',14);
        axis([0 maxTime -inf inf])
        grid on
        
        hold on
        
        colorMin=min(M2);
        colorMax=max(M2);
        
        if cycles==0
        else
            for n=1:cycles
                k=n*4-3;
                i=k+1;
                j=k+2;
                l=k+3;
                v2(k,1)=colorTime(n,1);
                v2(k,2)=colorMin;
                v2(i,1)=colorTime(n,2);
                v2(i,2)=colorMin;
                v2(j,1)=colorTime(n,2);
                v2(j,2)=colorMax;
                v2(l,1)=colorTime(n,1);
                v2(l,2)=colorMax;
            end
            
            for n=1:cycles
                k=n*4-3;
                i=k+1;
                j=k+2;
                l=k+3;
                f2(n,1)=k;
                f2(n,2)=i;
                f2(n,3)=j;
                f2(n,4)=l;
            end
            
            patch('Faces',f2,'Vertices',v2,'FaceColor',[0.9 1 0.9],'FaceAlpha',.5,'EdgeColor','none')
            
        end
        hold off
        
        subplot(24,3,[37 38 40 41 43 44]);
        plot(xx,M3)
        title('M3 - Skewness', 'Fontsize',14);
        axis([0 maxTime -inf inf])
        grid on
        
        hold on
        
        colorMin=min(M3);
        colorMax=max(M3);
        if cycles==0
        else
            for n=1:cycles
                k=n*4-3;
                i=k+1;
                j=k+2;
                l=k+3;
                v2(k,1)=colorTime(n,1);
                v2(k,2)=colorMin;
                v2(i,1)=colorTime(n,2);
                v2(i,2)=colorMin;
                v2(j,1)=colorTime(n,2);
                v2(j,2)=colorMax;
                v2(l,1)=colorTime(n,1);
                v2(l,2)=colorMax;
            end
            
            for n=1:cycles
                k=n*4-3;
                i=k+1;
                j=k+2;
                l=k+3;
                f2(n,1)=k;
                f2(n,2)=i;
                f2(n,3)=j;
                f2(n,4)=l;
            end
            
            patch('Faces',f2,'Vertices',v2,'FaceColor',[0.9 1 0.9],'FaceAlpha',.5,'EdgeColor','none')
            
        end
        hold off
        
        subplot(24,3,[49 50 52 53 55 56]);
        plot(xx,M4)
        title('M4 - Kurtosis', 'Fontsize',14);
        axis([0 maxTime -inf inf])
        xlabel('Time (s)', 'Fontsize',14);
        grid on
        
        hold on
        
        colorMin=min(M4);
        colorMax=max(M4);
        if cycles==0
        else
            for n=1:cycles
                k=n*4-3;
                i=k+1;
                j=k+2;
                l=k+3;
                v2(k,1)=colorTime(n,1);
                v2(k,2)=colorMin;
                v2(i,1)=colorTime(n,2);
                v2(i,2)=colorMin;
                v2(j,1)=colorTime(n,2);
                v2(j,2)=colorMax;
                v2(l,1)=colorTime(n,1);
                v2(l,2)=colorMax;
            end
            
            for n=1:cycles
                k=n*4-3;
                i=k+1;
                j=k+2;
                l=k+3;
                f2(n,1)=k;
                f2(n,2)=i;
                f2(n,3)=j;
                f2(n,4)=l;
            end
            
            patch('Faces',f2,'Vertices',v2,'FaceColor',[0.9 1 0.9],'FaceAlpha',.5,'EdgeColor','none')
            
        end
        hold off
        
        if pptt ==3
            %Temperature Variation
            pTemperature=subplot(24,3,[61 62 64 65 67 68]);
            yyaxis left
            plot(pTemperature,tempTime,tempData,'r')
            axis([0 maxTime -inf inf])
            ylabel('T (ºC)', 'Fontsize',12)
            grid on
            
            hold on
            %Pressure Variation
            pPressure=subplot(24,3,[61 62 64 65 67 68]);
            yyaxis right
            plot(pPressure,timepressure,pressure,'k')
            lgd=legend('Temperature','Pressure');
            set(lgd,'FontSize',6, 'color', 'none','Box','off');
            axis([0 maxTime -inf inf])
            xlabel('Time (s)', 'Fontsize',12);
            ylabel('P (mbar)', 'Fontsize',12)
            grid on
            hold off
            
        elseif pptt ==1
            %Pressure Variation
            pPressure=subplot(24,3,[61 62 64 65 67 68]);
            plot(pPressure,timepressure,pressure,'k')
            lgd=legend('Pressure');
            set(lgd,'FontSize',6, 'color', 'none','Box','off');
            axis([0 maxTime -inf inf])
            xlabel('Time (s)', 'Fontsize',12);
            ylabel('P (mbar)', 'Fontsize',12)
            grid on
            
        elseif pptt ==2
            %Temperature Variation
            pTemperature=subplot(24,3,[61 62 64 65 67 68]);
            plot(pTemperature,tempTime,tempData,'r')
            lgd=legend('Temperature');
            set(lgd,'FontSize',6, 'color', 'none','Box','off');
            axis([0 maxTime -inf inf])
            ylabel('T (ºC)', 'Fontsize',12)
            grid on
        else
        end
                
        subplot(24,3,[39 42 45 48 51 54 57 60 63 66 69 72]);
        plot(wBandcutPeak1CMom,pseusopdf(:,2),'*')
        set(gca,'fontsize',14);
        axis([-inf inf -inf inf]);
        ylabel('Normalized (1-Transmittance)', 'Fontsize',14);
        xlabel('Normalized Wavelength', 'Fontsize',14);
        
        
        grid on
        
        subplot(24,3,[3 6 9 12 15 18 21 24 27 30]);
        plot(wBandcutPeak1,TBandcutPeak1(:,2),'*')
        
        set(gca,'fontsize',14);
        axis([WLminL WLminR 0 100]);
        title(systemname, 'Fontsize',16);
        ylabel('Transmittance (%)', 'Fontsize',14);
        xlabel('Wavelength (nm)', 'Fontsize',14);
        
        grid on
        
        set(gcf, 'Position', get(0, 'Screensize'));
        
        %%
        %%Date Limit
        if Date > DateLimit
            CreateStruct.Interpreter = 'tex';
            CreateStruct.WindowStyle = 'modal';
            msgbox('\fontsize{14} Demonstration period finished. Please contact Mr. M.S.Rodrigues - marcopsr@gmail.com','Finished',CreateStruct);
            return
        end
        %%
        
        savefigname=fullfile(path,'results',strcat(basefoldername,'_CentralMoments_',smethod,num2str(smooth),'.png'));
        saveas(gcf,savefigname)
        savefigname=fullfile(path,'results',strcat(basefoldername,'_CentralMoments_',smethod,num2str(smooth),'.svg'));
        print('-f4',savefigname,'-dsvg')
        %         savefigname=fullfile('VariationFigures',strcat(basefoldername,'_CentralMoments_',smethod,num2str(smooth),'.png'));
        %         saveas(gcf,savefigname)
        savefigname=fullfile(path,'results',strcat(basefoldername,'_CentralMoments_',smethod,num2str(smooth),'.fig'));
        saveas(gcf,savefigname)
    else
        
    end
    %%
    
    
    
    %% Raw data Figure overtime monitorization
    
    figure(5)
    
    subplot(5,3,[1 2]);
    plot(xx,deltaTwl(:,1:5))
    lgd=legend('450 nm', '500 nm', '550 nm', '600 nm', '650 nm');
    set(lgd,'FontSize',6, 'color', 'none','Box','off');
    title('Transmittance shift discrete wavelengths', 'Fontsize',14);
    axis([0 maxTime -inf inf])
    grid on
    
    hold on
    
    colorMin=min(min(deltaTwl(:,1:5)));
    colorMax=max(max(deltaTwl(:,1:5)));
    
    if cycles==0
    else
        for n=1:cycles
            k=n*4-3;
            i=k+1;
            j=k+2;
            l=k+3;
            v2(k,1)=colorTime(n,1);
            v2(k,2)=colorMin;
            v2(i,1)=colorTime(n,2);
            v2(i,2)=colorMin;
            v2(j,1)=colorTime(n,2);
            v2(j,2)=colorMax;
            v2(l,1)=colorTime(n,1);
            v2(l,2)=colorMax;
        end
        
        for n=1:cycles
            k=n*4-3;
            i=k+1;
            j=k+2;
            l=k+3;
            f2(n,1)=k;
            f2(n,2)=i;
            f2(n,3)=j;
            f2(n,4)=l;
        end
        
        patch('Faces',f2,'Vertices',v2,'FaceColor',[0.9 1 0.9],'FaceAlpha',.5,'EdgeColor','none')
        
    end
    hold off
    
    subplot(5,3,[4 5]);
    plot(xx,deltaTwl(:,6:9))
    lgd=legend('700 nm', '750 nm', '800 nm', '850 nm');
    set(lgd,'FontSize',6, 'color', 'none','Box','off');
    title('Transmittance shift discrete wavelengths', 'Fontsize',14);
    axis([0 maxTime -inf inf])
    grid on
    
    hold on
    
    colorMin=min(min(deltaTwl(:,6:9)));
    colorMax=max(max(deltaTwl(:,6:9)));
    
    if cycles==0
    else
        for n=1:cycles
            k=n*4-3;
            i=k+1;
            j=k+2;
            l=k+3;
            v2(k,1)=colorTime(n,1);
            v2(k,2)=colorMin;
            v2(i,1)=colorTime(n,2);
            v2(i,2)=colorMin;
            v2(j,1)=colorTime(n,2);
            v2(j,2)=colorMax;
            v2(l,1)=colorTime(n,1);
            v2(l,2)=colorMax;
        end
        
        for n=1:cycles
            k=n*4-3;
            i=k+1;
            j=k+2;
            l=k+3;
            f2(n,1)=k;
            f2(n,2)=i;
            f2(n,3)=j;
            f2(n,4)=l;
        end
        
        patch('Faces',f2,'Vertices',v2,'FaceColor',[0.9 1 0.9],'FaceAlpha',.5,'EdgeColor','none')
        
    end
    
    hold off
    
    %LSPR band integral plot
    
    subplot(5,3,[7 8]);
    plot(xx,realDataIntegral)
    title('LSPR band Integral', 'Fontsize',14);
    axis([0 maxTime -inf inf])
    grid on
    
    hold on
    
    colorMin=min(realDataIntegral);
    colorMax=max(realDataIntegral);
    
    if cycles==0
    else
        for n=1:cycles
            k=n*4-3;
            i=k+1;
            j=k+2;
            l=k+3;
            v2(k,1)=colorTime(n,1);
            v2(k,2)=colorMin;
            v2(i,1)=colorTime(n,2);
            v2(i,2)=colorMin;
            v2(j,1)=colorTime(n,2);
            v2(j,2)=colorMax;
            v2(l,1)=colorTime(n,1);
            v2(l,2)=colorMax;
        end
        
        for n=1:cycles
            k=n*4-3;
            i=k+1;
            j=k+2;
            l=k+3;
            f2(n,1)=k;
            f2(n,2)=i;
            f2(n,3)=j;
            f2(n,4)=l;
        end
        
        patch('Faces',f2,'Vertices',v2,'FaceColor',[0.9 1 0.9],'FaceAlpha',.5,'EdgeColor','none')
        
    end
    hold off
    
    
    %Left LSPR band integral plot
    subplot(5,3,[10 11]);
    yyaxis left
    plot(xx,integralLeft)
    axis([0 maxTime -inf inf])
    hold on
    
    %Right LSPR band integral plot
    
    yyaxis right
    plot(xx,integralRight)
    title('LSPR band Integral L&R', 'Fontsize',14);
    axis([0 maxTime -inf inf])
    xlabel('Time (s)', 'Fontsize',14);
    hold off
    grid on
    
    hold on
    
    colorMin=min(integralRight);
    colorMax=max(integralRight);
    
    if cycles==0
    else
        for n=1:cycles
            k=n*4-3;
            i=k+1;
            j=k+2;
            l=k+3;
            v2(k,1)=colorTime(n,1);
            v2(k,2)=colorMin;
            v2(i,1)=colorTime(n,2);
            v2(i,2)=colorMin;
            v2(j,1)=colorTime(n,2);
            v2(j,2)=colorMax;
            v2(l,1)=colorTime(n,1);
            v2(l,2)=colorMax;
        end
        
        for n=1:cycles
            k=n*4-3;
            i=k+1;
            j=k+2;
            l=k+3;
            f2(n,1)=k;
            f2(n,2)=i;
            f2(n,3)=j;
            f2(n,4)=l;
        end
        
        patch('Faces',f2,'Vertices',v2,'FaceColor',[0.9 1 0.9],'FaceAlpha',.5,'EdgeColor','none')
        
    end
    hold off
    
    if pptt ==3
        %Temperature Variation
        pTemperature=subplot(5,3,[13 14]);
        yyaxis left
        plot(pTemperature,tempTime,tempData,'r')
        axis([0 maxTime -inf inf])
        ylabel('T (ºC)', 'Fontsize',12)
        grid on
        
        hold on
        %Pressure Variation
        pPressure=subplot(5,3,[13 14]);
        yyaxis right
        plot(pPressure,timepressure,pressure,'k')
        lgd=legend('Temperature','Pressure');
        set(lgd,'FontSize',6, 'color', 'none','Box','off');
        axis([0 maxTime -inf inf])
        xlabel('Time (s)', 'Fontsize',12);
        ylabel('P (mbar)', 'Fontsize',12)
        grid on
        hold off
        
    elseif pptt ==1
        %Pressure Variation
        pPressure=subplot(5,3,[13 14]);
        plot(pPressure,timepressure,pressure,'k')
        lgd=legend('Pressure');
        set(lgd,'FontSize',6, 'color', 'none','Box','off');
        axis([0 maxTime -inf inf])
        xlabel('Time (s)', 'Fontsize',12);
        ylabel('P (mbar)', 'Fontsize',12)
        grid on
        
    elseif pptt ==2
        %Temperature Variation
        pTemperature=subplot(5,3,[13 14]);
        plot(pTemperature,tempTime,tempData,'r')
        lgd=legend('Temperature');
        set(lgd,'FontSize',6, 'color', 'none','Box','off');
        axis([0 maxTime -inf inf])
        ylabel('T (ºC)', 'Fontsize',12)
        grid on
    else
    end
    
    subplot(5,3,[3 6 9]);
    plot(wavelength,T(:,2),'r')
    set(gca,'fontsize',14);
    axis([350 850 0 100]);
    title(systemname, 'Fontsize',16);
    ylabel('Transmittance (%)', 'Fontsize',14);
    hold on
    plot(wcut,Tcut,'g')
    grid on
    hold off
    
    %OTC
    if cycles==0
    else
        subplot(5,3,[12 15]);
        plot(wavelength,OTC)
        set(gca,'fontsize',14);
        lgd=legend(legendOTC);
        set(lgd,'FontSize',10, 'color', 'none','Box','off');
        axis([350 850 -inf inf]);
        ylabel('OTC', 'Fontsize',14);
        xlabel('Wavelength (nm)', 'Fontsize',14);
        grid on
    end
    set(gcf, 'Position', get(0, 'Screensize'));
    
    savefigname=fullfile(path,'results',strcat(basefoldername,'_rawData_',smethod,num2str(smooth),'.png'));
    saveas(gcf,savefigname)
    savefigname=fullfile(path,'results',strcat(basefoldername,'_rawData_',smethod,num2str(smooth),'.svg'));
    print('-f5',savefigname,'-dsvg')
    
    savefigname=fullfile(path,'results',strcat(basefoldername,'_rawData_',smethod,num2str(smooth),'.fig'));
    saveas(gcf,savefigname)
    
    %% Bar plots
    if cycles==0
    else
        
        %% Raw Data bar plots
        
        figure(8)
        
        % 450 nm
        subplot(2,6,1);
        bar(1:cycles,shift450mean(:,1))
        hold on
        errorbar(1:cycles,shift450mean(:,1), shift450mean(:,2),'.','LineWidth',2)
        hold off
        title({strcat('Shift = ',num2str(meanshift450all(1,1),'%0.3f'),'\pm',num2str(meanshift450all(1,2),'%0.3f'),' pp'),strcat('SNR =',32,num2str(SNR450,'%0.1f'))},'FontSize',12,'units','normalized');
        
        ylabel('450 nm Shift (pp)', 'Fontsize',14);
        xlabel('Cycle number', 'Fontsize',14);
        grid on
        
        if sensitiv ==1
            %RIS : RISXwl    RISPeakwl
            subplot(2,6,7);
            bar(1:cycles,RISX450(:,1))
            hold on
            errorbar(1:cycles,RISX450(:,1), RISX450(:,2),'.','LineWidth',2)
            hold off
            
            ylabel('450 nm Shift (pp/RIU)', 'Fontsize',14);
            xlabel('Cycle number', 'Fontsize',14);
            grid on
            title(strcat('RIS = ',num2str(RISPeak450(1,1),'%0.1f'),'\pm',num2str(RISPeak450(1,2),'%0.1f'),' pp/RIU'),'FontSize',12,'units','normalized');%strcat('RIS = ',num2str(RISPeakwl,'%0.3f'),' nm/RIU'), 'Fontsize',14);
        end
        
        % 500 nm
        subplot(2,6,2);
        bar(1:cycles,shift500mean(:,1))
        hold on
        errorbar(1:cycles,shift500mean(:,1), shift500mean(:,2),'.','LineWidth',2)
        hold off
        title({strcat('Shift = ',num2str(meanshift500all(1,1),'%0.3f'),'\pm',num2str(meanshift500all(1,2),'%0.3f'),' pp'),strcat('SNR =',32,num2str(SNR500,'%0.1f'))},'FontSize',12,'units','normalized');
        
        ylabel('500 nm Shift (pp)', 'Fontsize',14);
        xlabel('Cycle number', 'Fontsize',14);
        grid on
        
        if sensitiv ==1
            %RIS : RISXwl    RISPeakwl
            subplot(2,6,8);
            bar(1:cycles,RISX500(:,1))
            hold on
            errorbar(1:cycles,RISX500(:,1), RISX500(:,2),'.','LineWidth',2)
            hold off
            
            ylabel('500 nm Shift (pp/RIU)', 'Fontsize',14);
            xlabel('Cycle number', 'Fontsize',14);
            grid on
            title(strcat('RIS = ',num2str(RISPeak500(1,1),'%0.1f'),'\pm',num2str(RISPeak500(1,2),'%0.1f'),' pp/RIU'),'FontSize',12,'units','normalized');%strcat('RIS = ',num2str(RISPeakwl,'%0.3f'),' nm/RIU'), 'Fontsize',14);
        end
        
        % 550 nm
        subplot(2,6,3);
        bar(1:cycles,shift550mean(:,1))
        hold on
        errorbar(1:cycles,shift550mean(:,1), shift550mean(:,2),'.','LineWidth',2)
        hold off
        title({strcat('Shift = ',num2str(meanshift550all(1,1),'%0.3f'),'\pm',num2str(meanshift550all(1,2),'%0.3f'),' pp'),strcat('SNR =',32,num2str(SNR550,'%0.1f'))},'FontSize',12,'units','normalized');
        
        ylabel('550 nm Shift (pp)', 'Fontsize',14);
        xlabel('Cycle number', 'Fontsize',14);
        grid on
        
        if sensitiv ==1
            %RIS : RISXwl    RISPeakwl
            subplot(2,6,9);
            bar(1:cycles,RISX550(:,1))
            hold on
            errorbar(1:cycles,RISX550(:,1), RISX550(:,2),'.','LineWidth',2)
            hold off
            
            ylabel('550 nm Shift (pp/RIU)', 'Fontsize',14);
            xlabel('Cycle number', 'Fontsize',14);
            grid on
            title(strcat('RIS = ',num2str(RISPeak550(1,1),'%0.1f'),'\pm',num2str(RISPeak550(1,2),'%0.1f'),' pp/RIU'),'FontSize',12,'units','normalized');%strcat('RIS = ',num2str(RISPeakwl,'%0.3f'),' nm/RIU'), 'Fontsize',14);
        end
        
        % 600 nm
        subplot(2,6,4);
        bar(1:cycles,shift600mean(:,1))
        hold on
        errorbar(1:cycles,shift600mean(:,1), shift600mean(:,2),'.','LineWidth',2)
        hold off
        title({strcat('Shift = ',num2str(meanshift600all(1,1),'%0.3f'),'\pm',num2str(meanshift600all(1,2),'%0.3f'),' pp'),strcat('SNR =',32,num2str(SNR600,'%0.1f'))},'FontSize',12,'units','normalized');
        
        ylabel('600 nm Shift (pp)', 'Fontsize',14);
        xlabel('Cycle number', 'Fontsize',14);
        grid on
        
        if sensitiv ==1
            %RIS : RISXwl    RISPeakwl
            subplot(2,6,10);
            bar(1:cycles,RISX600(:,1))
            hold on
            errorbar(1:cycles,RISX600(:,1), RISX600(:,2),'.','LineWidth',2)
            hold off
            
            ylabel('600 nm Shift (pp/RIU)', 'Fontsize',14);
            xlabel('Cycle number', 'Fontsize',14);
            grid on
            title(strcat('RIS = ',num2str(RISPeak600(1,1),'%0.1f'),'\pm',num2str(RISPeak600(1,2),'%0.1f'),' pp/RIU'),'FontSize',12,'units','normalized');%strcat('RIS = ',num2str(RISPeakwl,'%0.3f'),' nm/RIU'), 'Fontsize',14);
        end
        
        % 650 nm
        subplot(2,6,5);
        bar(1:cycles,shift650mean(:,1))
        hold on
        errorbar(1:cycles,shift650mean(:,1), shift650mean(:,2),'.','LineWidth',2)
        hold off
        
        title({strcat('Shift = ',num2str(meanshift650all(1,1),'%0.3f'),'\pm',...
            num2str(meanshift650all(1,2),'%0.3f'),' pp'),strcat('SNR =',32,num2str(SNR650,'%0.1f'))},'FontSize',12,'units','normalized');
        
        ylabel('650 nm Shift (pp)', 'Fontsize',14);
        xlabel('Cycle number', 'Fontsize',14);
        grid on
        
        if sensitiv ==1
            %RIS : RISXwl    RISPeakwl
            subplot(2,6,11);
            bar(1:cycles,RISX650(:,1))
            hold on
            errorbar(1:cycles,RISX650(:,1), RISX650(:,2),'.','LineWidth',2)
            hold off
            
            ylabel('650 nm Shift (pp/RIU)', 'Fontsize',14);
            xlabel('Cycle number', 'Fontsize',14);
            grid on
            title(strcat('RIS = ',num2str(RISPeak650(1,1),'%0.1f'),'\pm',num2str(RISPeak650(1,2),'%0.1f'),' pp/RIU'),'FontSize',12,'units','normalized');%strcat('RIS = ',num2str(RISPeakwl,'%0.3f'),' nm/RIU'), 'Fontsize',14);
        end
        
        % 700 nm
        subplot(2,6,6);
        bar(1:cycles,shift700mean(:,1))
        hold on
        errorbar(1:cycles,shift700mean(:,1), shift700mean(:,2),'.','LineWidth',2)
        hold off
        title({strcat('Shift = ',num2str(meanshift700all(1,1),'%0.3f'),'\pm',num2str(meanshift700all(1,2),'%0.3f'),' pp'),strcat('SNR =',32,num2str(SNR700,'%0.1f'))},'FontSize',12,'units','normalized');
        
        ylabel('700 nm Shift (pp)', 'Fontsize',14);
        xlabel('Cycle number', 'Fontsize',14);
        grid on
        
        if sensitiv ==1
            %RIS : RISXwl    RISPeakwl
            subplot(2,6,12);
            bar(1:cycles,RISX700(:,1))
            hold on
            errorbar(1:cycles,RISX700(:,1), RISX700(:,2),'.','LineWidth',2)
            hold off
            
            ylabel('700 nm Shift (pp/RIU)', 'Fontsize',14);
            xlabel('Cycle number', 'Fontsize',14);
            grid on
            title(strcat('RIS = ',num2str(RISPeak700(1,1),'%0.1f'),'\pm',num2str(RISPeak700(1,2),'%0.1f'),' pp/RIU'),'FontSize',12,'units','normalized');%strcat('RIS = ',num2str(RISPeakwl,'%0.3f'),' nm/RIU'), 'Fontsize',14);
        end
        
        set(gcf, 'Position', get(0, 'Screensize'));
        
        savefigname=fullfile(path,'results',strcat(basefoldername,'_RAW1ShiftBarError_',smethod,num2str(smooth),'.png'));
        saveas(gcf,savefigname)
        savefigname=fullfile(path,'results',strcat(basefoldername,'_RAW1hiftBarError_',smethod,num2str(smooth),'.svg'));
        print('-f8',savefigname,'-dsvg')
        
        savefigname=fullfile(path,'results',strcat(basefoldername,'_RAW1hiftBarError_',smethod,num2str(smooth),'.fig'));
        saveas(gcf,savefigname)
        
        
        figure(9)
        
        % 750 nm
        subplot(2,6,1);
        bar(1:cycles,shift750mean(:,1))
        hold on
        errorbar(1:cycles,shift750mean(:,1), shift750mean(:,2),'.','LineWidth',2)
        hold off
        title({strcat('Shift = ',num2str(meanshift750all(1,1),'%0.3f'),'\pm',num2str(meanshift750all(1,2),'%0.3f'),' pp'),strcat('SNR =',32,num2str(SNR750,'%0.1f'))},'FontSize',12,'units','normalized');
        
        ylabel('750 nm Shift (pp)', 'Fontsize',14);
        xlabel('Cycle number', 'Fontsize',14);
        grid on
        
        if sensitiv ==1
            %RIS : RISXwl    RISPeakwl
            subplot(2,6,7);
            bar(1:cycles,RISX750(:,1))
            hold on
            errorbar(1:cycles,RISX750(:,1), RISX750(:,2),'.','LineWidth',2)
            hold off
            
            ylabel('750 nm Shift (pp/RIU)', 'Fontsize',14);
            xlabel('Cycle number', 'Fontsize',14);
            grid on
            title({strcat('RIS = ',num2str(RISPeak750(1,1),'%0.1f'),'\pm',num2str(RISPeak750(1,2),'%0.1f'),' pp/RIU'),''},'FontSize',12,'units','normalized');%strcat('RIS = ',num2str(RISPeakwl,'%0.3f'),' nm/RIU'), 'Fontsize',14);
        end
        
        % 800 nm
        subplot(2,6,2);
        bar(1:cycles,shift800mean(:,1))
        hold on
        errorbar(1:cycles,shift800mean(:,1), shift800mean(:,2),'.','LineWidth',2)
        hold off
        title({strcat('Shift = ',num2str(meanshift800all(1,1),'%0.3f'),'\pm',num2str(meanshift800all(1,2),'%0.3f'),' pp'),strcat('SNR =',32,num2str(SNR800,'%0.1f'))},'FontSize',12,'units','normalized');
        
        ylabel('800 nm Shift (pp)', 'Fontsize',14);
        xlabel('Cycle number', 'Fontsize',14);
        grid on
        
        if sensitiv ==1
            %RIS : RISXwl    RISPeakwl
            subplot(2,6,8);
            bar(1:cycles,RISX800(:,1))
            hold on
            errorbar(1:cycles,RISX800(:,1), RISX800(:,2),'.','LineWidth',2)
            hold off
            
            ylabel('800 nm Shift (pp/RIU)', 'Fontsize',14);
            xlabel('Cycle number', 'Fontsize',14);
            grid on
            title({strcat('RIS = ',num2str(RISPeak800(1,1),'%0.1f'),'\pm',num2str(RISPeak800(1,2),'%0.1f'),' pp/RIU'),''},'FontSize',12,'units','normalized');%strcat('RIS = ',num2str(RISPeakwl,'%0.3f'),' nm/RIU'), 'Fontsize',14);
        end
        
        % 850 nm
        subplot(2,6,3);
        bar(1:cycles,shift850mean(:,1))
        hold on
        errorbar(1:cycles,shift850mean(:,1), shift850mean(:,2),'.','LineWidth',2)
        hold off
        title({strcat('Shift = ',num2str(meanshift850all(1,1),'%0.3f'),'\pm',num2str(meanshift850all(1,2),'%0.3f'),' pp'),strcat('SNR =',32,num2str(SNR850,'%0.1f'))},'FontSize',12,'units','normalized');
        
        ylabel('850 nm Shift (pp)', 'Fontsize',14);
        xlabel('Cycle number', 'Fontsize',14);
        grid on
        
        if sensitiv ==1
            %RIS : RISXwl    RISPeakwl
            subplot(2,6,9);
            bar(1:cycles,RISX850(:,1))
            hold on
            errorbar(1:cycles,RISX850(:,1), RISX850(:,2),'.','LineWidth',2)
            hold off
            
            ylabel('850 nm Shift (pp/RIU)', 'Fontsize',14);
            xlabel('Cycle number', 'Fontsize',14);
            grid on
            title({strcat('RIS = ',num2str(RISPeak850(1,1),'%0.1f'),'\pm',num2str(RISPeak850(1,2),'%0.1f'),' pp/RIU'),''},'FontSize',12,'units','normalized');%strcat('RIS = ',num2str(RISPeakwl,'%0.3f'),' nm/RIU'), 'Fontsize',14);
        end
        
        % realDataIntegral
        subplot(2,6,4);
        bar(1:cycles,shiftrealDataIntegralmean(:,1))
        hold on
        errorbar(1:cycles,shiftrealDataIntegralmean(:,1), shiftrealDataIntegralmean(:,2),'.','LineWidth',2)
        hold off
        title({strcat('Shift = ',num2str(meanshiftrealDataIntegralall(1,1),'%1.2e'),'\pm',num2str(meanshiftrealDataIntegralall(1,2),'%1.1e')),strcat('SNR =',32,num2str(SNRrealDataIntegral,'%0.1f'))},'FontSize',12,'units','normalized');
        
        ylabel('realDataIntegral Shift (pp.nm)', 'Fontsize',14);
        xlabel('Cycle number', 'Fontsize',14);
        grid on
        
        if sensitiv ==1
            %RIS : RISXwl    RISPeakwl
            subplot(2,6,10);
            bar(1:cycles,RISXrealDataIntegral(:,1))
            hold on
            errorbar(1:cycles,RISXrealDataIntegral(:,1), RISXrealDataIntegral(:,2),'.','LineWidth',2)
            hold off
            
            ylabel('realDataIntegral Shift (pp.nm.RIU^{-1})', 'Fontsize',14);
            xlabel('Cycle number', 'Fontsize',14);
            grid on
            title({strcat('RIS = ',num2str(RISPeakrealDataIntegral(1,1),'%1.2e'),'\pm',num2str(RISPeakrealDataIntegral(1,2),'%1.0e')),''},'FontSize',12,'units','normalized');%strcat('RIS = ',num2str(RISPeakwl,'%0.3f'),' nm/RIU'), 'Fontsize',14);
        end
        
        % integralLeft
        subplot(2,6,5);
        bar(1:cycles,shiftintegralLeftmean(:,1))
        hold on
        errorbar(1:cycles,shiftintegralLeftmean(:,1), shiftintegralLeftmean(:,2),'.','LineWidth',2)
        hold off
        title({strcat('Shift = ',num2str(meanshiftintegralLeftall(1,1),'%1.2e'),'\pm',num2str(meanshiftintegralLeftall(1,2),'%1.1e')),strcat('SNR =',32,num2str(SNRintegralLeft,'%0.1f'))},'FontSize',12,'units','normalized');
        
        ylabel('integralLeft Shift (pp.nm)', 'Fontsize',14);
        xlabel('Cycle number', 'Fontsize',14);
        grid on
        
        if sensitiv ==1
            %RIS : RISXwl    RISPeakwl
            subplot(2,6,11);
            bar(1:cycles,RISXintegralLeft(:,1))
            hold on
            errorbar(1:cycles,RISXintegralLeft(:,1), RISXintegralLeft(:,2),'.','LineWidth',2)
            hold off
            
            ylabel('integralLeft Shift (pp.nm.RIU^{-1})', 'Fontsize',14);
            xlabel('Cycle number', 'Fontsize',14);
            grid on
            title({strcat('RIS = ',num2str(RISPeakintegralLeft(1,1),'%1.2e'),'\pm',num2str(RISPeakintegralLeft(1,2),'%1.0e')),''},'FontSize',12,'units','normalized');%strcat('RIS = ',num2str(RISPeakwl,'%0.3f'),' nm/RIU'), 'Fontsize',14);
        end
        
        % integralRight
        subplot(2,6,6);
        bar(1:cycles,shiftintegralRightmean(:,1))
        hold on
        errorbar(1:cycles,shiftintegralRightmean(:,1), shiftintegralRightmean(:,2),'.','LineWidth',2)
        hold off
        title({strcat('Shift = ',num2str(meanshiftintegralRightall(1,1),'%1.2e'),'\pm',num2str(meanshiftintegralRightall(1,2),'%1.1e')),strcat('SNR =',32,num2str(SNRintegralRight,'%0.1f'))},'FontSize',12,'units','normalized');
        
        ylabel('integralRight Shift (pp.nm)', 'Fontsize',14);
        xlabel('Cycle number', 'Fontsize',14);
        grid on
        
        if sensitiv ==1
            %RIS : RISXwl    RISPeakwl
            subplot(2,6,12);
            bar(1:cycles,RISXintegralRight(:,1))
            hold on
            errorbar(1:cycles,RISXintegralRight(:,1), RISXintegralRight(:,2),'.','LineWidth',2)
            hold off
            
            ylabel('integralRight Shift (pp.nm.RIU^{-1})', 'Fontsize',14);
            xlabel('Cycle number', 'Fontsize',14);
            grid on
            title({strcat('RIS = ',num2str(RISPeakintegralRight(1,1),'%1.2e'),'\pm',num2str(RISPeakintegralRight(1,2),'%1.0e')),''},'FontSize',12,'units','normalized');%strcat('RIS = ',num2str(RISPeakwl,'%0.3f'),' nm/RIU'), 'Fontsize',14);
        end
        
        set(gcf, 'Position', get(0, 'Screensize'));
        
        
        savefigname=fullfile(path,'results',strcat(basefoldername,'_RAW2ShiftBarError_',smethod,num2str(smooth),'.png'));
        saveas(gcf,savefigname)
        savefigname=fullfile(path,'results',strcat(basefoldername,'_RAW2hiftBarError_',smethod,num2str(smooth),'.svg'));
        print('-f9',savefigname,'-dsvg')
        
        savefigname=fullfile(path,'results',strcat(basefoldername,'_RAW2hiftBarError_',smethod,num2str(smooth),'.fig'));
        saveas(gcf,savefigname)
        
        %% Polinomial fit data bar plots
        if Peak==0
        else
            figure(2)
            %shift
            subplot(2,5,1);
            bar(1:cycles,shiftWLmean(:,1))
            hold on
            errorbar(1:cycles,shiftWLmean(:,1), shiftWLmean(:,2),'.','LineWidth',2)
            hold off
            title({strcat('Shift = ',num2str(meanshiftWLall(1,1),'%0.3f'),'\pm',num2str(meanshiftWLall(1,2),'%0.3f'),' nm'),strcat('SNR =',32,num2str(SNRWL,'%0.1f'))},'FontSize',12,'units','normalized');
            
            ylabel('LSPR band Wavelength Shift (nm)', 'Fontsize',14);
            xlabel('Cycle number', 'Fontsize',14);
            grid on
            
            if sensitiv ==1
                %RIS : RISXwl    RISPeakwl
                subplot(2,5,6);
                bar(1:cycles,RISXwl(:,1))
                hold on
                errorbar(1:cycles,RISXwl(:,1), RISXwl(:,2),'.','LineWidth',2)
                hold off
                
                ylabel('LSPR band Wavelength Shift (nm/RIU)', 'Fontsize',14);
                xlabel('Cycle number', 'Fontsize',14);
                grid on
                title({strcat('RIS = ',num2str(RISPeakwl(1,1),'%0.3f'),'\pm',num2str(RISPeakwl(1,2),'%0.3f'),' nm/RIU ')},'FontSize',12,'units','normalized');%strcat('RIS = ',num2str(RISPeakwl,'%0.3f'),' nm/RIU'), 'Fontsize',14);
            end
            
            
            subplot(2,5,2);
            bar(1:cycles,shiftTmean(:,1))
            hold on
            errorbar(1:cycles,shiftTmean(:,1), shiftTmean(:,2),'.','LineWidth',2)
            hold off
            
            title({strcat('Shift = ',num2str(meanshiftTall(1,1),'%0.3f'),'\pm',num2str(meanshiftTall(1,2),'%0.3f'),' pp'),strcat('SNR =',32,num2str(SNRT,'%0.1f'))},'FontSize',12,'units','normalized');%systemname,'Fontsize',16);
            ylabel('LSPR band Transmittance Shift (pp)', 'Fontsize',14);
            xlabel('Cycle number', 'Fontsize',14);
            grid on
            
            if sensitiv ==1
                %RIS : RISXwl    RISPeakwl
                subplot(2,5,7);
                bar(1:cycles,RISXT(:,1))
                hold on
                errorbar(1:cycles,RISXT(:,1), RISXT(:,2),'.','LineWidth',2)
                hold off
                
                ylabel('LSPR band Transmittance Shift (pp/RIU)', 'Fontsize',14);
                xlabel('Cycle number', 'Fontsize',14);
                grid on
                title({strcat('RIS = ',num2str(RISPeakT(1,1),'%0.3f'),'\pm',num2str(RISPeakT(1,2),'%0.3f'),' pp/RIU '),''},'FontSize',12,'units','normalized');
            end
            
            
            subplot(2,5,4);
            bar(1:cycles,shiftLWHHmean(:,1))
            hold on
            errorbar(1:cycles,shiftLWHHmean(:,1), shiftLWHHmean(:,2),'.','LineWidth',2)
            hold off
            
            title({strcat('Shift = ',num2str(meanshiftLWHHall(1,1),'%0.3f'),'\pm',num2str(meanshiftLWHHall(1,2),'%0.3f'),' nm'),strcat('SNR =',32,num2str(SNRLWHH,'%0.1f'))},'FontSize',12,'units','normalized');
            ylabel('LWHH Shift (nm)', 'Fontsize',14);
            xlabel('Cycle number', 'Fontsize',14);
            grid on
            
            if sensitiv ==1
                %RIS : RISXwl    RISPeakwl
                subplot(2,5,9);
                bar(1:cycles,RISXLWHH(:,1))
                hold on
                errorbar(1:cycles,RISXLWHH(:,1), RISXLWHH(:,2),'.','LineWidth',2)
                hold off
                
                ylabel('LWHH Shift (nm/RIU)', 'Fontsize',14);
                xlabel('Cycle number', 'Fontsize',14);
                grid on
                title({strcat('RIS = ',num2str(RISPeakLWHH(1,1),'%0.3f'),'\pm',num2str(RISPeakLWHH(1,2),'%0.3f'),' nm/RIU '),''},'FontSize',12,'units','normalized');%strcat('RIS = ',num2str(RISPeakwl,'%0.3f'),' nm/RIU'), 'Fontsize',14);
            end
            
            subplot(2,5,3);
            bar(1:cycles,shiftFWHHmean(:,1))
            hold on
            errorbar(1:cycles,shiftFWHHmean(:,1), shiftFWHHmean(:,2),'.','LineWidth',2)
            hold off
            
            title({strcat('Shift = ',num2str(meanshiftFWHHall(1,1),'%0.3f'),'\pm',num2str(meanshiftFWHHall(1,2),'%0.3f'),' nm'),strcat('SNR =',32,num2str(SNRFWHH,'%0.1f'))},'FontSize',12,'units','normalized');
            
            ylabel('FWHH Shift (nm)', 'Fontsize',14);
            xlabel('Cycle number', 'Fontsize',14);
            grid on
            
            if sensitiv ==1
                %RIS : RISXwl    RISPeakwl
                subplot(2,5,8);
                bar(1:cycles,RISXFWHH(:,1))
                hold on
                errorbar(1:cycles,RISXFWHH(:,1), RISXFWHH(:,2),'.','LineWidth',2)
                hold off
                
                ylabel('FWHH Shift (nm/RIU)', 'Fontsize',14);
                xlabel('Cycle number', 'Fontsize',14);
                grid on
                title({strcat('RIS = ',num2str(RISPeakFWHH(1,1),'%0.3f'),'\pm',num2str(RISPeakFWHH(1,2),'%0.3f'),' nm/RIU '),''},'FontSize',12,'units','normalized');%strcat('RIS = ',num2str(RISPeakwl,'%0.3f'),' nm/RIU'), 'Fontsize',14);
            end
            
            subplot(2,5,5);
            bar(1:cycles,shiftRWHHmean(:,1))
            hold on
            errorbar(1:cycles,shiftRWHHmean(:,1), shiftRWHHmean(:,2),'.','LineWidth',2)
            hold off
            
            title({strcat('Shift = ',num2str(meanshiftRWHHall(1,1),'%0.3f'),'\pm',num2str(meanshiftRWHHall(1,2),'%0.3f'),' nm'),strcat('SNR =',32,num2str(SNRRWHH,'%0.1f'))},'FontSize',12,'units','normalized');
            
            ylabel('RWHH Shift (nm)', 'Fontsize',14);
            xlabel('Cycle number', 'Fontsize',14);
            grid on
            
            if sensitiv ==1
                %RIS : RISXwl    RISPeakwl
                subplot(2,5,10);
                bar(1:cycles,RISXRWHH(:,1))
                hold on
                errorbar(1:cycles,RISXRWHH(:,1), RISXRWHH(:,2),'.','LineWidth',2)
                hold off
                
                ylabel('RWHH Shift (nm/RIU)', 'Fontsize',14);
                xlabel('Cycle number', 'Fontsize',14);
                grid on
                title({strcat('RIS = ',num2str(RISPeakRWHH(1,1),'%0.3f'),'\pm',num2str(RISPeakRWHH(1,2),'%0.3f'),' nm/RIU '),''},'FontSize',12,'units','normalized');%strcat('RIS = ',num2str(RISPeakwl,'%0.3f'),' nm/RIU'), 'Fontsize',14);
            end
            
            
            set(gcf, 'Position', get(0, 'Screensize'));
            
            
            
            savefigname=fullfile(path,'results',strcat(basefoldername,'_LSPRShiftBarError_',smethod,num2str(smooth),'.png'));
            saveas(gcf,savefigname)
            savefigname=fullfile(path,'results',strcat(basefoldername,'_LSPRShiftBarError_',smethod,num2str(smooth),'.svg'));
            print('-f2',savefigname,'-dsvg')
            
            savefigname=fullfile(path,'results',strcat(basefoldername,'_LSPRShiftBarError_',smethod,num2str(smooth),'.fig'));
            saveas(gcf,savefigname)
        end
        
        %% Central moments data bar plots
        if CentMom==1 & Peak ==1
            figure(7)
            subplot(2,5,1);
            bar(1:cycles,shiftMzeromean(:,1))
            hold on
            errorbar(1:cycles,shiftMzeromean(:,1), shiftMzeromean(:,2),'.','LineWidth',2)
            hold off
            title({strcat('Shift = ',num2str(meanshiftMzeroall(1,1),'%1.1e'),'\pm',num2str(meanshiftMzeroall(1,2),'%1.0e')),strcat('SNR =',32,num2str(SNRMzero,'%0.1f'))},'FontSize',12,'units','normalized');
            
            ylabel('Mzero Shift ', 'Fontsize',14);
            xlabel('Cycle number', 'Fontsize',14);
            grid on
            
            if sensitiv ==1
                %RIS : RISXwl    RISPeakwl
                subplot(2,5,6);
                bar(1:cycles,RISXMzero(:,1))
                hold on
                errorbar(1:cycles,RISXMzero(:,1), RISXMzero(:,2),'.','LineWidth',2)
                hold off
                
                ylabel('Mzero Shift (RIU^{-1})', 'Fontsize',14);
                xlabel('Cycle number', 'Fontsize',14);
                grid on
                title({strcat('RIS = ',num2str(RISPeakMzero(1,1),'%0.3f'),'\pm',num2str(RISPeakMzero(1,2),'%0.3f'),' RIU^{-1}'),''},'FontSize',12,'units','normalized');%strcat('RIS = ',num2str(RISPeakwl,'%0.3f'),' nm/RIU'), 'Fontsize',14);
            end
            
            subplot(2,5,2);
            bar(1:cycles,shiftM1mean(:,1))
            hold on
            errorbar(1:cycles,shiftM1mean(:,1), shiftM1mean(:,2),'.','LineWidth',2)
            hold off
            title({strcat('Shift = ',num2str(meanshiftM1all(1,1),'%1.1e'),'\pm',num2str(meanshiftM1all(1,2),'%1.0e')),strcat('SNR =',32,num2str(SNRM1,'%0.1f'))},'FontSize',12,'units','normalized');
            
            ylabel('M1 Shift ', 'Fontsize',14);
            xlabel('Cycle number', 'Fontsize',14);
            grid on
            
            if sensitiv ==1
                %RIS : RISXwl    RISPeakwl
                subplot(2,5,7);
                bar(1:cycles,RISXM1(:,1))
                hold on
                errorbar(1:cycles,RISXM1(:,1), RISXM1(:,2),'.','LineWidth',2)
                hold off
                
                ylabel('M1 Shift (RIU^{-1})', 'Fontsize',14);
                xlabel('Cycle number', 'Fontsize',14);
                grid on
                title({strcat('RIS = ',num2str(RISPeakM1(1,1),'%0.3f'),'\pm',num2str(RISPeakM1(1,2),'%0.3f'),' RIU^{-1}'),''},'FontSize',12,'units','normalized');%strcat('RIS = ',num2str(RISPeakwl,'%0.3f'),' nm/RIU'), 'Fontsize',14);
            end
            
            subplot(2,5,3);
            bar(1:cycles,shiftM2mean(:,1))
            hold on
            errorbar(1:cycles,shiftM2mean(:,1), shiftM2mean(:,2),'.','LineWidth',2)
            hold off
            title({strcat('Shift = ',num2str(meanshiftM2all(1,1),'%1.1e'),'\pm',num2str(meanshiftM2all(1,2),'%1.0e')),strcat('SNR =',32,num2str(SNRM2,'%0.1f'))},'FontSize',12,'units','normalized');
            
            ylabel('M2 Shift ', 'Fontsize',14);
            xlabel('Cycle number', 'Fontsize',14);
            grid on
            
            if sensitiv ==1
                %RIS : RISXwl    RISPeakwl
                subplot(2,5,8);
                bar(1:cycles,RISXM2(:,1))
                hold on
                errorbar(1:cycles,RISXM2(:,1), RISXM2(:,2),'.','LineWidth',2)
                hold off
                
                ylabel('M2 Shift (RIU^{-1})', 'Fontsize',14);
                xlabel('Cycle number', 'Fontsize',14);
                grid on
                title({strcat('RIS = ',num2str(RISPeakM2(1,1),'%0.3f'),'\pm',num2str(RISPeakM2(1,2),'%0.3f'),' RIU^{-1}'),''},'FontSize',12,'units','normalized');%strcat('RIS = ',num2str(RISPeakwl,'%0.3f'),' nm/RIU'), 'Fontsize',14);
            end
            
            subplot(2,5,4);
            bar(1:cycles,shiftM3mean(:,1))
            hold on
            errorbar(1:cycles,shiftM3mean(:,1), shiftM3mean(:,2),'.','LineWidth',2)
            hold off
            title({strcat('Shift = ',num2str(meanshiftM3all(1,1),'%1.1e'),'\pm',num2str(meanshiftM3all(1,2),'%1.0e')),strcat('SNR =',32,num2str(SNRM3,'%0.1f'))},'FontSize',12,'units','normalized');
            
            ylabel('M3 Shift ', 'Fontsize',14);
            xlabel('Cycle number', 'Fontsize',14);
            grid on
            
            if sensitiv ==1
                %RIS : RISXwl    RISPeakwl
                subplot(2,5,9);
                bar(1:cycles,RISXM3(:,1))
                hold on
                errorbar(1:cycles,RISXM3(:,1), RISXM3(:,2),'.','LineWidth',2)
                hold off
                
                ylabel('M3 Shift (RIU^{-1})', 'Fontsize',14);
                xlabel('Cycle number', 'Fontsize',14);
                grid on
                title({strcat('RIS = ',num2str(RISPeakM3(1,1),'%0.3f'),'\pm',num2str(RISPeakM3(1,2),'%0.3f'),' RIU^{-1}'),''},'FontSize',12,'units','normalized');%strcat('RIS = ',num2str(RISPeakwl,'%0.3f'),' nm/RIU'), 'Fontsize',14);
            end
            
            subplot(2,5,5);
            bar(1:cycles,shiftM4mean(:,1))
            hold on
            errorbar(1:cycles,shiftM4mean(:,1), shiftM4mean(:,2),'.','LineWidth',2)
            hold off
            title({strcat('Shift = ',num2str(meanshiftM4all(1,1),'%1.1e'),'\pm',num2str(meanshiftM4all(1,2),'%1.0e')),strcat('SNR =',32,num2str(SNRM4,'%0.1f'))},'FontSize',12,'units','normalized');
            
            ylabel('M4 Shift ', 'Fontsize',14);
            xlabel('Cycle number', 'Fontsize',14);
            grid on
            
            if sensitiv ==1
                %RIS : RISXwl    RISPeakwl
                subplot(2,5,10);
                bar(1:cycles,RISXM4(:,1))
                hold on
                errorbar(1:cycles,RISXM4(:,1), RISXM4(:,2),'.','LineWidth',2)
                hold off
                
                ylabel('M4 Shift (RIU^{-1})', 'Fontsize',14);
                xlabel('Cycle number', 'Fontsize',14);
                grid on
                title({strcat('RIS = ',num2str(RISPeakM4(1,1),'%0.3f'),'\pm',...
                    num2str(RISPeakM4(1,2),'%0.3f'),' RIU^{-1}'),''},'FontSize',12,'units','normalized');
            end
            
            set(gcf, 'Position', get(0, 'Screensize'));
            
            
            %%
            savefigname=fullfile(path,'results',strcat(basefoldername,'_CMShiftBarError_',smethod,num2str(smooth),'.png'));
            saveas(gcf,savefigname)
            savefigname=fullfile(path,'results',strcat(basefoldername,'_CMShiftBarError_',smethod,num2str(smooth),'.svg'));
            print('-f7',savefigname,'-dsvg')
            
            savefigname=fullfile(path,'results',strcat(basefoldername,'_CMShiftBarError_',smethod,num2str(smooth),'.fig'));
            saveas(gcf,savefigname)
        end
    end
    %% Make Table with all LSPR Band parameters
    if cycles==0
    else
        OTCTable = table(wavelength,OTC);
        OTCTable.Properties.VariableNames(1) = {'Wavelength'};
        
        %Write Table with OTC
        linefilename = fullfile(path,'results',strcat(basefoldername,'_OTC','.txt'));
        fileID = fopen(linefilename,'at');
        writetable(OTCTable,linefilename,'Delimiter','\t');
        fclose(fileID);
    end
    
    %%
    MAC=1;
    if MAC==1
        %Make Table with all LSPR Band parameters
        if Peak==0
            
        else
            %Write Table with all LSPR Band parameters
            LSPRparTable = table(xx,minwlfitn,minTfitn);
            LSPRparTable.Properties.VariableNames = {'Time','LSPRBandWL','LSPRBandTr'};
            
            linefilename = fullfile(path,'results',strcat(basefoldername,'_LSPRBandandParameters','.txt'));
            fileID = fopen(linefilename,'at');
            writetable(LSPRparTable,linefilename,'Delimiter','\t');
            fclose(fileID);
            
            %Write Table with all LSPR Band parameters
            CMparTable = table(xx,Mzero.',M1.',M2.',M3.',M4.');
            CMparTable.Properties.VariableNames = {'Time','Mzero','M1','M2','M3','M4'};
            
            linefilename = fullfile(path,'results',strcat(basefoldername,'_CMomentsParameters','.txt'));
            fileID = fopen(linefilename,'at');
            writetable(CMparTable,linefilename,'Delimiter','\t');
            fclose(fileID);
            
            %Write Table with all LSPR Band FIT coeficients
            FITCoefTable = table(xx,fliplr(coeffitPeakn));
            FITCoefTable.Properties.VariableNames = {'Time','Coefx0to9'};
            
            linefilename = fullfile(path,'results',strcat(basefoldername,'_FIT_coeficients_','.txt'));
            fileID = fopen(linefilename,'at');
            writetable(FITCoefTable,linefilename,'Delimiter','\t');
            fclose(fileID);
            
            %Write Table with all analysis parameters
            Algpar=[{'valey','savefig','cycleFrom','cycleTo',...
                'sensitiv','RI1','RI2','drift','LSPRBand','range',...
                'minWL','maxWL','scantime','cycles','system0name','lic','path','file'} ; ...
                {valey,savefig,cycleFrom,cycleTo,sensitiv,RI1,RI2,drift,Peak,range,...
                minWL,maxWL,scantime,cycles,system0name,lic,path,file}];
            
            AlgparTable = table(Algpar.');
            
            linefilename = fullfile(path,'results',strcat(basefoldername,'_Algorithm_parameters_','.txt'));
            fileID = fopen(linefilename,'at');
            writetable(AlgparTable,linefilename,'Delimiter','\t');
            fclose(fileID);
            
        end
    else
        
       
        
        %file open and LSPR FIT peak funtion save
        
        coeffitPeakn=coeffitPeakn.';
        
        LSPRfitfilename = fullfile(path,'results',strcat(basefoldername,'_FIT_coeficients_',smethod,num2str(smooth),'.txt'));
        fileID = fopen(LSPRfitfilename,'at');
        fprintf(fileID,'%10E %10E %10E %10E %10E %10E %10E %10E %10E %10E \n', coeffitPeakn);
        fclose(fileID);
        
        
        
        
        LSPRBandCondfilename = fullfile(path,'results',strcat(basefoldername,'_Conditions_',smethod,num2str(smooth),'.txt'));
        fileID = fopen(LSPRBandCondfilename,'at');
        fprintf(fileID,'%14s \t %14s \n','System Name = ',systemname2);
        fprintf(fileID,'%14s \t %5E \n','wlminnm = ',wlminnm);
        fprintf(fileID,'%14s \t %5E \n','wlmaxnm = ',wlmaxnm);
        fclose(fileID);
        
        
        
        
    end
    
    %%
    %     %Save Workspace to file
    %     filename = fullfile(path,'results',strcat(basefoldername,'.mat'));
    %     save(filename)
    %%
    
    CreateStruct.Interpreter = 'tex';
    CreateStruct.WindowStyle = 'modal';
    deltadate=DateLimit-Date;
    if deltadate < 30
        message = strcat('\fontsize{14}',sprintf(strcat('All Operations Completed!\nTrial will end on:\n',DateString,'\n','Elapsed time: ',num2str(toc),'s')));
    else
        message = strcat('\fontsize{14}',sprintf(strcat('All Operations Completed!','\n','Elapsed time: ',num2str(toc),'s')));
    end
    f=msgbox(message,'Success',CreateStruct);
    close all
    winopen(fullfile(path,'results'));
    
else
    %a='Folder does not exist or could not be created'
    msgbox('Folder does not exist or could not be created. Operations Canceled', 'Canceled','error');
end
