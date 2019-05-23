%% clear and add libs to path
clear all; clc; close all;

%% Add mex functions for ECC simulation
addpath(fullfile('..','Matlab'))


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                          PATH SETTINGS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inputFolderName  = '\\\\NAS-skope\\SoftwareDevelopment\\ReconCustomerData\\2019-04-12_ECC_Validation\\phantom\\siemensData';
outputFolderName = fullfile(pwd,'ECC_Validation');

if ~isfolder(outputFolderName)
   mkdir(outputFolderName) 
end

scanWithGradients    = 'meas_MID155_epi_withGradients';
scanWithoutGradients = 'meas_MID156_epi_withoutGradients';
dspFile              = 'meas_MID155_epi_withGradients';


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                          CONVERSION TO MRD
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Style sheet for conversion
stylesheetFullFilePath      = fullfile('..','..','parameter_maps','StyleSheet_EPI_3D.xsl');

%% Data with gradients
siemensFileNameWithGradient = fullfile(inputFolderName,[scanWithGradients '.dat']);
mrdFileNameWithGradient     = fullfile(outputFolderName,[scanWithGradients '.h5']);

% Run system command
if system(['"../../build/Release/siemens_to_ismrmrd.exe"' ' -f ' siemensFileNameWithGradient ' -o ' mrdFileNameWithGradient ' --user-stylesheet ' stylesheetFullFilePath ])
   error('Siemens to ISMRMRD converter failed.'); 
end

%% Data without gradients
siemensFileNameWithoutGradient              = fullfile(inputFolderName,[scanWithoutGradients '.dat']);
mrdFileNameWithoutGradient                  = fullfile(outputFolderName,[scanWithoutGradients '.h5']);

% Run system command
if system(['"../../build/Release/siemens_to_ismrmrd.exe"' ' -f ' siemensFileNameWithoutGradient ' -o ' mrdFileNameWithoutGradient ' --user-stylesheet ' stylesheetFullFilePath ])
   error('Siemens to ISMRMRD converter failed.'); 
end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                          MEASURED ECC PHASE          
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Open both datasets
dsetGradOn  = ISMRMRDBuffer(mrdFileNameWithGradient);
dsetGradOff = ISMRMRDBuffer(mrdFileNameWithoutGradient);

header = dsetGradOff.GetXMLHeader();

% Allocate memory
aqInGradOn  = dsetGradOn.GetAcquisition(1);
ECCCompMeasured = zeros(size(aqInGradOn.data{1},1),dsetGradOn.GetNumberOfAcquisitions);

% Loop through file
for indAqIn = 1:dsetGradOn.GetNumberOfAcquisitions
    
    % Print progress
    if mod(indAqIn,1000)==0
       disp(indAqIn) 
    end
    
    % Get acquisition
    aqInGradOn  = dsetGradOn.GetAcquisition(indAqIn);
    aqInGradOff = dsetGradOff.GetAcquisition(indAqIn);
        
    % Get data phase difference btw GradOn and GradOff
    dataGradOn  = aqInGradOn.data{1}(:,1);
    dataGradOff = aqInGradOff.data{1}(:,1);
    
    % Phase difference
    phase = angle(dataGradOn.*conj(dataGradOff));
    
    ECCCompMeasured(:,indAqIn) = phase;
     
end

% Calls destructor of ISMRMRDBuffer class and closes file
clear dsetGradOff dsetGradOn


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                          ECC COEFFICIENTS
%                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=0:2
    for para = 1:numel(header.userParameters.userParameterDouble)
        if strcmpi(header.userParameters.userParameterDouble(para).name,['sB0CompensationX_Amp_' num2str(i)])
            CoeffsX.Amp(i+1) = header.userParameters.userParameterDouble(para).value;
        end
        if strcmpi(header.userParameters.userParameterDouble(para).name,['sB0CompensationX_Tau_' num2str(i)])
            CoeffsX.Tau(i+1) = header.userParameters.userParameterDouble(para).value;
        end
    end
    for para = 1:numel(header.userParameters.userParameterDouble)
        if strcmpi(header.userParameters.userParameterDouble(para).name,['sB0CompensationY_Amp_' num2str(i)])
            CoeffsY.Amp(i+1) = header.userParameters.userParameterDouble(para).value;
        end
        if strcmpi(header.userParameters.userParameterDouble(para).name,['sB0CompensationY_Tau_' num2str(i)])
            CoeffsY.Tau(i+1) = header.userParameters.userParameterDouble(para).value;
        end
    end
     for para = 1:numel(header.userParameters.userParameterDouble)
        if strcmpi(header.userParameters.userParameterDouble(para).name,['sB0CompensationZ_Amp_' num2str(i)])
            CoeffsZ.Amp(i+1) = header.userParameters.userParameterDouble(para).value;
        end
        if strcmpi(header.userParameters.userParameterDouble(para).name,['sB0CompensationZ_Tau_' num2str(i)])
            CoeffsZ.Tau(i+1) = header.userParameters.userParameterDouble(para).value;
        end
     end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                         SIMULATED ECC PHASE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
interpolateToRX = true;
[RXSamples, RXTime, TXTime, TrigTime, ECCCompSimulated] = SeqSim_mex(fullfile(inputFolderName,[dspFile,'.xml']),'eddyPhase', interpolateToRX, CoeffsX,CoeffsY,CoeffsZ);      
ECCCompSimulated = reshape(ECCCompSimulated,RXSamples(1),[]);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                        COMPARE BOTH PHASES          
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Undo drift of measured data
ECCCompMeasuredWoDrift = ECCCompMeasured(:) + 0.054*double(RXTime);

plot(unwrap(ECCCompMeasuredWoDrift),'LineWidth',2)
hold on
plot(ECCCompSimulated(:)+2,'LineWidth',2)

xlabel('ADC samples')
ylabel('ECC phase [rad]')
title('ECC simulation vs measurement')
legend({'Measured','Simulated'})



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                      ECC IN CHECK CONVERTER    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Data with gradients
% Run system command
mrdFileNameWithGradientECCUndone     = fullfile(outputFolderName,[scanWithGradients '2.h5']);
if system(['"../../build/Release/siemens_to_ismrmrd.exe"' ' -f ' siemensFileNameWithGradient ' -o ' mrdFileNameWithGradientECCUndone ' --user-stylesheet ' stylesheetFullFilePath ' -d ' fullfile(inputFolderName,[dspFile '.xml'])])
   error('Siemens to ISMRMRD converter failed.'); 
end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                          MEASURED ECC PHASE          
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Open both datasets
dsetGradOn  = ISMRMRDBuffer(mrdFileNameWithGradientECCUndone);
dsetGradOff = ISMRMRDBuffer(mrdFileNameWithoutGradient);

header = dsetGradOff.GetXMLHeader();

% Allocate memory
aqInGradOn  = dsetGradOn.GetAcquisition(1);
ECCCompMeasured = zeros(size(aqInGradOn.data{1},1),dsetGradOn.GetNumberOfAcquisitions);

% Loop through file
for indAqIn = 1:dsetGradOn.GetNumberOfAcquisitions
    
    % Print progress
    if mod(indAqIn,1000)==0
       disp(indAqIn) 
    end
    
    % Get acquisition
    aqInGradOn  = dsetGradOn.GetAcquisition(indAqIn);
    aqInGradOff = dsetGradOff.GetAcquisition(indAqIn);
        
    % Get data phase difference btw GradOn and GradOff
    dataGradOn  = aqInGradOn.data{1}(:,1);
    dataGradOff = aqInGradOff.data{1}(:,1);
    
    % Phase difference
    phase = angle(dataGradOn.*conj(dataGradOff));
    
    ECCCompMeasured(:,indAqIn) = phase;
     
end

% Calls destructor of ISMRMRDBuffer class and closes file
clear dsetGradOff dsetGradOn

%% Plot difference to check that ECC has been undon
% The linear drift is also removed
plot(ECCCompMeasured(:)+ 0.054*double(RXTime))







