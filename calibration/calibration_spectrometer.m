function calibration_spectrometer(screenNumber)

Nsteps = 10;

% setup
% 4 fields for specifing equipment,nominal sync frequency (1=adaptive),
% integration time (ms), measurements to average, metric units(1/0), enter
setup = 'S,,,,,100,3,1\n';


% measurement
% response code 5 = spectral data
measurement = 'M5';

try
    close all;
    clear all;
    sca
%% start screen
    fprintf('Press a key to start the psychtoolbox screen')
    KbStrokeWait;
    PsychImaging('PrepareConfiguration');
    PsychDefaultSetup(0)
    PsychImaging('AddTask', 'General', 'UseFastOffscreenWindows')
    
    
    white = WhiteIndex(screenNumber);
    black = BlackIndex(screenNumber);
    
    
    if screenNumber ==0
        [windowNr,rect]=PsychImaging('OpenWindow',screenNumber,black+0.5*(white-black),[100 100 400 400]);
    else
        [windowNr,rect]=PsychImaging('OpenWindow',screenNumber,black+0.5*(white-black));
    end
    
%% start spectrometer
    fprintf('Press a key to start the spectrometer')
    KbStrokeWait;
    spectrometer = serial('COM1','BaudRate',9600);
    spectrometer.Timeout = 30;
    
    fopen(spectrometer);
    fprintf(spectrometer,'%s',setup);
    
%% measurements 
    %first a measurement of the gray screen to see that all settings work
    fprintf('Press a key to start the test measurement')
    KbStrokeWait;
    Screen('Flip',windowNr)
    fprintf(spectrometer,'%s',measurement);
    mtest = fscanf(spectrometer);
    
    % now do the whole measurement (Nsteps per gun + grayscale to check
    % additivity)
    fprintf('As all configurations seem to work, your next keypress will start the full measurement.\n You may want to leave, as this takes time')
    KbStrokeWait;
catch e
    sca;
    rethrow(e)
end