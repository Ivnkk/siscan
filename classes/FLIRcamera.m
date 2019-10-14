classdef FLIRcamera < Camera
    %FLIRCAMERA Class to control FLIR cameras
    %   Detailed explanation goes here
    
    properties
        cam = [];
        exposure = 0.3e-4; % default: 0.3ms exposure time
        binit = false;
        bits = [];
        serial = [];
        img = [];
        fps = 26; %26 Hz default framerate
        live = false;
        colormode = [];
    end
    
    methods
        function this = FLIRcamera(serial_) % constructor method for FLIR camera class
            
            this.serial = serial_; %Update serial number property of instatiated ob
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

