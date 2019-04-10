classdef DummyCamera < Camera
    %APOGEEALTA Class to control the apogee alta camera
    % Install 64 Bit driver from the ApogeeUnifiedInstaller.exe. You can
    % get it from the apogee support. Link that worked:
    % http://www.andor.com/ftp/21429/ApogeeUnifiedInstaller.zip
    % First version: Jan Vogelsang, August 2018
    
    properties
        
        cam = [];
        exposure = 1; % default: 1s exposure time
        binit = false;
        bits = 16;
        
        roi_w = [];
        roi_h = [];
        
    end
    
    methods
        
        function this = DummyCamera()
        end
        
        % init
        function ret = init(this)
            
                        
            
            this.width = 800;
            this.height = 600;
            
            fprintf('Dummy camera initialized.\n')
            this.binit = true;
            ret = 0;
            return
        end
        
        function setROICam(this,~,v_px,~,h_px,~,~)
            
            this.roi_w = h_px;
            this.roi_h = v_px;

        end
        
        function getImageAsyncStartCam(this) %#ok
            
        end
        
        function data = getImageAsyncEndCam(this)
            
            pause(this.exposure);
            if ~isempty(this.roi_w) && ~isempty(this.roi_h)
                data = uint16(700*rand(this.roi_h,this.roi_w));
            else
                data = uint16(700*rand(this.height,this.width));
            end
            
        end
        
        % set the exposure time in seconds
        function setExposure(this,expTime) % in s
            
            this.exposure = expTime;
        end
        
        % get the exposure time in seconds
        function expTime = getExposure(this) % in s
            expTime = this.exposure;
        end
        
        % clean up
        function delete(this)
            this.binit = false;
        end
        
    end
    
end
