classdef IDSCamera < Camera
    %IDSCAMERA Class to control IDS cameras
    % https://de.ids-imaging.com/manuals/uEye_SDK/EN/uEye_DotNET_Manual/index.html
    % First version: Jan Vogelsang, February 2019
    
    properties
        
        cam = [];
        exposure = 0.3e-4; % default: 0.3ms exposure time
        binit = false;
        bits = [];
        serial = [];
        img = [];
        fps = 14; %14 Hz default framerate
        live = false;
    end
    
    methods
        
        function this = IDSCamera(serial_)
            this.serial = serial_;
        end
        
        % init
        function ret = init(this)
            
            %   Add NET assembly if it does not exist
            asm = System.AppDomain.CurrentDomain.GetAssemblies;
            if ~any(arrayfun(@(n) strncmpi(char(asm.Get(n-1).FullName), ...
                    'uEyeDotNet', length('uEyeDotNet')), 1:asm.Length))
                NET.addAssembly(...
                    'C:\lab programs\siscan\dll\uEyeDotNet.dll');
            end
            
            %  Create camera object handle
            this.cam = uEye.Camera;
            device_id = 0; % starts from 0
            ret = char(this.cam.Init(device_id));
            if ~strcmpi(ret, 'SUCCESS')
                error('Could not initialize camera');
            end
            [~,caminfo] = this.cam.Information.GetCameraInfo;
            if ~strcmp(char(caminfo.SerialNumber),this.serial)
                error('Used the first device and got serial %s, but we wanted %s!',char(caminfo.SerialNumber),this.serial)
            end
            
            %   Set display mode to bitmap (DiB)
            if ~strcmpi(char(this.cam.Display.Mode.Set(uEye.Defines.DisplayMode.DiB)), ...
                    'SUCCESS')
                error('Could not set display mode');
            end
            %   Set colormode 
            if ~strcmpi(char(this.cam.PixelFormat.Set(uEye.Defines.ColorMode.Mono8)), ...
                    'SUCCESS')
                error('Could not set pixel format.');
            end
%             %   Set trigger mode to software (single image acquisition)
%             if ~strcmpi(char(this.cam.Trigger.Set(uEye.Defines.TriggerMode.Software)), 'SUCCESS')
%                 error('Could not set trigger format to software.');
%             end
            
            % get /set framerate
            setFramerate(this,this.fps);
            
            %   Allocate memory for a single image
            [ErrChk, this.img.ID] = this.cam.Memory.Allocate(true);
            if ~strcmpi(char(ErrChk), 'SUCCESS')
                error('Could not allocate memory');
            end
            
            % set active image id. not needed since allocate gets the "true" parameter.
            this.cam.Memory.SetActive(this.img.ID);
            
            %   Obtain image information
            [ErrChk, this.img.Width, this.img.Height, this.img.Bits, this.img.Pitch] ...
                = this.cam.Memory.Inquire(this.img.ID);
            if ~strcmpi(char(ErrChk), 'SUCCESS')
                error('Could not get image information');
            end

            this.width = this.img.Width;
            this.height = this.img.Height;
            this.bits = this.img.Bits;
            
            fprintf('IDS camera initialized.\n')
            
            this.binit = true;
            ret = 0;
            return
        end
        
        function setROICam(this,v_start,v_px,h_start,h_px,bv,bh) %#ok
            error('Method not implemented yet.')
        end
        
        % get single image
        function getImageAsyncStartCam(this)
                        
            % start recording an image with the uEye
            if ~strcmpi(char(this.cam.Acquisition.Freeze(false)), 'SUCCESS') % true would block
                error('Could not acquire image');
            end
            
            [~,has_started] = this.cam.Acquisition.HasStarted;
            if ~has_started
                error('error starting image acquisition on the uEye camera.')
            end
            
        end
        
        function data = getImageAsyncEndCam(this)
            
            % wait until ueye is finished
            pause_time = 0.05; % ms
            cnt_max = max(2*this.exposure/pause_time,500e-3/pause_time); % wait twice the exposure time, but at least 500ms
            cnt = 0;
            while 1
                [~,ueye_finished] = this.cam.Acquisition.IsFinished;
                if ueye_finished
                    break
                end
                cnt = cnt+1;
                if cnt > cnt_max
                    error('uEye does not finish.')
                end
                pause(pause_time);
            end
            
            %   Extract uEye image
            [ErrChk, tmp] = this.cam.Memory.CopyToArray(this.img.ID,uEye.Defines.ColorMode.Mono8);
            if ~strcmpi(char(ErrChk), 'SUCCESS')
                error('Could not obtain image data');
            end
            
            %   Reshape image
            data = reshape(uint8(tmp), [this.img.Width, this.img.Height, this.img.Bits/this.cam.PixelFormat.GetBitsPerPixel]).';
            
        end
        
        % go into LiveView and display the image
        
        function startLiveView(this)
            this.cam.Acquisition.Capture; 
            
            if ~strcmpi(char(this.cam.Acquisition.Capture), 'CAPTURE_RUNNING')
                error('Could not start the live video');
            else
                this.live = true;
            end 
%            this.img.MemorySequence.InitImageQueue() %initialize image queue
        end
            
        function data = extractImage(this)    
                %   Extract uEye image
            [ErrChk, tmp] = this.cam.Memory.CopyToArray(this.img.ID,uEye.Defines.ColorMode.Mono8);
            if ~strcmpi(char(ErrChk), 'SUCCESS')
                error('Could not obtain image data');
            end
            
            %   Reshape image
            data = reshape(uint8(tmp), [this.img.Width, this.img.Height, this.img.Bits/this.cam.PixelFormat.GetBitsPerPixel]).';
         
        end
        
        
        
        % stop live acquisition
        
        function stopLiveView(this)
            if this.live
%                 this.img.MemorySequence.ExitImageQueue()
                this.cam.Acquisition.Stop
                fprintf('Acquisition stopped.\n')
            else
                error('Camera is not live');
            end
            this.live = false;
        end
        
        
        %uEye.GainHardwareScaled.SetMaster(int s32Value)  
      
        % gain setting for later implementation
        
        % set the exposure time in seconds
        function setExposure(this,expTime) % in s
            
            exposure_time = expTime*1e3; % ms
            
            [~,exposure_range] = this.cam.Timing.Exposure.GetRange;
                        
            if exposure_time >= exposure_range.Maximum
                error('Maximum exposure is %d ms.\n',round(exposure_range.Maximum));
            end
            
            if exposure_time <= exposure_range.Minimum
                error('Minimum exposure is %d ms.\n',exposure_range.Minimum);
            end
            
            this.cam.Timing.Exposure.Set(exposure_time);
            [~,val] = this.cam.Timing.Exposure.Get;
            fprintf('Exposure time of uEye set to %.1fms.\n',val)
            
            this.exposure = expTime;
        end
        
        % get the exposure time in seconds
        function expTime = getExposure(this) % in s
            expTime = this.exposure;
        end
        
        % set the frame rate
        function fps_out = setFramerate(this,fps_in) % in 1/s
            % get /set framerate
            this.cam.Timing.Framerate.Set(fps_in);
            [~,fps_out] = this.cam.Timing.Framerate.Get;
            fprintf('Framerate set to %.1f Hz.\n',fps_out)
            
        end
        
        % clean up
        function delete(this)
            if this.binit
               
                this.cam.Exit
                clear cam
                
                fprintf('Camera closed.\n')
            end
            this.binit = false;
        end
        
    end
    
end
