classdef Camera < handle
    %CAMERA superclass
    %   An instance of this class cannot be created
    
    properties
        background = [];
        last_image = [];
        width = 2;
        height = 2;
        
        roi_v_px = 1;
        roi_h_px = 1;
        roi_v_start = 1;
        roi_h_start = 1;
        subroi_v_px = 1;
        subroi_h_px = 1;
        subroi_v_start = 1;
        subroi_h_start = 1;
        roi_v_binning = 1;
        roi_h_binning = 1;
        
        bSubROI = false;
    end
    
    methods (Abstract)
         ret = init(this); % returns 0 if no error occured
         setExposure(this,exposure_time); % s
         setROICam(this,v_start,v_px,h_start,h_px,bv,bh); % vertical start (unbinned), vertical pixels (binned), horizontal start, horizontal pixels (binned), binning vertical/horizontal
         getImageAsyncStartCam(this);
         data = getImageAsyncEndCam(this);
    end
    
    methods
        
        function ret_image = getImage(this)
            this.getImageAsyncStart();
            ret_image = this.getImageAsyncEnd();
        end
        
        function getImageAsyncStart(this)
            this.getImageAsyncStartCam;
        end
                
        function ret_image = getImageAsyncEnd(this)
                   
%             if this.bSubROI
%                 data = zeros(this.subroi_v_px,this.subroi_h_px,'uint16');
%             else
%                 data = zeros(this.roi_v_px,this.roi_h_px,'uint16');
%             end
            data = this.getImageAsyncEndCam;
            
            if this.bSubROI
                ret_image = this.last_image;
                ret_image( this.subroi_v_start : this.subroi_v_start+this.subroi_v_px-1,...
                    this.subroi_h_start : this.subroi_h_start+this.subroi_h_px-1) = data;
            else
                ret_image = data;
            end
            
            this.last_image = ret_image;
            
        end
        
        function setROI(this,v_start,v_px,h_start,h_px,bv,bh) % vertical start (unbinned), vertical pixels (binned), horizontal start, horizontal pixels (binned), binning vertical/horizontal
            
            % check size
            if (h_start-1) + h_px*bh > this.width || (v_start-1) + v_px*bv > this.height
                error('ROI too large.')
            end
            
            % clear last_image and background if size changed
            if v_px~=this.roi_v_px || h_px~=this.roi_h_px || bv~=this.roi_v_binning || bh~=this.roi_h_binning
                this.last_image = zeros(v_px,h_px,'uint16');
                this.background = zeros(v_px,h_px,'uint16');
            end
            
            % set values in class
            this.roi_v_px = v_px;
            this.roi_h_px = h_px;
            this.roi_v_start = v_start;
            this.roi_h_start = h_start;
            this.roi_v_binning = bv;
            this.roi_h_binning = bh;
            
            this.setROICam(this.roi_v_start, this.roi_v_px,this.roi_h_start,this.roi_h_px,this.roi_v_binning,this.roi_h_binning);
        end
        
        function setSubROI(this,v_start,v_px,h_start,h_px) % all binned units
            
            % check size
            if (h_start-1) + h_px > this.roi_h_px || (v_start-1) + v_px > this.roi_v_px
                error('SubROI too large.')
            end
            
            this.subroi_v_px = v_px;
            this.subroi_h_px = h_px;
            this.subroi_v_start = v_start;
            this.subroi_h_start = h_start;
            this.bSubROI = true;
            this.setROICam(this.roi_v_start + (v_start-1)*this.roi_v_binning,...
                v_px,...
                this.roi_h_start + (h_start-1)*this.roi_h_binning,...
                h_px,...
                this.roi_v_binning,this.roi_h_binning);
        end
        
        function clearSubROI(this)
            this.bSubROI = false;
            this.setROICam(this.roi_v_start, this.roi_v_px,this.roi_h_start,this.roi_h_px,this.roi_v_binning,this.roi_h_binning);
        end
        
    end
    
end

