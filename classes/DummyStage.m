classdef DummyStage < Stage
    %DUMMYSTAGE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        binit = false;
    end
    
    methods
        
        function this = DummyStage()
            this.bPiezo = false;
        end
        
        function ret = init(this)
            
            fprintf('Dummy stage initialized.\n')
            ret = 0;
            this.binit = true;
        end
        
        
        function remoteControl(this,bRemote) % boolean input

        end
            
        function setPosition(this,pos_) % m
           
        end
        
        function pos = getPosition(this)
            pos = 0;
        end
        
        function home(this)
            
        end
        
        function delete(this)

            this.binit = false;
        end
        
    end
    
end

