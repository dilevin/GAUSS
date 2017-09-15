%CLASS_INTERFACE Example MATLAB class wrapper to an underlying C++ class
classdef WorldFEM < handle
    properties (SetAccess = private, Hidden = true)
        objectHandle; % Handle to the underlying C++ class instance
    end
    methods
        %% Constructor - Create a new C++ class instance 
        function this = WorldFEM(varargin)
            this.objectHandle = WorldFEM_Interface('new', varargin{:});
        end
        
        %% Destructor - Destroy the C++ class instance
        function delete(this)
            WorldFEM_Interface('delete', this.objectHandle);
        end

        %% Train - an example class method call
        function varargout = state(this, varargin)
            [varargout{1:nargout}] = WorldFEM_Interface('state', this.objectHandle, varargin{:});
        end

        %% Test - another example class method call
        function varargout = mass(this, varargin)
            [varargout{1:nargout}] = WorldFEM_Interface('M', this.objectHandle, varargin{:});
        end
        
        function varargout = stiffness(this, varargin)
            [varargout{1:nargout}] = WorldFEM_Interface('K', this.objectHandle, varargin{:});
        end
        
        function varargout = force(this, varargin)
            [varargout{1:nargout}] = WorldFEM_Interface('f', this.objectHandle, varargin{:});
        end
    end
end
