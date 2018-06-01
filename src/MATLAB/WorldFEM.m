%CLASS_INTERFACE Example MATLAB class wrapper to an underlying C++ class
classdef WorldFEM < handle
    properties (SetAccess = private, Hidden = true)
        objectHandle; % Handle to the underlying C++ class instance
    end
    methods
        %% Constructor - Create a new C++ class instance
        function this = WorldFEM(varargin)
            varargin{3} = int32(varargin{3});
            this.objectHandle = WorldFEM_Interface('new', varargin{:});
        end

        %% Destructor - Destroy the C++ class instance
        function delete(this)
            WorldFEM_Interface('delete', this.objectHandle);
        end

        %% setMeshParameters - set mesh parameters
        function varargout = setMeshParameters(this, varargin)
            [varargout{1:nargout}] = WorldFEM_Interface('setMeshParameters', this.objectHandle, varargin{:});
        end

        %% State - get the current state of the world
        function varargout = state(this, varargin)
            [varargout{1:nargout}] = WorldFEM_Interface('state', this.objectHandle, varargin{:});
        end

        %% setState - set the current state of the world
        function varargout = setState(this, varargin)
            [varargout{1:nargout}] = WorldFEM_Interface('setState', this.objectHandle, varargin{:});
        end
        %% setQ - set the current positional state of the world
        function varargout = setQ(this, varargin)
            [varargout{1:nargout}] = WorldFEM_Interface('setQ', this.objectHandle, varargin{:});
        end
        %% setState - set the current state of the world
        function varargout = setQDot(this, varargin)
            [varargout{1:nargout}] = WorldFEM_Interface('setQDot', this.objectHandle, varargin{:});
        end
        %% Mass - get the world mass matrix
        function varargout = mass(this, varargin)
            [varargout{1:nargout}] = WorldFEM_Interface('M', this.objectHandle, varargin{:});
        end

        %% Stiffness - get the world stiffness matrix
        function varargout = stiffness(this, varargin)
            [varargout{1:nargout}] = WorldFEM_Interface('K', this.objectHandle, varargin{:});
        end

        %% strainEnergy - get the strain energy
         function varargout = strainEnergy(this, varargin)
            [varargout{1:nargout}] = WorldFEM_Interface('strener', this.objectHandle, varargin{:});
         end

          function varargout = strainEnergyPerElement(this, varargin)
            [varargout{1:nargout}] = WorldFEM_Interface('strenertet', this.objectHandle, varargin{:});
          end          
          %% bodyForceEnergy - get the body force energy
         function varargout = bodyForceEnergy(this, varargin)
            [varargout{1:nargout}] = WorldFEM_Interface('bdfener', this.objectHandle, varargin{:});
         end
          %% Force - get the world foce
        function varargout = force(this, varargin)
            [varargout{1:nargout}] = WorldFEM_Interface('f', this.objectHandle, varargin{:});
        end
        %% Internal Force - get the world internal (strain) force
        function varargout = internalForce(this, varargin)
            [varargout{1:nargout}] = WorldFEM_Interface('if', this.objectHandle, varargin{:});
        end
         %% Stress - get the stress on each element
        function varargout = stress(this, varargin)
            [varargout{1:nargout}] = WorldFEM_Interface('stress', this.objectHandle, varargin{:});
        end

        %% convenience functions
        function c = strainEnergyFromQ(this, q)
            setQ(this,q);
            c = strainEnergy(this);
        end

         %% convenience functions
         
        function c = bodyForceEnergyFromQ(this, q)
            setQ(this,q);
            c = bodyForceEnergy(this);
        end
        
        %% convenience functions
        function c = strainEnergyPerElementFromQ(this, q)
            setQ(this,q);
            c = strainEnergyPerElement(this);
        end
    end
end
