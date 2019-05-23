% ISMRMRDBuffer   Buffer class for fast ISMRMRD file access. Abstracts ISMRMRD
% read/write operations into blocks for faster file read/write and easy handling.  

% Copyright 2017 Skope Magnetic Resonance Technologies AG
% Authors: Benjamin Dietrich, 2017-11-17

classdef ISMRMRDBuffer < handle
    
    properties (Constant, Access = private)
        maxNumberOfBlocks = 5;
    end
    
    properties (Access = private)
        dset            % ISMRMRD dataset object
        blockSize
        bufferMap
        blockAge = []
        ageInd = 1
        age = 0
        writeBuffer
        writeLevel
    end
    
    methods (Access = public)
        
        function this = ISMRMRDBuffer(mrxFile)
            % Constructor. Input parameter: ISMRMRD dataset.
            
            % assign data set and determine buffer block size
            this.dset = ismrmrd.Dataset(mrxFile);
            
            % determine read/buffer block size
            try
                % read first acquisitions
                aq = this.dset.readAcquisition(1);
                this.blockSize = round(1E8/(numel(aq.data{1})*8)); % magic ;)
            catch
                this.blockSize = 5000;
            end
            
            % create read buffer map
            this.bufferMap = containers.Map('KeyType','int32','ValueType','any');
            
            % initialize write buffer if it does not exist
            this.writeBuffer = ismrmrd.Acquisition(this.blockSize);
            this.writeLevel = 0;
        end
        
        function delete(this)
            % Destructor => close ISMRMRD file automatically
            if ~isempty(this.dset)
                
                % write remaining acquisitions in the write buffer to the file
                if this.writeLevel > 0
                    tempBuffer = ismrmrd.Acquisition(this.writeLevel);
                    tempBuffer.head = this.writeBuffer.head.select(1:this.writeLevel);
                    tempBuffer.traj = this.writeBuffer.traj(1:this.writeLevel);
                    tempBuffer.data = this.writeBuffer.data(1:this.writeLevel);
                    
                    this.dset.appendAcquisition(tempBuffer);
                end
                
                % close the file
                this.dset.close();
            end
        end
        
        function N = GetNumberOfAcquisitions(this)
            % Returns the total number of acquisitions within the dataset.
            N = this.dset.getNumberOfAcquisitions();
        end
        
        function aq = GetAcquisition(this,ind)
            % Returns the requested acquisition with index ind.
            
            % check if requested acquisition is buffered
            blockNr = ceil(ind/this.blockSize);
            blockIndex = mod(ind-1,this.blockSize)+1;
            
            try
                % requested acquisition is in the buffer => retreive
                block = this.bufferMap(blockNr);
                aq = block.select(blockIndex);
            catch
                % requested acquisition is not in the buffer => load missing block
                indBlockStart = (blockNr-1)*this.blockSize + 1;
                indBlockEnd = min(indBlockStart - 1 + this.blockSize, ...
                                  this.dset.getNumberOfAcquisitions());
                block = this.dset.readAcquisition(indBlockStart,indBlockEnd);
                
                % get requested acquisition
                aq = block.select(blockIndex);
                
                % add the new block to the buffer
                this.bufferMap(blockNr) = block;
                
                % update block age tracker
                this.blockAge = [this.blockAge,blockNr];
                
                % remove oldest block in the buffer if to many
                if this.bufferMap.Count > this.maxNumberOfBlocks
                    this.bufferMap.remove(this.blockAge(1));
                    this.blockAge = this.blockAge(2:end);
                end
            end
        end
        
        function AppendAcquisition(this,aq)
            % Append acquisition to dataset. Write to file if buffer is full.
            
%             % check the input data type
%             if ~isa(aq,'ismrmrd.Acquisition')
%                 error('Invalid input type!')
%             end
            
            % add acquisition to the buffer
            this.writeLevel = this.writeLevel + 1;
            this.writeBuffer.head.set(this.writeLevel,aq.head);
            this.writeBuffer.traj(this.writeLevel) = aq.traj;
            this.writeBuffer.data(this.writeLevel) = aq.data;
            
            % check if buffer is full and adjust write level
            if this.writeLevel >= this.blockSize
                % write buffer to file
                this.dset.appendAcquisition(this.writeBuffer);
                
                % reset buffer level
                this.writeLevel = 0;
            end
        end
        
        function header = GetXMLHeader(this)
            % Read xml header.
            header = ismrmrd.xml.deserialize(this.dset.readxml());
        end
        
        function SetXMLHeader(this,header)
            % Write xml header.
            this.dset.writexml(ismrmrd.xml.serialize(header)); 
        end
    end
end

