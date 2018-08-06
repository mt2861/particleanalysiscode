classdef Imseq
    %IMSEQ Summary of this class goes here
    %   Detailed explanation goes here    
    properties
        directory
        dimensions
        no_frames
        zero_index
        extension
        pad_len
        step
        flip_axis
    end
    
    methods
        function obj = Imseq(directory, step, flip_axis)
            obj.step = step;
            obj.flip_axis = flip_axis
            obj.directory = directory;
            frames = dir(directory);
            % Find the maximum frame
            max_frame = 0;
            for i = 1:length(frames)
                frame = strsplit(frames(i).name, '.');
                if length(frame{1}) < 1
                    continue
                end
                [frame_no, status] = str2num(frame{1});
                % Is frame{1} is a number?
                if status == 0
                    continue
                end
                % If we find a zero frame, then our imseq is zero indexed
                if frame_no == 0
                    obj.zero_index = 1;
                end
                if frame_no > max_frame
                    max_frame = frame_no;
                    obj.extension = frame{2};
                    obj.pad_len = length(frame{1});
                end
            end
            if obj.zero_index == 1
                obj.no_frames = floor((max_frame + 1)/step);
            else
                obj.zero_index = 0;
                obj.no_frames = floor(max_frame/step);
            end
            obj.dimensions = size(obj.get(1));
        end
        function image = get(obj, index)
            image_filename = strcat('%0',num2str(obj.pad_len, '%d'),'i','.',obj.extension);
            if obj.zero_index
                index = index*obj.step-1;
            else
                index = index*obj.step;
            end
            image_path = fullfile(obj.directory, num2str(index, image_filename));
            if obj.extension == 'tif'
                tiff = Tiff(image_path);
                image = read(tiff);
            else
                image = imread(image_path);
            end
            if obj.flip_axis
                image = flip(image, 1);
            end
        end
    end
    
end

