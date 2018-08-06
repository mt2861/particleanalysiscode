function [ output_args ] = saveMovie( dir_to_save_data, filename_movie, mov )
%SAVEMOVIE Summary of this function goes here
%   Detailed explanation goes here
       writerObj2 = VideoWriter([dir_to_save_data '/' filename_movie '.avi']);
       writerObj2.Quality = 25;
       writerObj2.FrameRate = 5;
       open(writerObj2)
       for t = 1:length(mov)
               if size(mov(t).cdata,1)>0 &&  size(mov(t).cdata,2)>0
                writeVideo(writerObj2,mov(t))
               end
       end
       close(writerObj2)

end

