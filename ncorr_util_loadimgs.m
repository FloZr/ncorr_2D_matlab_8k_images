function [imgs,outstate] = ncorr_util_loadimgs(lazyparam)
% This function loads image(s) from a GUI and returns the corresponding 
% ncorr_class_img(s). If a single image is selected, it can be named anything. 
% If more than one image is selected, it must have the format of
% 'name_#.ext'. These image(s) have type of either 'file' or 'lazy'.
%
% Inputs -----------------------------------------------------------------%
%   lazyparam - logical; if true, then use lazy loading.
%
% Outputs ----------------------------------------------------------------%
%   imgs - ncorr_class_img; the image(s) selected through the GUI.
%   outstate - integer; returns either out.cancelled, out.failed, or out.success.
    
    % Initialize output
    outstate = out.cancelled;
    imgs = ncorr_class_img.empty;
    
    % Images set through this utility will either be 'lazy' or 'file'
    if (lazyparam)
        type = 'lazy';
    else
        type = 'file';
    end    
    
    % Select images
    [filename,pathname] = uigetfile({'*.jpg;*.jpeg;*.tif;*.tiff;*.png;*.bmp;','Image files (.jpg,.tif,.png,.bmp)'},'Select Image(s) (must be .jpg, .tif, .png, or .bmp)','MultiSelect', 'on');     
    
    % Check to see if operation was cancelled 
    if (isequal(filename,0) || isequal(pathname,0))  
        return;
    end
    
    % See how many images were selected (char is 1 image, while a cell
    % array is multiple)
    if (ischar(filename)) 
        % Only 1 image selected 
        if (ncorr_util_properimgfmt(fullfile(pathname,filename),'Image') == out.success)
            imgs = ncorr_class_img;
            if (strcmp(type,'lazy'))
                imgs.set_img(type,struct('name',filename,'path',pathname)); 
            else
                imgs.set_img(type,struct('img',imread(fullfile(pathname,filename)),'name',filename,'path',pathname)); 
            end          
                        
            outstate = out.success;
        else
            h_error = errordlg('Loading image failed. The data was determined to have incorrect formatting.','Error','modal');
            uiwait(h_error);

            outstate = out.failed;      
        end
    else
        % Put in try block because Ncorr can run out of memory here,
        % especially if the user is not using lazy loading.
        try
            % Multiple images selected
            % Files must have format of "name_#.ext"                   
            % Split at each filename's last underscore. R2026a requires
            % scalar colon operands; strfind can return several positions.
            [~,stem,ext] = fileparts(filename{1});
            separator = find(stem == '_',1,'last');
            name = '';
            if (~isempty(separator))
                name = stem(1:separator-1);
            end

            if (~isempty(name) && ~isempty(ext))
                list_num = zeros(length(filename),1);
                imgs_prelim = ncorr_class_img.empty;
                incorrectformat = false;
                for i = 0:length(filename)-1
                    [~,stem_buffer,ext_buffer] = fileparts(filename{i+1});
                    separator_buffer = find(stem_buffer == '_',1,'last');
                    name_buffer = '';
                    num_buffer = NaN;
                    if (~isempty(separator_buffer))
                        name_buffer = stem_buffer(1:separator_buffer-1);
                        num_buffer = str2double(stem_buffer(separator_buffer+1:end));
                    end
                    % A missing separator/extension cannot match the first
                    % filename's nonempty prefix/extension.
                    if (strcmp(name_buffer,name) && ...
                        strcmp(ext_buffer,ext) && ...
                        ncorr_util_isintbb(num_buffer,0,inf,'Image number') == out.success && ...
                        ncorr_util_properimgfmt(fullfile(pathname,filename{i+1}),'Image') == out.success)   
                        % Store img and number
                        list_num(i+1) = num_buffer;
                        imgs_prelim(i+1) = ncorr_class_img;
                        if (strcmp(type,'lazy'))
                            imgs_prelim(i+1).set_img(type,struct('name',filename{i+1},'path',pathname)); 
                        else
                            imgs_prelim(i+1).set_img(type,struct('img',imread(fullfile(pathname,filename{i+1})),'name',filename{i+1},'path',pathname)); 
                        end
                    else
                        incorrectformat = true;
                        break;
                    end
                end      

                if (~incorrectformat)
                    % Now sort list_number
                    [val_sorted,idx_sorted] = sort(list_num); %#ok<ASGLU>

                    % At this point all images are correct; return the
                    % sorted imgs list
                    imgs = imgs_prelim(idx_sorted);  
                    
                    outstate = out.success;
                else
                    h_error = errordlg('Loading images failed. They were determined to have incorrect formatting. Make sure images are in the form of "name_#.ext" and that data has the correct format.','Error','modal');
                    uiwait(h_error);

                    outstate = out.failed;    
                end
            else
                h_error = errordlg('Loading images failed. Name or extension could not be obtained. Make sure images are in the form of "name_#.ext".','Error','modal');
                uiwait(h_error);

                outstate = out.failed;    
            end
        catch ME
            h_error = errordlg(sprintf('Loading images failed.\n%s: %s',ME.identifier,ME.message),'Error','modal');
            uiwait(h_error);

            outstate = out.failed;    
        end
    end
end
