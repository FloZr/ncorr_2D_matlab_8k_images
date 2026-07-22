function [seedinfo,threaddiagram,outstate] = ncorr_gui_setseeds(reference,current,roi,num_region,spacing,radius,cutoff_diffnorm,cutoff_iteration,total_threads,enabled_stepanalysis,subsettrunc,num_img,total_imgs,pos_parent)
% This is a GUI for setting the rgdic seeds. The number of seeds placed in
% each region depends on the number of threads. 
%
% Inputs -----------------------------------------------------------------%
%   reference - ncorr_class_img; used for displaying the background
%   image and calculations
%   current - ncorr_class_img(s); used for displaying the background image
%   and calculations
%   roi - ncorr_class_roi; roi corresponding to the reference image.
%   num_region - integer; number corresponding to the region; 
%   spacing - integer; spacing between subsets
%   radius - integer; radius of subsets
%   cutoff_diffnorm - double; cutoff of norm of the difference vector
%   cutoff_iteration - integer; cutoff for the number of iterations
%   total_threads - integer; total number of threads
%   enabled_stepanalysis - logical; if true, then process as many images as
%   possible. If false, process all images.
%   subsettrunc - logical; if true, then subset truncation is enabled
%   num_img - integer; current number of reference image
%   total_imgs - integer; total number of images being processed
%   pos_parent - integer array; this is the position of the parent figure
%   which determines where to position this figure
%
% Outputs ----------------------------------------------------------------%
%   seedinfo - struct; contains
%   struct('paramvector',{},'num_region',{},'num_thread',{},'computepoints'
%   ,{}) 
%   threaddiagram - integer array; array the same size as the reduced
%   ROI. The elements indicate which points will be computed by which
%   threads
%   outstate - integer; returns either out.cancelled, out.failed, or
%   out.success.
%
% Note that if step analysis is enabled, it is gauranteed that at least one 
% image will be succesfully seeded if outstate returns true. If step 
% analysis is disabled, it is gauranteed that all images will be seeded or 
% if outstate returns true.

    % Data ---------------------------------------------------------------%     
    % Initialize Outputs 
    outstate = out.cancelled;
    threaddiagram = [];
    seedinfo = struct('paramvector',{},'num_region',{},'num_thread',{},'computepoints',{}); % paramvector = [x y u v du/dx du/dy dv/dx dv/dy corrcoef]
    % Get GUI handles    
    handles_gui = init_gui();
    % Run c-tor
    feval(ncorr_util_wrapcallbacktrycatch(@constructor,handles_gui.figure));
    
    % Callbacks and functions --------------------------------------------%
    function constructor()                        
        % Initialize buffers
        % Initialize reduced ref
        ref_reduced = reference.reduce(spacing);
        % Initialize roi_reduced
        roi_reduced = roi.reduce(spacing);
        % Get region mask
        regionmask = roi_reduced.get_regionmask(num_region);   
        % Initialize threaddiagram_prelim and threaddiagram preview - these
        % are modified in-place
        threaddiagram_prelim = -ones(size(roi_reduced.mask));
        preview_threaddiagram = zeros(size(roi_reduced.mask));
        
        % Subset Preview Array
        preview_subset = zeros(2*radius+1);
        
        % Store data   
        setappdata(handles_gui.figure,'ref_reduced',ref_reduced);  
        setappdata(handles_gui.figure,'roi_reduced',roi_reduced);     
        setappdata(handles_gui.figure,'handle_preview',[]);    
        setappdata(handles_gui.figure,'handle_points',impoint.empty);          
        setappdata(handles_gui.figure,'regionmask',regionmask);        
        setappdata(handles_gui.figure,'threaddiagram_prelim',threaddiagram_prelim);       
        setappdata(handles_gui.figure,'preview_threaddiagram',preview_threaddiagram);   
        setappdata(handles_gui.figure,'preview_subset',preview_subset);
        setappdata(handles_gui.figure,'handle_subset',[]);
        setappdata(handles_gui.figure,'pos_preview',[]); % Store the x, y of the last active point
        setappdata(handles_gui.figure,'active_seed',[]); % One-based index of the seed selected for coordinate editing
        setappdata(handles_gui.figure,'placing_with_mouse',false);
        
        % Update thread diagram and preview - must do this after setting the data
        get_threaddiagram([]);        
        
        % Update
        update_axes('set');
        update_sidemenu();
        
        % Set visible        
        set(handles_gui.figure,'Visible','on');         
    end    

    function callback_button_setseeds(hObject,eventdata) %#ok<INUSD>
        freeze_menu();
        
        % Get data
        roi_reduced = getappdata(handles_gui.figure,'roi_reduced');
        handle_points = getappdata(handles_gui.figure,'handle_points');
        num_existing = length(handle_points);
        
        % Set seeds ------------------------------------------------------%        
        % Continue after any seeds which were entered manually.
        if (isempty(handle_points))
            handle_points = impoint.empty;
        end
        for i = num_existing:total_threads-1
            % Place seed and make sure it's within region
            seedfinished = false;                
            while (~seedfinished)   
                handle_point_buffer = impoint(handles_gui.axes_setseeds);  
                
                % impoint will be empty if esc is pressed
                if (~isempty(handle_point_buffer))   
                    % Get seed position - convert to zero based indexing
                    % and round
                    pos_seed = round(getPosition(handle_point_buffer))-1; 
                    
                    % Make sure point is within the correct region
                    if (roi_reduced.get_num_region(pos_seed(1),pos_seed(2),false(1,length(roi_reduced.region))) == num_region)       
                        % Store and format impoint
                        handle_points(i+1) = handle_point_buffer;
                        setColor(handle_points(i+1),'g');
                        set(handle_points(i+1),'UIContextMenu','');                    

                        % Set constraint function
                        api_point = iptgetapi(handle_points(i+1));
                        constrainfcn = ncorr_util_formregionconstraint(roi_reduced.region(num_region+1));
                        api_point.setPositionConstraintFcn(constrainfcn);    

                        % Assign Position callback
                        addNewPositionCallback(handle_points(i+1),ncorr_util_wrapcallbacktrycatch(@(pos)callback_impoint(pos,i+1),handles_gui.figure));
                        
                        % Update threaddiagram and threaddiagram preview           
                        get_threaddiagram(handle_points);

                        % Set data
                        setappdata(handles_gui.figure,'handle_points',handle_points);
                        setappdata(handles_gui.figure,'pos_preview', pos_seed);
                        setappdata(handles_gui.figure,'active_seed',i+1);

                        % Update here to ensure preview is updated as the seeds
                        % are being placed.
                        update_sidemenu();
                        update_axes('update');                        
                        
                        % Exit
                        seedfinished = true;  
                    else
                        delete(handle_point_buffer);
                        
                        h_error = errordlg('Seed point not within region, try again.','Error','modal');     
                        uiwait(h_error);
                    end 
                else
                    % Escape was pressed - remove only points added during
                    % this mouse-placement pass. Preserve manually entered
                    % seeds which existed before the pass began.
                    if (length(handle_points) > num_existing)
                        delete(handle_points(num_existing+1:end));
                    end
                    if (num_existing == 0)
                        handle_points = impoint.empty;
                        active_seed = [];
                        pos_preview = [];
                    else
                        handle_points = handle_points(1:num_existing);
                        active_seed = num_existing;
                        pos_preview = round(getPosition(handle_points(active_seed)))-1;
                    end
                    
                    % Update threaddiagram and threaddiagram preview           
                    get_threaddiagram(handle_points);

                    % Set data
                    setappdata(handles_gui.figure,'handle_points',handle_points);
                    setappdata(handles_gui.figure,'pos_preview',pos_preview);
                    setappdata(handles_gui.figure,'active_seed',active_seed);
                    
                    % Update
                    update_sidemenu();
                    update_axes('update');
                    
                    % Return since analysis was cancelled
                    unfreeze_menu();
                    return;
                end
            end
        end
                                
        % Update
        update_sidemenu();
        
        unfreeze_menu();
    end    

    function callback_button_addseed(hObject,eventdata) %#ok<INUSD>
        handle_points = getappdata(handles_gui.figure,'handle_points');
        if (length(handle_points) >= total_threads)
            show_coordinate_error('All required seeds are already set. Select a seed and use "Apply to Active" to move it.');
            return;
        end

        [pos_reduced,isvalid] = read_coordinate_fields();
        if (~isvalid)
            return;
        end
        if (position_is_duplicate(pos_reduced,handle_points,[]))
            show_coordinate_error('Another seed already occupies these coordinates.');
            return;
        end

        seed_number = length(handle_points)+1;
        handle_point_buffer = impoint(handles_gui.axes_setseeds,pos_reduced(1)+1,pos_reduced(2)+1);
        handle_points(seed_number) = handle_point_buffer;
        configure_seed_point(handle_points(seed_number),seed_number);

        setappdata(handles_gui.figure,'handle_points',handle_points);
        setappdata(handles_gui.figure,'active_seed',seed_number);
        setappdata(handles_gui.figure,'pos_preview',pos_reduced);
        get_threaddiagram(handle_points);
        update_sidemenu();
        update_axes('update');
    end

    function callback_button_applycoords(hObject,eventdata) %#ok<INUSD>
        handle_points = getappdata(handles_gui.figure,'handle_points');
        active_seed = getappdata(handles_gui.figure,'active_seed');
        if (isempty(active_seed) || active_seed < 1 || active_seed > length(handle_points))
            show_coordinate_error('Select an active seed before applying coordinates.');
            return;
        end

        [pos_reduced,isvalid] = read_coordinate_fields();
        if (~isvalid)
            update_sidemenu();
            return;
        end
        if (position_is_duplicate(pos_reduced,handle_points,active_seed))
            show_coordinate_error('Another seed already occupies these coordinates.');
            update_sidemenu();
            return;
        end

        setPosition(handle_points(active_seed),pos_reduced+1);
        setappdata(handles_gui.figure,'pos_preview',pos_reduced);
        get_threaddiagram(handle_points);
        update_sidemenu();
        update_axes('update');
    end

    function callback_popup_activeseed(hObject,eventdata) %#ok<INUSD>
        handle_points = getappdata(handles_gui.figure,'handle_points');
        if (isempty(handle_points))
            return;
        end

        active_seed = get(handles_gui.popup_activeseed,'Value');
        setappdata(handles_gui.figure,'active_seed',active_seed);
        setappdata(handles_gui.figure,'pos_preview',round(getPosition(handle_points(active_seed)))-1);
        update_sidemenu();
        update_axes('update');
    end

    function callback_button_finish(hObject,eventdata) %#ok<INUSD>
        % Get data
        handle_points = getappdata(handles_gui.figure,'handle_points'); 
        threaddiagram_prelim = getappdata(handles_gui.figure,'threaddiagram_prelim');
        
        % Get position - and make sure no two seeds overlap
        pos_seed = [];
        for i = 0:length(handle_points)-1
            % Convert to zero based indexing and round
            pos_buffer = round(getPosition(handle_points(i+1)))-1; 
            pos_seed = vertcat(pos_seed,pos_buffer); %#ok<AGROW>
        end
        if (size(unique(pos_seed,'rows'),1) == size(pos_seed,1))
            % Convert pos_seed to regular coordinates since they are WRT to
            % the reduced ROI at this point
            pos_seed = pos_seed*(spacing+1);
            
            % Get seedinfo and preview of seeds. Note that if step analysis
            % is enabled, then ncorr_gui_seedpreview will return at least one
            % seed if it returns success. If step analysis is disabled, 
            % then ncorr_gui_seedpreview will provide seeds for 
            % all images if it returns true. 
            [seedinfo_prelim,outstate_preview] = ncorr_gui_seedpreview(reference, ...
                                                                       current, ...
                                                                       roi, ...
                                                                       num_region, ...
                                                                       pos_seed, ...
                                                                       radius, ...
                                                                       cutoff_diffnorm, ...
                                                                       cutoff_iteration, ...
                                                                       enabled_stepanalysis, ...
                                                                       subsettrunc, ...
                                                                       num_img, ...
                                                                       total_imgs, ...
                                                                       get(handles_gui.figure,'OuterPosition'));
            if (outstate_preview == out.success)
                % Set outputs
                seedinfo = seedinfo_prelim;
                % Determine number of compute points and append this to
                % seedinfo_buffer
                for i = 0:size(seedinfo,3)-1
                    for j = 0:size(seedinfo,1)-1
                        seedinfo(j+1,1,i+1).computepoints = length(find(threaddiagram_prelim == j));
                    end
                end
                % Assign thread diagram
                threaddiagram = threaddiagram_prelim;
                outstate = out.success;
                
                % Exit                
                close(handles_gui.figure);
            end
        else
            % Two seeds occupy same location, tell user to move
            % generators until they occupy unique locations
            h_error = errordlg('Two or more seeds occupy the same location, please move all seeds to unique locations.','Error','modal');
            uiwait(h_error);
        end
    end

    function callback_button_cancel(hObject,eventdata) %#ok<INUSD>
        close(handles_gui.figure);
    end

    function callback_impoint(pos,seed_number) %#ok<INUSD>
        % Get data
        handle_points = getappdata(handles_gui.figure,'handle_points');
        
        % Update threaddiagram and threaddiagram preview           
        get_threaddiagram(handle_points);
        
        % Update pos preview tracking variable
        setappdata(handles_gui.figure,'pos_preview', round(pos)-1);
        setappdata(handles_gui.figure,'active_seed',seed_number);
        
        % Update
        update_sidemenu();
        update_axes('update');
    end

    function configure_seed_point(handle_point,seed_number)
        setColor(handle_point,'g');
        set(handle_point,'UIContextMenu','');

        roi_reduced = getappdata(handles_gui.figure,'roi_reduced');
        api_point = iptgetapi(handle_point);
        constrainfcn = ncorr_util_formregionconstraint(roi_reduced.region(num_region+1));
        api_point.setPositionConstraintFcn(constrainfcn);
        addNewPositionCallback(handle_point,ncorr_util_wrapcallbacktrycatch(@(pos)callback_impoint(pos,seed_number),handles_gui.figure));
    end

    function [pos_reduced,isvalid] = read_coordinate_fields()
        pos_reduced = [];
        isvalid = false;
        x = str2double(get(handles_gui.edit_seedx,'String'));
        y = str2double(get(handles_gui.edit_seedy,'String'));
        grid_step = spacing+1;

        if (~isfinite(x) || ~isfinite(y))
            show_coordinate_error('X and Y must be finite numeric pixel coordinates.');
            return;
        end
        if (x < 0 || x >= reference.width || y < 0 || y >= reference.height)
            show_coordinate_error(['Coordinates must satisfy 0 <= X <= ' num2str(reference.width-1) ...
                                   ' and 0 <= Y <= ' num2str(reference.height-1) '.']);
            return;
        end

        % Snap typed coordinates to the nearest in-bounds DIC grid point.
        % This avoids requiring the user to calculate multiples of the grid
        % step manually. Write the snapped full-resolution coordinates back
        % to the controls so the final seed location is always explicit.
        max_grid_x = floor((reference.width-1)/grid_step)*grid_step;
        max_grid_y = floor((reference.height-1)/grid_step)*grid_step;
        x = min(max(round(x/grid_step)*grid_step,0),max_grid_x);
        y = min(max(round(y/grid_step)*grid_step,0),max_grid_y);
        set(handles_gui.edit_seedx,'String',num2str(x));
        set(handles_gui.edit_seedy,'String',num2str(y));

        pos_reduced = [x y]/grid_step;
        roi_reduced = getappdata(handles_gui.figure,'roi_reduced');
        if (roi_reduced.get_num_region(pos_reduced(1),pos_reduced(2),false(1,length(roi_reduced.region))) ~= num_region)
            show_coordinate_error('The specified seed coordinate is not within the selected ROI region.');
            pos_reduced = [];
            return;
        end
        isvalid = true;
    end

    function duplicate = position_is_duplicate(pos_reduced,handle_points,ignore_seed)
        duplicate = false;
        for seed_number = 1:length(handle_points)
            if (~isempty(ignore_seed) && seed_number == ignore_seed)
                continue;
            end
            if (isequal(round(getPosition(handle_points(seed_number)))-1,pos_reduced))
                duplicate = true;
                return;
            end
        end
    end

    function show_coordinate_error(message)
        h_error = errordlg(message,'Invalid Seed Coordinates','modal');
        uiwait(h_error);
    end

    function freeze_menu(hObject,eventdata) %#ok<INUSD>
        setappdata(handles_gui.figure,'placing_with_mouse',true);
        update_sidemenu();
    end

    function unfreeze_menu(hObject,eventdata) %#ok<INUSD>
        setappdata(handles_gui.figure,'placing_with_mouse',false);
        update_sidemenu();
    end

    function update_sidemenu()
        % Get data
        handle_points = getappdata(handles_gui.figure,'handle_points');
        active_seed = getappdata(handles_gui.figure,'active_seed');
        placing_with_mouse = getappdata(handles_gui.figure,'placing_with_mouse');
        
        % Update text
        set(handles_gui.text_seedcount,'String',['Seed(s) Set: ' num2str(length(handle_points)) ' of ' num2str(total_threads)]) 

        % Update the active seed selector and full-resolution coordinate
        % fields. Ncorr stores seed coordinates using zero-based indexing.
        if (isempty(handle_points))
            setappdata(handles_gui.figure,'active_seed',[]);
            set(handles_gui.popup_activeseed,'String',{'None'},'Value',1,'Enable','off');
            set(handles_gui.edit_seedx,'String','');
            set(handles_gui.edit_seedy,'String','');
        else
            if (isempty(active_seed) || active_seed > length(handle_points))
                active_seed = length(handle_points);
                setappdata(handles_gui.figure,'active_seed',active_seed);
            end
            seed_labels = arrayfun(@(n)['Seed ' num2str(n)],1:length(handle_points),'UniformOutput',false);
            set(handles_gui.popup_activeseed,'String',seed_labels,'Value',active_seed,'Enable','on');
            pos_reduced = round(getPosition(handle_points(active_seed)))-1;
            pos_normal = pos_reduced*(spacing+1);
            set(handles_gui.edit_seedx,'String',num2str(pos_normal(1)));
            set(handles_gui.edit_seedy,'String',num2str(pos_normal(2)));
        end
        
        % Set buttons
        if (length(handle_points) == total_threads)
            % Seeds are finished
            set(handles_gui.button_setseeds,'Enable','off');
            set(handles_gui.button_addseed,'Enable','off');
            set(handles_gui.button_finish,'Enable','on'); 
        else
            set(handles_gui.button_setseeds,'Enable','on');
            set(handles_gui.button_addseed,'Enable','on');
            set(handles_gui.button_finish,'Enable','off'); 
        end
        if (isempty(handle_points))
            set(handles_gui.button_applycoords,'Enable','off');
        else
            set(handles_gui.button_applycoords,'Enable','on');
        end
        set(handles_gui.edit_seedx,'Enable','on');
        set(handles_gui.edit_seedy,'Enable','on');

        % Avoid re-entrant edits while impoint is waiting for mouse input.
        if (placing_with_mouse)
            set(handles_gui.button_setseeds,'Enable','off');
            set(handles_gui.button_addseed,'Enable','off');
            set(handles_gui.button_applycoords,'Enable','off');
            set(handles_gui.popup_activeseed,'Enable','off');
            set(handles_gui.edit_seedx,'Enable','off');
            set(handles_gui.edit_seedy,'Enable','off');
        end
    end

    function update_axes(action)
        % Get data    
        handle_preview = getappdata(handles_gui.figure,'handle_preview');  
        preview_threaddiagram = getappdata(handles_gui.figure,'preview_threaddiagram');    
        handle_subset = getappdata(handles_gui.figure,'handle_subset');
        preview_subset = getappdata(handles_gui.figure,'preview_subset');
        pos_preview = getappdata(handles_gui.figure,'pos_preview');
        
        if (strcmp(action,'set'))
            handle_preview = imshow(preview_threaddiagram,[reference.min_gs 2*reference.max_gs],'Parent',handles_gui.axes_setseeds);
            set(handles_gui.axes_setseeds,'Visible','off');
            set(handles_gui.text_name,'String',['Name: ' reference.name(1:min(end,22))]);
            
            handle_subset = imshow(preview_subset,[reference.min_gs 1.25*reference.max_gs],'Parent',handles_gui.axes_subset);
            set(handles_gui.axes_subset,'Visible','off');
        end
        
        if (strcmp(action,'set') || strcmp(action,'update'))
            % Update the overlay image
            set(handle_preview,'CData',preview_threaddiagram)
            
            % Update Subset Preview if a seed is active
            if ~isempty(pos_preview)
                preview_subset(:) = 0;
                gs_buffer = reference.get_gs();
                
                % Note: pos_preview is w.r.t the reduced ROI. Convert to normal coords.
                pos_normal = pos_preview * (spacing + 1);
                
                % Get cirroi at normal scale
                cirroi = roi.get_cirroi(pos_normal(1), pos_normal(2), radius, subsettrunc);
                
                for i = 0:size(cirroi.region.noderange,1)-1
                    x = i+(cirroi.x-cirroi.radius);
                    for j = 0:2:cirroi.region.noderange(i+1)-1
                        y_vec = cirroi.region.nodelist(i+1,j+1):cirroi.region.nodelist(i+1,j+2);
                        preview_subset(y_vec-(cirroi.y-cirroi.radius)+1,i+1) = gs_buffer(y_vec+1,x+1) + 0.25*reference.max_gs;                            
                    end     
                end    

                % Paint points over subset to illustrate subset spacing
                counter_column = mod(cirroi.x-cirroi.radius,spacing+1);
                if (counter_column == 0)
                    counter_column = spacing; 
                else
                    counter_column = counter_column-1;
                end
                
                for i = 0:size(cirroi.region.noderange,1)-1
                    counter_column = counter_column+1;
                    if (counter_column == spacing+1)
                        counter_column = 0; 
                    end
                    
                    x = i+(cirroi.x-cirroi.radius);
                    for j = 0:2:cirroi.region.noderange(i+1)-1
                        if (counter_column == 0)
                            counter_row = mod(cirroi.region.nodelist(i+1,j+1),spacing+1);
                            if (counter_row == 0)
                                counter_row = spacing;
                            else
                                counter_row = counter_row-1;
                            end
                        end
                        
                        for k = cirroi.region.nodelist(i+1,j+1):cirroi.region.nodelist(i+1,j+2) 
                            y = k;
                            if (counter_column == 0)
                                counter_row = counter_row + 1;
                                if (counter_row == spacing+1)  
                                    preview_subset(y-(cirroi.y-cirroi.radius)+1,i+1) = gs_buffer(y+1,x+1) + 0.5*reference.max_gs;
                                    counter_row = 0;
                                end
                            end
                        end            
                    end  
                end 
                set(handle_subset,'CData',preview_subset);
            end
        end
        
        % Set data 
        setappdata(handles_gui.figure,'handle_preview',handle_preview);   
        setappdata(handles_gui.figure,'handle_subset',handle_subset);
    end

    function get_threaddiagram(handle_points)
        % Get data
        ref_reduced = getappdata(handles_gui.figure,'ref_reduced');     
        regionmask = getappdata(handles_gui.figure,'regionmask');      
        threaddiagram_prelim = getappdata(handles_gui.figure,'threaddiagram_prelim');    
        preview_threaddiagram = getappdata(handles_gui.figure,'preview_threaddiagram');    
        
        % Draw thread diagram - get generators first
        generators = [];
        for i = 0:length(handle_points)-1
            % Convert to zero based indexing and round
            pos_buffer = round(getPosition(handle_points(i+1)))-1; 
            generators = vertcat(generators,pos_buffer); %#ok<AGROW>
        end
        
        % This function modifies the arrays in-place
        ncorr_alg_formthreaddiagram(threaddiagram_prelim,preview_threaddiagram,int32(generators),regionmask,ref_reduced.formatted());
    end

    function handles_gui = init_gui()
    % GUI controls -------------------------------------------------------%
        % Figure
        handles_gui.figure = figure( ...
            'Tag', 'figure', ...
            'Units', 'characters', ...
            'Position', ncorr_util_figpos(pos_parent,[35 175]), ...
            'Name', 'Set Seeds', ...
            'MenuBar', 'none', ...
            'NumberTitle', 'off', ...
            'Color', get(0,'DefaultUicontrolBackgroundColor'), ...
            'handlevisibility','off', ...
            'DockControls','off', ...
            'Resize','off', ...
            'WindowStyle','modal', ...
            'Visible','off', ...
            'IntegerHandle','off', ...
            'Interruptible','off');

        % Panels
        handles_gui.group_seedoptions = uipanel( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'group_seedoptions', ...
            'Units', 'characters', ...
            'Position', [2 19.0 35 15.5], ...
            'Title', 'Seed Options', ...
            'Interruptible','off');
        
        handles_gui.group_menu = uipanel( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'group_menu', ...
            'Units', 'characters', ...
            'Position', [2 11.7 35 6.4], ...
            'Title', 'Menu', ...
            'Interruptible','off');
        
        handles_gui.group_setseeds = uipanel( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'group_setseeds', ...
            'Units', 'characters', ...
            'Position', [39.0 0.75 92 33.75], ...
            'Title', 'Select Region', ...
            'Interruptible','off');
            
        handles_gui.group_preview = uipanel( ...
            'Parent', handles_gui.figure, ...
            'Tag', 'group_preview', ...
            'Units', 'characters', ...
            'Position', [133 0.75 40 33.75], ...
            'Title', 'Subset Preview', ...
            'Interruptible','off');

        % Axes
        handles_gui.axes_setseeds = axes( ...
            'Parent', handles_gui.group_setseeds, ...
            'Tag', 'axes_setseeds', ...
            'Units', 'characters', ...
            'Position', [2.4 3 86.2 28.8], ...
            'Visible', 'off', ...
            'Interruptible','off');
            
        handles_gui.axes_subset = axes( ...
            'Parent', handles_gui.group_preview, ...
            'Tag', 'axes_subset', ...
            'Units', 'characters', ...
            'Position', [2 2 36 29], ...
            'Visible','off', ...
            'Interruptible','off');

        % Static Texts
        handles_gui.text_seedcount = uicontrol( ...
            'Parent', handles_gui.group_seedoptions, ...
            'Tag', 'text_seedcount', ...
            'Style', 'text', ...
            'Units', 'characters', ...
            'Position', [2.2 12.0 29.7 1.5], ...
            'String', 'Seeds Set: 0 of ', ...
            'HorizontalAlignment', 'left', ...
            'Interruptible','off');

        handles_gui.text_activeseed = uicontrol( ...
            'Parent', handles_gui.group_seedoptions, ...
            'Tag', 'text_activeseed', ...
            'Style', 'text', ...
            'Units', 'characters', ...
            'Position', [2.2 9.8 12.5 1.5], ...
            'String', 'Active seed:', ...
            'HorizontalAlignment', 'left', ...
            'Interruptible','off');

        handles_gui.popup_activeseed = uicontrol( ...
            'Parent', handles_gui.group_seedoptions, ...
            'Tag', 'popup_activeseed', ...
            'Style', 'popupmenu', ...
            'Units', 'characters', ...
            'Position', [15.0 10.0 16.8 1.5], ...
            'String', {'None'}, ...
            'Enable', 'off', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_popup_activeseed,handles_gui.figure), ...
            'Interruptible','off');

        handles_gui.text_seedx = uicontrol( ...
            'Parent', handles_gui.group_seedoptions, ...
            'Tag', 'text_seedx', ...
            'Style', 'text', ...
            'Units', 'characters', ...
            'Position', [2.2 7.5 4.0 1.5], ...
            'String', 'X:', ...
            'HorizontalAlignment', 'left', ...
            'Interruptible','off');

        handles_gui.edit_seedx = uicontrol( ...
            'Parent', handles_gui.group_seedoptions, ...
            'Tag', 'edit_seedx', ...
            'Style', 'edit', ...
            'Units', 'characters', ...
            'Position', [5.2 7.8 10.0 1.5], ...
            'String', '', ...
            'HorizontalAlignment', 'left', ...
            'Interruptible','off');

        handles_gui.text_seedy = uicontrol( ...
            'Parent', handles_gui.group_seedoptions, ...
            'Tag', 'text_seedy', ...
            'Style', 'text', ...
            'Units', 'characters', ...
            'Position', [17.2 7.5 4.0 1.5], ...
            'String', 'Y:', ...
            'HorizontalAlignment', 'left', ...
            'Interruptible','off');

        handles_gui.edit_seedy = uicontrol( ...
            'Parent', handles_gui.group_seedoptions, ...
            'Tag', 'edit_seedy', ...
            'Style', 'edit', ...
            'Units', 'characters', ...
            'Position', [20.2 7.8 11.6 1.5], ...
            'String', '', ...
            'HorizontalAlignment', 'left', ...
            'Interruptible','off');

        handles_gui.text_coordinatehelp = uicontrol( ...
            'Parent', handles_gui.group_seedoptions, ...
            'Tag', 'text_coordinatehelp', ...
            'Style', 'text', ...
            'Units', 'characters', ...
            'Position', [2.2 5.5 29.7 1.5], ...
            'String', ['Full-res, 0-based; grid step ' num2str(spacing+1) ' px'], ...
            'HorizontalAlignment', 'left', ...
            'Interruptible','off');
        
        handles_gui.text_name = uicontrol( ...
            'Parent', handles_gui.group_setseeds, ...
            'Tag', 'text_name', ...
            'Style', 'text', ...
            'Units', 'characters', ...
            'Position', [3 0.7 56.9 1.3], ...
            'String', 'Name: ', ...
            'HorizontalAlignment', 'left', ...
            'Interruptible','off');    

        % Pushbuttons
        handles_gui.button_setseeds = uicontrol( ...
            'Parent', handles_gui.group_seedoptions, ...
            'Tag', 'button_setseeds', ...
            'Style', 'pushbutton', ...
            'Units', 'characters', ...
            'Position', [2.1 0.8 29.7 1.7], ...
            'String', 'Place Remaining with Mouse', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_button_setseeds,handles_gui.figure), ...
            'Interruptible','off');

        handles_gui.button_addseed = uicontrol( ...
            'Parent', handles_gui.group_seedoptions, ...
            'Tag', 'button_addseed', ...
            'Style', 'pushbutton', ...
            'Units', 'characters', ...
            'Position', [2.1 3.0 14.3 1.7], ...
            'String', 'Add Seed', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_button_addseed,handles_gui.figure), ...
            'Interruptible','off');

        handles_gui.button_applycoords = uicontrol( ...
            'Parent', handles_gui.group_seedoptions, ...
            'Tag', 'button_applycoords', ...
            'Style', 'pushbutton', ...
            'Units', 'characters', ...
            'Position', [17.5 3.0 14.3 1.7], ...
            'String', 'Apply to Active', ...
            'Enable', 'off', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_button_applycoords,handles_gui.figure), ...
            'Interruptible','off');
        
        handles_gui.button_finish = uicontrol( ...
            'Parent', handles_gui.group_menu, ...
            'Tag', 'button_finish', ...
            'Style', 'pushbutton', ...
            'Units', 'characters', ...
            'Position', [2.1 3 29.7 1.7], ...
            'String', 'Finish', ...
            'Enable', 'off', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_button_finish,handles_gui.figure), ...
            'Interruptible','off');

        handles_gui.button_cancel = uicontrol( ...
            'Parent', handles_gui.group_menu, ...
            'Tag', 'button_cancel', ...
            'Style', 'pushbutton', ...
            'Units', 'characters', ...
            'Position', [2.1 1 29.7 1.7], ...
            'String', 'Cancel', ...
            'Callback', ncorr_util_wrapcallbacktrycatch(@callback_button_cancel,handles_gui.figure), ...
            'Interruptible','off');
    end
        
    % Pause until figure is closed ---------------------------------------%
    waitfor(handles_gui.figure);    
end
