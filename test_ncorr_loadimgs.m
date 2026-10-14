function test_ncorr_loadimgs()
%TEST_NCORR_LOADIMGS Multi-underscore names retain numeric image order.
% Substitute only the file picker; load real TIFFs through the production code.
work = tempname;
mkdir(work);
oldpath = path;
cleanup = onCleanup(@() restore_environment(work,oldpath));
names = {'DIC_B_010.tif','DIC_B_002.tif','DIC_B_001.tif'};
values = [30 20 10];
for k = 1:numel(names)
    imwrite(repmat(uint16(values(k)),32,32),fullfile(work,names{k}));
end
% Keep file selection deterministic, and fail immediately on an unexpected
% dialog rather than hanging a noninteractive regression run.
picker = sprintf(['function [names,folder] = uigetfile(varargin)\n' ...
    'names = {''DIC_B_010.tif'',''DIC_B_002.tif'',''DIC_B_001.tif''};\n' ...
    'folder = ''%s'';\nend\n'],strrep(work,'''',''''''));
write_helper(fullfile(work,'uigetfile.m'),picker);
write_helper(fullfile(work,'errordlg.m'),sprintf(['function h = errordlg(varargin)\n' ...
    'error(''Ncorr:test:imageLoadDialog'',''%%s'',varargin{1});\nend\n']));
addpath(work,'-begin');
clear uigetfile errordlg
for lazy = [false true]
    [imgs,state] = ncorr_util_loadimgs(lazy);
    assert(state == out.success && numel(imgs) == 3, ...
        'Ncorr:test:imageLoad','Multi-image loading did not succeed (lazy=%d).',lazy);
    expected = [10 20 30]/65535;
    for k = 1:3
        gs = imgs(k).get_gs();
        assert(all(abs(gs(:)-expected(k)) <= eps), ...
            'Ncorr:test:imageOrder','Image content is not in numeric suffix order (lazy=%d, index=%d).',lazy,k);
    end
end
fprintf('PASS: eager/lazy multi-underscore TIFF selection loads in numeric order.\n');
end

function write_helper(filename,text)
fid = fopen(filename,'w');
assert(fid ~= -1,'Ncorr:test:temporaryFile','Cannot create %s.',filename);
close_file = onCleanup(@() fclose(fid));
fprintf(fid,'%s',text);
end

function restore_environment(work,oldpath)
path(oldpath);
clear uigetfile errordlg
rmdir(work,'s');
end
