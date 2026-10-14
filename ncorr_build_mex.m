function build_info = ncorr_build_mex(support_openmp,total_cores)
%NCORR_BUILD_MEX Build this checkout with the separate-complex MATLAB API.
if ~(isscalar(support_openmp) && (islogical(support_openmp) || isnumeric(support_openmp)) && isreal(support_openmp) && any(support_openmp == [0 1])) || ...
        ~(isnumeric(total_cores) && isscalar(total_cores) && isreal(total_cores) && isfinite(total_cores) && total_cores >= 1 && total_cores == fix(total_cores))
    error('Ncorr:build:invalidOptions','Expected OpenMP 0 or 1 and a finite positive integer core count.');
end
support_openmp = logical(support_openmp);
total_cores = double(total_cores);
root = fileparts(mfilename('fullpath'));
stamp = fullfile(root,'ncorr_buildinfo.mat');
if exist(stamp,'file'), delete(stamp); end
compiler = mex.getCompilerConfigurations('C++','Selected');
if isempty(compiler)
    error('Ncorr:build:noCompiler','Run mex -setup C++ in this MATLAB release.');
end
flags = {};
openmp_flags = {};
if contains(compiler.Manufacturer,'GNU','IgnoreCase',true)
    flags = {'CXXFLAGS=$CXXFLAGS -std=c++11'};
    if support_openmp
        openmp_flags = {'CXXFLAGS=$CXXFLAGS -std=c++11 -fopenmp -DNCORR_OPENMP','CXXLIBS=$CXXLIBS -lgomp'};
    end
elseif contains(compiler.Manufacturer,'Microsoft','IgnoreCase',true)
    if support_openmp
        openmp_flags = {'COMPFLAGS=$COMPFLAGS /openmp /DNCORR_OPENMP'};
    end
elseif support_openmp
    error('Ncorr:build:unsupportedOpenMP','OpenMP flags are not mapped for compiler %s.',compiler.Name);
end
if ~support_openmp, openmp_flags = flags; end
libs = {'standard_datatypes','ncorr_datatypes','ncorr_lib'};
funcs = {'ncorr_alg_formmask','ncorr_alg_formregions','ncorr_alg_formboundary', ...
    'ncorr_alg_formthreaddiagram','ncorr_alg_formunion','ncorr_alg_extrapdata', ...
    'ncorr_alg_adddisp','ncorr_alg_convert','ncorr_alg_dispgrad', ...
    'ncorr_alg_testopenmp','ncorr_alg_calcseeds','ncorr_alg_rgdic'};
for k = 1:numel(funcs), clear(funcs{k}); end
if ispc, objext = '.obj'; else, objext = '.o'; end
objects = cellfun(@(name) fullfile(root,[name objext]),libs,'UniformOutput',false);
common = {'-R2017b',['-I' root],'-outdir',root};
for k = 1:numel(libs)
    source = fullfile(root,[libs{k} '.cpp']);
    require_file(source);
    require_file(fullfile(root,[libs{k} '.h']));
    fprintf('Installing %s...\n',libs{k});
    mex(common{:},flags{:},'-c',source);
end
for k = 1:numel(objects), require_file(objects{k}); end
for k = 1:numel(funcs)
    source = fullfile(root,[funcs{k} '.cpp']);
    require_file(source);
    unit_flags = flags;
    if k > 9, unit_flags = openmp_flags; end
    fprintf('Installing %s...\n',funcs{k});
    mex(common{:},unit_flags{:},source,objects{:},'-output',funcs{k});
end
% Resolve the probe from this root even when the caller is elsewhere.
oldpath = path;
restore_path = onCleanup(@() path(oldpath)); %#ok<NASGU>
addpath(root,'-begin');
if support_openmp && ~ncorr_alg_testopenmp()
    error('Ncorr:build:openMPUnavailable','Compiled OpenMP probe returned false.');
end
sources = [dir(fullfile(root,'*.cpp')); dir(fullfile(root,'*.h')); ...
    dir(fullfile(root,'ncorr.m')); dir(fullfile(root,'ncorr_build_mex.m'))];
[~,order] = sort({sources.name});
sources = sources(order);
build_info = struct('matlab_release',version('-release'),'mex_extension',mexext, ...
    'compiler_name',compiler.Name,'compiler_version',compiler.Version, ...
    'compiler_location',compiler.Location,'api','-R2017b', ...
    'support_openmp',support_openmp,'total_cores',total_cores, ...
    'source_files',{{sources.name}},'source_bytes',[sources.bytes], ...
    'source_datenums',[sources.datenum],'openmp_flags',{openmp_flags});
% Publish the success stamp last, so any failure forces a fresh installation.
dlmwrite(fullfile(root,'ncorr_installinfo.txt'),[support_openmp total_cores]);
save(stamp,'build_info');
end

function require_file(filename)
if ~isfile(filename)
    error('Ncorr:build:missingFile','Required build input is missing: %s',filename);
end
end
