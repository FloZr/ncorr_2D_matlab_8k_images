function test_ncorr_rgdic_spacing()
%TEST_NCORR_RGDIC_SPACING Prove RG-DIC spacing is independent of call order.
% Run against freshly built ncorr_alg_rgdic. The inherited speedup source
% (git show 8c5462f:ncorr_alg_rgdic.cpp) should fail the warm 0 -> 5 mask check.

previous_rng = rng;
restore_rng = onCleanup(@() rng(previous_rng)); %#ok<NASGU>
rng(7,'twister');
reference_image = imgaussfilt(rand(384),0.7);
[x,y] = meshgrid(0:383,0:383);
current_image = interp2(x,y,reference_image,192+(x-192)/1.25,y,'linear',0);
reference = ncorr_class_img();
reference.set_img('load',struct('img',reference_image,'name','','path',''));
current = ncorr_class_img();
current.set_img('load',struct('img',current_image,'name','','path',''));
mask = false(384);
mask(65:320,65:320) = true;
roi = ncorr_class_roi();
roi.set_roi('load',struct('mask',mask,'cutoff',0));
radius = 21;
[seed,convergence,state] = ncorr_alg_calcseeds(reference.formatted(), ...
    current.formatted(),roi.formatted(),int32(0),int32([192 192]), ...
    int32(radius),1e-6,int32(100),false,false);
assert(state == out.success,'Ncorr:test:spacingSeed', ...
    'Seed failed: state=%g, diffnorm=%.17g, iterations=%g.', ...
    state,convergence.diffnorm,convergence.num_iterations);
seed.num_region = int32(0);
seed.num_thread = int32(0);
fprintf('Spacing fixture seed: diffnorm=%.6g, iterations=%g.\n', ...
    convergence.diffnorm,convergence.num_iterations);

% Clearing between cold calls resets the MEX-local static state. Never clear
% inside a warm pair: the second call must inherit the first call's state.
clear ncorr_alg_rgdic
cold0 = run_spacing(reference,current,roi,seed,radius,0,'cold 0');
clear ncorr_alg_rgdic
cold5 = run_spacing(reference,current,roi,seed,radius,5,'cold 5');
check_analytic(cold5,roi,radius);

clear ncorr_alg_rgdic
warm0 = run_spacing(reference,current,roi,seed,radius,0,'0 -> 5, first');
warm5 = run_spacing(reference,current,roi,seed,radius,5,'0 -> 5, second');
compare_result(cold0,warm0,0,'0 -> 5');
compare_result(cold5,warm5,5,'0 -> 5');

clear ncorr_alg_rgdic
warm5 = run_spacing(reference,current,roi,seed,radius,5,'5 -> 0, first');
warm0 = run_spacing(reference,current,roi,seed,radius,0,'5 -> 0, second');
compare_result(cold5,warm5,5,'5 -> 0');
compare_result(cold0,warm0,0,'5 -> 0');
clear ncorr_alg_rgdic
fprintf('PASS: mixed-spacing masks match and displacement differences <= 1e-10 px.\n');
end

function raw = run_spacing(reference,current,roi,seed,radius,spacing,label)
reduced_roi = roi.reduce(spacing);
regionmask = reduced_roi.get_regionmask(0);
% One region and one thread: assign thread 0 only to its reduced-grid points.
diagram = -ones(size(reduced_roi.mask),'int32');
diagram(regionmask) = int32(0);
seed.computepoints = int32(nnz(regionmask));
[raw,state] = ncorr_alg_rgdic(reference.formatted(),current.formatted(), ...
    roi.formatted(),seed,diagram,int32(radius),int32(spacing),1e-6, ...
    int32(100),false,int32(0),int32(1));
assert(state == out.success,'Ncorr:test:spacingAnalysis', ...
    '%s (spacing %d) failed: state=%g.',label,spacing,state);
assert(any(raw.plot_validpoints(:)),'Ncorr:test:spacingEmpty', ...
    '%s (spacing %d) returned no valid points.',label,spacing);
fprintf('%s: spacing=%d, valid=%d/%d reduced-region points.\n', ...
    label,spacing,nnz(raw.plot_validpoints),nnz(regionmask));
end

function check_analytic(raw,roi,radius)
spacing = 5;
reduced_roi = roi.reduce(spacing);
[x,y] = meshgrid((0:size(raw.plot_u,2)-1)*(spacing+1), ...
    (0:size(raw.plot_u,1)-1)*(spacing+1));
% The full-resolution mask bounds are [64,319] in zero-based coordinates.
interior = reduced_roi.get_regionmask(0) & ...
    x >= 64+radius & x <= 319-radius & ...
    y >= 64+radius & y <= 319-radius;
valid = logical(raw.plot_validpoints) & interior;
coverage = nnz(valid)/nnz(interior);
assert(coverage >= 0.8,'Ncorr:test:spacingCoverage', ...
    'Cold spacing 5 interior coverage %.6f (%d/%d) is below 0.8.', ...
    coverage,nnz(valid),nnz(interior));
u_error = abs(raw.plot_u(valid)-0.25*(x(valid)-192));
v_error = abs(raw.plot_v(valid));
assert(all(isfinite(u_error)) && all(isfinite(v_error)), ...
    'Ncorr:test:spacingFinite','Cold spacing 5 has nonfinite valid displacement.');
median_u_error = median(u_error);
median_v_error = median(v_error);
fprintf('Cold spacing 5: interior coverage=%.6f, median |u error|=%.6g px, |v error|=%.6g px.\n', ...
    coverage,median_u_error,median_v_error);
assert(median_u_error <= 0.15 && median_v_error <= 0.05, ...
    'Ncorr:test:spacingAccuracy', ...
    'Cold spacing 5 median errors exceed limits: u=%.17g (0.15), v=%.17g (0.05) px.', ...
    median_u_error,median_v_error);
end

function compare_result(cold,warm,spacing,order)
assert(isequal(cold.plot_validpoints,warm.plot_validpoints), ...
    'Ncorr:test:spacingMask', ...
    '%s spacing %d mask differs from cold: cold=%d valid, warm=%d valid.', ...
    order,spacing,nnz(cold.plot_validpoints),nnz(warm.plot_validpoints));
valid = logical(cold.plot_validpoints);
u_delta = abs(cold.plot_u(valid)-warm.plot_u(valid));
v_delta = abs(cold.plot_v(valid)-warm.plot_v(valid));
assert(all(isfinite(u_delta)) && all(isfinite(v_delta)), ...
    'Ncorr:test:spacingFinite', ...
    '%s spacing %d has nonfinite valid displacement differences.',order,spacing);
u_difference = max(u_delta);
v_difference = max(v_delta);
fprintf('%s spacing %d: max |delta u|=%.17g px, |delta v|=%.17g px.\n', ...
    order,spacing,u_difference,v_difference);
assert(isfinite(u_difference) && isfinite(v_difference) && ...
    u_difference <= 1e-10 && v_difference <= 1e-10, ...
    'Ncorr:test:spacingDisplacement', ...
    '%s spacing %d differs from cold: u=%.17g, v=%.17g px (limit 1e-10).', ...
    order,spacing,u_difference,v_difference);
end
