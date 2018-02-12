%% compute_hull: compute hull by 2d decomposition
function [C] = compute_hull(f, v_ival, sanity_check)

    if nargin < 2
        sanity_check = 0
    end

    Htot = zeros(0,length(f));
    htot = zeros(0,0);

    for idx = nchoosek(1:length(f), 2)'
        f1 = f(idx(1));
        f2 = f(idx(2));
        [H12, h12] = bounding_polytope_2d(f1, f2, v_ival, 2);
        H_new_blk = zeros(size(H12, 1), length(f));
        H_new_blk(:, idx) = H12;

        Htot = [Htot; H_new_blk];
        htot = [htot; h12];
    end

    C = Polyhedron(Htot, htot);
    C.minHRep;
    C.minVRep;

    if sanity_check
        % are points along curve contained in polytope?
        for val = v_ival(1):range(v_ival)/10:v_ival(2)
            pt = double(subs(f', val));
            assert(C.contains(pt));
        end
    end

%% bounding_polytope: compute hull around monotone 2D curve
function [H,h] = bounding_polytope_2d(f1, f2, ival, num_pts)

    % TODO: add verification that curve is monotone!

    % endpoint values
    f11 = double(subs(f1, ival(1)));
    f12 = double(subs(f1, ival(2)));
    f21 = double(subs(f2, ival(1)));
    f22 = double(subs(f2, ival(2)));

    % start with hyperbox
    H = [eye(2); -eye(2)];
    h = [max(f11, f12); max(f21, f22); -min(f11, f12); -min(f21, f22)];

    if range(ival) == 0
        return
    end

    % add diagonal constraint
    sol = null([f11 f21 -1; f12 f22 -1]);
    c = sol(1:2)'; % normal vector s.t. c x = d
    d = sol(3);    % offset

    % If positive, curving right, if negative, curving left
    direction = sign(double(subs([diff(f1, 2), diff(f2, 2)], (ival(2) + ival(1))/2)) * c');

    H = [H; direction*sol(1:2)'];
    h = [h; direction*sol(3)];

    % add gradient constraints
    for pt = ival(1) : range(ival)/num_pts : ival(2);
        grad = double(subs([diff(f1, 1), diff(f2, 1)], pt));
        slope = null(grad)';
        direction = sign(slope * double(subs([diff(f1, 2); diff(f2, 2)], pt)));
        int = slope * double([subs(f1, pt); subs(f2, pt)]);
        H = [H; -direction*slope];
        h = [h; -direction*int];
    end

    % sanity check
    for val = ival(1) : range(ival)/(2*num_pts) : ival(2)
        assert( all(H*double(subs([f1; f2], val)) <= h + 1e-6) )
    end

