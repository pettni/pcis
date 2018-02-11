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