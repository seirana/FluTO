function [model, flipped] = canonicalizeReactionDirections(model, tolerance)
%CANONICALIZEREACTIONDIRECTIONS Flip non-positive irreversible reactions.
%
%   [MODEL, FLIPPED] = fluto.canonicalizeReactionDirections(MODEL, TOL)
%   changes reactions constrained to non-positive flux into an equivalent
%   non-negative orientation. Stoichiometric columns and bounds are flipped
%   together so the feasible flux space is preserved under the sign change.

    if nargin < 2 || isempty(tolerance)
        tolerance = 1e-9;
    end
    validateattributes(tolerance, {"numeric"}, {"scalar", "real", "nonnegative", "finite"});

    model = fluto.validateModel(model);
    flipped = model.lb < -tolerance & model.ub <= tolerance;

    if any(flipped)
        oldLower = model.lb(flipped);
        oldUpper = model.ub(flipped);

        model.lb(flipped) = -oldUpper;
        model.ub(flipped) = -oldLower;
        model.S(:, flipped) = -model.S(:, flipped);
    end
end
