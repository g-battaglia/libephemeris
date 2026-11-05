//! Search implementation for the parent module.
//! Public entry points and types are exposed through the parent API.

use super::*;

impl SearchConfig {
    /// Checked constructor; invalid fields are [`Error`], never silent
    /// clamping.
    pub fn new(tolerance: f64, min_step: f64, max_iter: usize) -> Result<Self> {
        if !tolerance.is_finite() || tolerance <= 0.0 {
            return Err(Error::InvalidSearchTolerance { value: tolerance });
        }
        if !min_step.is_finite() || min_step <= 0.0 {
            return Err(Error::InvalidSearchStep { value: min_step });
        }
        if max_iter == 0 {
            return Err(Error::InvalidSearchIterations { max_iter });
        }
        Ok(Self {
            tolerance,
            min_step,
            max_iter,
        })
    }
}

impl NewtonConfig {
    /// Checked constructor.
    pub fn new(tolerance: f64, deriv_floor: f64, max_iter: usize) -> Result<Self> {
        if !tolerance.is_finite() || tolerance <= 0.0 {
            return Err(Error::InvalidSearchTolerance { value: tolerance });
        }
        if !deriv_floor.is_finite() || deriv_floor <= 0.0 {
            return Err(Error::InvalidSearchStep { value: deriv_floor });
        }
        if max_iter == 0 {
            return Err(Error::InvalidSearchIterations { max_iter });
        }
        Ok(Self {
            tolerance,
            deriv_floor,
            max_iter,
        })
    }
}

pub(super) fn check_interval(a: f64, b: f64) -> Result<()> {
    if !a.is_finite() || !b.is_finite() || a >= b {
        return Err(Error::InvalidSearchInterval { a, b });
    }
    Ok(())
}

pub(super) fn eval_finite(f: &mut impl FnMut(f64) -> f64, x: f64) -> Result<f64> {
    let fx = f(x);
    if !fx.is_finite() {
        return Err(Error::NonFiniteSearchValue { x });
    }
    Ok(fx)
}

/// Validate the bracket `[a, b]` for `f` under `cfg`.
///
/// Evaluates both endpoints exactly once and reports [`BracketOutcome`].
/// Invalid intervals or non-finite endpoint values are `Err`.
pub fn check_bracket(
    f: &mut impl FnMut(f64) -> f64,
    a: f64,
    b: f64,
    cfg: SearchConfig,
) -> Result<BracketOutcome> {
    check_interval(a, b)?;
    let fa = eval_finite(f, a)?;
    let fb = eval_finite(f, b)?;
    if fa.abs() <= cfg.tolerance {
        return Ok(BracketOutcome::EndpointRoot { x: a, fx: fa });
    }
    if fb.abs() <= cfg.tolerance {
        return Ok(BracketOutcome::EndpointRoot { x: b, fx: fb });
    }
    if fa * fb < 0.0 {
        Ok(BracketOutcome::Bracketed { a, b, fa, fb })
    } else {
        Ok(BracketOutcome::NoBracket { fa, fb })
    }
}

/// Bisect the bracket `[a, b]` to the function tolerance.
///
/// Requires a strict sign change (see [`check_bracket`]); otherwise the
/// `NoBracket` (or endpoint) outcome is returned without iterating.
pub fn bisection(
    f: &mut impl FnMut(f64) -> f64,
    a: f64,
    b: f64,
    cfg: SearchConfig,
) -> Result<RootOutcome> {
    let (mut lo, mut hi, mut flo, mut fhi) = match check_bracket(f, a, b, cfg)? {
        BracketOutcome::Bracketed { a, b, fa, fb } => (a, b, fa, fb),
        BracketOutcome::EndpointRoot { x, fx } => {
            return Ok(RootOutcome::EndpointRoot { x, fx });
        }
        BracketOutcome::NoBracket { fa, fb } => return Ok(RootOutcome::NoBracket { fa, fb }),
    };
    for iter in 1..=cfg.max_iter {
        let mid = 0.5 * (lo + hi);
        let fmid = eval_finite(f, mid)?;
        if fmid.abs() <= cfg.tolerance {
            return Ok(RootOutcome::Converged {
                x: mid,
                fx: fmid,
                iterations: iter,
            });
        }
        if flo * fmid < 0.0 {
            hi = mid;
            fhi = fmid;
        } else {
            lo = mid;
            flo = fmid;
        }
        if hi - lo <= cfg.min_step {
            // Root localized to `min_step`; report the closer endpoint.
            let (x, fx) = if flo.abs() <= fhi.abs() {
                (lo, flo)
            } else {
                (hi, fhi)
            };
            return Ok(RootOutcome::Converged {
                x,
                fx,
                iterations: iter,
            });
        }
    }
    let (best_x, best_fx) = if flo.abs() <= fhi.abs() {
        (lo, flo)
    } else {
        (hi, fhi)
    };
    Ok(RootOutcome::NonConverged {
        best_x,
        best_fx,
        iterations: cfg.max_iter,
    })
}

/// Refine the bracket `[a, b]` with Brent's method (Brent–Dekker).
///
/// Same bracket gate as [`bisection`]; convergence is declared on
/// `|f| <= tolerance` only. On cap exhaustion the best estimate is
/// reported as [`RootOutcome::NonConverged`], never as a root.
pub fn brent(
    f: &mut impl FnMut(f64) -> f64,
    a: f64,
    b: f64,
    cfg: SearchConfig,
) -> Result<RootOutcome> {
    let (mut xa, mut xb, mut fa, mut fb) = match check_bracket(f, a, b, cfg)? {
        BracketOutcome::Bracketed { a, b, fa, fb } => (a, b, fa, fb),
        BracketOutcome::EndpointRoot { x, fx } => {
            return Ok(RootOutcome::EndpointRoot { x, fx });
        }
        BracketOutcome::NoBracket { fa, fb } => return Ok(RootOutcome::NoBracket { fa, fb }),
    };
    // `xb` holds the better estimate: |fa| >= |fb|.
    if fa.abs() < fb.abs() {
        std::mem::swap(&mut xa, &mut xb);
        std::mem::swap(&mut fa, &mut fb);
    }
    let mut c = xa;
    let mut fc = fa;
    let mut mflag = true;
    let mut d = 0.0;
    for iter in 1..=cfg.max_iter {
        if fb.abs() <= cfg.tolerance {
            return Ok(RootOutcome::Converged {
                x: xb,
                fx: fb,
                iterations: iter - 1,
            });
        }
        let s = if fa != fc && fb != fc {
            // Inverse quadratic interpolation.
            xa * fb * fc / ((fa - fb) * (fa - fc))
                + xb * fa * fc / ((fb - fa) * (fb - fc))
                + c * fa * fb / ((fc - fa) * (fc - fb))
        } else {
            // Secant step.
            xb - fb * (xb - xa) / (fb - fa)
        };
        // Bisection fallback guards.
        let in_cone = (3.0 * xa + xb) / 4.0 < s && s < xb || xb < s && s < (3.0 * xa + xb) / 4.0;
        let cond1 = !in_cone;
        let cond2 = mflag && (s - xb).abs() >= (xb - c).abs() / 2.0;
        let cond3 = !mflag && (s - xb).abs() >= (c - d).abs() / 2.0;
        let cond4 = mflag && (xb - c).abs() < cfg.min_step;
        let cond5 = !mflag && (c - d).abs() < cfg.min_step;
        let s = if cond1 || cond2 || cond3 || cond4 || cond5 {
            mflag = true;
            0.5 * (xa + xb)
        } else {
            mflag = false;
            s
        };
        if !s.is_finite() {
            return Ok(RootOutcome::NonConverged {
                best_x: xb,
                best_fx: fb,
                iterations: iter,
            });
        }
        let fs = eval_finite(f, s)?;
        d = c;
        c = xb;
        fc = fb;
        if fa * fs < 0.0 {
            xb = s;
            fb = fs;
        } else {
            xa = s;
            fa = fs;
        }
        if fa.abs() < fb.abs() {
            std::mem::swap(&mut xa, &mut xb);
            std::mem::swap(&mut fa, &mut fb);
        }
    }
    Ok(RootOutcome::NonConverged {
        best_x: xb,
        best_fx: fb,
        iterations: cfg.max_iter,
    })
}

/// Newton refinement from `x0` with a derivative floor.
///
/// `fdf` returns `(f(x), f'(x))`. Convergence is `|f| <= tolerance`. A
/// step from `|df| < deriv_floor`, or a non-finite Newton step, ends the
/// search as [`RootOutcome::NonConverged`] with the last finite point —
/// the caller switches to a bracketed method there (station guard).
/// Non-finite `f`/`df` evaluations are [`Error::NonFiniteSearchValue`].
pub fn newton(
    fdf: &mut impl FnMut(f64) -> (f64, f64),
    x0: f64,
    cfg: NewtonConfig,
) -> Result<RootOutcome> {
    if !x0.is_finite() {
        return Err(Error::NonFiniteSearchValue { x: x0 });
    }
    let mut x = x0;
    let mut last = x0;
    let mut last_f = f64::INFINITY;
    for iter in 1..=cfg.max_iter {
        let (fx, dfx) = fdf(x);
        if !fx.is_finite() {
            return Err(Error::NonFiniteSearchValue { x });
        }
        if !dfx.is_finite() {
            return Err(Error::NonFiniteSearchValue { x });
        }
        last = x;
        last_f = fx;
        if fx.abs() <= cfg.tolerance {
            return Ok(RootOutcome::Converged {
                x,
                fx,
                iterations: iter - 1,
            });
        }
        if dfx.abs() < cfg.deriv_floor {
            return Ok(RootOutcome::NonConverged {
                best_x: x,
                best_fx: fx,
                iterations: iter - 1,
            });
        }
        let step = fx / dfx;
        let next = x - step;
        if !next.is_finite() {
            return Ok(RootOutcome::NonConverged {
                best_x: x,
                best_fx: fx,
                iterations: iter,
            });
        }
        x = next;
    }
    Ok(RootOutcome::NonConverged {
        best_x: last,
        best_fx: last_f,
        iterations: cfg.max_iter,
    })
}

/// Scan `[lo, hi]` uniformly for the first strict sign change.
///
/// `count` sub-intervals (`count + 1` samples, endpoints exact). A sample
/// within `tolerance` is [`ScanOutcome::EndpointRoot`]; the first adjacent
/// sign change is [`ScanOutcome::Bracket`]; otherwise [`ScanOutcome::NoRoot`].
pub fn scan_first_bracket(
    f: &mut impl FnMut(f64) -> f64,
    lo: f64,
    hi: f64,
    count: usize,
    tolerance: f64,
) -> Result<ScanOutcome> {
    check_interval(lo, hi)?;
    if count == 0 {
        return Err(Error::InvalidScanSamples { count });
    }
    if !tolerance.is_finite() || tolerance <= 0.0 {
        return Err(Error::InvalidSearchTolerance { value: tolerance });
    }
    let step = (hi - lo) / count as f64;
    let mut prev_x = lo;
    let mut prev_f = eval_finite(f, lo)?;
    if prev_f.abs() <= tolerance {
        return Ok(ScanOutcome::EndpointRoot { x: lo, fx: prev_f });
    }
    for i in 1..=count {
        let cur_x = if i == count { hi } else { lo + i as f64 * step };
        let cur_f = eval_finite(f, cur_x)?;
        if cur_f.abs() <= tolerance {
            return Ok(ScanOutcome::EndpointRoot {
                x: cur_x,
                fx: cur_f,
            });
        }
        if prev_f * cur_f < 0.0 {
            return Ok(ScanOutcome::Bracket {
                a: prev_x,
                b: cur_x,
                fa: prev_f,
                fb: cur_f,
            });
        }
        prev_x = cur_x;
        prev_f = cur_f;
    }
    Ok(ScanOutcome::NoRoot)
}
