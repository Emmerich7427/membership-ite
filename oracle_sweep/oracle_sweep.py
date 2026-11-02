#!/usr/bin/env python3
"""Population calculation for Section 5.5, with unrestricted cross-world coupling.

The potential outcomes share a symmetric Gaussian-mixture marginal. Calibration
uses P(C-Y subset L_tau), evaluated from Gaussian probabilities. Closed level
sets include score ties; probability grids serve only as diagnostic checks.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.special import ndtr, ndtri

ALPHA_CF = 0.025
ALPHA_INC = 0.025
N_ROOT_GRID = 12001
OUTPUT_DIR = Path(__file__).resolve().parent / 'oracle_sweep_output'
ROOT_TOL = 2e-13


def phi(x):
    return np.exp(-0.5 * np.asarray(x) ** 2) / np.sqrt(2 * np.pi)


def mix_pdf(x, mu):
    return 0.5 * (phi(np.asarray(x) - mu) + phi(np.asarray(x) + mu))


def normal_mass(lo, hi):
    """Stable normal interval probabilities, including in the right tail."""
    lo, hi = np.asarray(lo), np.asarray(hi)
    return np.where(lo > 0, ndtr(-lo) - ndtr(-hi), ndtr(hi) - ndtr(lo))


def mix_mass(lo, hi, mu):
    return 0.5 * (normal_mass(np.asarray(lo) - mu, np.asarray(hi) - mu)
                  + normal_mass(np.asarray(lo) + mu, np.asarray(hi) + mu))


def mix_cdf(x, mu):
    return 0.5 * (ndtr(np.asarray(x) - mu) + ndtr(np.asarray(x) + mu))


def mix_ppf(p, mu):
    if not 0 < p < 1:
        raise ValueError('The probability must lie strictly between zero and one.')
    bound = mu + max(12.0, abs(float(ndtri(min(p, 1 - p)))) + 2)
    return brentq(lambda x: float(mix_cdf(x, mu) - p), -bound, bound, xtol=ROOT_TOL)


def merge_intervals(intervals, split_at=()):
    """Merge intervals, retaining deleted points as barriers when requested.

    Strict level-set components may have touching closures. Those components
    must stay separate for whole-set containment. No tolerance fills gaps.
    """
    result = []
    barriers = set(split_at)
    for lo, hi in sorted((float(lo), float(hi)) for lo, hi in intervals if hi >= lo):
        if result and (lo < result[-1][1]
                       or (lo == result[-1][1] and lo not in barriers)):
            result[-1][1] = max(result[-1][1], hi)
        else:
            result.append([lo, hi])
    return result


def intersect_intervals(first, second):
    return merge_intervals([(max(a, c), min(b, d)) for a, b in first
                            for c, d in second if max(a, c) <= min(b, d)])


def interval_length(intervals):
    return float(sum(b - a for a, b in intervals))


def mixture_mass(intervals, mu):
    return float(sum(mix_mass(a, b, mu) for a, b in intervals))


def oracle_hpd_set(mu, alpha_cf=ALPHA_CF):
    """Use the known stationary points of the symmetric mixture density."""
    mode = (brentq(lambda x: x - mu * np.tanh(mu * x), 1e-8, mu,
                   xtol=ROOT_TOL) if mu > 1 + 1e-8 else 0.0)
    upper = mu + 14.0
    peak = float(mix_pdf(mode, mu))

    def endpoints(level):
        outer = brentq(lambda x: float(mix_pdf(x, mu) - level), mode, upper, xtol=ROOT_TOL)
        inner = (brentq(lambda x: float(mix_pdf(x, mu) - level), 0, mode, xtol=ROOT_TOL)
                 if level > mix_pdf(0, mu) else 0.0)
        return inner, outer

    def excess_mass(level):
        u, v = endpoints(level)
        return 2 * float(mix_mass(u, v, mu)) - (1 - alpha_cf)

    level = brentq(excess_mass, 1e-16, peak * (1 - 1e-13), xtol=2e-15)
    u, v = endpoints(level)
    return ([[-v, -u], [u, v]] if u > 0 else [[-v, v]]), level


def membership(d, mu, C):
    return sum(mix_mass(a - np.asarray(d), b - np.asarray(d), mu) for a, b in C)


def membership_derivative(d, mu, C):
    return sum(mix_pdf(a - np.asarray(d), mu) - mix_pdf(b - np.asarray(d), mu) for a, b in C)


def lei_first_stage_halfwidth(mu, alpha_cf=ALPHA_CF):
    return mix_ppf(1 - alpha_cf / 2, mu)


def lei_exante_length(mu, alpha_cf=ALPHA_CF, alpha_inc=ALPHA_INC):
    return 2 * (mix_ppf(1 - alpha_cf / 2, mu) + mix_ppf(1 - alpha_inc / 2, mu))


class OracleProblem:
    def __init__(self, mu, alpha_cf=ALPHA_CF, alpha_inc=ALPHA_INC, n_root_grid=N_ROOT_GRID):
        if mu <= 0 or not (1e-8 < alpha_cf < 1 - 1e-8 and 1e-8 < alpha_inc < 1 - 1e-8):
            raise ValueError('Require mu > 0 and error rates in (1e-8, 1-1e-8).')
        self.mu, self.alpha_cf, self.alpha_inc = mu, alpha_cf, alpha_inc
        self.C, self.hpd_level = oracle_hpd_set(mu, alpha_cf)
        # Beyond max(C)+mu every component interval probability is decreasing.
        tail_start = self.C[-1][1] + mu
        self.bound = tail_start + 14.0
        grid = np.linspace(1e-7, tail_start, n_root_grid)
        values = self.dm(grid)
        indices = np.flatnonzero(np.signbit(values[:-1]) != np.signbit(values[1:]))
        self.critical = sorted(set(brentq(self.dm, grid[i], grid[i + 1], xtol=ROOT_TOL)
                                   for i in indices))
        self.valleys = [x for x in self.critical if self.dm(x - 1e-5) < 0 < self.dm(x + 1e-5)]

    def m(self, d):
        # Enforce the analytic symmetry so opposite stationary points share
        # the identical floating-point level when the score has an atom.
        return membership(np.abs(d), self.mu, self.C)

    def dm(self, d):
        return membership_derivative(d, self.mu, self.C)

    def level_set(self, level, strict=False):
        if level <= 0:
            return [[-np.inf, np.inf]]
        if level > self.m(0):
            return []
        bound = self.bound
        while self.m(bound) >= level:
            bound += 8
        nodes = [0.0] + self.critical + [bound]
        values = [float(self.m(x) - level) for x in nodes]
        roots = [x for x, v in zip(nodes, values) if v == 0]
        for left, right, vl, vr in zip(nodes[:-1], nodes[1:], values[:-1], values[1:]):
            if (vl < 0 < vr) or (vr < 0 < vl):
                roots.append(brentq(lambda d: float(self.m(d) - level), left, right, xtol=ROOT_TOL))
        points = sorted(set([0.0, bound] + roots))
        positive = [[lo, hi] for lo, hi in zip(points[:-1], points[1:])
                    if self.m((lo + hi) / 2) > level]
        if not strict:
            positive += [[x, x] for x in roots]
        pieces = positive + [[-hi, -lo] for lo, hi in positive]
        barriers = roots + [-x for x in roots] if strict else []
        return merge_intervals(pieces, split_at=barriers)

    def containment(self, level, strict=False):
        """P(Q >= level), or its right limit P(Q > level).

        C_j=[a,b] shifted by -y fits inside [l,r] iff y is in [b-r,a-l].
        Intersect the resulting unions over j and integrate the marginal law.
        Openness affects only probability-zero endpoints, but a deleted point
        separating components must not be filled in the strict calculation.
        """
        L = self.level_set(level, strict=strict)
        feasible = [[-np.inf, np.inf]]
        for a, b in self.C:
            shifts = merge_intervals([[b - r, a - l] for l, r in L if r - l >= b - a])
            feasible = intersect_intervals(feasible, shifts)
        return mixture_mass(feasible, self.mu), L, feasible

    def threshold(self):
        target = 1 - self.alpha_inc
        for level in sorted(set(float(self.m(x)) for x in self.critical)):
            closed = self.containment(level)[0]
            above = self.containment(level, strict=True)[0]
            if closed >= target and above <= target:
                return level, above < closed - 1e-10
        cutoff = brentq(lambda t: self.containment(t)[0] - target,
                        0.0, float(self.m(0)), xtol=2e-15, rtol=1e-13)
        return cutoff, False

    def whole_score(self, y):
        y = np.asarray(y)
        score = np.ones_like(y, dtype=float)
        stationary = self.critical + [-x for x in self.critical] + [0.0]
        for a, b in self.C:
            score = np.minimum(score, np.minimum(self.m(a - y), self.m(b - y)))
            for x in stationary:
                inside = (a - y <= x) & (x <= b - y)
                score = np.minimum(score, np.where(inside, self.m(x), np.inf))
        return score

    def result(self):
        tau, atom = self.threshold()
        probability, L, feasible = self.containment(tau)
        return dict(mu=self.mu, alpha_cf=self.alpha_cf, alpha_inc=self.alpha_inc,
                    hpd_level=self.hpd_level, tau_star=tau, threshold_has_atom=atom,
                    inclusion_probability=probability,
                    strict_inclusion_probability=self.containment(tau, strict=True)[0],
                    ours_exante_length=interval_length(L),
                    lei_exante_length=lei_exante_length(self.mu, self.alpha_cf, self.alpha_inc),
                    n_components=len(L), hpd_stage1_length=interval_length(self.C),
                    lei_stage1_length=2 * lei_first_stage_halfwidth(self.mu, self.alpha_cf),
                    C_intervals=self.C, L_tau_intervals=L, inclusion_Y_intervals=feasible,
                    valleys=[dict(location=x, level=float(self.m(x)),
                                  closed_probability=self.containment(float(self.m(x)))[0],
                                  strict_probability=self.containment(float(self.m(x)), True)[0])
                             for x in self.valleys])


def compute_one_mu(mu, alpha_cf=ALPHA_CF, alpha_inc=ALPHA_INC, n_root_grid=N_ROOT_GRID):
    return OracleProblem(mu, alpha_cf, alpha_inc, n_root_grid).result()


def finite_transitions(df, alpha_cf=ALPHA_CF, alpha_inc=ALPHA_INC, n_root_grid=N_ROOT_GRID):
    """Refine the first transitions bracketed by the sweep; no uniqueness claim."""
    mus = df.mu.to_numpy()
    difference = (df.ours_exante_length - df.lei_exante_length).to_numpy()
    components = df.n_components.to_numpy()
    # Ignore exact agreement before the HPD set separates. The manuscript's
    # length boundary is the crossing after the final set has disconnected.
    indices = np.flatnonzero((difference[:-1] > 1e-9) & (difference[1:] < -1e-9)
                            & (components[:-1] == 3))
    results = {}
    if len(indices):
        i = indices[0]
        def length_gap(mu):
            row = compute_one_mu(mu, alpha_cf, alpha_inc, n_root_grid)
            return row['ours_exante_length'] - row['lei_exante_length']
        results['mu_len'] = brentq(length_gap, mus[i], mus[i + 1], xtol=1e-10)
    indices = np.flatnonzero((components[:-1] == 1) & (components[1:] == 3))
    if len(indices):
        i = indices[0]
        def topology_gap(mu):
            problem = OracleProblem(mu, alpha_cf, alpha_inc, n_root_grid)
            if len(problem.valleys) != 1:
                raise RuntimeError('Topology refinement requires one positive local minimum.')
            level = float(problem.m(problem.valleys[0]))
            return problem.containment(level, strict=True)[0] - (1 - alpha_inc)
        results['mu_disc'] = brentq(topology_gap, mus[i], mus[i + 1], xtol=1e-10)
        problem = OracleProblem(results['mu_disc'], alpha_cf, alpha_inc, n_root_grid)
        level = float(problem.m(problem.valleys[0]))
        results['topology_boundary'] = dict(level=level,
            closed_probability=problem.containment(level)[0],
            strict_probability=problem.containment(level, strict=True)[0])
    return results


def run_sweep(alpha_cf=ALPHA_CF, alpha_inc=ALPHA_INC, mu_step=0.01):
    rows = []
    mus = np.linspace(2, 8, round(6 / mu_step) + 1)
    for index, mu in enumerate(mus, start=1):
        rows.append(compute_one_mu(float(mu), alpha_cf, alpha_inc))
        if index % 100 == 0 or index == len(mus):
            print(f'Computed {index}/{len(mus)} separations', flush=True)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--alpha-cf', type=float, default=ALPHA_CF)
    parser.add_argument('--alpha-inc', type=float, default=ALPHA_INC)
    parser.add_argument('--mu-step', type=float, default=0.01)
    parser.add_argument('--output-dir', type=Path, default=OUTPUT_DIR)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows = run_sweep(args.alpha_cf, args.alpha_inc, args.mu_step)
    df = pd.DataFrame([{k: v for k, v in row.items() if not isinstance(v, list)} for row in rows])
    transitions = finite_transitions(df, args.alpha_cf, args.alpha_inc)
    df.to_csv(args.output_dir / 'oracle_sweep_results.csv', index=False)
    (args.output_dir / 'oracle_sweep_details.json').write_text(json.dumps(rows, indent=2), encoding='utf8')
    (args.output_dir / 'oracle_transitions.json').write_text(json.dumps(transitions, indent=2), encoding='utf8')
    from plot_oracle_sweep import make_figure
    make_figure(df, args.output_dir)
    print(json.dumps(transitions, indent=2), flush=True)


if __name__ == '__main__':
    main()
