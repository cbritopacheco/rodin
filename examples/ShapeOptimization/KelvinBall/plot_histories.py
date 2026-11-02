#!/usr/bin/env python3
"""Generate a parameter-first LaTeX/PGFPlots PDF from KelvinBall's CSV.

Requires Python's standard library and pdflatex with PGFPlots installed.
Reads a snapshot; ignores only an incomplete trailing CSV record.
Nonfinite values remain gaps, never zeros. All columns are accounted for.
"""
import argparse
import csv
import datetime
import io
import math
from pathlib import Path
import shutil
import subprocess
import tempfile
from zoneinfo import ZoneInfo


PARAMETERS = (
    'outer_radius initial_grid_spacing hmin_factor hmax_factor requested_hmin '
    'requested_hmax dt advection_quadrature_order regularization_length '
    'normal_smoothing_length nitsche_penalty stabilization_factor assembly_backend '
    'direct_solver level_set_penalty thickness_min thickness_guard_distance '
    'target_volume remesh_hausdorff remesh_hmin remesh_hmax').split()

LABELS = {
    'rho': r'$\rho$', 'k': '$k$', 'c': '$c$', 'q': '$q$',
    'h': 'Mean Edge $h$', 'background_h': 'Background Mean Edge',
    'dt': r'Advection Step $\Delta t$',
    'initial_grid_spacing': 'Reference Grid Spacing $h_0$',
    'requested_hmin': r'Requested $h_{\min}$',
    'requested_hmax': r'Requested $h_{\max}$',
    'thickness_min': 'Minimum Thickness $D_0$',
    'thickness_guard_distance': 'Guard Distance $D_g$',
    'thickness_penalty': r'Guard $P_{D_g}$',
    'thickness_nominal_penalty': r'Nominal $P_{D_0}$',
    'thickness_minimum_exit': r'$\min\{D_g,\min_x t(x)\}$',
    'actual_delta_rho': 'Actual', 'predicted_delta_rho': 'Predicted',
    'actual_delta_volume': 'Actual', 'predicted_delta_volume': 'Predicted',
    'actual_delta_thickness': 'Actual', 'predicted_delta_thickness': 'Predicted',
    'mesh_quality_min': 'Minimum', 'mesh_quality_mean': 'Mean',
    'mesh_quality_max': 'Maximum', 'cells': 'Current Mesh',
    'remesh_cells_after': 'After Reconstruction',
    'remesh_cells_before': 'Before Reconstruction',
    'obstacle_cells': 'Solid', 'fluid_cells': 'Fluid',
    'rho_gradient_residual': r'$\rho$ Gradient',
    'volume_gradient_residual': 'Volume Gradient',
    'rho_gradient_jump': r'$\rho$ Gradient',
    'volume_gradient_jump': 'Volume Gradient',
    'eikonal_rotated_jump': 'Eikonal', 'projected_rotated_jump': 'Projected',
    'advected_rotated_jump': 'Advected',
}

PARAMETER_DEFINITIONS = {
    'outer_radius': (r'$R$', 'Outer chamber radius'),
    'initial_grid_spacing': (r'$h_0$', r'$R/(n-1)$; reference grid spacing'),
    'hmin_factor': (r'$a_{\min}$', r'$h_{\min}/h_0$'),
    'hmax_factor': (r'$a_{\max}$', r'$h_{\max}/h_0$'),
    'requested_hmin': (r'$h_{\min}$', r'$a_{\min}h_0$; minimum target size'),
    'requested_hmax': (r'$h_{\max}$', r'$a_{\max}h_0$; maximum target size'),
    'dt': (r'$\Delta t$', r'$s h_0$; advection step'),
    'advection_quadrature_order': (r'$p_{\mathrm{quad}}$', 'Advection quadrature order'),
    'regularization_length': (r'$\ell_\theta$', r'$a_\theta h_0$; Hilbert extension length'),
    'normal_smoothing_length': (r'$\ell_n$', r'$a_n h_0$; normal smoothing length'),
    'nitsche_penalty': (r'$\gamma_u$', 'Stokes Nitsche penalty coefficient'),
    'stabilization_factor': (r'$\tau_0$', 'Pressure stabilisation factor'),
    'assembly_backend': ('---', 'Assembly backend'),
    'direct_solver': ('---', 'Linear solver'),
    'level_set_penalty': (r'$\gamma_\phi$', 'Scalar trace penalty coefficient'),
    'thickness_min': (r'$D_0$', 'Requested minimum thickness'),
    'thickness_guard_distance': (r'$D_g$', r'$D_0+2\Delta t$; guard distance'),
    'target_volume': (r'$V_0$', r'$|\Omega_{\mathrm{Ch}}^0|$; initial solid volume'),
    'remesh_hausdorff': (r'$\varepsilon_H$', r'$0.1h_0^2$; MMG geometric tolerance'),
    'remesh_hmin': (r'$h_{\min}^{\mathrm{MMG}}$', 'Minimum size passed to MMG'),
    'remesh_hmax': (r'$h_{\max}^{\mathrm{MMG}}$', 'Maximum size passed to MMG'),
}

LABELS.update({
    'd_rho_theta': r'$D\rho[\theta]$',
    'd_volume_theta': r'$DV[\theta]$',
    'required_d_volume_theta': r'Required $DV[\theta]$',
    'd_thickness_theta': r'$\widetilde{D}P[\theta]$',
    'xi_rho_inf_norm': r'$\|\xi_\rho\|_\infty$',
    'theta_inf_norm': r'$\|\theta\|_\infty$',
    'thickness_violating_rays': 'Guard-Active Rays',
    'thickness_active_feasible': 'Correction Target Feasible',
    **{f'z_{a}': f'$Z_{a}$' for a in 'xyz'},
    **{f'omega_{a}': rf'$\omega_{a}$' for a in 'xyz'},
})

# Importance order. The remaining changing columns are included afterward.
GROUPS = [
    ('Objective', r'$\rho=|c|/\sqrt{kq}$', ['rho']),
    ('Volume Constraint', r'$(V-V_0)/V_0\;(\%)$', ['volume_error']),
    ('First-Exit Thickness Penalties', '$P_D$', ['thickness_penalty', 'thickness_nominal_penalty']),
    ('Resistance Coupling', '$c$', ['c']),
    ('Mesh Quality', 'Quality', ['mesh_quality_min', 'mesh_quality_mean', 'mesh_quality_max']),
    ('Mesh Growth', 'Tetrahedra', ['cells', 'remesh_cells_after']),
    ('Translational Resistance', '$k$', ['k']),
    ('Rotational Resistance', '$q$', ['q']),
    ('First-Exit Distance', 'Length', ['thickness_minimum_exit', 'thickness_min']),
    ('Thickness Activity', 'Count', ['thickness_violating_rays']),
    ('Objective Increment', r'$\Delta\rho$', ['actual_delta_rho', 'predicted_delta_rho']),
    ('Volume Increment', r'$\Delta V$', ['actual_delta_volume', 'predicted_delta_volume']),
    ('Thickness Increment', r'$\Delta P$', ['actual_delta_thickness', 'predicted_delta_thickness']),
    ('Directional Objective Derivative', r'$D\rho[\theta]$', ['d_rho_theta']),
    ('Directional Volume Derivative', r'$DV[\theta]$', ['d_volume_theta', 'required_d_volume_theta']),
    ('Thickness Surrogate Differential', r'$\widetilde{D}P[\theta]$', ['d_thickness_theta']),
    ('Null-Space Multiplier', 'Multiplier', ['null_space_multiplier']),
    ('Thickness Correction', 'Multiplier', ['thickness_active_multiplier']),
    ('Direction Magnitude', 'Infinity Norm', ['xi_rho_inf_norm', 'theta_inf_norm']),
    ('Direction Tangential Fractions', 'Fraction', ['thickness_descent_tangent_fraction', 'theta_tangent_fraction']),
    ('Thickness Direction Components', 'RMS', ['thickness_descent_normal_rms', 'thickness_descent_tangent_rms']),
    ('Pair-Load Magnitude', 'Magnitude', ['thickness_maximum_pair_load']),
    ('Pair-Load Concentration', 'Top One Percent Share', ['thickness_top_one_percent_load_share']),
    ('Thickness Deficit', 'Length', ['thickness_maximum_deficit']),
    ('First-Exit Transversality', r'$m\cdot d$', ['thickness_minimum_transversality']),
    ('Normal Alignment', 'Dot Product', ['thickness_normal_alignment_min', 'thickness_normal_alignment_mean']),
    ('Smoothed Curvature', 'Curvature', ['thickness_curvature_min', 'thickness_curvature_max']),
    ('Projected Normal Magnitude', 'Magnitude', ['thickness_normal_magnitude_min', 'thickness_normal_magnitude_max']),
    ('Thickness Seam Jumps', 'Jump', ['thickness_geometric_seam_jump', 'thickness_normal_seam_jump', 'thickness_ray_seam_jump']),
    ('Gradient Solve Residuals', 'Relative Residual', ['rho_gradient_residual', 'volume_gradient_residual']),
    ('Gradient Trace Jumps', 'Jump', ['rho_gradient_jump', 'volume_gradient_jump']),
    ('Stokes Trace Jump', 'Jump', ['nitsche_jump']),
    ('Resistance Reciprocity', 'Residual', ['coupling_symmetry']),
    ('Scalar Trace Jumps', 'Jump', ['eikonal_rotated_jump', 'projected_rotated_jump', 'advected_rotated_jump']),
    ('Distance Projection', 'Correction', ['distance_correction', 'interface_shift_max']),
    ('Advection Increment', 'Infinity Norm', ['advection_increment']),
    ('Mesh Scales', 'Length', ['h', 'background_h']),
    ('Material Partition', 'Tetrahedra', ['obstacle_cells', 'fluid_cells']),
    ('Interface Resolution', 'Triangles', ['interface_triangles']),
    ('Mesh Vertices', 'Vertices', ['vertices']),
    ('Free Translation', 'Velocity', ['z_x', 'z_y', 'z_z']),
    ('Free Rotation', 'Angular Velocity', ['omega_x', 'omega_y', 'omega_z']),
    ('Mobility Coefficients', 'Mobility', ['translational_mobility', 'coupling_mobility', 'rotational_mobility']),
    ('Revolution Period', 'Period', ['revolution_period']),
    ('Helical Pitch', 'Pitch', ['pitch']),
    ('Resistance Length', r'$\sqrt{q/k}$', ['resistance_radius']),
    ('Principal Stage Timings', 'Time (s)', ['stage_3_seconds', 'stage_6_seconds', 'stage_9_seconds']),
    ('Transport And Output Timings', 'Time (s)', ['stage_7_seconds', 'stage_8_seconds']),
    ('Remaining Stage Timings', 'Time (s)', ['stage_2_seconds', 'stage_4_seconds', 'stage_5_seconds']),
]

# Direction of preference concerns the named diagnostic, not every component.
PREFERENCES = {
    'Objective': r'Higher $\rho$ is better.',
    'Volume Constraint': 'Closer to zero is better; the sign distinguishes excess and deficit.',
    'First-Exit Thickness Penalties': 'Lower penalties are better; zero means no sampled active deficit.',
    'Resistance Coupling': r'Larger $|c|$ increases coupling at fixed $kq$; the sign records handedness.',
    'Mesh Quality': 'Higher quality is better, particularly the minimum.',
    'Mesh Growth': 'Fewer cells reduce cost at comparable accuracy; more cells are not intrinsically better.',
    'Translational Resistance': r'Lower $k$ increases $\rho$ only at fixed $c$ and $q$.',
    'Rotational Resistance': r'Lower $q$ increases $\rho$ only at fixed $c$ and $k$.',
    'First-Exit Distance': r'Higher exit distances are preferable up to $D_0$; the reported minimum is capped at $D_g$.',
    'Thickness Activity': 'Fewer guard-active rays are preferable at comparable sampling density.',
    'Objective Increment': 'Positive actual increments are preferable; agreement with the prediction tests the update.',
    'Volume Increment': 'Actual and predicted increments closer to zero are preferable for volume preservation.',
    'Thickness Increment': 'Negative actual increments reduce the penalty; agreement with the surrogate prediction is preferable.',
    'Directional Objective Derivative': 'Positive values predict objective increase; magnitude depends on direction normalisation.',
    'Directional Volume Derivative': 'Closer agreement with the required derivative is better, not a larger signed value.',
    'Thickness Surrogate Differential': 'Negative values predict reduction of the surrogate thickness penalty.',
    'Null-Space Multiplier': 'Neither higher nor lower is intrinsically better; this is a constraint correction coefficient.',
    'Thickness Correction': 'Neither higher nor lower is intrinsically better; this is the selected repair multiplier.',
    'Direction Magnitude': r'$\|\theta\|_\infty=1$ is the normalisation target; gradient magnitude has no monotonic preference.',
    'Direction Tangential Fractions': 'Lower is preferable for predominantly normal interface motion, not a convergence certificate.',
    'Thickness Direction Components': 'Normal and tangential RMS values describe the direction; larger magnitudes are not intrinsically better.',
    'Pair-Load Magnitude': 'Lower maxima indicate less concentrated forcing; they do not alone establish thickness feasibility.',
    'Pair-Load Concentration': 'Lower shares indicate less concentrated forcing.',
    'Thickness Deficit': 'Lower deficits are better; zero satisfies the sampled guard threshold.',
    'First-Exit Transversality': 'Higher values indicate less grazing; this is a ray-conditioning diagnostic.',
    'Normal Alignment': 'Values closer to one indicate better alignment.',
    'Smoothed Curvature': 'Neither sign nor magnitude is intrinsically better; curvature describes local bending.',
    'Projected Normal Magnitude': 'Values closer to one indicate better unit-normal recovery.',
    'Thickness Seam Jumps': 'Lower jumps are better for rotational compatibility.',
    'Gradient Solve Residuals': 'Lower relative residuals are better.',
    'Gradient Trace Jumps': 'Lower jumps are better.',
    'Stokes Trace Jump': 'Lower jumps are better for weak periodic compatibility.',
    'Resistance Reciprocity': 'Closer to zero is better.',
    'Scalar Trace Jumps': 'Lower jumps are better.',
    'Distance Projection': 'Smaller corrections indicate less trace repair; zero interface shift is preferable.',
    'Advection Increment': 'Neither higher nor lower is intrinsically better; the increment depends on the step and velocity.',
    'Mesh Scales': 'Smaller edges improve resolution but increase cost; comparisons require comparable accuracy.',
    'Material Partition': 'Neither count is intrinsically better; these record the solid and fluid partitions.',
    'Interface Resolution': 'More triangles can resolve finer geometry but increase cost; count alone does not measure accuracy.',
    'Mesh Vertices': 'Fewer vertices reduce cost at comparable accuracy.',
    'Free Translation': 'Neither signed component is intrinsically better for the coupling objective.',
    'Free Rotation': 'Neither signed component is intrinsically better; the sign records rotation sense.',
    'Mobility Coefficients': 'Neither signed value is intrinsically better for the normalised coupling objective.',
    'Revolution Period': 'Shorter periods mean faster rotation under the prescribed force, not necessarily higher coupling.',
    'Helical Pitch': 'Neither larger nor smaller pitch is intrinsically better for the coupling objective.',
    'Resistance Length': 'Neither higher nor lower is intrinsically better; this is a resistance length scale.',
    'Principal Stage Timings': 'Lower times are better at comparable accuracy and mesh size.',
    'Transport And Output Timings': 'Lower times are better at comparable accuracy and mesh size.',
    'Remaining Stage Timings': 'Lower times are better at comparable accuracy and mesh size.',
}

SUPPLEMENT_PREFERENCES = {
    'volume': 'Closer to the target volume is better.',
    'remesh_cells_before': 'Fewer cells reduce cost at comparable accuracy.',
    'thickness_active_feasible': 'One indicates a feasible correction target in the tested span; zero indicates infeasibility.',
    'thickness_seam_samples': 'Neither higher nor lower is intrinsically better; this records diagnostic sampling density.',
    'thickness_seam_unmatched': 'Fewer unmatched samples are better; zero is the target.',
}


def escape(text):
    return ''.join({'_': r'\_', '%': r'\%', '&': r'\&', '#': r'\#',
                    '$': r'\$', '{': r'\{', '}': r'\}'}.get(c, c) for c in str(text))


class HistoryReport:
    def __init__(self, csv_path):
        self.path = csv_path
        text = csv_path.read_text()
        reader = csv.DictReader(io.StringIO(text))
        self.keys = reader.fieldnames
        records = list(reader)
        self.rows = []
        for index, row in enumerate(records):
            if None in row or any(v is None or v == '' for v in row.values()):
                if index == len(records) - 1 and not text.endswith('\n'):
                    print('Ignored An Incomplete Trailing Record.')
                    continue
                raise ValueError(f'Malformed CSV Record {index + 2}')
            self.rows.append(row)
        if not self.rows:
            raise ValueError('No Complete Design Records')
        self.x = self.numeric('iteration')
        if any(b <= a for a, b in zip(self.x, self.x[1:])):
            raise ValueError('Design Iterations Must Be Strictly Increasing')
        self.used = {'iteration'}

    def numeric(self, key):
        return [float(row[key]) for row in self.rows]

    def changing(self, key):
        try:
            data = self.numeric(key)
        except ValueError:
            return False
        finite = [v for v in data if math.isfinite(v)]
        return len(set(finite)) > 1

    def label(self, key):
        if key.startswith('stage_'):
            return 'Stage ' + key.split('_')[1]
        return LABELS.get(key, escape(key.replace('_', ' ').title().replace('Rms', 'RMS').replace('Bbox', 'Bounding Box')))

    def table(self, keys, title):
        result = ''
        for start in range(0, len(keys), 26):
            result += r'\newpage' if result else ''
            result += ('\\subsection*{' if start else '\\section{') + title.capitalize() + (' (continued)' if start else '') + '}\n'
            parameters = title == 'Parameters'
            result += (r'\small\begin{tabular}{@{}p{2.1cm}p{9.1cm}p{4.1cm}@{}}\toprule Symbol & Definition & Value \\ \midrule'
                       if parameters else r'\small\begin{tabular}{@{}p{10cm}p{5.3cm}@{}}\toprule Quantity & Recorded Value \\ \midrule') + '\n'
            for key in keys[start:start+26]:
                self.used.add(key)
                raw = [r[key] for r in self.rows]
                try:
                    data = self.numeric(key)
                    finite = [v for v in data if math.isfinite(v)]
                    value = (f'{min(finite):.6g}' if len(set(finite)) == 1
                             else f'{min(finite):.6g} To {max(finite):.6g}') if finite else 'Unavailable'
                    missing = sum(math.isnan(v) for v in data)
                    infinite = sum(math.isinf(v) for v in data)
                    if missing: value += f'; {missing} NaN Records'
                    if infinite: value += f'; {infinite} Infinite Records'
                except ValueError:
                    value = ', '.join(dict.fromkeys(raw))
                if parameters:
                    symbol, definition = PARAMETER_DEFINITIONS[key]
                    label = symbol + ' & ' + definition
                elif key.startswith(('actual_delta_', 'predicted_delta_')):
                    kind, _, quantity = key.split('_', 2)
                    symbol = {'rho': r'\rho', 'volume': 'V', 'thickness': 'P'}[quantity]
                    label = kind.title() + rf' Increment $\Delta {symbol}$'
                else:
                    label = self.label(key)
                result += label + ' & ' + escape(value) + r' \\[4pt]' + '\n'
            result += r'\bottomrule\end{tabular}\par\normalsize' + '\n'
        return result

    def plot(self, title, ylabel, keys, options=''):
        keys = [k for k in keys if k in self.keys]
        self.used.update(keys)
        result = f'\\nextgroupplot[title={{{title}}},ylabel={{{ylabel}}}{options}]\n'
        for key in keys:
            data = self.numeric(key)
            if key == 'volume_error':
                data = [100*v/float(r['target_volume']) for v, r in zip(data, self.rows)]
            points = ' '.join(f'({x:g},{v:.12g})' if math.isfinite(v) else f'({x:g},nan)' for x,v in zip(self.x, data))
            result += '\\addplot coordinates {' + points + '};\n'
            result += '\\addlegendentry{' + self.label(key) + '}\n'
        return result

    def generate(self):
        preamble = r'''\documentclass[10pt,reqno]{amsart}
\usepackage[letterpaper,margin=0.9in]{geometry}
\usepackage[T1]{fontenc}
\usepackage{lmodern,amsmath,pgfplots,booktabs,microtype,float}
\usepackage[hidelinks]{hyperref}
\makeatletter
\renewcommand{\section}{\@startsection{section}{1}{\z@}
{1.5\baselineskip}{0.7\baselineskip}{\normalfont\normalsize\bfseries}}
\makeatother
\usepgfplotslibrary{groupplots}
\pgfplotsset{compat=1.18}
\definecolor{CurveBlue}{HTML}{245C85}
\definecolor{CurveOrange}{HTML}{B65A32}
\definecolor{CurveGreen}{HTML}{437A5A}
\pagestyle{plain}\setlength{\parindent}{0pt}
\pgfplotsset{history/.style={width=7.15cm,height=5.6cm,
axis lines=left,axis line style={black!65,thin},tick align=outside,
grid=major,major grid style={black!10,thin},xlabel={Design Iteration},
tick label style={font=\scriptsize},title style={font=\small,yshift=6pt},
legend style={font=\scriptsize,draw=none,fill=none,at={(0.5,-0.25)},anchor=north,
legend columns=1},scaled ticks=true,unbounded coords=jump,enlarge x limits=0.025,
every axis plot/.append style={line width=0.9pt,mark=*,mark size=1.4pt},
cycle list={{CurveBlue,mark=*},{CurveOrange,mark=square*},{CurveGreen,mark=triangle*}}}}
\begin{document}
'''
        timestamp = datetime.datetime.now(ZoneInfo('Europe/Paris')).strftime('%d %B %Y, %H:%M %Z')
        result = preamble + '\\begin{center}{\\large\\bfseries Shape optimisation of the Kelvin ball}\\par\\medskip Numerical histories\\par\\medskip ' + timestamp + '\\end{center}\n'
        result += '\\bigskip\n'
        result += self.table([k for k in PARAMETERS if k in self.keys], 'Parameters')
        result += r'''\clearpage
\tableofcontents
\bigskip
\textbf{Thickness Definition}\par
\[
P_D=24\sum_{T\subset\Gamma_{\mathrm{Ch}}}\frac{|T|}{3}
\sum_{j=1}^3(D-t(x_{T,j}))_+^2,\quad (s)_+=\max(s,0),\quad D_g=D_0+2\Delta t.
\]
Here $t(x)$ is the first outward exit distance along the inward geometric normal
through the full rotated solid, and $x_{T,j}$ are the three quadrature points
on each chamber interface triangle $T$. The factor $24$ reconstructs the full
solid penalty. The guard $P_{D_g}$ anticipates approach, whereas $P_{D_0}$ measures
sampled violations of $D_0$. The update uses a paired-normal surrogate differential,
not the exact shape derivative of this measured diagnostic. Missing values are
gaps, not zeros. The symbols $s$, $a_\theta$ and $a_n$ denote the dimensionless
step, regularisation and normal-smoothing factors. Only parameters recorded in
the CSV are tabulated.
'''
        result += '\\clearpage\n'
        panels = []
        for title, ylabel, keys in GROUPS:
            available = [k for k in keys if k in self.keys]
            if any(self.changing(k) for k in available):
                panels.append((title, ylabel, available))
                self.used.update(available)
        for key in self.keys:
            if key not in self.used and self.changing(key):
                panels.append((self.label(key), 'Recorded Value', [key]))
                self.used.add(key)
        span = self.x[-1]-self.x[0]
        step = max(1, math.ceil(span/8))
        ticks = ','.join(str(i) for i in range(int(self.x[0]), int(self.x[-1])+1, step))
        sections = [
            ('Objective and constraints', 0, 13),
            ('Shape direction and thickness correction', 13, 24),
            ('Normal fields and numerical consistency', 24, 36),
            ('Mesh resolution and free motion', 36, 46),
            ('Computational cost and supplementary diagnostics', 46, len(GROUPS)),
        ]
        for section, lo, hi in sections:
            titles = {g[0] for g in GROUPS[lo:hi]}
            selected = [p for p in panels if p[0] in titles or
                        (hi == len(GROUPS) and p[0] not in {g[0] for g in GROUPS})]
            if not selected:
                continue
            result += '\\section{' + section + '}\n'
            for start in range(0, len(selected), 4):
                batch = selected[start:start+4]
                columns = 1 if len(batch) == 1 else 2
                rows = math.ceil(len(batch) / columns)
                result += '\\vspace*{\\fill}\n\\begin{figure}[H]\\centering\n\\makebox[\\textwidth][c]{\\begin{tikzpicture}\n\\begin{groupplot}[history,group style={'
                result += f'group size={columns} by {rows},horizontal sep=1.5cm,vertical sep=3.9cm}},'
                result += f'xmin={self.x[0]:g},xmax={max(self.x[-1],self.x[0]+1):g},xtick={{{ticks}}}]\n'
                caption = []
                for i, (title, ylabel, keys) in enumerate(batch):
                    letter = chr(97+i)
                    options = ',xshift=4.325cm' if len(batch) == 3 and i == 2 else ''
                    result += self.plot(f'({letter}) {title}', ylabel, keys, options)
                    preference = PREFERENCES.get(title, SUPPLEMENT_PREFERENCES.get(keys[0],
                        'Neither higher nor lower is intrinsically better; this is a recorded diagnostic.'))
                    caption.append(f'({letter}) {title} ({ylabel}). {preference}')
                result += '\\end{groupplot}\n\\end{tikzpicture}}\n'
                result += ('\\caption{' + ' '.join(caption) +
                    ' Curves are identified by the panel legends; markers denote recorded designs, joined by straight segments. Axes are linear; powers of ten indicate tick scaling. Nonfinite values are omitted.}\\end{figure}\n\\vspace*{\\fill}\\clearpage\n')
        remaining = [k for k in self.keys if k not in self.used]
        if remaining:
            result += '\\newpage\n' + self.table(remaining, 'Constant Or Unavailable Diagnostics')
        nonfinite = []
        for key in self.keys:
            try:
                if any(not math.isfinite(v) for v in self.numeric(key)) and key not in remaining and key not in PARAMETERS:
                    nonfinite.append(key)
            except ValueError:
                pass
        if nonfinite:
            result += '\\newpage\n' + self.table(nonfinite, 'Nonfinite Data Summary')
        assert set(self.keys) <= self.used
        return result + '\\end{document}\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('csv', type=Path)
    parser.add_argument('--output', type=Path, default=Path('kelvin-scientific-histories.pdf'))
    parser.add_argument('--tex', type=Path, help='LaTeX Source; Defaults To The PDF Basename')
    parser.add_argument('--latex', default='pdflatex', help='LaTeX Executable')
    parser.add_argument('--source-only', action='store_true', help='Generate LaTeX Without Compiling')
    args = parser.parse_args()
    report = HistoryReport(args.csv.resolve())
    output = args.output.resolve()
    tex = args.tex.resolve() if args.tex else output.with_suffix('.tex')
    tex.parent.mkdir(parents=True, exist_ok=True)
    tex.write_text(report.generate())
    if not args.source_only:
        output.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix='kelvin-history-latex-') as scratch:
            for _ in range(3):
                run = subprocess.run([args.latex, '-interaction=nonstopmode', '-halt-on-error',
                                      '-output-directory='+scratch, str(tex)],
                                     capture_output=True, text=True)
                if run.returncode:
                    raise RuntimeError(run.stdout[-5000:] + run.stderr[-1000:])
                if 'Overfull' in run.stdout:
                    raise RuntimeError('LaTeX Reported An Overfull Box; Adjust The Layout Before Delivery.')
            shutil.copyfile(Path(scratch)/tex.with_suffix('.pdf').name, output)
        print('PDF:', output)
    print('LaTeX:', tex)
    print(f'Designs {report.x[0]:g}--{report.x[-1]:g}; {len(report.keys)} Columns Accounted For.')


if __name__ == '__main__':
    main()
