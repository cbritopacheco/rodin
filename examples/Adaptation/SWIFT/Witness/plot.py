#!/usr/bin/env python3
# Copyright (c) 2026 Rodin contributors.
# Distributed under the Boost Software License, Version 1.0.
# See accompanying file LICENSE_1_0.txt or https://www.boost.org/LICENSE_1_0.txt.
"""Plot a Witness JSON result as an SVG without a plotting dependency."""

import argparse
import html
import json
from pathlib import Path

import numpy as np

class Reference:
    """Reference drawing data supplied by the C++ generator, without inference."""

    def __init__(self, result):
        self.name = result['geometry']
        self.dimension = result['dimension']
        self.vertices = np.array(result['reference_vertices'])
        self.edges = result['reference_edges']


class Plot:
    """Initial and final sets share one scale; 3D views are projections only."""

    def __init__(self, result):
        self.result = result
        self.geometry = Reference(result)

    def project(self, points, view=None):
        """Orthonormal view axes avoid collapsing the cube diagonal."""
        if self.geometry.dimension < 2:
            return np.column_stack((points[:, 0], np.zeros(len(points))))
        if self.geometry.dimension == 2:
            return points
        if view is not None and view != 'XYZ':
            return points[:, {'XY': (0, 1), 'YZ': (1, 2), 'XZ': (0, 2)}[view]]
        return np.column_stack((.8*points[:, 0]-.6*points[:, 1],
                                np.sqrt(.75)*points[:, 2]-.3*points[:, 0]-.4*points[:, 1]))

    def edges(self):
        """Use reference edges written by C++; no hull reconstruction."""
        return self.geometry.edges

    def save(self, path):
        geometry = self.geometry
        dimension = geometry.dimension
        views = ('XYZ', 'XY', 'YZ', 'XZ') if dimension == 3 else (None,)
        height = 680+540*(len(views)-1)
        svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="{height}" viewBox="0 0 1200 {height}">',
               f'<rect width="1200" height="{height}" fill="white"/>',
               '<g font-family="sans-serif" fill="#25352f">',
               f'<text x="600" y="32" text-anchor="middle" font-size="23">{geometry.name.title()}: {self.result["count"]} witnesses</text>']

        projections = np.concatenate([self.project(geometry.vertices, view) for view in views])
        lower, upper = projections.min(axis=0), projections.max(axis=0)
        if dimension == 3:
            radius = max(self.result['initial_radius'], self.result['covering_radius'])
            lower -= radius
            upper += radius
        extent = max(float((upper-lower).max()), 1.)
        scale = 400/extent
        initialization = 'Lattice-derived initialization'
        panels = [(row, view, panel, name, key)
                  for row, view in enumerate(views)
                  for panel, (name, key) in enumerate(((initialization, 'initial_points'),
                                                       ('Numerical covering candidate', 'points')))]
        for row, view, panel, name, key in panels:
            offset = row*540
            projected = self.project(geometry.vertices, view)
            lower = projected.min(axis=0)
            if dimension == 3:
                lower -= radius
            points = np.array(self.result[key])
            initial = key == 'initial_points'
            evaluation = {
                'radius': self.result['initial_radius' if initial else 'covering_radius'],
                'holes': np.array(self.result['initial_worst_locations' if initial else 'worst_locations']),
                'cells': [np.array(cell) for cell in self.result['initial_cells' if initial else 'cells']]}

            def screen(points):
                origin = [300+panel*600, 350] if dimension == 0 else [80+panel*600, 350 if dimension == 1 else 550]
                origin[1] += offset
                return np.array(origin) + (self.project(points, view)-lower)*[scale, -scale]

            def polygon(points):
                return ' '.join(f'{x:.5f},{y:.5f}' for x, y in screen(points))

            label = name if view is None else f'{name}: {view}'
            svg.append(f'<text x="{70+panel*600}" y="{76+offset}" font-size="21">{label}</text>')
            svg.append(f'<text x="{70+panel*600}" y="{102+offset}" font-size="18">Covering radius: {evaluation["radius"]:.6f}</text>')
            if view is not None and view != 'XYZ':
                svg.append(f'<text x="{280+panel*600}" y="{580+offset}" font-size="17">{view[0]}</text>')
                svg.append(f'<text x="{45+panel*600}" y="{350+offset}" font-size="17">{view[1]}</text>')
            if dimension == 2:
                boundary = geometry.vertices
                outline = polygon(boundary)
                svg.append(f'<defs><clipPath id="domain{panel}"><polygon points="{outline}"/></clipPath></defs>')
                svg.append(f'<polygon points="{outline}" fill="#f2f7f5"/>')
                for cell in evaluation['cells']:
                    svg.append(f'<polygon points="{polygon(cell)}" fill="none" stroke="#91a9a0"/>')
                for x, y in screen(evaluation['holes']):
                    svg.append(f'<circle cx="{x}" cy="{y}" r="{evaluation["radius"]*scale}" fill="#d75448" fill-opacity=".13" stroke="#bf342b" stroke-dasharray="6 5" clip-path="url(#domain{panel})"/>')
                svg.append(f'<polygon points="{outline}" fill="none" stroke="#203e35" stroke-width="2"/>')
            elif dimension == 3:
                for x, y in screen(evaluation['holes']):
                    svg.append(f'<circle class="uncovered-ball" cx="{x}" cy="{y}" r="{evaluation["radius"]*scale}" fill="#d75448" fill-opacity=".06" stroke="#bf342b" stroke-dasharray="6 5"/>')
                for i, j in self.edges():
                    a, b = screen(geometry.vertices[[i, j]])
                    svg.append(f'<path d="M{a[0]},{a[1]} L{b[0]},{b[1]}" stroke="#91a9a0" fill="none"/>')
            elif dimension == 1:
                a, b = screen(geometry.vertices)
                svg.append(f'<path d="M{a[0]},{a[1]} L{b[0]},{b[1]}" stroke="#91a9a0" fill="none"/>')
            for x, y in screen(evaluation['holes']):
                svg.append(f'<path d="M{x-5},{y-5} l10,10 M{x-5},{y+5} l10,-10" stroke="#bf342b" stroke-width="2"/>')
            for i, (x, y) in enumerate(screen(points)):
                svg.append(f'<circle cx="{x}" cy="{y}" r="5" fill="#12664f" stroke="white"/>')
                if len(points) <= 30:
                    svg.append(f'<text x="{x+8}" y="{y-8}" font-size="12">{i+1}</text>')
        svg.append(f'<text x="600" y="{height-60}" text-anchor="middle" font-size="16">Green: witnesses. Red crosses: globally worst-covered locations for each set.</text>')
        note = ('Red shading: largest empty disks clipped to the domain. Gray: Voronoi cells.'
                if dimension == 2 else
                'XYZ is oblique. Circles project uncovered 3D balls; projected witnesses may overlap them.'
                if dimension == 3 else 'Distances and radii are computed in the reference geometry.')
        svg.append(f'<text x="600" y="{height-32}" text-anchor="middle" font-size="15">{html.escape(note)}</text>')
        svg.append('</g></svg>')
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text('\n'.join(svg)+'\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("result", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    Plot(json.loads(args.result.read_text())).save(args.output)
