# Copyright (c) 2026 Rodin contributors.
# Distributed under the Boost Software License, Version 1.0.
"""Independent SciPy/Qhull oracle for tests, never used by the visualizers."""

import numpy as np
from scipy.spatial import ConvexHull, HalfspaceIntersection, cKDTree


class Reference:
    VERTICES = {
        'point': [[0]], 'segment': [[0], [1]],
        'triangle': [[0, 0], [1, 0], [0, 1]],
        'quadrilateral': [[0, 0], [1, 0], [1, 1], [0, 1]],
        'tetrahedron': [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]],
        'pyramid': [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1]],
        'hexahedron': [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                       [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]],
        'wedge': [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [0, 1, 1]]}

    def __init__(self, name):
        self.vertices = np.array(self.VERTICES[name], dtype=float)
        self.dimension = 0 if name == 'point' else self.vertices.shape[1]

    def radius(self, points):
        points = np.array(points)
        if not self.dimension:
            return 0.
        if self.dimension == 1:
            ordered = np.sort(points[:, 0])
            return max(ordered[0], 1-ordered[-1],
                       np.diff(ordered).max()/2 if len(points) > 1 else 0)
        domain = np.unique(ConvexHull(self.vertices).equations, axis=0)
        center = self.vertices.mean(axis=0)
        cells = []
        for i, site in enumerate(points):
            differences = np.delete(points, i, axis=0)-site
            lengths = np.linalg.norm(differences, axis=1)
            halfspaces = domain.copy()
            halfspaces[:, -1] += halfspaces[:, :-1] @ site
            halfspaces = np.vstack((halfspaces, np.column_stack((differences/lengths[:, None], -lengths/2))))
            direction = center-site
            slopes = halfspaces[:, :-1] @ direction
            positive = slopes > 0
            limit = min(1., np.min(-halfspaces[positive, -1]/slopes[positive]) if positive.any() else 1.)
            cells.append(HalfspaceIntersection(halfspaces, .5*limit*direction).intersections+site)
        return float(cKDTree(points).query(np.vstack(cells))[0].max())
