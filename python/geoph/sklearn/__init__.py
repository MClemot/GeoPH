from .. import delaunayRipsPersistenceDiagram

from joblib import Parallel, delayed
from sklearn.base import BaseEstimator, TransformerMixin

__all__ = ["DelaunayRipsPersistence"]

class DelaunayRipsPersistence(BaseEstimator, TransformerMixin):
    """
    Computes Delaunay-Rips persistence diagrams from a collection of point clouds.

    Parameters
    ----------
    homology_dimensions : list
        The returned persistence diagrams dimensions.

    reduced_homology : bool, default=True
        If True, discards the 0-dimensional essential pair (0.0, inf).

    n_jobs : int, default=-1
        Number of parallel workers.

    Input
    -----
    X : iterable of ndarray
        Each element is a point cloud of shape (n_points, n_dimensions).

    Returns
    -------
    diagrams : list
        One persistence diagram per input point cloud.
    """

    def __init__(self, homology_dimensions, reduced_homology=True, n_jobs=-1):
        self.homology_dimensions = homology_dimensions
        self.reduced_homology = reduced_homology
        self.n_jobs = n_jobs

    def fit(self, X, y=None):
        return self

    def __transform(self, x):
        diags = delaunayRipsPersistenceDiagram(x, "contiguous", self.reduced_homology)
        return [diags[d] if d < len(diags) else np.empty((0, 2)) for d in self.homology_dimensions]

    def transform(self, X):
        return Parallel(n_jobs=self.n_jobs, prefer="threads")(delayed(self.__transform)(x) for x in X)