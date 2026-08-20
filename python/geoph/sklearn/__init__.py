from .. import *

from joblib import Parallel, delayed
from sklearn.base import BaseEstimator, TransformerMixin

class DelaunayRipsPersistence(BaseEstimator, TransformerMixin):
    """
    Computes Delaunay-Rips persistence diagrams from a collection of point clouds.

    Parameters
    ----------
    homology_dimensions : list
        The returned persistence diagrams dimensions.

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

    def __init__(self, homology_dimensions, n_jobs=-1):
        self.homology_dimensions = homology_dimensions
        self.n_jobs = n_jobs

    def fit(self, X, y=None):
        return self

    def transform(self, X):
        def func(x):
            diags = delaunayRipsPersistenceDiagram(x, "contiguous")
            return [diags[d] for d in self.homology_dimensions]

        return Parallel(n_jobs=self.n_jobs)(delayed(func)(x) for x in X)