from __future__ import annotations
from .geoph_impl import *

__version__ = "0.1.2"

def delaunayRipsPersistenceDiagram(X, format:str='default', reduced_homology:bool=True) -> list[list[tuple[float, float]]] | list[tuple[int, tuple[float, float]]]:
    """
    Computes Delaunay-Rips persistence diagram.

    Args:
        X: an array with shape (*, *)
        format: 'default' or 'gudhi'
        reduced_homology: if True, discards the 0-dimensional essential pair (0.0, inf)

    Returns:
        A list of list of points (birth,death) when format is 'default'

        A list of tuples (dim,(birth,death)) when format is 'gudhi'
    """
    return geoph_impl.delaunayRipsPersistenceDiagram(X, format, reduced_homology)

def delaunayRipsPersistencePairs(X, reduced_homology:bool=True) -> list[list[tuple[list[int], list[int]]]]:
    """
    Computes Delaunay-Rips persistence pairs.

    Args:
        X: an array with shape (*, *)
        reduced_homology: if True, discards the 0-dimensional essential pair (0.0, inf)

    Returns:
        A list of list of simplex pairs (birthSimplex, deathSimplex)
    """
    return geoph_impl.delaunayRipsPersistencePairs(X, reduced_homology)

def delaunayRipsPersistenceGenerators2(X) -> list[tuple[list[tuple[int, int]], tuple[float, float]]]:
    """
    Computes Delaunay-Rips 1-dimensional persistent generators for a point cloud in R^2.
    Args:
        X: an array with shape (*, 2)
    Returns:
        A list of 1-dimensional persistent generators.
        Each generators consists in a list of edge, and a pair of values (birth, death)
    """
    return geoph_impl.delaunayRipsPersistenceGenerators2(X)

def delaunayRipsPersistenceGenerators3(X) -> tuple[list[tuple[list[tuple[int, int]], tuple[float, float]]], list[tuple[list[list[int]], tuple[float, float]]]]:
    """
    Computes Delaunay-Rips 1- and 2-dimensional persistent generators for a point cloud in R^3.
    Args:
        X: an array with shape (*, 3)
    Returns:
        A list of 1-dimensional persistent generators, and a list of 2-dimensional persistent generators.
    """
    return geoph_impl.delaunayRipsPersistenceGenerators3(X)

def ripsPersistenceGenerators2(X) -> list[tuple[list[tuple[int, int]], tuple[float, float]]]:
    """
    Computes Rips 1-dimensional persistent generators for a point cloud in R^2.
    Args:
        X: an array with shape (*, 2)
    Returns:
        A list of 1-dimensional persistent generators.
        Each generators consists in a list of edge, and a pair of values (birth, death)
    """
    return geoph_impl.ripsPersistenceGenerators2(X)

def ripsPersistenceDiagram2(X, format:str='default', reduced_homology:bool=True) -> list[list[tuple[float, float]]] | list[tuple[int, tuple[float, float]]]:
    """
    Computes Rips persistence diagram for a point cloud in R^2.
    Args:
        X: an array with shape (*, *)
        format: 'default' or 'gudhi'
        reduced_homology: if True, discards the 0-dimensional essential pair (0.0, inf)
    Returns:
        A list of list of points (birth,death) when format is 'default'

        A list of tuples (dim,(birth,death)) when format is 'gudhi'
    """
    return geoph_impl.ripsPersistenceDiagram2(X, format, reduced_homology)

def ripsPersistencePairs2(X, reduced_homology:bool=True) -> list[list[tuple[list[int], list[int]]]]:
    """
    Computes Rips persistence pairs for a point cloud in R^2.
    Args:
        X: an array with shape (*, *)
        reduced_homology: if True, discards the 0-dimensional essential pair (0.0, inf)
    Returns:
        A list of list of simplex pairs (birthSimplex, deathSimplex)
    """
    return geoph_impl.ripsPersistencePairs2(X, reduced_homology)