from __future__ import absolute_import
from structrans.general_imports import *

logger = logging.getLogger(__name__)


DFS_CACHE = {}
global DFS_CACHE_COUNT
DFS_CACHE_COUNT = 0

logger.setLevel(logging.CRITICAL)


def dfs(start, cache=DFS_CACHE):

    flat = start.elem.flatten()
    key = tuple(flat)
    # find local minimum
    logger.debug("DFS search started from {:s}".format(str(flat)))

    if key in cache:
        cache['counter'] += 1
        logger.info("retrieve result from DFS_CACHE; counter = {:d}".format(cache['counter']))
        return cache[key]

    minnode = start

    cachenodes = [start]

    nextmin = minnode.steepest_des()
    minlociter = 0  # debug
    while nextmin is not None:
        minnode = nextmin
        cachenodes.append(minnode)
        nextmin = minnode.steepest_des()
        minlociter += 1 # debug
    logger.debug("Found local min at {:g} after {:d} iterations".format(minnode.elem_dist, minlociter))

    # applying cache
    for node in cachenodes:
        cache[tuple(node.elem.flatten())] = minnode
    logger.debug("DFS_CACHE has {:d} results now".format(len(cache)))

    return minnode