import collections
import pdb
import sys

#Thanks to Jochen Ritzel
# https://stackoverflow.com/a/2912455
class keydefaultdict(collections.defaultdict):
    def __missing__(self, key):
        if self.default_factory is None:
            raise KeyError( key )
        else:
            ret = self[key] = self.default_factory(key)
            return ret

# Thanks to Romuald Brunet
# https://stackoverflow.com/a/23654936
class ForkedPdb(pdb.Pdb):
    """A Pdb subclass that may be used
    from a forked multiprocessing child
   
    Use 'ForkedPdb().set_trace()' as a replacement for 'breakpoint()'
    """
    def interaction(self, *args, **kwargs):
        _stdin = sys.stdin
        try:
            sys.stdin = open('/dev/stdin')
            pdb.Pdb.interaction(self, *args, **kwargs)
        finally:
            sys.stdin = _stdin
mpbreakpoint = ForkedPdb().set_trace

