import collections
import pdb
import sys


def nacast(var, nullval='.'):
    """
    Cast var to str, but use a specific character for instances of None.
    For the purpose of consistently logging null values.

    :param var: The variable to consider
    :param nullval: What to display for a None value of var
    """
    return str(var) if var is not None else nullval

#Thanks to Jochen Ritzel
# https://stackoverflow.com/a/2912455
class keydefaultdict(collections.defaultdict):
    def __missing__(self, key):
        if self.default_factory is None:
            raise KeyError( key )
        else:
            ret = self[key] = self.default_factory(key, self)
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

