"Tools for benchmarking functions."

import cProfile
import pstats
import io


def profile(func):
    "Decorator for profiling the performance of an intensive function."
    def wrapper(*args, **kwargs):
        pr = cProfile.Profile()
        try:
            pr.enable()
            # PROFILE CODE
            retval = func(*args, **kwargs)
            # PROFILE CODE
            pr.disable()
            return retval
        finally:
            s = io.StringIO()
            sortby = 'time'
            ps = pstats.Stats(pr, stream=s).strip_dirs().sort_stats(sortby)
            ps.print_stats(20)
            print(s.getvalue())
            # END PROFILE
    return wrapper
