import functools
import time

## decorator for running time
def timerun(func):
    """ Calculate the execution time of a method and return it back"""
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        duration = time.time() - start

        print(f"(debugging : Duration of {func.__name__} function was {duration}.)")
        return result
    return wrapper
