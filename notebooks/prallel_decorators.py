"""
A few useful function for parallel and asynchronous processing. 

BUGGY
"""
from threading import Thread
from functools import wraps


class ThreadWithReturnValue(Thread):
    def __init__(self, group=None, target=None, name=None,
                 args=(), kwargs={}, Verbose=None):
        Thread.__init__(self, group, target, name, args, kwargs)
        self._return = None

    def run(self):
        if self._target is not None:
            self._return = self._target(*self._args, **self._kwargs)

    def result(self, *args):
        Thread.join(self, *args)
        return self._return
    
    def __repr__(self):
        while self.is_alive():
            return 'process is running silently in thread'
        return 'process is done, get results with Thread[:]'
    def __getitem__(self, index):
        return self.result()

def asynk(func):
    """Function decorator, intended to make "func" run in a separate
    thread  (asynchronously).
    Example:
    .. code-block:: python
        @async
        def task1():
            do_something
        @async
        def task2():
            do_something_too
        t1 = task1()
        t2 = task2()
        t1.join()
        t2.join()
        """
    @wraps(func)
    def _inner(*args, **kwargs):
        func_th = ThreadWithReturnValue(target=func, args=args, kwargs=kwargs)
        func_th.start()
        return func_th
    return _inner


def synk(func):
    """Function decorator, intended to make "func" wait for threading async
    finished.
    Example:
    .. code-block:: python
        @async
        def task1():
            do_something
        @async
        def task2():
            do_something_too
        @sync
        def task3():
            do_something_when_task1_and_task2_finished()
        t1 = task1()
        t2 = task2()
        t3 = task3() # Only runs when task1 and task2 finished.
    """
    @wraps(func)
    def _inner(*args, **kwargs):
        import threading
        import time
        while threading.activeCount() > 1:
            time.sleep(1)
        return func(*args, **kwargs)
    return _inner


def parallel(func, n_jobs=8, verbose=True, threading=False):
    """
    Parallel implementation for any function
    
    It's quick, it's dirty, it might fail, but it's beautiful when it works.
    This wrapper uses joblib in the backend to run scripts in parallel. 
    """
    from functools import wraps
    
    @wraps(func)
    def run_parallel(*args, n_jobs=n_jobs, verbose=verbose, **kwargs):
        """Runs the function through joblib. limited funcionality
        """
        from joblib import Parallel, delayed, parallel_backend
        from collections.abc import Iterable

        def isiter(v):
            not_string = not isinstance(v, str)
            is_iter = isinstance(v, Iterable)
            return not_string and is_iter
        
        if not all([isiter(a) for a in args]):
            raise ValueError(
                'Note that this function has been parallelised. You can thus '
                'only pass arguements that are iterable to the function. '
                'These can be lists, tuples, etc. but not strings. '
                'All the iterable items must be the same length, if not, then '
                'an error will be raised. ')

        lengths_args = set([len(a) for a in args if isiter(a)])
        not_iters = all([not isiter(v) for v in kwargs.values()])
        
        assert len(lengths_args) <= 1, 'All parallel inputs must have the same length on the 1st dimension'
        assert not_iters, 'keyword arguements cannot be iterable'
        
        len_arg = list(lengths_args)[0]
        if len_arg < n_jobs:
            n_jobs = len_arg
            
        function = delayed(func)
        parallel = Parallel(
            verbose=verbose,
            prefer='threading' if threading else 'processes',
            n_jobs=n_jobs)

        delayed_calls = []
        for arg in zip(*args):
            delayed_calls += function(*arg, **kwargs),

        return parallel(delayed_calls)
    return run_parallel