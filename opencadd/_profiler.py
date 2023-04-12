"""
Tools for code profiling.
"""

from typing import Callable, Generator, Sequence, Union, NoReturn, List
from collections.abc import Iterable
import timeit

import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt


__author__ = "Armin Ariamajd"


class FunctionProfiler:
    """
    Profile (and compare) one or several functions.
    """
    def __init__(
            self,
            funcs: Union[Callable, Sequence[Callable]],
            arg_gens: Union[
                Callable[[Sequence[int]], Iterable[Sequence]],
                Sequence[Callable[[Sequence[int]], Iterable[Sequence]]]
            ],
    ):
        """
        Parameters
        ----------
        funcs : callable or sequence of callable
            Function(s) to be profiled.
        arg_gens : callable or sequence of callable
            Generators for each function in `funcs`, which receive a sequence of numbers
            indicating different input sizes for the respective function's arguments, and yield
            input arguments for that function in form of a tuple.
            That is, each function `f` in `funcs` should be callable
            with the corresponding generator `g` in `arg_gens`, for example as
            `[f(*args) for args in g([1, 10, 100])]`.
        """
        self._funcs: List[Callable] = funcs if isinstance(funcs, Sequence) else [funcs]
        self._arg_gens: List[
            Callable[[Sequence[int]], Iterable[Sequence]]
        ] = arg_gens if isinstance(arg_gens, Sequence) else [arg_gens]
        if len(self._funcs) != len(self._arg_gens):
            raise ValueError(
                "Parameters `funcs` and `arg_gens` expect inputs with identical lengths, "
                f"but `funcs` had a length of {len(self._funcs)}, "
                f"while `arg_gens` had  {len(self._arg_gens)}."
            )
        self._arg_sizes: List[int] = None
        self._runs: int = None
        self._loops_per_run: int = None
        self._results: np.ndarray = None
        return

    def profile(self, arg_sizes: Sequence[int], runs: int = 100, loops_per_run: int = 1) -> NoReturn:
        """
        Profile all functions.

        Parameters
        ----------
        arg_sizes : sequence of int
            Different argument sizes to profile the functions with.
        runs : int, default: 100
            Number of times the profiling is repeated. Higher values provide more accurate results.
            The shortest duration between all runs will be selected,
            as it represents the most accurate duration.
        loops_per_run : int, default: 1
            Number of times each run is repeated. For each run, the average of all loops will be selected.

        Returns
        -------
        None
        """
        self._arg_sizes = arg_sizes
        self._runs = runs
        self._loops_per_run = loops_per_run
        self._results = np.array(
            [
                [
                    np.min(
                        timeit.repeat(
                            lambda: func(*args),
                            repeat=self._runs,
                            number=self._loops_per_run
                        )
                    ) / self._loops_per_run for args in arg_gen(arg_sizes)
                ] for func, arg_gen in zip(self._funcs, self._arg_gens)
            ]
        )
        return

    def plot(self, show: bool = True):
        if self._results is None:
            raise ValueError("No profiling has been performed yet; call `FunctionProfiler.profile` first.")
        fig, ax = plt.subplots()
        artists = []
        for func, result in zip(self._funcs, self._results):
            line, = ax.plot(
                self._arg_sizes, result,
                marker=".",
                label=f"\u200b{func.__module__}.{func.__qualname__}"
            )
            artists.append(line)
        ax.legend(handles=artists, loc="best")
        plt.xlabel("Input size")
        plt.ylabel("Time [s]")
        plt.xscale("log")
        plt.yscale("log")
        if show:
            plt.show()
            return
