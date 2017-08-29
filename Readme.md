# Calculate drift velocity

Calculates the drift velocity as in the `diskevolution` code.

# Installation

Either use `pip`: go to the directory where `setup.py` is located and run

    pip install -e .

Or compile the fortran code directly with

```
cd calc_drift # where the *.f90 files are
make
```

This should compile the module. you need to have the path of the module in your `PYTHONPATH` to import it, or be in the same directory.

# Calling the routine

See the `jupyter` notebook [`notebook/call_routine.ipynb`](notebook/call_routine.ipynb) for how to call it.

# Uninstall

You can just use `make clobber` if you used `make` for compilation or `pip uninstall calc_drift` if you used `pip`. You can also clean up everything in the directory with `git clean -xdf`, but **beware**: this will remove everything that is in the repository but not under version control (i.e. files you created, run `git clean -xdn` to check what will be deleted).
