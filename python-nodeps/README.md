# No-Dependency Package

This is a minimal Frechet distance package which has been created due to one scientist
having problems using boost python and numpy. It is not tested, yet it should do roughly as it works.
Its error reporting is agressive (c++ exceptions just killing the python3 process).


# Install:
```
martin@TULRBGD-PROF01:/mnt/c/Users/ge98siw/Documents/git/frechetrange/python-nodeps$ pip3 install -e .
Obtaining file:///mnt/c/Users/ge98siw/Documents/git/frechetrange/python-nodeps
Installing collected packages: frechet
  Found existing installation: frechet 0.1.0
    Uninstalling frechet-0.1.0:
      Successfully uninstalled frechet-0.1.0
  Running setup.py develop for frechet
Successfully installed frechet
martin@TULRBGD-PROF01:/mnt/c/Users/ge98siw/Documents/git/frechetrange/python-nodeps$
```

# Usage
```
martin@TULRBGD-PROF01:/mnt/c/Users/ge98siw/Documents/git/frechetrange/python-nodeps$ python3
Python 3.7.3 (default, Dec 20 2019, 18:57:59)
[GCC 8.3.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> from frechet import frechet_distance
>>> a = frechet_distance([[1.1,1.2],[1.3,1.4]],[[2.2,2.3],[2.4,2.5]])
>>> print(a)
1.555633544921875
>>>
```

# TODO

Before using this, make sure that it works on a wide range of test cases. Don't blindly trust the numbers. Enable the debug output (commented) and check that it works.

The algorithm is not based on enumerating critical points (as we should), but starts with a hard-coded rnage (low,high), replaces this range with (high, 2*high) as long as a test on high is positive. Afterwards, it takes m = (low+high)/2 and halves intervals until some convergence. To me, this looks sound, yet it should (!) be test-driven. There are lots of things one can do wrong here.

Read and carefully check the source code! It is your duty as a scientist. We are providing this just as a jump-start to having a working environment without boost python and boost numpy.

Look at numpy-test.py how you can use numpy anyway

