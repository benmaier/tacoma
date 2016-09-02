# cFlockwork

For fast simulations of epidemic spreading on flockworks, random rewiring networks and static networks.

## Install

For all systems, first clone this repository.

### Matlab (Mac OSX)

You need to have the current XCode version installed (free in AppStore). Open Matlab and change into the directory of the repository. At first, there's two files you need to change.

```matlab
>>> cd /path/to/repository
>>> edit ([matlabroot '/bin/maci64/mexopts/clang++_maci64.xml'])
>>> edit ([matlabroot '/bin/maci64/mexopts/clang_maci64.xml'])
```

In both files, copy lines matching occurences of `MacOSX10.x.sdk` and change `MacOSX10.x.sdk` to `MacOSX10.11.sdk`(or whichever current version of XCode you're using).

Now, run


```matlab
>>> setup
>>> cd sandbox
>>> flockworktest
```

### Python

    $ sudo pip install ./EpiFlockwork

## Example

### Python

    $ python sandbox/meanfieldtest.py

### Matlab

    >>> cd sandbox
    >>> flockworktest
