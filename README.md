# EpiFlockworks

For a fast epidemic simulation on flockworks.

## Install

For all systems, first clone this repository.

### Matlab (Mac OSX)

You need to have the current XCode version installed (free in AppStore). Open Matlab and change into the directory of the repository. Type

```matlab
>>> edit ([matlabroot '/bin/maci64/mexopts/clang++_maci64.xml'])
>>> edit ([matlabroot '/bin/maci64/mexopts/clang_maci64.xml'])
```

In both files, copy lines matching occurences of `MacOSX10.x.sdk` and change `MacOSX10.x.sdk` to `MacOSX10.11.sdk`(or whichever current version of XCode you're using).

Now, run


```matlab
>>> setup
>>> epitest
```

### Python

    $ sudo pip ./EpiFlockwork install

## Example

    $ python sandbox/eq_test.py
