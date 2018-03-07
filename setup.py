#This file is forked from https://github.com/pybind/pbtest, original author: Sylvain Corlay

from setuptools import setup, Extension
import setuptools
from setuptools.command.build_ext import build_ext
import setuptools
import os, sys
from pip import locations

class get_pybind_include(object):

    def __init__(self,user=False):
        self.user = user

    def __str__(self):
        pybind_include = os.path.dirname(locations.distutils_scheme('pybind11',self.user)['headers'])
        return pybind_include

ext_modules = [
    Extension(
        '_tacoma',
        [ 
            '_tacoma/Utilities.cpp', 
            '_tacoma/Events.cpp', 
            '_tacoma/SIS.cpp', 
            '_tacoma/SIR.cpp',             
            '_tacoma/SI.cpp', 
            '_tacoma/measurements.cpp', 
            '_tacoma/resampling.cpp', 
            '_tacoma/social_trajectories.cpp', 
            '_tacoma/ZSBB_model.cpp', 
            '_tacoma/dyn_RGG.cpp', 
            '_tacoma/FW_P_varying.cpp', 
            '_tacoma/Flockwork.cpp', 
            '_tacoma/_tacoma.cpp', 
        ],
        include_dirs=[
            get_pybind_include(),
            get_pybind_include(user=True),
            "./_tacoma/"
        ],
        language='c++',
    ),
]

def has_flag(compiler, flagname):
    """Return a boolean indicating whether a flag name is supported on
    the specified compiler.
    """
    import tempfile
    fd, fname = tempfile.mkstemp('.cpp', 'main', text=True)
    with os.fdopen(fd, 'w') as f:
        f.write('int main (int argc, char **argv) { return 0; }')
    try:
        compiler.compile([fname], extra_postargs=[flagname])
    except setuptools.distutils.errors.CompileError:
        return False
    return True

def cpp_flag(compiler):
    """Return the -std=c++[11/14] compiler flag.
    The c++14 is preferred over c++11 (when it is available).
    """
    if has_flag(compiler, '-std=c++14'):
        return '-std=c++14'
    elif has_flag(compiler, '-std=c++11'):
        return '-std=c++11'
    else:
        raise RuntimeError('Unsupported compiler -- at least C++11 support is needed!')


class BuildExt(build_ext):
    """A custom build extension for adding compiler-specific options."""
    c_opts = {
        'msvc': ['/EHsc'],
        'unix': [],
    }

    if sys.platform == 'darwin':
        c_opts['unix'] += ['-stdlib=libc++', '-mmacosx-version-min=10.7','-ftemplate-depth=1024']

    def build_extensions(self):
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        if ct == 'unix':
            opts.append(cpp_flag(self.compiler))
            if has_flag(self.compiler, '-fvisibility=hidden'):
                opts.append('-fvisibility=hidden')
        for ext in self.extensions:
            ext.extra_compile_args = opts
        build_ext.build_extensions(self)

setup(
    name = 'tacoma',
    version = '0.0.17',
    author = 'Benjamin F. Maier',
    author_email = 'bfmaier@physik.hu-berlin.de',
    url = 'https://github.com/benmaier/tacoma',
    license = 'MIT',
    description = 'A package to both analyze real-world temporal networks as well as to simulate a variety of temporal network models and spreading processes on them. Analyze contact durations, group sizes, group durations, social trajectories, aggregated social network, etc.',
    long_description = '',
    packages = setuptools.find_packages(),
    ext_modules = ext_modules,
    install_requires = ['pybind11'],
    cmdclass = {'build_ext': BuildExt},
    zip_safe = False,
)
