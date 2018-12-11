#This file is forked from https://github.com/pybind/pbtest, original author: Sylvain Corlay

from setuptools import setup, Extension
import setuptools
from setuptools.command.build_ext import build_ext
import os, sys

# get __version__, __author__, and __email__
exec(open("./tacoma/metadata.py").read())

class get_pybind_include(object):
    """Helper class to determine the pybind11 include path
    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked. """

    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        import pybind11
        return pybind11.get_include(self.user)

ext_modules = [
    Extension(
        '_tacoma',
        [ 
            '_tacoma/Utilities.cpp', 
            '_tacoma/Events.cpp', 
            '_tacoma/slice.cpp', 
            '_tacoma/eSIS.cpp', 
            '_tacoma/EdgeActivityModel.cpp',
            '_tacoma/activity_model.cpp',
            '_tacoma/verify_formats.cpp', 
            '_tacoma/measurements.cpp', 
            '_tacoma/social_trajectories.cpp', 
            '_tacoma/edge_trajectories.cpp', 
            '_tacoma/concatenation.cpp', 
            '_tacoma/FW_P_varying.cpp', 
            '_tacoma/FW_P_varying_alpha_beta.cpp', 
            '_tacoma/resampling.cpp', 
            '_tacoma/flockwork_parameter_estimation.cpp', 
            '_tacoma/conversion.cpp', 
            '_tacoma/SIS.cpp', 
            '_tacoma/SIS_node_based.cpp', 
            '_tacoma/SIR.cpp',             
            '_tacoma/SIRS.cpp',
            '_tacoma/SI.cpp', 
            '_tacoma/ZSBB_model.cpp', 
            '_tacoma/dyn_RGG.cpp', 
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
    version = __version__,
    author = __author__,
    author_email = __email__,
    url = 'https://github.com/benmaier/tacoma',
    license = __license__,
    description = 'A package to both analyze real-world temporal networks as well as to simulate a variety of temporal network models and spreading processes on them. Analyze contact durations, group sizes, group durations, social trajectories, aggregated social network, etc.',
    long_description = '',
    packages = setuptools.find_packages(),
    ext_modules = ext_modules,
    setup_requires = [
            'pybind11>=2.0.0'
            ],
    install_requires = [
            'pybind11>=2.0.0',
            'wget',
            'numpy',
            'scipy',
            'lmfit',
            ],
    include_package_data = True,
    cmdclass = {'build_ext': BuildExt},
    zip_safe = False,
)
