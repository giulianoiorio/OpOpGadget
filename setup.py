from setuptools import setup
from distutils.core import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext
import numpy
import shutil
import os
import sysconfig


def get_ext_filename_without_platform_suffix(filename):
    name, ext = os.path.splitext(filename)
    ext_suffix = sysconfig.get_config_var('EXT_SUFFIX')

    if ext_suffix == ext:
        return filename

    ext_suffix = ext_suffix.replace(ext, '')
    idx = name.find(ext_suffix)

    if idx == -1:
        return filename
    else:
        return name[:idx] + ext


class BuildExtWithoutPlatformSuffix(build_ext):
    def get_ext_filename(self, ext_name):
        filename = super().get_ext_filename(ext_name)
        return get_ext_filename_without_platform_suffix(filename)


# Parse options; current options
# --no-openmp: compile without OpenMP support

extra_compile_args = ['-std=c99']

option_list = ['-CC', '--O1', '--O2', '--O0', '--no-openmp']
command_list = ['build', 'install', 'develop', 'registrate']

# optional params
# CC
'''
try:
    cc_pos = sys.argv.index('--CC')
    cc_val=sys.argv[cc_pos+1]
    CC=cc_val
    if (cc_val in option_list) or (cc_val in command_list):
        print(cc_val)
        raise IOError()
except ValueError:
    CC='gcc'
except IOError:
    CC='gcc'
    del sys.argv[cc_pos]
else:
    del sys.argv[cc_pos]
    del sys.argv[cc_pos]
os.environ["CC"] = CC
#optimization
#no-openmp
try:
    openmp_pos = sys.argv.index('--no-openmp')
except ValueError:
    extra_compile_args.append("-fopenmp")
else:
    del sys.argv[openmp_pos]
#optimazione level
try:
    o_lev = sys.argv.index('--O0')
except ValueError:
    pass
else:
    extra_compile_args.append("-O0")
    del sys.argv[o_lev]
#optimazione level
try:
    o_lev = sys.argv.index('--O1')
except ValueError:
    pass
else:
    extra_compile_args.append("-O1")
    del sys.argv[o_lev]
	#optimazione level
try:
    o_lev = sys.argv.index('--O2')
except ValueError:
    pass
else:
    extra_compile_args.append("-O2")
    del sys.argv[o_lev]
'''

# df C extension
df_c_src = ['OpOp/src/df_src/df_c_ext/spherical.c']
df_c_ext = Extension('OpOp/src/df_src/df_c_ext/df_spherical',
                     sources=df_c_src,
                     extra_compile_args=extra_compile_args
                     )

# Model C extension
model_c_src = ['OpOp/src/model_src/model_c_ext/GeneralModel.c']
model_c_ext = Extension('OpOp/src/model_src/model_c_ext/GeneralModel',
                        sources=model_c_src,
                        extra_compile_args=extra_compile_args
                        )

# Generate Model C extension
genmod_c_src = ['OpOp/src/model_src/model_c_ext/GenerateModel.c', 'OpOp/src/model_src/model_c_ext/MT_random.c']
genmod_c_ext = Extension('OpOp/src/model_src/model_c_ext/GenerateModel',
                         sources=genmod_c_src,
                         extra_compile_args=extra_compile_args)

jsolv_c_src = ['OpOp/src/jsolver_src/jsolver_c_ext/CJsolver.pyx']
jsolv_c_ext = Extension('OpOp/src/jsolver_src/jsolver_c_ext/CJsolver', sources=jsolv_c_src,
                        extra_compile_args=['-fopenmp', ], extra_link_args=['-fopenmp', ])

#io ext read
io_c_src_read = ['OpOp/src/io_src/io_c_ext/cread_fvfps.pyx']
io_c_ext_read = Extension('OpOp/src/io_src/io_c_ext/cread_fvfps', sources=io_c_src_read)

#io ext write
io_c_src_write = ['OpOp/src/io_src/io_c_ext/cwrite_fvfps.pyx']
io_c_ext_write = Extension('OpOp/src/io_src/io_c_ext/cwrite_fvfps', sources=io_c_src_write)

ext_modules = [df_c_ext, model_c_ext, genmod_c_ext] + cythonize(jsolv_c_ext) + cythonize(io_c_ext_read) + cythonize(io_c_ext_write)

setup(
    name='OpOpGadget',
    version='1.8',
    author='Giuliano Iorio',
    author_email='giuliano.iorio@unibo.it',
    url='http://github.com/iogiul/OpOp',
    cmdclass={'build_ext': BuildExtWithoutPlatformSuffix},
    package_dir={'OpOp/src/': ''},
    packages=['OpOp', 'OpOp/src/df_src', 'OpOp/src/grid_src', 'OpOp/src/model_src', 'OpOp/src/particle_src',
              'OpOp/src/analysis_src', 'OpOp/src/io_src', 'OpOp/src/utility_src', 'OpOp/src/jsolver_src',
              'OpOp/src/densityprofile_src'],
    install_requires=['numpy>=1.9', 'scipy>=0.16', 'matplotlib', 'astropy>=1', 'fermi'],
    ext_modules=ext_modules,
    scripts=['OpOp/script/analyse_sims.py'],
    include_dirs=[numpy.get_include()]
)

shutil.rmtree('build')
shutil.rmtree('dist')
shutil.rmtree('OpOpGadget.egg-info')