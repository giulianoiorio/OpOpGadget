#from setuptools import setup
from distutils.core import setup, Extension

print('Testing Cython installation:')
try:
    from Cython.Build import cythonize
    from Cython.Distutils import build_ext
except ImportError:
    import pip
    pip.main(['install', 'Cython'])
    from Cython.Build import cythonize
    from Cython.Distutils import build_ext

	
	
import numpy
import shutil
import os
import sysconfig
import sys

if sys.version_info[0]==2:
    #time.sleep(5)
    cmdclass_option = {}
    print('You are using Python2, what a shame!')
    #raise ValueError('You are using Python2, what a shame! Download Python3 to use this module. \n If you are using anaconda you can install a python3 virtual env just typing:\n "conda create -n yourenvname python=3.6 anaconda". \n Then you can activate the env with the bash command  "source activate yourenvname"')

elif sys.version_info[0]==3:
    print('You are using Python3, you are a wise person!')

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
    cmdclass_option = {'build_ext': BuildExtWithoutPlatformSuffix}

else:
    raise ValueError('You are not using neither Python2 nor Python3, probably you are a time traveller from the Future or from the Past')

'''
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
'''


# Parse options; current options
# --no-openmp: compile without OpenMP support

extra_compile_args = ['-std=c99','-w']

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
df_c_ext = Extension('OpOp.src.df_src.df_c_ext.df_spherical',
                     sources=df_c_src,
                     extra_compile_args=extra_compile_args
                     )

# Model C extension
model_c_src = ['OpOp/src/model_src/model_c_ext/GeneralModel.c']
model_c_ext = Extension('OpOp.src.model_src.model_c_ext.GeneralModel',
                        sources=model_c_src,
                        extra_compile_args=extra_compile_args
                        )

# Generate Model C extension
genmod_c_src = ['OpOp/src/model_src/model_c_ext/GenerateModel.c', 'OpOp/src/model_src/model_c_ext/MT_random.c']
genmod_c_ext = Extension('OpOp.src.model_src.model_c_ext.GenerateModel',
                         sources=genmod_c_src,
                         extra_compile_args=extra_compile_args)

jsolv_c_src = ['OpOp/src/jsolver_src/jsolver_c_ext/CJsolver.pyx']
jsolv_c_ext = Extension('OpOp.src.jsolver_src.jsolver_c_ext.CJsolver', sources=jsolv_c_src,
                        extra_compile_args=extra_compile_args)

#io ext read
io_c_src_read = ['OpOp/src/io_src/io_c_ext/cread_fvfps.pyx']
io_c_ext_read = Extension('OpOp.src.io_src.io_c_ext.cread_fvfps', sources=io_c_src_read,extra_compile_args=extra_compile_args)

#io ext write
io_c_src_write = ['OpOp/src/io_src/io_c_ext/cwrite_fvfps.pyx']
io_c_ext_write = Extension('OpOp.src.io_src.io_c_ext.cwrite_fvfps', sources=io_c_src_write,extra_compile_args=extra_compile_args)

ext_modules = [df_c_ext, model_c_ext, genmod_c_ext] + cythonize(jsolv_c_ext) + cythonize(io_c_ext_read) + cythonize(io_c_ext_write)

setup(
    name='OpOpGadget',
    version='1.8.dev0',
    author='Giuliano Iorio',
    author_email='giuliano.iorio@unibo.it',
    url='http://github.com/iogiul/OpOp',
    cmdclass=cmdclass_option,
    packages=['OpOp', 'OpOp/src','OpOp/src/df_src', 'OpOp/src/df_src/df_c_ext', 'OpOp/src/grid_src', 'OpOp/src/model_src', 'OpOp/src/model_src/model_c_ext', 'OpOp/src/particle_src',
              'OpOp/src/analysis_src', 'OpOp/src/io_src', 'OpOp/src/io_src/io_c_ext', 'OpOp/src/utility_src', 'OpOp/src/jsolver_src', 'OpOp/src/jsolver_src/jsolver_c_ext',
              'OpOp/src/densityprofile_src'],
    install_requires=['numpy>=1.9', 'scipy>=0.19', 'matplotlib', 'astropy>=1', 'Cython'],
    ext_modules=ext_modules,
    scripts=['OpOp/script/analyse_sims.py'],
    include_dirs=[numpy.get_include()],
    zip_save=False
)

try:
    #shutil.rmtree('build')
    shutil.rmtree('dist')
    shutil.rmtree('OpOpGadget.egg-info')
except:
    pass