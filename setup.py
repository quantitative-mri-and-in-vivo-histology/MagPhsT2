from distutils.core import setup
from Cython.Build import cythonize
from Cython.Distutils import build_ext
import numpy

# Avoid a gcc warning below:

class BuildExt(build_ext):
    def build_extensions(self):
        #self.compiler.compiler_so.remove('-Wstrict-prototypes')
        super(BuildExt, self).build_extensions()

setup(
    name='PyEPG',
    version='0.1',
    cmdclass={'build_ext': BuildExt},
    ext_modules=cythonize('pyepg.pyx', compiler_directives={'embedsignature': True}),
    include_dirs = [numpy.get_include()]
)
