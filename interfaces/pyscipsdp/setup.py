import glob
import os
import platform
import subprocess
import sys

import numpy as np
from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext

import pyscipopt

# PySCIPSDP links against SCIP-SDP and the libscip bundled with the PySCIPOpt wheel, so that there is
# only one SCIP in the process. SCIP's headers come from SCIPOPTDIR.

# look for environment variable that specifies path to SCIP (only its headers are used)
scipoptdir = os.environ.get("SCIPOPTDIR", "").strip('"')
if not scipoptdir or not os.path.exists(os.path.join(scipoptdir, "include", "scip", "scip.h")):
    sys.exit("Set SCIPOPTDIR to a SCIP installation (for its headers) of the version used by PySCIPOpt.")
scip_includedir = os.path.abspath(os.path.join(scipoptdir, "include"))

# the headers must be of the same SCIP version as the libscip in the PySCIPOpt wheel
version = {}
with open(os.path.join(scip_includedir, "scip", "config.h")) as f:
    for line in f:
        if line.startswith("#define SCIP_VERSION_"):
            _, key, value = line.split()[:3]
            version[key] = value
model = pyscipopt.Model()
if (version.get("SCIP_VERSION_MAJOR"), version.get("SCIP_VERSION_MINOR")) != (str(model.getMajorVersion()), str(model.getMinorVersion())):
    sys.exit("SCIPOPTDIR has SCIP %s.%s, but PySCIPOpt uses SCIP %d.%d."
             % (version.get("SCIP_VERSION_MAJOR"), version.get("SCIP_VERSION_MINOR"), model.getMajorVersion(), model.getMinorVersion()))

# the libscip bundled with the PySCIPOpt wheel
pyscipoptdir = os.path.dirname(os.path.abspath(pyscipopt.__file__))
if platform.system() == "Darwin":
    libscipdir = os.path.join(pyscipoptdir, ".dylibs")
    libscip = glob.glob(os.path.join(libscipdir, "libscip.*.dylib"))
elif platform.system() == "Linux":
    libscipdir = os.path.join(os.path.dirname(pyscipoptdir), "pyscipopt.libs")
    libscip = glob.glob(os.path.join(libscipdir, "libscip*.so*"))
else:
    sys.exit("PySCIPSDP currently supports Linux and macOS.")
if len(libscip) != 1:
    sys.exit("Could not find the libscip of the PySCIPOpt wheel in %s." % libscipdir)
libscip = libscip[0]

# look for environment variable that specifies path to SCIP-SDP: an installation (cmake --install),
# or by default this repository, which is then built with CMake into build/ (SDP solver from SDPS)
scipsdpdir = os.environ.get("SCIPSDPDIR", "").strip('"')
if scipsdpdir:
    scipsdp_includedirs = [os.path.join(scipsdpdir, "include", "scip"), os.path.join(scipsdpdir, "include")]
    scipsdp_libdir = os.path.join(scipsdpdir, "lib64" if os.path.exists(os.path.join(scipsdpdir, "lib64")) else "lib")
else:
    repodir = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
    builddir = os.path.join(repodir, "build")
    if not os.path.exists(os.path.join(builddir, "CMakeCache.txt")) or "SDPS" in os.environ:
        # no build rpath: libscipsdp must find libscip through our runtime search path (PySCIPOpt's first)
        subprocess.check_call(["cmake", "-S", repodir, "-B", builddir, "-DCMAKE_BUILD_TYPE=Release",
                               "-DCMAKE_PREFIX_PATH=" + scipoptdir, "-DCMAKE_SKIP_BUILD_RPATH=ON",
                               "-DSDPS=" + os.environ.get("SDPS", "none")])
    subprocess.check_call(["cmake", "--build", builddir, "--target", "libscipsdp", "--parallel"])
    scipsdp_includedirs = [os.path.join(repodir, "src", "scipsdp"), os.path.join(repodir, "src")]
    scipsdp_libdir = os.path.join(builddir, "lib")
scipsdp_includedirs = [os.path.abspath(d) for d in scipsdp_includedirs]
scipsdp_libdir = os.path.abspath(scipsdp_libdir)
if not glob.glob(os.path.join(scipsdp_libdir, "libscipsdp.*")):
    sys.exit("Could not find libscipsdp in %s; build SCIP-SDP with CMake or set SCIPSDPDIR." % scipsdp_libdir)

# libscipsdp refers to libscip through @rpath/$ORIGIN, so with the PySCIPOpt wheel's directory first
# in our runtime search path, it uses the same libscip as PySCIPOpt
if platform.system() == "Darwin":
    extra_objects = [libscip]
    extra_link_args = ["-Wl,-rpath,@loader_path/../pyscipopt/.dylibs"]
else:
    extra_objects = []
    extra_link_args = ["-L" + libscipdir, "-l:" + os.path.basename(libscip), "-Wl,-rpath,$ORIGIN/../pyscipopt.libs"]
extra_link_args.append("-Wl,-rpath," + scipsdp_libdir)

# enable debug mode if requested
extra_compile_args = []
if "--debug" in sys.argv:
    extra_compile_args.append("-UNDEBUG")
    sys.argv.remove("--debug")


class BuildExt(build_ext):
    """Refers to the libscip of the PySCIPOpt wheel through @rpath, as the wheel's copy has a placeholder install name."""

    def build_extension(self, ext):
        super().build_extension(ext)
        if platform.system() == "Darwin":
            installname = subprocess.check_output(["otool", "-D", libscip], text=True).split()[-1]
            subprocess.check_call(["install_name_tool", "-change", installname,
                                   "@rpath/" + os.path.basename(libscip), self.get_ext_fullpath(ext.name)])


from Cython.Build import cythonize

packagedir = os.path.join("src", "pyscipsdp")

sources = [os.path.join(packagedir, "scipsdp.pyx")]

extensions = [
    Extension(
        "pyscipsdp.scipsdp",
        sources,
        include_dirs=[scip_includedir] + scipsdp_includedirs + [np.get_include()],
        library_dirs=[scipsdp_libdir],
        libraries=["scipsdp"],
        extra_objects=extra_objects,
        extra_compile_args=extra_compile_args,
        extra_link_args=extra_link_args,
    )
]

# the .pxd files of the installed PySCIPOpt are cimported
extensions = cythonize(extensions, include_path=[os.path.dirname(pyscipoptdir)], compiler_directives={"language_level": 3})

setup(
    ext_modules=extensions,
    packages=["pyscipsdp"],
    package_dir={"pyscipsdp": packagedir},
    package_data={"pyscipsdp": ["scipsdp.pyx", "scipsdp.pxd"]},
    cmdclass={"build_ext": BuildExt},
)
