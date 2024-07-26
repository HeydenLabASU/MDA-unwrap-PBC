# print("Running setup.py...")

from setuptools import setup, find_packages
from setuptools.command.install import install
from setuptools.command.build_clib import build_clib
import subprocess
import os
# import some_build_toolkit

project_name = 'MDA-unwrap-PBC'

project_folder = 'MDA_unwrap_PBC'
clib_folder = 'ctypes_lib'

# def get_version():
#     version = some_build_toolkit.compute_version()
#     return version
class CustomBuildClibCommand(build_clib):
    """Custom build command class to compile the required C-library."""
    def run(self):
        # Directory containing the compile script
        # print("Custom build_clib command running...")
        ctypes_lib_dir = os.path.join(os.path.dirname(__file__), project_folder, clib_folder)
        compile_filename = 'compile.sh'
        # Changing file permission
        subprocess.check_call(
            ['chmod', '+x', compile_filename],
            cwd=ctypes_lib_dir
        )
        # Compiling the C-library
        subprocess.check_call(
            ['sh', './' + compile_filename],
            cwd=ctypes_lib_dir
        )
        # print("Compilation finished")
        # Call the standard build_clib command
        super().run()

class CustomInstallCommand(install):
    def run(self):
        # print("Custom install command running...")
        self.run_command('build_clib')  # Ensures build_clib is run
        # Call the standard install command
        super().run()

setup(
    # name=project_name,
    # # version='0.1.0',  # Replace with accuarte version
    # # description='A brief description of the MDA_unwrap_PBC package',
    # # author='Matthias Heyden',
    # # author_email='mheyden1@asu.edu',
    # # url='https://github.com/HeydenLabASU/MDA-unwrap-PBC/',  # Project's URL
    # packages=find_packages(include=[project_folder, project_folder + '.*']),
    # install_requires=[
    #     # Package dependencies
    #     # e.g. 'numpy>=1.18.1', leaving empty for now
    # ],
    # extras_require={
    #     'dev': [
    #         # Additional dependencies for development
    #         # e.g. 'pytest', 'sphinx'
    #         # Base depends
    #         'python',
    #         'pip',

    #         # MDAnalysis
    #         'MDAnalysis',

    #         # Testing
    #         'pytest',
    #         'pytest-cov',
    #         'pytest-xdist',
    #         'codecov'

    #         # Pip-only installs
    #         #pip:
    #         #  codecov
    #         ],
    #     'docs': [
    #         # Dependencies for building documentation
    #         'python',
    #         'pip',

    #         'mdanalysis-sphinx-theme >=1.0.1',
    #         # Other documentation-related packages
    #     ],
    # },
    # python_requires='>=3.7',  # Python versions it supports
    cmdclass={
        'build_clib': CustomBuildClibCommand,  # Custom command class for building requisite C libraries
        'install': CustomInstallCommand,
    },
)
