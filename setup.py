import os
from setuptools import setup, find_packages
from setuptools.command.install import install as _install


class install(_install):
    def run(self):
        self._fetch_lfs_objects()
        _install.run(self)

    def _fetch_lfs_objects(self):
        os.system("git lfs pull")


with open('requirements.txt') as f:
    required = f.read().splitlines()


setup(
    name='ppiref',
    author='Anton Bushuiev',
    version='1.2.0',
    packages=find_packages(),
    install_requires=required,
    include_package_data=True,
    cmdclass={
        'install': install,
    }
)
