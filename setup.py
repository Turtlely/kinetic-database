#!/usr/bin/env python
#test
from distutils.core import setup
from distutils.command.build_py import build_py

packages = ['kdb']
package_dir = {'kdb':'kdb'}
scripts = ['bin/kdb_local_client.py','bin/kdb_remote_client.py']
package_data = {'kdb' : ['pymysql/*']}

setup(name='kdb',
    version='1.0',
    description='kinetic database',
    author='Henkelman Group',
    url='http://www.henkelmanlab.org/',
    packages=packages,
    scripts=scripts,
    package_dir=package_dir,
    package_data=package_data,
    cmdclass={'build_py' : build_py}
    )
