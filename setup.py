#! /usr/bin/env python
'''
$Id: setup_qubicpack.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Mon 07 Aug 2017 08:35:24 CEST
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

setup.py for qubicpack only.  
Use this to install qubicpack without pystudio
'''
from __future__ import division, print_function
import os,sys,subprocess
from setuptools import setup

DISTNAME         = 'qubicpack'
DESCRIPTION      = 'Utilities for QUBIC detector data visualization'
AUTHOR           = 'Steve Torchinsky'
AUTHOR_EMAIL     = 'satorchi@apc.in2p3.fr'
MAINTAINER       = 'Steve Torchinsky'
MAINTAINER_EMAIL = 'satorchi@apc.in2p3.fr'
URL              = 'https://github.com/satorchi/qubicpack'
LICENSE          = 'GPL'
DOWNLOAD_URL     = 'https://github.com/satorchi/qubicpack'
VERSION          = '3.0'

with open('README.md') as f:
    long_description = f.read()


setup(install_requires=['numpy'],
      name=DISTNAME,
      version=VERSION,
      packages=[DISTNAME],
      zip_safe=False,
      package_data={DISTNAME: ['data/*']},
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      maintainer=MAINTAINER,
      maintainer_email=MAINTAINER_EMAIL,
      description=DESCRIPTION,
      license=LICENSE,
      url=URL,
      download_url=DOWNLOAD_URL,
      long_description=long_description,
      classifiers=[
          'Development Status :: 3 - Alpha',
          'Environment :: Console',
          'Intended Audience :: Developers',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License (GPL)',
          'Topic :: Scientific/Engineering'],
)

# install the executable scripts
exec_dir = '/usr/local/bin'
scripts = ['scripts/quicklook.py']
if len(sys.argv)>1 and sys.argv[1]=='install':
    print('installing executable scripts...')
    for F in scripts:
        basename = os.path.basename(F)
        cmd = 'rm -f %s/%s; cp -puv %s %s;chmod +x %s/%s' % (exec_dir,basename,F,exec_dir,exec_dir,basename)
        proc=subprocess.Popen(cmd,stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out,err=proc.communicate()
        if out:print(out.decode().strip())
        if err:print(err.decode().strip())

