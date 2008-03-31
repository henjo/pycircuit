from distutils.core import setup

setup(name='Pycircuit',
      version='0.0',
      description='Python circuit design tools',
      author='Henrik Johansson, Joacim Olsson, Andreas Drejfert',
      author_email='henjo2006@gmail.com',
      url='http://rigel.johome.net/svn/pycircuit'
      packages=['pycircuit', 'pycircuit.circuit', 'pycircuit.cds'],
      package_dir={'pycircuit': 'src'}
     )
