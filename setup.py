from distutils.core import setup

setup(name='Pycircuit',
      version='0.0',
      description='Python circuit design tools',
      author='Henrik Johansson, Joacim Olsson, Andreas Drejfert',
      author_email='henjo2006@gmail.com',
      url='http://rigel.johome.net/svn/pycircuit',
      packages=['pycircuit', 'pycircuit.circuit', 'pycircuit.post', 'pycircuit.utilities',
                'pycircuit.sim', 'pycircuit.post.cds', 'pycircuit.post.jwdb'],
     )
