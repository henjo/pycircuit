from distutils.core import setup, Extension

gnucap_ext = Extension('gnucap',
                       sources = ['gnucap.cc'],
                       libraries = ['gnucap']
                       )

setup (name = 'gnucap-ext',
       version = '1.0',
       description = 'Gnucap python extension',
       ext_modules = [gnucap_ext])
