from distutils.core import setup
import numpy
import py2exe
import glob
import matplotlib
import pytz
import pytz.tzinfo

# # cd Z:\Software Projects\siFi_2015
# inno

data_files = matplotlib.get_py2exe_datafiles()

setup(windows=[{"script": "main.py",
                }],

      options={"py2exe":
                {"includes": ['sip', 'decimal', 'PyQt4.QtSql', "PyQt4.QtCore", "PyQt4.QtGui", "PyQt4.QtNetwork",
                              "matplotlib.backends.backend_tkagg", r'scipy.sparse.csgraph._validation',
                                r'scipy.special._ufuncs_cxx', 'zmq.backend.cython'],
                  #"bundle_files":1,
                  #"optimize": 2,
                  "dll_excludes": ["mswsock.dll", "powrprof.dll", "MSVCP90.dll", 'libzmq.pyd', 'libgdk-win32-2.0-0.dll',
                                   'libgobject-2.0-0.dll', "libgdk_pixbuf-2.0-0.dll"],
                 "packages": ['matplotlib','pytz']
                }},
      data_files=data_files)





#options = {"build_exe" : {"includes" : "atexit" }}

# from distutils.core import setup
# import py2exe, sys, os
#
# import matplotlib as mpl
# mpl.use('Qt4Agg') #I use this is GUI.py
#
# sys.argv.append('py2exe')
#
# includes = ['sip', 'PyQt4', 'PyQt4.QtGui', 'PyQt4.QtCore', 'matplotlib.backends']
#
# excludes = ['_gtkagg', '_tkagg', 'bsddb', 'curses', 'pywin.debugger',
#             'pywin.debugger.dbgcon', 'pywin.dialogs', 'tcl',
#             'Tkconstants', 'pydoc', 'doctest', 'test', 'wx']
#
# packages = ['matplotlib'] #should my custom modules go in here?
#
# dll_excludes = ["MSVCP90.dll", "MSWSOCK.dll", "mswsock.dll", "powrprof.dll",
#             'libgdk-win32-2.0-0.dll', 'libgobject-2.0-0.dll', 'tcl84.dll',
#             "libgdk_pixbuf-2.0-0.dll"]
#
# data_files = mpl.get_py2exe_datafiles()
#
#
# setup(name = "main",
#       windows = [{"script":"main.py"}],
#       options = {'py2exe': {
#
#           'includes': includes,
#           'excludes': excludes,
#           'packages': packages,
#           'dll_excludes': dll_excludes,
#           'bundle_files': 1}},
#       zipfile = None,
#       data_files = data_files
# )