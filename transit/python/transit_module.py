# This file was automatically generated by SWIG (http://www.swig.org).
# Version 3.0.5
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.





from sys import version_info
if version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_transit_module', [dirname(__file__)])
        except ImportError:
            import _transit_module
            return _transit_module
        if fp is not None:
            try:
                _mod = imp.load_module('_transit_module', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _transit_module = swig_import_helper()
    del swig_import_helper
else:
    import _transit_module
del version_info
try:
    _swig_property = property
except NameError:
    pass  # Python < 2.2 doesn't have 'property'.


def _swig_setattr_nondynamic(self, class_type, name, value, static=1):
    if (name == "thisown"):
        return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name, None)
    if method:
        return method(self, value)
    if (not static):
        if _newclass:
            object.__setattr__(self, name, value)
        else:
            self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)


def _swig_setattr(self, class_type, name, value):
    return _swig_setattr_nondynamic(self, class_type, name, value, 0)


def _swig_getattr_nondynamic(self, class_type, name, static=1):
    if (name == "thisown"):
        return self.this.own()
    method = class_type.__swig_getmethods__.get(name, None)
    if method:
        return method(self)
    if (not static):
        return object.__getattr__(self, name)
    else:
        raise AttributeError(name)

def _swig_getattr(self, class_type, name):
    return _swig_getattr_nondynamic(self, class_type, name, 0)


def _swig_repr(self):
    try:
        strthis = "proxy of " + self.this.__repr__()
    except:
        strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except AttributeError:
    class _object:
        pass
    _newclass = 0



def transit_init(argc, argv):
    return _transit_module.transit_init(argc, argv)
transit_init = _transit_module.transit_init

def get_no_samples():
    return _transit_module.get_no_samples()
get_no_samples = _transit_module.get_no_samples

def get_waveno_arr(waveno_arr):
    return _transit_module.get_waveno_arr(waveno_arr)
get_waveno_arr = _transit_module.get_waveno_arr

def run_transit(re_input, transit_out):
    return _transit_module.run_transit(re_input, transit_out)
run_transit = _transit_module.run_transit

def free_memory():
    return _transit_module.free_memory()
free_memory = _transit_module.free_memory
# This file is compatible with both classic and new-style classes.

