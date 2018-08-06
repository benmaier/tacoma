import _tacoma as _tc
import api as tc

original_methods = [d for d in dir(_tc) if not d.startswith('_') ]
api_methods = [d for d in dir(tc) if not d.startswith('_') ]

print(original_methods, api_methods)

mport 
