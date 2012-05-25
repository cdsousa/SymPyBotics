import pickle
import hashlib
import os
import time
import re



def memoized( func, extra_deps='', cache_folder=None, hash_args_by_str=False, debug=False ):

  encode_pkl = lambda x: pickle.dumps( x )
  #encode_str = lambda x: str( x ).encode()
  encode_str = lambda x: re.sub( ' at 0x[0-9a-f]+', ' at 0xXXXXXXXX', str(x) ).encode()
  
  if hash_args_by_str:
    encode_args = encode_str
  else:
    encode_args = encode_pkl

  def decorated_function( *args, **kwargs ):

    funcstr = encode_pkl(func)
    argsstr = encode_args(args)
    kwargsstr = encode_args(kwargs)
    extrastr = encode_args(extra_deps)

    if debug:
      print('funcstr  ',hashlib.sha1(funcstr).hexdigest())
      print('argsstr  ',hashlib.sha1(argsstr).hexdigest())
      print('kwargsstr',hashlib.sha1(kwargsstr).hexdigest())
      print('extrastr ',hashlib.sha1(extrastr).hexdigest())
    
    bytestr = funcstr + argsstr + kwargsstr + extrastr

    sha1 = hashlib.sha1(bytestr).hexdigest()
    
    if debug:
      print('all      ',sha1)

    filename = func.__name__ + '-' + sha1 + '.pkl'

    pathfile = os.path.join( cache_folder, filename )

    try:
      with open( pathfile, 'rb' ) as file:
        cached_return = pickle.load( file )

    except:
      print('(memoization cache miss ...)')

      if cache_folder and not os.path.exists( cache_folder ):
        os.makedirs( cache_folder )

      t0 = time.time()

      call_return = func( *args, **kwargs )

      print( '(... memoization call wall time: %.3f)' % (time.time()-t0) )

      if debug:
        print('call     ',hashlib.sha1(str(call_return).encode()).hexdigest())

      with open( pathfile, 'wb' ) as file:
          pickle.dump( call_return, file )

      return call_return
      
    print('(memoization cache hit)')
      
    if debug:
      print('cache    ',hashlib.sha1(str(cached_return).encode()).hexdigest())
        
    return cached_return

  return decorated_function
