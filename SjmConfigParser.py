
import ConfigParser
# import os.path
from collections import OrderedDict

# Support multiple entries for same key:
class MultiOrderedDict( OrderedDict ):

    def __setitem__( self, key, value ):
        if key in self and isinstance( value, list ):
            self[key].extend( value )
        else:
            super( MultiOrderedDict, self ).__setitem__( key, value )

def readConfig( name="test.cfg" ):
    config= ConfigParser.ConfigParser( dict_type=MultiOrderedDict )
    config.read( name )    
    return config

