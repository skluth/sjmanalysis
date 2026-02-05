

import configparser
from collections import OrderedDict

# Support multiple entries for same key:
# From the net somewhere ...
# https://stackoverflow.com/questions/15848674/how-to-configparse-a-file-keeping-multiple-values-for-identical-keys/15848928#15848928

class MultiOrderedDict( OrderedDict ):

    def __setitem__( self, key, value ):
        if key in self and isinstance( value, list ):
            self[key].extend( value )
        else:
            super().__setitem__( key, value )

def readConfig( name="test.cfg" ):
    config= configparser.RawConfigParser( dict_type=MultiOrderedDict, strict=False )
    config.read( name )    
    return config

