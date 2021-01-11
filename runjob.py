
import SjmConfigParser

import os.path

def checkFiles( runconfig="sjmconfig_91_96.cfg" ):
    configfiles= [ runconfig, "sjmGeneralOptions.cfg" ]
    config= SjmConfigParser.readConfig( configfiles )
    files= config.get( "Data", "files" ).split()
    files+= config.get( "Signal", "files" ).split()
    files+= config.get( "AltSignal", "files" ).split()
    for section in [ "BkgWWllqq", "BkgWWqqqq", "BkgWWeeqq" ]:
        if config.has_section( section ):
            files+= config.get( section, "files" ).split()
    path= config.get( "General", "path" )
    for fname in files:
        filepath= os.path.join( path, fname )
        if os.path.isfile( filepath ):
            print fname, "exists", "in", path
        else:
            download( fname, path, config )
    return

def download( fname, path, config ):
    import wget
    import ssl
    url= config.get( "General", "url" )+"/"+fname
    print "Downloading", url, "into", path
    ssl._https_verify_certificates( 0 )
    wget.download( url, path )
    print
    return

