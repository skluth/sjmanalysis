
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
    print files
    for fname in files:
        if os.path.isfile( fname ):
            print fname, "exists"
        else:
            download( fname, config )
    return

def download( fname, config ):
    import wget
    import ssl
    url= config.get( "General", "url" )+"/"+fname
    print "Downloading", url
    ssl._https_verify_certificates( 0 )
    wget.download( url, fname )
    print
    return

