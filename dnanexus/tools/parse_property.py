#!/usr/bin/env python2.7
# parse_property.py  Reads property string and prses a requested value.
#                    Write request to stdout and verbose info to stderr.  This allows easy use in dx app scripts.

# imports needed for Settings class:
import os, sys, string, argparse, json
import dxpy

def env_get_current_project_id():
    ''' Returns the current project name for the command-line environment '''
    err, proj_name = commands.getstatusoutput('cat ~/.dnanexus_config/DX_PROJECT_CONTEXT_NAME')
    if err != 0:
        return None
    proj = dxencode.get_project(proj_name)
    return proj.get_id()
    
    return proj_name

def file_get_property(filePath,key,subkey,return_json=False,verbose=False):
    '''Returns dx file's property matching 'key'.'''
    dxfile = None
    try:
        dxfile = dxpy.get_handler(filePath)
    except:
        try:
            dxlink = dxpy.dxlink(filePath)
            dxfile = dxpy.get_handler(dxlink)
        except:
            try:
                proj_id = env_get_current_project_id()
                dxfile = dxpy.DXFile(filePath,project=proj_id)
            except:
                sys.stderr.write('ERROR: unable to find file "' + filePath + '": \n')
                sys.exit(0)  # Do not error on tool run in dx script 
                
    if dxfile == None:
        sys.stderr.write('ERROR: unable to find file "' + filePath + '": \n')
        sys.exit(0)  # Do not error on tool run in dx script 
    
    props = dxfile.get_properties()
    if not props:
        sys.stderr.write('ERROR: unable to find properties for file "' + filePath + '": \n') 
        sys.exit(0)  # Do not error on tool run in dx script 
    
    if key not in props:
        sys.stderr.write('ERROR: unable to find "'+key+'" in properties for file "' + filePath + '": \n') 
        sys.exit(0)  # Do not error on tool run in dx script
    props = props[key]
         
    if return_json or subkey != None:
        try:
            props = json.loads(props)
        except:
            try:
                props = json.loads("{"+props+"}")
            except:
                sys.stderr.write('Failure parsing "'+props+'" as json.\n') 
                sys.exit(0)  # Do not error on tool run in dx script

    if subkey != None:
        if subkey not in props:
            sys.stderr.write('ERROR: unable to find "'+subkey+'" in properties for file "' + filePath + '": \n') 
            sys.exit(0)  # Do not error on tool run in dx script
        props = props[subkey]
        
    if verbose:
        sys.stderr.write(props + '\n')
    
    return props

def main():
    parser = argparse.ArgumentParser(description =  "Creates a json string of qc_metrics for a given applet. " + \
                                                    "Returns string to stdout and formatted json to stderr.")
    parser.add_argument('-f', '--file',
                        help='DX id, link or path to file.',
                        required=True)
    parser.add_argument('-p','--property',
                        help="Property name.",
                        default='QC',
                        required=False)
    parser.add_argument('-s','--subproperty',
                        help="Property name.",
                        default=None,
                        required=False)
    parser.add_argument('-k', '--key',
                        help='Prints just the value for this key.',
                        default=None,
                        required=False)
    parser.add_argument('--keypair',
                        help='Prints the key: value pair for this key.',
                        default=None,
                        required=False)
    parser.add_argument('-j', '--json', action="store_true", required=False, default=True, 
                        help="Expect json.")
    parser.add_argument('-q', '--quiet', action="store_true", required=False, default=False, 
                        help="Suppress non-error stderr messages.")
    parser.add_argument('-v', '--verbose', action="store_true", required=False, default=False, 
                        help="Make some noise.")

    args = parser.parse_args(sys.argv[1:])
    if len(sys.argv) < 2:
        parser.print_usage()
        return
        
    properties = file_get_property(args.file,args.property,args.subproperty,return_json=args.json,verbose=args.verbose)
    
    # Print out the properties
    if args.key != None:
        if args.key in properties:
            print json.dumps(properties[args.key])
            if not args.quiet:
                sys.stderr.write(json.dumps(properties[args.key],indent=4) + '\n')
        else:
            print ''   
            if not args.quiet:
                sys.stderr.write('(not found)\n')
    elif args.keypair != None:
        if args.keypair in properties:
            print '"' + args.keypair + '": ' + json.dumps(properties[args.keypair])
            if not args.quiet:
                sys.stderr.write('"' + args.keypair + '": ' + json.dumps(properties[args.keypair],indent=4) + '\n')
        else:
            print '"' + args.keypair + '": '
            if not args.quiet:
                sys.stderr.write('"' + args.keypair + '": \n')
    else: 
        print json.dumps(properties)
        if not args.quiet:
            sys.stderr.write(json.dumps(properties,indent=4) + '\n')
    
if __name__ == '__main__':
    main()

