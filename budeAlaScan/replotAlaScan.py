#!/usr/bin/env python3
# encoding: utf-8
'''
replotAlaScan -- Re-plot the results from the budeAlaScan by reading the json file created by the app.

@author:     Amaurys Ávila Ibarra

@copyright:  2017 University of Bristol. All rights reserved.

@license:    license

@contact:    user_email
@deffield    updated: Updated
'''

from argparse import ArgumentParser, RawDescriptionHelpFormatter
import sys, os

from budeAlaScan.plots.replot_ala_scan import replot_ala_scan

__all__ = []
__version__ = 0.1
__date__ = '2017-07-21'
__updated__ = '2017-07-21'


def main(argv=None):
    '''Command line options.'''

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
    program_shortdesc = __import__('__main__').__doc__.split("\n")[1]

    program_license = '''%s

  Created on %s.

  Licensed under the Apache License 2.0
  http://www.apache.org/licenses/LICENSE-2.0

  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE
''' % (program_shortdesc, str(__date__))

    try:
        # Setup argument parser
        parser = ArgumentParser(description=program_license, formatter_class=RawDescriptionHelpFormatter)
        parser.add_argument('-V', '--version', action='version', version=program_version_message)
        mandatory_options = parser.add_argument_group('mandatory_args')
        mandatory_options.add_argument("-j", "--json-file", dest="json_fname", help="Name of the json file created by budeAlaScan [default: %(default)s]", metavar="fname_scan.json", required=True)

        # Process arguments
        args = parser.parse_args()

    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0
    except Exception as e:
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help")
        return 2

    replot_ala_scan(args.json_fname)


if __name__ == "__main__":
    # If no arguments is given show help.
    if len(sys.argv) is 1:
        sys.argv.append("-h")
    sys.exit(main())
