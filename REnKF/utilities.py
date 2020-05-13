""" utility to read input file """

# standard library imports
import sys


def read_input_data(input_file):
    """ Parse an input file into a dictionary.
    """
    param_dict = {}
    try:
        file = open(input_file, "r")
    except IOError:
        print('Cannot open file: {}'.format(input_file))
        sys.exit(1)
    else:
        for line in file:
            # remove comments and trailing/leading whitespace
            line = line.split('#', 1)[0]
            line = line.strip()
            # ignore empty or comment lines
            if not line:
                continue
            # allow ":" and "=" in input file
            line = line.replace(':', ' ').replace('=', ' ')
            # parse into dictionary
            param, value = line.split(None, 1)
            param_dict[param] = value
        file.close()
    return param_dict

