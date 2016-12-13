import re
import fileinput
import sys

if __name__ == '__main__':
    for f in sys.argv[1:]:
        for line in fileinput.input(f, inplace=1):
            re.sub(r"_\d", "", line)
