import sys, re

input = open(sys.argv[1],'rb').read()
output = re.compile('(?P<a>[\d\.]+)(/[\d\.]*)*').sub('\g<a>',input)
open(sys.argv[2],'wb').write(output)
