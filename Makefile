doc:
	python tools/implementation.py cylon/cylon.cpp  >docs/implementation.tex
	cd docs; pdflatex cylon.tex
	cd docs; pdflatex cylon.tex
	cd docs; pdflatex cylon.tex	
code:
	#python tools/make_cylon.py docs/cylon.tex > cylon.cpp
	#g++ cylon.cpp
bin:
	cd cylon; g++ -I../include -framework GLUT -framework OpenGL -framework Cocoa cylon.cpp -o ../cylon.bin 