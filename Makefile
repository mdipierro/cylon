doc:
	python tools/implementation.py cylon/cylon.cpp  >docs/implementation.tex
	cd docs; pdflatex cylon.tex
	cd docs; pdflatex cylon.tex
	cd docs; pdflatex cylon.tex	
bin:
	cd cylon; g++ -I../include -framework GLUT -framework OpenGL -framework Cocoa cylon.cpp -o ../cylon.bin 
publish:
	make doc
	make bin
	cp docs/cylon.pdf ~/Dropbox/Public/DePaul/GAM450/game_physics_nutshell.pdf 
	cd ..; zip -r ~/Dropbox/Public/DePaul/GAM450/cylon.zip cylon/*
	git push
