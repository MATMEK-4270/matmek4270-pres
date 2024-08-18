render:
	quarto render

publish:
    ghp-import -n -p -f _build    

pdf:
	decktape reveal https://matmek-4270.github.io/matmek4270-pres/intro.html#/title-slide intro.pdf
	decktape reveal https://matmek-4270.github.io/matmek4270-pres/analysis.html#/title-slide analysis.pdf	
