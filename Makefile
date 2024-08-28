render:
	quarto render

publish:
	ghp-import -n -p -f _build    

pdfs:
	decktape reveal https://matmek-4270.github.io/matmek4270-pres/intro.html#/title-slide pdfs/intro.pdf
	decktape reveal https://matmek-4270.github.io/matmek4270-pres/analysis.html#/title-slide pdfs/analysis.pdf	
	decktape reveal https://matmek-4270.github.io/matmek4270-pres/vibration.html#/title-slide pdfs/vibration.pdf
	decktape reveal https://matmek-4270.github.io/matmek4270-pres/finitedifferences.html#/title-slide pdfs/finitedifferences.pdf