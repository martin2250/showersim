CPPFLAGS=-O3

.PHONY: plot
plot: plot.py result.txt
	./$<

result.txt: showersim
	./showersim > $@