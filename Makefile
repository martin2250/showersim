CPPFLAGS=-O3 -std=c++17 -Wall -Wextra -Wno-cpp -Werror=reorder -Werror=unused-variable 

.PHONY: plot
plot: plot.py result.txt
	./$<

result.txt: showersim
	./showersim > $@