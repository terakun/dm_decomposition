CXX = g++
FLAGS = -std=c++11 -O0 -g
# FLAGS = -std=c++11 -O2 -march=native
# FLAGS = -std=c++11 -O2
INCLUDES = -I/usr/local/include -I/usr/local/include/eigen3/
LIBS = -lm -L/usr/local/lib

TARGET = main
SRCS = main.cc dm_decomp.cc

OBJS = $(patsubst %.cc,%.o,$(SRCS))

.cc.o:
	$(CXX) $(FLAGS) $(INCLUDES) -c $< -o $@

$(TARGET): $(OBJS)
		$(CXX) $(FLAGS) $(FLAG) -o $@ $(OBJS) $(LIBS)

clean:
	rm $(TARGET) $(OBJS)
