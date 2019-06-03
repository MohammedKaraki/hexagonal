EXECUTABLE = x
SOURCES = main.cpp
OBJECTS = $(patsubst %.cpp, %.o, $(SOURCES))

CXX = g++
CXX_GENERAL_FLAGS = -g -O2 -std=c++1z -Wall -Wpedantic -Wextra
CXX_LIBRARY_FLAGS = -lncurses

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(CXX_LIBRARY_FLAGS) -o $(EXECUTABLE) $(OBJECTS)

$(OBJECTS): %.o: %.cpp
	$(CXX) $(CXX_GENERAL_FLAGS)  -c $< -o $@

.dep: $(SOURCES)
	rm -f ./.dep
	$(CXX) $(CXX_GENERAL_FLAGS) -MM $^ > .dep


-include .dep

.PHONY: clean

clean:
	rm -rf *.o
	rm -rf .dep
