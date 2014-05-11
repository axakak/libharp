
OUT = tqc
TARGET = bin/$(OUT)
FLAGS = -std=c++11 -ansi -Wall -g -Wno-c++11-extensions
OBJS = tmp/$(OUT).o \
       tmp/traceData.o

$(TARGET): $(OBJS)
	g++ -o $(TARGET) $(OBJS) -lyaml-cpp

tmp/$(OUT).o: src/$(OUT).cpp src/traceData.h
	g++ -c src/$(OUT).cpp -o tmp/$(OUT).o $(FLAGS)

tmp/traceData.o: src/traceData.cpp src/traceData.h
	g++ -c src/traceData.cpp -o tmp/traceData.o $(FLAGS)

clean:
	rm -f $(TARGET) $(OBJS)

