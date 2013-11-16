
OUT = tqc
TARGET = bin/$(OUT)
FLAGS = -ansi -Wall -g
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

