
OUT = tqc
TARGET = bin/$(OUT)
FLAGS = -ansi -Wall -g
OBJS = tmp/$(OUT).o


$(TARGET): $(OBJS)
	g++ -o $(TARGET) $(OBJS) -lyaml-cpp

tmp/$(OUT).o: src/$(OUT).cpp src/traceData.h
	g++ -c src/$(OUT).cpp -o tmp/$(OUT).o $(FLAGS)

clean:
	rm -f $(TARGET) $(OBJS)

