
OUT = tqc
TARGET = bin/$(OUT)
FLAGS = -ansi -Wall
OBJS = tmp/$(OUT).o

YAML_CPP_DIR = /usr/local/Cellar/yaml-cpp/0.5.1
YAML_CPP_LIB = $(YAML_CPP_DIR)/lib/ -lyaml-cpp
YAML_CPP_INC = $(YAML_CPP_DIR)/include/

$(TARGET): $(OBJS)
	g++ -o $(TARGET) $(OBJS) -L$(YAML_CPP_LIB) -I$(YAML_CPP_INC)

tmp/$(OUT).o: src/$(OUT).cpp src/traceData.h
	g++ -c src/$(OUT).cpp -o tmp/$(OUT).o $(FLAGS) -I$(YAML_CPP_INC)

clean:
	rm -f $(TARGET) $(OBJS)

