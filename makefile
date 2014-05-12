
FLAGS = -std=c++11 -ansi -Wall -g -Wno-c++11-extensions
BINDIR = bin
OBJDIR = tmp
OBJS = $(addprefix $(OBJDIR)/, traceDataTester.o traceData.o)

bin/traceDataTester: $(OBJS)
	g++ -o bin/traceDataTester $(OBJS) -lyaml-cpp

$(OBJDIR)/traceDataTester.o: $(addprefix src/, traceDataTester.cpp traceData.h)
	g++ -c src/traceDataTester.cpp -o $(OBJDIR)/traceDataTester.o $(FLAGS)

$(OBJDIR)/traceData.o: $(addprefix src/, traceData.cpp traceData.h)
	g++ -c src/traceData.cpp -o tmp/traceData.o $(FLAGS)

$(OBJS): | $(OBJDIR) $(BINDIR)

$(OBJDIR):
	mkdir $(OBJDIR)

$(BINDIR):
	mkdir $(BINDIR)

clean:
	rm -f $(BINDIR)/* $(OBJDIR)/*

