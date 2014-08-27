
FLAGS = -std=c++11 -stdlib=libc++ -Wall -g -Wno-c++11-extensions
BINDIR = bin
OBJDIR = tmp
OBJS = $(addprefix $(OBJDIR)/, traceDataTester.o traceData.o c-rbf.o)

bin/traceDataTester: $(OBJS)
	g++ -o bin/traceDataTester $(OBJS) -lyaml-cpp

$(OBJDIR)/traceDataTester.o: $(addprefix src/, traceDataTester.cpp traceData.h)
	g++ -c src/traceDataTester.cpp -o $(OBJDIR)/traceDataTester.o $(FLAGS)

$(OBJDIR)/traceData.o: $(addprefix src/, traceData.cpp traceData.h)
	g++ -c src/traceData.cpp -o tmp/traceData.o $(FLAGS)

$(OBJDIR)/c-rbf.o: $(addprefix src/, c-rbf.cpp c-rbf.h)
	g++ -c src/c-rbf.cpp -o tmp/c-rbf.o $(FLAGS)

$(OBJS): | $(OBJDIR) $(BINDIR)

$(OBJDIR):
	mkdir $(OBJDIR)

$(BINDIR):
	mkdir $(BINDIR)

clean:
	rm -f $(BINDIR)/* $(OBJDIR)/*

