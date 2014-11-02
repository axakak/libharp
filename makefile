FLAGS = -std=c++11 -stdlib=libc++ -Wall -g -Wno-c++11-extensions
BINDIR = bin
OBJDIR = tmp
#OBJS = $(addprefix $(OBJDIR)/, crbfTester.o traceDataTester.o traceData.o c-rbf.o)

all: dir $(addprefix bin/, traceDataTester crbfTester)

$(BINDIR)/traceDataTester: $(addprefix $(OBJDIR)/, traceDataTester.o traceData.o)
	g++ -o bin/traceDataTester $(addprefix $(OBJDIR)/, traceDataTester.o traceData.o) -lyaml-cpp

$(BINDIR)/crbfTester: $(addprefix $(OBJDIR)/, crbfTester.o c-rbf.o traceData.o)
	g++ -o bin/crbfTester $(addprefix $(OBJDIR)/, crbfTester.o c-rbf.o traceData.o) -lyaml-cpp

$(OBJDIR)/traceDataTester.o: $(addprefix src/, traceDataTester.cpp traceData.h)
	g++ -c src/traceDataTester.cpp -o $(OBJDIR)/traceDataTester.o $(FLAGS)

$(OBJDIR)/crbfTester.o: $(addprefix src/, crbfTester.cpp c-rbf.h)
	g++ -c src/crbfTester.cpp -o $(OBJDIR)/crbfTester.o $(FLAGS)

$(OBJDIR)/traceData.o: $(addprefix src/, traceData.cpp traceData.h)
	g++ -c src/traceData.cpp -o tmp/traceData.o $(FLAGS)

$(OBJDIR)/c-rbf.o: $(addprefix src/, c-rbf.cpp c-rbf.h traceData.h)
	g++ -c src/c-rbf.cpp -o tmp/c-rbf.o $(FLAGS)

dir: | $(OBJDIR) $(BINDIR)

$(OBJDIR):
	mkdir $(OBJDIR)

$(BINDIR):
	mkdir $(BINDIR)

clean:
	rm -f $(BINDIR)/* $(OBJDIR)/*

