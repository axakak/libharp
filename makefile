FLAGS = -std=c++11 -stdlib=libc++ -Wall -g -Wno-c++11-extensions
BINDIR = bin
OBJDIR = tmp

all: dir $(addprefix bin/, traceDataTester crbfTrainer)

$(BINDIR)/traceDataTester: $(addprefix $(OBJDIR)/, traceDataTester.o traceData.o)
	g++ -o bin/traceDataTester $(addprefix $(OBJDIR)/, traceDataTester.o traceData.o) -lyaml-cpp

$(BINDIR)/crbfTrainer: $(addprefix $(OBJDIR)/, crbfTrainer.o c-rbf.o traceData.o)
	g++ -o bin/crbfTrainer $(addprefix $(OBJDIR)/, crbfTrainer.o c-rbf.o traceData.o) -lyaml-cpp

$(OBJDIR)/traceDataTester.o: $(addprefix src/, traceDataTester.cpp traceData.h)
	g++ -c src/traceDataTester.cpp -o $(OBJDIR)/traceDataTester.o $(FLAGS)

$(OBJDIR)/crbfTrainer.o: $(addprefix src/, crbfTrainer.cpp c-rbf.h)
	g++ -c src/crbfTrainer.cpp -o $(OBJDIR)/crbfTrainer.o $(FLAGS)

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

