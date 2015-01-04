CXXFLAGS = -std=c++11 -stdlib=libc++ -Wall -g -Wno-c++11-extensions
LDLIBS = -lyaml-cpp
BINDIR = bin
OBJDIR = tmp
SRCDIR = src
SUBDIRS =  $(BINDIR) $(OBJDIR)

all: $(SUBDIRS) $(addprefix $(BINDIR)/, traceDataTester crbfTrainer)

$(BINDIR)/traceDataTester: $(addprefix $(OBJDIR)/, traceDataTester.o traceData.o)
	$(CXX) -o $@ $^ $(LDLIBS)

$(BINDIR)/crbfTrainer: $(addprefix $(OBJDIR)/, crbfTrainer.o c-rbf.o traceData.o)
	$(CXX) -o $@ $^ $(LDLIBS)

$(OBJDIR)/traceDataTester.o: $(addprefix $(SRCDIR)/, traceDataTester.cpp traceData.h)
	$(CXX) -c $< -o $@ $(CXXFLAGS)

$(OBJDIR)/crbfTrainer.o: $(addprefix $(SRCDIR)/, crbfTrainer.cpp c-rbf.h)
	$(CXX) -c $< -o $@ $(CXXFLAGS)

$(OBJDIR)/traceData.o: $(addprefix $(SRCDIR)/, traceData.cpp traceData.h)
	$(CXX) -c $< -o $@ $(CXXFLAGS)

$(OBJDIR)/c-rbf.o: $(addprefix $(SRCDIR)/, c-rbf.cpp c-rbf.h traceData.h)
	$(CXX) -c $< -o $@ $(CXXFLAGS)

$(SUBDIRS):
	mkdir $@

clean:
	rm -rfv $(BINDIR) $(OBJDIR)

