CXX ?= g++
CXXFLAGS = -O2 -Wall -std=c++17 -pthread
PREFIX ?= /usr/local

TARGET = gffsub
SRCS = src/gffsub.cpp src/annotation_index.cpp src/gff3_parser.cpp src/region.cpp src/annotation_filter.cpp src/isoform_filter.cpp src/annotation_output.cpp
HDRS = src/gff3.hpp

.PHONY: all clean test install uninstall

all: $(TARGET)

$(TARGET): $(SRCS) $(HDRS)
	$(CXX) $(CXXFLAGS) -o $@ $(SRCS)

clean:
	rm -f $(TARGET) src/*.o

test: $(TARGET)
	./$(TARGET) ../SoyL04.gene.gff -r Chr01:1-35000 -f gene
	./$(TARGET) ../SoyL04.gene.gff --longest | head -10
	./$(TARGET) ../SoyL04.gene.gff -r Chr01:1-35000 -t gtf3 | head -10

install: $(TARGET)
	install -d $(PREFIX)/bin
	install -m 755 $(TARGET) $(PREFIX)/bin/

uninstall:
	rm -f $(PREFIX)/bin/$(TARGET)
