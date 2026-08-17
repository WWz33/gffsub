CXX ?= g++
CXXFLAGS = -O2 -Wall -std=c++17 -pthread
AR ?= ar
PREFIX ?= /usr/local

TARGET = gffsub
LIB = libgffsub_core.a

LIB_SRCS = src/feature_types.cpp \
       src/string_utils.cpp \
       src/selector_filter.cpp \
       src/gtf_parser.cpp \
       src/query_summary.cpp \
       src/attributes.cpp \
       src/annotation_index.cpp \
       src/gff3_parser.cpp \
       src/region.cpp \
       src/annotation_filter.cpp \
       src/isoform_filter.cpp \
       src/annotation_output.cpp \
       src/query.cpp \
       src/window.cpp \
       src/subset.cpp

CLI_SRCS = src/cli.cpp \
       src/cli_usage.cpp \
       src/gffsub.cpp

SRCS = $(LIB_SRCS) $(CLI_SRCS)

HDRS = src/annotation.hpp \
       src/cli.hpp \
       src/cli_usage.hpp \
       src/expr_parser.hpp \
       src/feature_types.hpp \
       src/filter.hpp \
       src/gff3.hpp \
       src/gtf_parser.hpp \
       src/output.hpp \
       src/parser.hpp \
       src/query.hpp \
       src/query_summary.hpp \
       src/record.hpp \
       src/region.hpp \
       src/selector_filter.hpp \
       src/string_utils.hpp \
       src/subset.hpp \
       src/version.hpp \
       src/window.hpp

LIB_OBJS = $(LIB_SRCS:.cpp=.o)
CLI_OBJS = $(CLI_SRCS:.cpp=.o)
OBJS = $(LIB_OBJS) $(CLI_OBJS)

ANNOTATION_INDEX_SMOKE = annotation_index_smoke
CLI_OUTPUT_ATTRS_SMOKE = cli_output_attrs_smoke
CLI_SELECTOR_SMOKE = cli_selector_smoke
REGRESSION_SMOKE = regression_smoke
FEATURE_TYPES_SMOKE = feature_types_smoke

.PHONY: all clean test install uninstall

all: $(TARGET)

$(LIB): $(LIB_OBJS)
	$(AR) rcs $@ $(LIB_OBJS)

$(TARGET): $(CLI_OBJS) $(LIB)
	$(CXX) $(CXXFLAGS) -o $@ $(CLI_OBJS) $(LIB)

%.o: %.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(ANNOTATION_INDEX_SMOKE): tests/annotation_index_smoke.cpp src/feature_types.cpp src/string_utils.cpp src/attributes.cpp src/annotation_index.cpp src/gff3_parser.cpp src/gtf_parser.cpp src/region.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -Isrc -o $@ tests/annotation_index_smoke.cpp src/feature_types.cpp src/string_utils.cpp src/attributes.cpp src/annotation_index.cpp src/gff3_parser.cpp src/gtf_parser.cpp src/region.cpp

$(CLI_OUTPUT_ATTRS_SMOKE): tests/cli_output_attrs_smoke.cpp tests/test_utils.hpp
	$(CXX) $(CXXFLAGS) -o $@ tests/cli_output_attrs_smoke.cpp

$(CLI_SELECTOR_SMOKE): tests/cli_selector_smoke.cpp tests/test_utils.hpp
	$(CXX) $(CXXFLAGS) -o $@ tests/cli_selector_smoke.cpp

$(REGRESSION_SMOKE): tests/regression_smoke.cpp tests/test_utils.hpp
	$(CXX) $(CXXFLAGS) -o $@ tests/regression_smoke.cpp

$(FEATURE_TYPES_SMOKE): tests/feature_types_smoke.cpp src/feature_types.cpp src/string_utils.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -Isrc -o $@ tests/feature_types_smoke.cpp src/feature_types.cpp src/string_utils.cpp

clean:
	rm -f $(TARGET) $(LIB) $(OBJS) $(ANNOTATION_INDEX_SMOKE) $(CLI_OUTPUT_ATTRS_SMOKE) $(CLI_SELECTOR_SMOKE) $(REGRESSION_SMOKE) $(FEATURE_TYPES_SMOKE)

test: $(TARGET) $(ANNOTATION_INDEX_SMOKE) $(CLI_OUTPUT_ATTRS_SMOKE) $(CLI_SELECTOR_SMOKE) $(REGRESSION_SMOKE) $(FEATURE_TYPES_SMOKE)
	./$(ANNOTATION_INDEX_SMOKE)
	./$(CLI_OUTPUT_ATTRS_SMOKE) ./$(TARGET)
	./$(CLI_SELECTOR_SMOKE) ./$(TARGET)
	./$(REGRESSION_SMOKE) ./$(TARGET)
	./$(FEATURE_TYPES_SMOKE)

install: $(TARGET)
	install -d $(PREFIX)/bin
	install -m 755 $(TARGET) $(PREFIX)/bin/

uninstall:
	rm -f $(PREFIX)/bin/$(TARGET)
