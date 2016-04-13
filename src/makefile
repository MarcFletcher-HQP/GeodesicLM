R = R CMD SHLIB

RFLAGS = --clean

CSRCS = $(wildcard *.c)
FSRCS = $(wildcard *.f)
F95SRCS = $(wildcard *.f90)

TARGET = geodesicLM.dll
all: $(TARGET)

$(TARGET): $(CSRCS) $(FSRCS) $(F95SRCS)
	$(R) $(RFLAGS) $^ -o $@
