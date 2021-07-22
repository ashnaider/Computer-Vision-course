src_dir = ../src

TARGET = test
DEBUG=0
CFLAGS = -c -DDEBUG=$(DEBUG) -MD $(addprefix -I,$(src_dir))

source_dirs = $(src_dir)

src_search_wildcards = $(addsuffix /*.c,$(src_dir))

all: $(TARGET)

$(TARGET): $(notdir $(patsubst %.c,%.o,$(wildcard $(src_search_wildcards))))
	$(CC) $^ -o $@

vpath %.h $(src_dir)
vpath %.c $(src_dir)

%.o: %.c
	$(CC) $(CFLAGS) $<

include $(wildcard *.d)


.PHONY: clean

clean:
	rm $(TARGET) *.d *.o
