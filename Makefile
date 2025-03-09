CC = gcc
CFLAGS = -Wall -O3

TARGET = bin/main

SRC_DIR = src
INCLUDE_DIR = include
BUILD_DIR = build

SRCS = $(shell find $(SRC_DIR) -name '*.c') # Automatically find all .c files in src and its subdirectories
OBJS = $(SRCS:$(SRC_DIR)/%.c=$(BUILD_DIR)/%.o) # Convert .c files to .o files in build directory

HEADER_DIRS = $(shell find $(SRC_DIR) $(INCLUDE_DIR) -type d -print)
INCLUDES = $(addprefix -I, $(HEADER_DIRS))

all: $(TARGET)

# ===== Main target rule
$(TARGET): $(OBJS)
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(TARGET) $(OBJS) -L$(BUILD_DIR)

# Generic rule for .c to .o
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@ $(DEFS)

clean:
	$(RM) -r bin
	$(RM) -r $(BUILD_DIR) $(TARGET)

