# Compiler and flags
CXX      := g++
CXXFLAGS := -std=c++17 -Wall -Iinclude -I/usr/include/eigen3/ -O3

# Directories
SRC_DIR  := src
BUILD_DIR := build
TARGET  := main

# Sources and objects
SRCS    := $(wildcard $(SRC_DIR)/*.cpp)
OBJS    := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SRCS))

# Default target
all: $(TARGET)

# Link object files into executable
$(TARGET): $(OBJS)
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Compile each .cpp into .o
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean build artifacts
clean:
	rm -rf $(BUILD_DIR) $(TARGET)

# Rebuild everything from scratch
rebuild: clean all

.PHONY: all clean rebuild

