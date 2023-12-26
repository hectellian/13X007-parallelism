# Common settings
CXX = g++ # Change to your preferred C++ compiler
CXXFLAGS = -Wall -O2 -std=c++17
TARGET_EXEC := heat

BUILD_DIR := ./build
SRC_DIRS := ./src
INC_DIRS := ./include

SRCS := $(shell find $(SRC_DIRS) -name '*.cpp')
OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)
DEPS := $(OBJS:.o=.d)

# LDFLAGS := -ltbb
INC_FLAGS := $(addprefix -I,$(INC_DIRS))
CPPFLAGS := $(INC_FLAGS) -MMD -MP

# CPU compilation settings
CPU_CXX = $(CXX)
CPU_CXXFLAGS = $(CXXFLAGS)

# GPU compilation settings (Adjust as per your GPU compiler and flags)
GPU_CXX = nvcc
GPU_CXXFLAGS = --compiler-options "$(CXXFLAGS)"

# Default rule
all: cpu

# CPU compilation
cpu: $(BUILD_DIR)/$(TARGET_EXEC)_cpu

$(BUILD_DIR)/$(TARGET_EXEC)_cpu: $(OBJS)
	$(CPU_CXX) $(OBJS) -o $@ $(LDFLAGS)

# GPU compilation
gpu: $(BUILD_DIR)/$(TARGET_EXEC)_gpu

$(BUILD_DIR)/$(TARGET_EXEC)_gpu: $(OBJS)
	$(GPU_CXX) $(OBJS) -o $@ $(LDFLAGS)

# Object file rules
$(BUILD_DIR)/%.c.o: %.c
	mkdir -p $(dir $@)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/%.cpp.o: %.cpp
	mkdir -p $(dir $@)
	$(CPU_CXX) $(CPPFLAGS) $(CPU_CXXFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -r $(BUILD_DIR)

-include $(DEPS)
