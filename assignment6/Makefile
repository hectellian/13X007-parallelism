CXX = g++
CXXFLAGS = -Wall -O2 -fopenmp
TARGET_EXEC1 := fractal_static
TARGET_EXEC2 := fractal_dynamic

BUILD_DIR := ./build
SRC_DIRS := ./src
INC_DIRS := ./include

# Find source files, exclude main_static.cpp and main_dynamic.cpp
SRCS := $(filter-out $(SRC_DIRS)/main_static.cpp $(SRC_DIRS)/main_dynamic.cpp, $(shell find $(SRC_DIRS) -name '*.cpp'))
OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)
DEPS := $(OBJS:.o=.d)

# Object files for the individual main files
MAIN_OBJ1 := $(BUILD_DIR)/$(SRC_DIRS)/main_static.cpp.o
MAIN_OBJ2 := $(BUILD_DIR)/$(SRC_DIRS)/main_dynamic.cpp.o

INC_FLAGS := $(addprefix -I,$(INC_DIRS))

CPPFLAGS := $(INC_FLAGS) -MMD -MP

# Default target to build both executables
all: $(BUILD_DIR)/$(TARGET_EXEC1) $(BUILD_DIR)/$(TARGET_EXEC2)

# Build rules for main_static and main_dynamic
$(BUILD_DIR)/$(TARGET_EXEC1): $(OBJS) $(MAIN_OBJ1)
	$(CXX) $(OBJS) $(MAIN_OBJ1) -o $@ $(CXXFLAGS) $(LDFLAGS)

$(BUILD_DIR)/$(TARGET_EXEC2): $(OBJS) $(MAIN_OBJ2)
	$(CXX) $(OBJS) $(MAIN_OBJ2) -o $@ $(CXXFLAGS) $(LDFLAGS)

# Pattern rule for .cpp file compilation
$(BUILD_DIR)/%.cpp.o: %.cpp
	mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -r $(BUILD_DIR)

-include $(DEPS)
