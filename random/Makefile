APPLICATION = random_stream

SRC_DIR := ./src
CODELETS_SRC_DIR := ./src/codelets
OBJ_DIR := ./build


CXXFLAGS ?= -std=c++17 -Wall -Wextra -MMD -g $(EXTRA_FLAGS)
LDFLAGS ?= -lpoplar -lgcl_ct -lpvti
POPC ?= popc
POPXXFLAGS_LIB=-O2
POPXXFLAGS ?= $(POPXXFLAGS_LIB) -target=ipu2 -DNDEBUG $(EXTRA_FLAGS)

# Files that generate .gp (target IPU).
CODELET_FILES = $(wildcard $(CODELETS_SRC_DIR)/*.cpp)
CODELET_TARGETS = $(patsubst $(CODELETS_SRC_DIR)/%.cpp,%.gp,$(CODELET_FILES))

# Files that generate .o (target host).
SRC_FILES = $(wildcard $(SRC_DIR)/*.cpp)
SRC_FILES += $(CODELET_FILES)
OBJ_FILES = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(patsubst $(SRC_DIR)/%.c,$(OBJ_DIR)/%.o,$(SRC_FILES)))

ALL = $(APPLICATION) $(CODELET_TARGETS)

all: $(ALL)

run: all
	./$(APPLICATION) $(ARGS)

info:
	@echo "SRC_FILES"
	@echo $(SRC_FILES)
	@echo "\nOBJ_FILES"
	@echo $(OBJ_FILES)
	@echo "\nCODELET_FILES"
	@echo $(CODELET_FILES)
	@echo "\nCODELET_TARGETS"
	@echo $(CODELET_TARGETS)
	@echo "\nAPPLICATION"
	@echo $(APPLICATION)

$(APPLICATION): $(OBJ_FILES)
	@echo $(OBJ_FILES)
	$(CXX) -o $@ $^ $(LDFLAGS)

$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp
	mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(OBJ_DIR)/%.o : $(SRC_DIR)/%.c
	mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

%.gp: $(CODELETS_SRC_DIR)/%.cpp
	$(POPC) -I $(CODELETS_SRC_DIR) $(POPXXFLAGS) $+ -o $@

.PHONY: clean

SUBDIRS = . codelets
clean:
	$(RM) $(ALL)
	$(foreach SD,$(SUBDIRS),$(RM) $(OBJ_DIR)/$(SD)/*.o;$(RM) $(OBJ_DIR)/$(SD)/*.d;)

include .makefileok
.makefileok: Makefile
	@touch $@
	$(MAKE) clean

-include $(OBJ_FILES:.o=.d)
