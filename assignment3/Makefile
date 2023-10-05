#*****************************************************************              
# File:   Makefile                                                              
# Author:
# Author:
# Date:   
# Coms:       
#*****************************************************************              

#---------------------------------------------------------
# macros
RM = /bin/rm -rf
CC = g++
#---------------------------------------------------------
# vars
EXEC = main

SRC_DIR = src
BUILD_DIR = build
TARGET = bin

SOURCES = ${EXEC}.o
OBJECTS = $(patsubst %.o,$(BUILD_DIR)/%.o,$(SOURCES))

#---------------------------------------------------------
# flags
CPPFLAGS = -std=gnu++17 -O3 -fopenmp
#CPPFLAGS = -std=gnu++17 -O0 -g -fsanitize=address -DDEBUG

all: $(EXEC)
#---------------------------------------------------------
$(EXEC): $(OBJECTS)
	@mkdir -p $(TARGET)
	${CC} $^ -o $(TARGET)/$@ ${CPPFLAGS}
#---------------------------------------------------------
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(BUILD_DIR)
	$(CC) -c $(CPPFLAGS) $< -o $@
#---------------------------------------------------------
process: ${SRC_DIR}/process.cpp
	@mkdir -p ${TARGET}
	${CC} $^ -o $(TARGET)/$@ ${CPPFLAGS}
#---------------------------------------------------------
.PHONY: clean
#---------------------------------------------------------
clean:
	$(RM) $(TARGET) $(BUILD_DIR) ${GEN_DOCS}
