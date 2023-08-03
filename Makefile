################################################################################
#                                   LICENSE                                    #
################################################################################
#   This file is part of error_function.                                       #
#                                                                              #
#   error_function is free software: you can redistribute it and/or modify     #
#   it under the terms of the GNU General Public License as published by       #
#   the Free Software Foundation, either version 3 of the License, or          #
#   (at your option) any later version.                                        #
#                                                                              #
#   error_function is distributed in the hope that it will be useful,          #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of             #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              #
#   GNU General Public License for more details.                               #
#                                                                              #
#   You should have received a copy of the GNU General Public License          #
#   along with error_function.  If not, see <https://www.gnu.org/licenses/>.   #
################################################################################
TARGET_LIB=liberrorfunction.so
BUILD_DIR=./build
SRC_DIRS=./src
CC=cc
CFLAGS=-Wall -Wextra -Wpedantic $(EXTRA_FLAGS) -O3 -fPIC -flto -DNDEBUG -c
LFLAGS=-O3 -flto -shared -ltmpl
SRCS=$(shell find $(SRC_DIRS) -name "*.c" -and -not -name "setup.c")
OBJS=$(SRCS:%=$(BUILD_DIR)/%.o)

.PHONY: clean all

all: $(BUILD_DIR) $(TARGET_LIB)

$(TARGET_LIB): $(OBJS)
	$(CC) $(OBJS) $(LFLAGS) -o $@

$(BUILD_DIR)/%.c.o: %.c
	$(CC) $(CFLAGS) $< -o $@

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)/src/

clean:
	rm -rf $(BUILD_DIR)
	rm -f $(TARGET_LIB)
