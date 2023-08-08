# Compiler settings - Can be customized.
CC = g++
CXXFLAGS = -std=c++17 -Wall -O3 -I src
LDFLAGS = -static -static-libgcc -static-libstdc++

# Makefile settings - Can be customized.
APPNAME = quarto_solver
EXT = .cpp
SRCDIR = .
OBJDIR = obj
DEPDIR = dep

############## Do not change anything from here downwards! #############
rwildcard = $(wildcard $1$2) $(foreach d,$(wildcard $1*),$(call rwildcard,$d/,$2))
SRC = $(call rwildcard,$(SRCDIR)/,*$(EXT))
OBJ = $(SRC:$(SRCDIR)/%$(EXT)=$(OBJDIR)/%.o)
DEP = $(OBJ:$(OBJDIR)/%.o=$(DEPDIR)/%.d)
# UNIX-based OS variables & settings
RM = rm
DELOBJ = $(OBJ)
# Windows OS variables & settings
DEL = rm
EXE = .exe
WDELOBJ = $(SRC:$(SRCDIR)/%$(EXT)=$(OBJDIR)/%.o)

########################################################################
####################### Targets beginning here #########################
########################################################################

all: $(APPNAME)
	$(DEL) $(DEP)

# Builds the app
$(APPNAME): $(OBJ)
	@mkdir -p $(@D)
	$(CC) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# Creates the dependecy rules
$(DEPDIR)/%.d: $(SRCDIR)/%$(EXT)
	@mkdir -p $(@D)
	@$(CC) $(CXXFLAGS) $< -MM -MT $(@:$(DEPDIR)/%.d=$(OBJDIR)/%.o) >$@

# Includes all .h files
-include $(DEP)

# Building rule for .o files and its .c/.cpp in combination with all .h
$(OBJDIR)/%.o: $(SRCDIR)/%$(EXT)
	@mkdir -p $(@D)
	$(CC) $(CXXFLAGS) -o $@ -c $<
