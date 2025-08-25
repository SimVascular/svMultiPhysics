#
# Copyright (c) Stanford University, The Regents of the University of
#               California, and others.
#
# All Rights Reserved.
#
# See Copyright-SimVascular.txt for additional details.
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject
# to the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
# IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
# OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#--------------------------------------------------------------------
#
#	This is the Makefile to build GenBC.
#
#--------------------------------------------------------------------

# HERE COMES THE DEFINITION

include Makefile.in

INCLUDES = -I./include

MYFUN = Modules.f \
        USER.f \
        GenBC.f

GenBC_EXE = genBC.exe

DSYM_DIR =
ifeq ($(debug),1)
   ifeq ($(OS),Darwin)
      DSYM_DIR = $(GenBC_EXE).dSYM
   endif
endif
#############################################################################
# AND HERE ARE THE RULES

SRC = $(patsubst %,$(SRC_DIR)/%,$(MYFUN))
OBJ = $(patsubst %.f,$(OBJ_DIR)/%.o,$(MYFUN))

#.PHONY: $(TEST)
#$(TEST): $(TEST:.exe=.f) $(GenBC_EXE)
#	$(FORTRAN) $< $(FFLAGS) $(INCLUDES)  -o $@

.PHONY: $(GenBC_EXE)
$(GenBC_EXE): $(OBJ)
	$(FORTRAN) $(OBJ) $(FFLAGS) -o $@

$(OBJ): | $(OBJ_DIR)

$(OBJ_DIR):
	mkdir -p $@

$(OBJ_DIR)/%.o : $(SRC_DIR)/%.f
	$(FORTRAN) $(FFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -r -f $(OBJ_DIR) $(GenBC_EXE) $(TEST) $(DSYM_DIR)
