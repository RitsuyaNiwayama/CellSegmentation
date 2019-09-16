PROG := test.out
SRCS := levelset_memseg030414sf.cpp tiffio_030414sf.cpp levelset_innercell3d_func030414sf.cpp MPIio030414sf.cpp energy3cRN.cpp llistRN.cpp lsops3cRN.cpp sfm_local_chanvese_mexRN.cpp utility_sf.cpp parameters-tools_sf.cpp
OBJS := $(SRCS:%.cpp=%.o)
DEPS := $(SRCS:%.cpp=%.d)

CXX := g++
MPICC := mpiCC
INC  =  -I/usr/include/
LIB  = -L/usr/local/lib/ -L/usr/lib64/  -ltiff
all: $(PROG)
-include $(DEPS)
$(PROG): $(OBJS)
	$(MPICC)  -o $@ $^ $(LIB) -O3
%.o: %.cpp
	$(MPICC) -c -MMD -MP $< $(INC) -O3
clean:
	rm -f $(PROG) $(OBJS) $(DEPS)
