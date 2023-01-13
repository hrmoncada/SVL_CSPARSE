# >> gcc  -Wall   -o      out    C_Code_SVL_Mar_11_2015_CSparse.c -lcsparse -lfftw3 -lm 
#   (CC) (CFLAGS) -o (EXECUTABLE)           (SOURCE)                   (LFLAGS)
#
# Compiler (CC): gcc for C program,
# Compiler flags (CFLAGS):
#      -g    adds debugging information to the executable file
#      -Wall turns on most, but not all, compiler warnings
# Define any libraries to link into executable (LFLAGS):
# Define the C source files (SOURCE):
#       SOURCE = myprog
# Build OUTPUT executable (EXECUTABLE):
#       EXECUTABLE = output name

CC         := gcc  
CFLAGS     := -Wall
LFLAGS     := -lcsparse -lfftw3 -lm
TRASH	   :=  *.*~  *~  data* *.o
SOURCE     := SVL_main.c SVL_unit_cell.c SVL_creat_matrix.c SVL_cartesian_to_polar.c SVL_polar_to_cartesian.c SVL_print_1d_to_2d_array.c SVL_print_2d_to_2d_array.c SVL_orientation_function.c SVL_Swap_quadrants.c  SVL_fftw.c SVL_ifftw.c SVL_creat_matrix_complex.c SVL_truncate_fftw_spatial_harmonic.c SVL_grading_vector.c SVL_fdder.c SVL_loop.c 

OBJECTS    := $(SOURCE:.c=.o)   #C_Code_SVL_Mar_11_2015_CSparse.o
EXECUTABLE := out

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
		$(CC) $(CFLAGS) $(OBJECTS) -o $@ $(LFLAGS)
clean:
	rm -rf $(EXECUTABLE) $(OBJECTS) $(TRASH)
