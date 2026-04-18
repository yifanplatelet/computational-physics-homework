#include <stdio.h>
#include <float.h>

int main() {
    printf("FLT_EPSILON  = %.20e\n", FLT_EPSILON);
    printf("DBL_EPSILON  = %.20e\n", DBL_EPSILON);
    printf("LDBL_EPSILON = %.20Le\n", LDBL_EPSILON);
    return 0;
}