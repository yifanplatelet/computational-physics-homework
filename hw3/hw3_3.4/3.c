#include <stdio.h>

int main() {
    unsigned int x = 1;
    char *p = (char *)&x;

    if (*p == 1) {
        printf("当前系统是小端序（Little Endian）\n");
    } else {
        printf("当前系统是大端序（Big Endian）\n");
    }

    return 0;
}