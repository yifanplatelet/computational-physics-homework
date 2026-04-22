#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <png.h>

static void skip_comments(FILE *fp) {
    int c;
    while ((c = fgetc(fp)) == '#') {
        while ((c = fgetc(fp)) != '\n' && c != EOF) {
        }
    }
    if (c != EOF) {
        ungetc(c, fp);
    }
}

static int read_ppm_header(FILE *fp, int *width, int *height, int *maxval) {
    char magic[3] = {0};

    if (fscanf(fp, "%2s", magic) != 1) return 0;
    if (strcmp(magic, "P6") != 0) return 0;

    skip_comments(fp);
    if (fscanf(fp, "%d", width) != 1) return 0;

    skip_comments(fp);
    if (fscanf(fp, "%d", height) != 1) return 0;

    skip_comments(fp);
    if (fscanf(fp, "%d", maxval) != 1) return 0;

    if (*width <= 0 || *height <= 0) return 0;
    if (*maxval != 255) return 0;

    fgetc(fp); /* 吃掉 header 后面的单个空白字符 */
    return 1;
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s input.ppm output.png\n", argv[0]);
        return 1;
    }

    const char *input_path = argv[1];
    const char *output_path = argv[2];

    FILE *fp_in = fopen(input_path, "rb");
    if (!fp_in) {
        perror("fopen input");
        return 1;
    }

    int width, height, maxval;
    if (!read_ppm_header(fp_in, &width, &height, &maxval)) {
        fprintf(stderr, "Invalid or unsupported PPM file: %s\n", input_path);
        fclose(fp_in);
        return 1;
    }

    size_t row_bytes = (size_t)width * 3;
    size_t image_bytes = row_bytes * (size_t)height;

    unsigned char *image = (unsigned char *)malloc(image_bytes);
    if (!image) {
        fprintf(stderr, "Failed to allocate memory for image\n");
        fclose(fp_in);
        return 1;
    }

    size_t read_bytes = fread(image, 1, image_bytes, fp_in);
    fclose(fp_in);

    if (read_bytes != image_bytes) {
        fprintf(stderr, "Failed to read complete PPM pixel data\n");
        free(image);
        return 1;
    }

    FILE *fp_out = fopen(output_path, "wb");
    if (!fp_out) {
        perror("fopen output");
        free(image);
        return 1;
    }

    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png_ptr) {
        fprintf(stderr, "png_create_write_struct failed\n");
        fclose(fp_out);
        free(image);
        return 1;
    }

    png_infop info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr) {
        fprintf(stderr, "png_create_info_struct failed\n");
        png_destroy_write_struct(&png_ptr, NULL);
        fclose(fp_out);
        free(image);
        return 1;
    }

    if (setjmp(png_jmpbuf(png_ptr))) {
        fprintf(stderr, "Error during PNG creation\n");
        png_destroy_write_struct(&png_ptr, &info_ptr);
        fclose(fp_out);
        free(image);
        return 1;
    }

    png_init_io(png_ptr, fp_out);

    png_set_IHDR(
        png_ptr,
        info_ptr,
        (png_uint_32)width,
        (png_uint_32)height,
        8,
        PNG_COLOR_TYPE_RGB,
        PNG_INTERLACE_NONE,
        PNG_COMPRESSION_TYPE_DEFAULT,
        PNG_FILTER_TYPE_DEFAULT
    );

    png_write_info(png_ptr, info_ptr);

    png_bytep *rows = (png_bytep *)malloc((size_t)height * sizeof(png_bytep));
    if (!rows) {
        fprintf(stderr, "Failed to allocate row pointers\n");
        png_destroy_write_struct(&png_ptr, &info_ptr);
        fclose(fp_out);
        free(image);
        return 1;
    }

    for (int y = 0; y < height; ++y) {
        rows[y] = image + (size_t)y * row_bytes;
    }

    png_write_image(png_ptr, rows);
    png_write_end(png_ptr, NULL);

    free(rows);
    png_destroy_write_struct(&png_ptr, &info_ptr);
    fclose(fp_out);
    free(image);

    fprintf(stderr, "Converted: %s -> %s\n", input_path, output_path);
    return 0;
}