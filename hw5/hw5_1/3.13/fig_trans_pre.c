#define _POSIX_C_SOURCE 200809L
#define _FILE_OFFSET_BITS 64

#include <png.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void skip_ws_and_comments(FILE *fp) {
    int c;
    while ((c = fgetc(fp)) != EOF) {
        if (c == '#') {
            while ((c = fgetc(fp)) != '\n' && c != EOF) {
            }
        } else if (c == ' ' || c == '\n' || c == '\r' || c == '\t') {
            continue;
        } else {
            ungetc(c, fp);
            break;
        }
    }
}

static int read_ppm_header(FILE *fp, int *width, int *height, int *maxval) {
    char magic[3] = {0};

    if (fscanf(fp, "%2s", magic) != 1) return 0;
    if (strcmp(magic, "P6") != 0) return 0;

    skip_ws_and_comments(fp);
    if (fscanf(fp, "%d", width) != 1) return 0;

    skip_ws_and_comments(fp);
    if (fscanf(fp, "%d", height) != 1) return 0;

    skip_ws_and_comments(fp);
    if (fscanf(fp, "%d", maxval) != 1) return 0;

    if (*width <= 0 || *height <= 0 || *maxval != 255) return 0;

    /* 吃掉 header 后的一个空白字符 */
    fgetc(fp);
    return 1;
}

int main(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr,
                "Usage: %s input.ppm output.png scale\n"
                "Example: %s secant.ppm secant_preview.png 20\n",
                argv[0], argv[0]);
        return 1;
    }

    const char *input_path = argv[1];
    const char *output_path = argv[2];
    int scale = atoi(argv[3]);

    if (scale <= 0) {
        fprintf(stderr, "scale must be a positive integer\n");
        return 1;
    }

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

    int out_w = (width + scale - 1) / scale;
    int out_h = (height + scale - 1) / scale;

    fprintf(stderr, "Input : %dx%d\n", width, height);
    fprintf(stderr, "Output: %dx%d (scale=%d)\n", out_w, out_h, scale);

    FILE *fp_out = fopen(output_path, "wb");
    if (!fp_out) {
        perror("fopen output");
        fclose(fp_in);
        return 1;
    }

    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png_ptr) {
        fprintf(stderr, "png_create_write_struct failed\n");
        fclose(fp_in);
        fclose(fp_out);
        return 1;
    }

    png_infop info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr) {
        fprintf(stderr, "png_create_info_struct failed\n");
        png_destroy_write_struct(&png_ptr, NULL);
        fclose(fp_in);
        fclose(fp_out);
        return 1;
    }

    if (setjmp(png_jmpbuf(png_ptr))) {
        fprintf(stderr, "Error during PNG creation\n");
        png_destroy_write_struct(&png_ptr, &info_ptr);
        fclose(fp_in);
        fclose(fp_out);
        return 1;
    }

    png_init_io(png_ptr, fp_out);
    png_set_IHDR(png_ptr, info_ptr,
                 (png_uint_32)out_w,
                 (png_uint_32)out_h,
                 8,
                 PNG_COLOR_TYPE_RGB,
                 PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_DEFAULT,
                 PNG_FILTER_TYPE_DEFAULT);

    png_write_info(png_ptr, info_ptr);

    size_t in_row_bytes = (size_t)width * 3;
    unsigned char *in_row = (unsigned char *)malloc(in_row_bytes);
    unsigned char *out_row = (unsigned char *)malloc((size_t)out_w * 3);

    if (!in_row || !out_row) {
        fprintf(stderr, "Memory allocation failed\n");
        free(in_row);
        free(out_row);
        png_destroy_write_struct(&png_ptr, &info_ptr);
        fclose(fp_in);
        fclose(fp_out);
        return 1;
    }

    for (int y = 0, oy = 0; y < height; ++y) {
        size_t got = fread(in_row, 1, in_row_bytes, fp_in);
        if (got != in_row_bytes) {
            fprintf(stderr, "Failed to read row %d from PPM\n", y);
            free(in_row);
            free(out_row);
            png_destroy_write_struct(&png_ptr, &info_ptr);
            fclose(fp_in);
            fclose(fp_out);
            return 1;
        }

        if (y % scale != 0) continue;

        for (int ox = 0, x = 0; x < width; x += scale, ++ox) {
            out_row[3 * ox + 0] = in_row[3 * x + 0];
            out_row[3 * ox + 1] = in_row[3 * x + 1];
            out_row[3 * ox + 2] = in_row[3 * x + 2];
        }

        png_write_row(png_ptr, out_row);
        ++oy;
    }

    png_write_end(png_ptr, NULL);

    free(in_row);
    free(out_row);
    png_destroy_write_struct(&png_ptr, &info_ptr);
    fclose(fp_in);
    fclose(fp_out);

    fprintf(stderr, "Wrote preview PNG: %s\n", output_path);
    return 0;
}