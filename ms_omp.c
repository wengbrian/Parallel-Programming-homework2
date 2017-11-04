#define PNG_NO_SETJMP

#include <assert.h>
#include <png.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <unistd.h>

double left;
double right;
double lower;
double upper;
int width;
int height;
int *image;

void cal_pixel(int i, int j){
    double y0 = j * ((upper - lower) / height) + lower;
    double x0 = i * ((right - left) / width) + left;
    int repeats = 0;
    double x = 0;
    double y = 0;
    double length_squared = 0;
    while (repeats < 100000 && length_squared < 4) {
        double temp = x * x - y * y + x0;
        y = 2 * x * y + y0;
        x = temp;
        length_squared = x * x + y * y;
        ++repeats;
    }
    image[j * width + i] = repeats;
}
void write_png(const char* filename, const int width, const int height, const int* buffer) {
    FILE* fp = fopen(filename, "wb");
    assert(fp);
    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    assert(png_ptr);
    png_infop info_ptr = png_create_info_struct(png_ptr);
    assert(info_ptr);
    png_init_io(png_ptr, fp);
    png_set_IHDR(png_ptr, info_ptr, width, height, 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
    png_write_info(png_ptr, info_ptr);
    size_t row_size = 3 * width * sizeof(png_byte);
    png_bytep row = (png_bytep)malloc(row_size);
    for (int y = 0; y < height; ++y) {
        memset(row, 0, row_size);
        //    printf("a");
        for (int x = 0; x < width; ++x) {
            while(buffer[(height - 1 - y) * width + x] == -1){
                // if pixel is not ok
                cal_pixel(x, y);
            }
            int p = buffer[(height - 1 - y) * width + x];
            row[x * 3] = ((p & 0xf) << 4);
            
        }
        png_write_row(png_ptr, row);
       // printf("write row %d\n", y);
    }
    free(row);
    png_write_end(png_ptr, NULL);
    png_destroy_write_struct(&png_ptr, &info_ptr);
    fclose(fp);
}

int main(int argc, char** argv) {
    /* argument parsing */
    assert(argc == 9);
    int num_threads = strtol(argv[1], 0, 10);
    left = strtod(argv[2], 0);
    right = strtod(argv[3], 0);
    lower = strtod(argv[4], 0);
    upper = strtod(argv[5], 0);
    width = strtol(argv[6], 0, 10);
    height = strtol(argv[7], 0, 10);
    const char* filename = argv[8];

    /* allocate memory for image */
    image = (int*)malloc(width * height * sizeof(int));
    memset(image, -1, sizeof(int) * width * height);
    //for(int i = 0; i < height*width; i++)
      //  printf("%d ", image[i]);
    //printf("%d\n", sizeof(int) * width * height);
    //fflush(stdout);
    assert(image);

    /* mandelbrot set */
  //  printf("start\n");
    #pragma omp parallel num_threads(num_threads) shared(image)
    {
        # pragma omp master
        {
    //        printf("%d write\n", omp_get_thread_num());
            write_png(filename, width, height, image);
        }
        
        # pragma omp for schedule(dynamic) nowait
        for (int j = height-1; j >= 0; --j) {
            for (int i = 0; i < width; ++i) {
                if(image[j * width + i] != -1)
                    continue;
                cal_pixel(i,j);
            }
          //  printf("%d: row %d is ok\n", omp_get_thread_num(), j);
        }
    }
    /* draw and cleanup */
    free(image);
}
