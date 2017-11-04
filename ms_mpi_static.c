#define PNG_NO_SETJMP

#include <assert.h>
#include <png.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

double cummTime = 0;
double IOTime = 0;
double totalTime = 0;
struct timespec diff(struct timespec start, struct timespec end) {
    struct timespec temp;
    if ((end.tv_nsec-start.tv_nsec)<0) {
        temp.tv_sec = end.tv_sec-start.tv_sec-1;
        temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
    } else {
        temp.tv_sec = end.tv_sec-start.tv_sec;
        temp.tv_nsec = end.tv_nsec-start.tv_nsec;
    }
    return temp;
}

int myMPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm, MPI_Request *request){
    struct timespec start, end, temp;
    clock_gettime(CLOCK_MONOTONIC, &start);
    MPI_Isend(buf, count, datatype, dest, tag, comm, request);
    clock_gettime(CLOCK_MONOTONIC, &end);
    temp = diff(start, end);
    double time_used = temp.tv_sec + (double) temp.tv_nsec / 1000000000.0;
    cummTime += time_used;
}

int myMPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm){
    struct timespec start, end, temp;
    clock_gettime(CLOCK_MONOTONIC, &start);
    MPI_Send(buf, count, datatype, dest, tag, comm);
    clock_gettime(CLOCK_MONOTONIC, &end);
    temp = diff(start, end);
    double time_used = temp.tv_sec + (double) temp.tv_nsec / 1000000000.0;
    cummTime += time_used;
}

int myMPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Request *request){
    struct timespec start, end, temp;
    clock_gettime(CLOCK_MONOTONIC, &start);
    MPI_Irecv(buf, count, datatype, source, tag, comm, request);
    clock_gettime(CLOCK_MONOTONIC, &end);
    temp = diff(start, end);
    double time_used = temp.tv_sec + (double) temp.tv_nsec / 1000000000.0;
    cummTime += time_used;
}

int myMPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status){
    struct timespec start, end, temp;
    clock_gettime(CLOCK_MONOTONIC, &start);
    MPI_Recv(buf, count, datatype, source, tag, comm, status);
    clock_gettime(CLOCK_MONOTONIC, &end);
    temp = diff(start, end);
    double time_used = temp.tv_sec + (double) temp.tv_nsec / 1000000000.0;
    cummTime += time_used;
}

int myMPI_Wait(MPI_Request *request, MPI_Status *status){
    struct timespec start, end, temp;
    clock_gettime(CLOCK_MONOTONIC, &start);
    MPI_Wait(request, status);
    clock_gettime(CLOCK_MONOTONIC, &end);
    temp = diff(start, end);
    double time_used = temp.tv_sec + (double) temp.tv_nsec / 1000000000.0;
    cummTime += time_used;
}

int myMPI_Waitall(int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[]){
    struct timespec start, end, temp;
    clock_gettime(CLOCK_MONOTONIC, &start);
    MPI_Waitall(count, array_of_requests, array_of_statuses);
    clock_gettime(CLOCK_MONOTONIC, &end);
    temp = diff(start, end);
    double time_used = temp.tv_sec + (double) temp.tv_nsec / 1000000000.0;
    cummTime += time_used;
}

void write_png(const char* filename, const int width, const int height, const int* buffer) {
    struct timespec start, end, temp;
    clock_gettime(CLOCK_MONOTONIC, &start);
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
        for (int x = 0; x < width; ++x) {
            int p = buffer[(height - 1 - y) * width + x];
            row[x * 3] = ((p & 0xf) << 4);
        }
        png_write_row(png_ptr, row);
    }
    free(row);
    png_write_end(png_ptr, NULL);
    png_destroy_write_struct(&png_ptr, &info_ptr);
    fclose(fp);
    clock_gettime(CLOCK_MONOTONIC, &end);
    temp = diff(start, end);
    double time_used = temp.tv_sec + (double) temp.tv_nsec / 1000000000.0;
    IOTime += time_used;
}

MPI_Comm comm = MPI_COMM_WORLD;
int rank, n;

int cal(double x0, double y0){
    int repeats = 0;
    double x = 0;
    double y = 0;
    double length_squared = 0;
    while (repeats < 100000 && length_squared < 4){
        double temp = x * x - y * y + x0;
        y = 2 * x * y + y0;
        x = temp;
        length_squared = x * x + y * y;
        ++repeats;
    }
    return repeats;
}


int main(int argc, char** argv) {
    struct timespec start, end, temp;
    clock_gettime(CLOCK_MONOTONIC, &start);
    /* argument parsing */
    assert(argc == 9);
    int num_threads = strtol(argv[1], 0, 10);
    double left = strtod(argv[2], 0);
    double right = strtod(argv[3], 0);
    double lower = strtod(argv[4], 0);
    double upper = strtod(argv[5], 0);
    int width = strtol(argv[6], 0, 10);
    int height = strtol(argv[7], 0, 10);
    const char* filename = argv[8];

    // init MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &n);

    // calcalate part of eac rank
    int imgSize = width*height;
    int re = imgSize % n;
    int startPos[n], endPos[n]; 
    startPos[0] = 0;
    endPos[n-1] = imgSize;
    for(int i = 1; i < n; i++){
        startPos[i] = imgSize/n + (i < re ? 1 : 0) + startPos[i-1];
        endPos[i-1] = startPos[i];
    }

    //for(int i = 0; i < n; i++)
    //    printf("[%d, %d]\n", startPos[i], endPos[i]);
    /* allocate memory for image */
    //printf("[rank %d]: %d\n", rank, endPos[rank]-startPos[rank]);
    int* image = (int*)malloc(imgSize * sizeof(int));
    assert(image);

    MPI_Status status[n-1];
    MPI_Request request[n-1];
    if(rank == 0){ // master
        for(int i = 1; i < n; i++){
            myMPI_Irecv(image+startPos[i], endPos[i]-startPos[i], MPI_INT, i, MPI_ANY_TAG, comm, &request[i-1]);
        }
    }
    /* mandelbrot set */
    for(int k = startPos[rank]; k < endPos[rank]; k++){
        // infer x,y from i
        int j = k / width, i = k % width; 
        double y0 = j * ((upper - lower) / height) + lower;
        double x0 = i * ((right - left) / width) + left;
        int repeats = cal(x0,y0);
        image[j * width + i] = repeats;
    }

    if(rank != 0){ //slave
        myMPI_Send(image+startPos[rank], endPos[rank]-startPos[rank], MPI_INT, 0, 0, comm);
    }
    if(rank == 0){ // wait
        myMPI_Waitall(n-1, request, status); 
        /* draw and cleanup */
        write_png(filename, width, height, image);
    }
    free(image);
	MPI_Finalize();
    clock_gettime(CLOCK_MONOTONIC, &end);
    temp = diff(start, end);
    double time_used = temp.tv_sec + (double) temp.tv_nsec / 1000000000.0;
    totalTime += time_used;
    totalTime -= cummTime;
    totalTime -= IOTime;
    printf("%f %f %f\n", rank, totalTime, cummTime, IOTime);
}
