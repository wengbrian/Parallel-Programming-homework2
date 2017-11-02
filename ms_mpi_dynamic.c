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
int num_threads;
double left;
double right;
double lower;
double upper;
int width;
int height;
char* filename;

void ms(int *image, int startPos, int endPos){
    /* mandelbrot set */
    for(int k = startPos; k < endPos; k++){
        // infer x,y from i
        int j = k / height, i = k % height; 
        double y0 = j * ((upper - lower) / height) + lower;
        double x0 = i * ((right - left) / width) + left;
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
        image[j * width + i] = repeats;
    }
}

int calAvg(int current, int imgSize, int lastAvg){
    // use this function to make split more dynamic
    int re = imgSize - current;
    if(re < lastAvg * n){
        if(re / n < 1)
            return 1;
        return re / n;
    }
    return lastAvg;
}

int procIfReqComp(MPI_Request *requests, int *image, int current, int lastAvg){
    int index;
    int flag;
    MPI_Status status;
    MPI_Testany(n-1, requests, &index, &flag, &status);
    int bound[2] = {current, current+lastAvg};
    if(index != MPI_UNDEFINED){
        // Send to who has complete
        myMPI_Send(bound, 2, MPI_INT, status.MPI_SOURCE, 1, comm);
        printf("master send [%d, %d] to slave%d\n", bound[0], bound[1], status.MPI_SOURCE);
        myMPI_Irecv(image+current, lastAvg, MPI_INT, status.MPI_SOURCE, MPI_ANY_TAG, comm, &requests[status.MPI_SOURCE-1]);
        printf("master assume receive [%d, %d] from slave%d\n", bound[0], bound[1], status.MPI_SOURCE);
        return 1;
    }
    return 0;
}
void master(int *image, int startPos[], int endPos[], int imgSize){
    MPI_Request requests[n-1];
    // init receive
    for(int i = 1; i < n; i++){
        myMPI_Irecv(image+startPos[i], endPos[i]-startPos[i], MPI_INT, i, MPI_ANY_TAG, comm, &requests[i-1]);
    }
    /* mandelbrot set */
    int lastAvg = imgSize * 0.5 * 0.5 / n;
    int counter = 0;
    for(int k = endPos[n-1]; k < imgSize; k++){
        // infer x,y from i
        int j = k / height, i = k % height; 
        double y0 = j * ((upper - lower) / height) + lower;
        double x0 = i * ((right - left) / width) + left;
        int repeats = 0;
        double x = 0;
        double y = 0;
        double length_squared = 0;
        while (repeats < 100000 && length_squared < 4){
            if((n > 1) && (repeats%100==0) && (k+1 < imgSize)){
                // for every 100 timestamp, check if any request is complete
                int result = procIfReqComp(requests, image, k+1, lastAvg);
                if(result > 0){
                    // increase k
                    k += lastAvg;
                    lastAvg = calAvg(k+1, imgSize, lastAvg);
                }
            }

            double temp = x * x - y * y + x0;
            y = 2 * x * y + y0;
            x = temp;
            length_squared = x * x + y * y;
            ++repeats;
        }
        image[j * width + i] = repeats;
    }
    MPI_Status statuses[n-1];
    MPI_Waitall(n-1, requests, statuses);
    for(int i = 1; i < n; i++){
        myMPI_Isend(image, 1, MPI_INT, i, 0, comm, requests);
    }
    MPI_Waitall(n-1, requests, statuses);
    // write file
    /* draw and cleanup */
    write_png(filename, width, height, image);
}

void slave(int *image, int startPos[], int endPos[]){
    ms(image, startPos[rank], endPos[rank]);
    // init send
    myMPI_Send(image+startPos[rank], endPos[rank]-startPos[rank], MPI_INT, 0, 0, comm);
    int bound[2];
    MPI_Status status;
    while(1){
        myMPI_Recv(bound, 2, MPI_INT, 0, MPI_ANY_TAG, comm, &status);
        printf("slave%d receive [%d, %d], tag: %d\n",rank, bound[0], bound[1], status.MPI_TAG);
        if(status.MPI_TAG == 0)
            break;
        ms(image, bound[0], bound[1]);
        myMPI_Send(image+bound[0], bound[1]-bound[0], MPI_INT, 0, 0, comm);
        printf("slave%d send [%d, %d]\n", rank, bound[0], bound[1]);
    }
}

int main(int argc, char** argv) {
    // calculate start time
    struct timespec start, end, temp;
    clock_gettime(CLOCK_MONOTONIC, &start);

    /* argument parsing */
    assert(argc == 9);
    num_threads = strtol(argv[1], 0, 10);
    left = strtod(argv[2], 0);
    right = strtod(argv[3], 0);
    lower = strtod(argv[4], 0);
    upper = strtod(argv[5], 0);
    width = strtol(argv[6], 0, 10);
    height = strtol(argv[7], 0, 10);
    filename = argv[8];

    // init MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(comm, &rank);
	MPI_Comm_size(comm, &n);

    int imgSize = width*height;
    int estiSize = imgSize*0.5;

    int startPos[n];
    int endPos[n];
    if(n > 1){
        startPos[1] = 0;
        endPos[1] = estiSize/n;
    }
    for(int i = 2; i < n; i++){
        startPos[i] = estiSize/n + startPos[i-1];
        endPos[i] = startPos[i] + estiSize/n;
    }

    int *image;
    if(rank == 0){ // master
        image = (int*)malloc(imgSize * sizeof(int));
        assert(image);
        master(image, startPos, endPos, imgSize);
    }else{ // slave
        image = (int*)malloc(estiSize/(n-1) * sizeof(int));
        assert(image);
        slave(image, startPos, endPos);
    }

    free(image);
    printf("rank%d finalize\n", rank);
	MPI_Finalize();

    // calculate end time
    clock_gettime(CLOCK_MONOTONIC, &end);
    temp = diff(start, end);
    double time_used = temp.tv_sec + (double) temp.tv_nsec / 1000000000.0;
    totalTime += time_used;
    totalTime -= cummTime;
    totalTime -= IOTime;
    printf("%f %f %f\n", rank, totalTime, cummTime, IOTime);
}
