#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <sys/time.h>

#include "graphics.h"
#include "file_operations.h"

int calcAccel(const int N, const double G, double* mass, double* xPos, double* yPos, double* xAccelNew, double* yAccelNew, int nThreads){
    /* Fills up xAccelNew and yAccelNew with magnitude and direction of the acceleration */

    const double epsilon = 0.001;
    double magnitude, xPartialAccel, yPartialAccel;
    int i, j;

    #pragma omp parallel for private(i, j, magnitude, xPartialAccel, yPartialAccel) num_threads(nThreads)
    for (int ij=0; ij<N*N; ij++){

        i = floor(ij/N);
        j = ij%N;

        if (i!=j){
            magnitude     = sqrt(pow(xPos[i] - xPos[j], 2) + pow(yPos[i] - yPos[j], 2)) + epsilon;
            xPartialAccel = -G*mass[j]*(xPos[i] - xPos[j])/pow(magnitude, 3);
            yPartialAccel = -G*mass[j]*(yPos[i] - yPos[j])/pow(magnitude, 3);

            #pragma omp critical
            {
                xAccelNew[i] += xPartialAccel;
                yAccelNew[i] += yPartialAccel;
            }
        }
    }
    return 0;
}

int main(int argc, char *argv[]){
    if(argc!=7){
        printf("ArgumentError: 6 input arguments are expected but %i received.\n", argc-1);
        return -1;
    }

    const int N      = atoi(argv[1]);
    char filename[50];
    strcpy(filename, argv[2]);
    int nsteps       = atoi(argv[3]);
    double delta_t   = atof(argv[4]);
    int withGraphics = atoi(argv[5]);
    int nThreads     = atoi(argv[6]);
    const double G   = 100.0/N;
    int nData        = N*6;
    double* rawData  = (double*)malloc(nData*sizeof(double));
    struct timeval begin, end, beginGlob, endGlob;
    double elapsed = 0.0, elapsedGlob = 0.0;
    
    gettimeofday(&beginGlob, NULL);

    /* reading initial conditions from a file*/
    int isRead = read_doubles_from_file(nData, rawData, filename);

    /* checking whether data has been imported correctly or not */
    if (isRead != 0){
        return -1;
    }
    else {
        printf("File \"%s\" imported successfully.\n", filename);
        printf("nPlanets: %i\n", N);
        printf("nsteps: %i\n", nsteps);
        printf("timestep: %f\n", delta_t);
        printf("withGraphics: %i\n\n", withGraphics);
    }

    /*  xPos[i] is the latest x position of planet i. */
    double* xPos   = (double*)malloc(N*sizeof(double));
    double* yPos   = (double*)malloc(N*sizeof(double));
    double* xSpeed = (double*)malloc(N*sizeof(double));
    double* ySpeed = (double*)malloc(N*sizeof(double));
    double* mass   = (double*)malloc(N*sizeof(double));
    double* bright = (double*)malloc(N*sizeof(double));

    /* inserting initial data into new arrays */
    #pragma omp parallel for num_threads(nThreads)
    for (int i=0; i<N; i++){
        xPos[i]   = rawData[i*6+0];
        yPos[i]   = rawData[i*6+1];
        mass[i]   = rawData[i*6+2];
        xSpeed[i] = rawData[i*6+3];
        ySpeed[i] = rawData[i*6+4];
        bright[i] = rawData[i*6+5];
    }
    free(rawData);

    /* Empty arrays to store acceleration information */
    double* xAccelNew = (double*)malloc(N*sizeof(double));
    double* yAccelNew = (double*)malloc(N*sizeof(double));

    /* Graphics initialization */
    const float circleRadius = 0.0025, circleColor = 0;
    const int windowWidth = 800;
    float L = 1, W = 1;
    if (withGraphics){
        InitializeGraphics(argv[0],windowWidth,windowWidth);
        SetCAxes(0,1);
    }
    
    /* SIMULATION */
    for (int k=0; k<nsteps; k++){
        #pragma omp parallel for num_threads(nThreads)
        for (int i=0; i<N; i++){ //reset accelerations to 0
            xAccelNew[i] = 0.0;
            yAccelNew[i] = 0.0;
        }

        gettimeofday(&begin, NULL);
        calcAccel(N, G, mass, xPos, yPos, xAccelNew, yAccelNew, nThreads); //calculates accelerations at this iteration
        gettimeofday(&end, NULL);
        elapsed += (end.tv_sec - begin.tv_sec) + ((end.tv_usec - begin.tv_usec)/1000000.0);

        #pragma omp parallel for num_threads(nThreads)
        for (int i=0; i<N; i++){ //calculates new positions and speeds
            xSpeed[i] += xAccelNew[i]*delta_t;
            ySpeed[i] += yAccelNew[i]*delta_t;

            xPos[i] += xSpeed[i]*delta_t;
            yPos[i] += ySpeed[i]*delta_t;
        }

        if (withGraphics){ //graphics routine
            ClearScreen();
            #pragma omp parallel for num_threads(nThreads)
            for (int i=0; i<N; i++){
                DrawCircle(xPos[i], yPos[i], L, W, circleRadius, circleColor);
            }
            Refresh();
            usleep(1500);
        }

        if (nsteps>15 && (k+1)%((int)(nsteps/15)) == 0){ //progress info
            printf("Progress: %3.1f %%.\n", ((double)k/nsteps)*100);
        }
    }

    if (withGraphics){
        FlushDisplay();
        CloseDisplay();
    }

    /* preparing data for saving */
    double* output = (double*)malloc(nData*sizeof(double));
    #pragma omp parallel for num_threads(nThreads)
    for (int i=0; i<N; i++){
        output[i*6+0] = xPos[i];
        output[i*6+1] = yPos[i];
        output[i*6+2] = mass[i];
        output[i*6+3] = xSpeed[i];
        output[i*6+4] = ySpeed[i];
        output[i*6+5] = bright[i];
    }

    /* SAVING */
    int isWritten = write_doubles_to_file(nData, output, "result.gal");
    if (isWritten ==0){
        printf("Result successfully written to \"result.gal\"\n");
    }
    else {
        printf("Writing failed\n.");
    }

    /* freeing memory */
    free(xPos);
    free(yPos);
    free(xSpeed);
    free(ySpeed);
    free(mass);
    free(bright);
    free(xAccelNew);
    free(yAccelNew);
    free(output);

    gettimeofday(&endGlob, NULL);
    elapsedGlob = (endGlob.tv_sec - beginGlob.tv_sec) + ((endGlob.tv_usec - beginGlob.tv_usec)/1000000.0);
    printf("Progress: COMPLETE (calcAccel function runtime is %f seconds which is %.6f%% of the total runtime)\n", elapsed, elapsed/elapsedGlob*100);

    return 0;
}
