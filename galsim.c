#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "graphics.h"
#include "file_operations.h"

int calcAccel(int N, const double G, double* mass, double* xPos, double* yPos, double* xAccelNew, double* yAccelNew){
    /* Fills up xAccelNew and yAccelNew with magnitude and direction of the acceleration.
       xPotentSum = x component of the sum of mass_j*rVec/r^3 values experienced by planet i from planet j.
       yPotentSum = y component of the sum of mass_j*rVec/r^3 values experienced by planet i from planet j. */

    const double epsilon = 0.001;
    double magnitude, xPotentSum, yPotentSum;

    for (int i=0; i<N; i++){
        xPotentSum = 0;
        yPotentSum = 0;
        for (int j=0; j<N; j++){
            if (i != j){
                magnitude   = sqrt(pow(xPos[i] - xPos[j], 2) + pow(yPos[i] - yPos[j], 2)) + epsilon;
                xPotentSum += mass[j]*(xPos[i] - xPos[j])/pow(magnitude, 3);
                yPotentSum += mass[j]*(yPos[i] - yPos[j])/pow(magnitude, 3);
            }
        }
        xAccelNew[i] = -G*xPotentSum;
        yAccelNew[i] = -G*yPotentSum;
    }
    return 0;
}

int main(int argc, char *argv[]){
    if(argc!=6){
        printf("ArgumentError: 5 input arguments are expected but %i received.\n", argc-1);
        return -1;
    }

    int N            = atoi(argv[1]);
    char filename[50];
    strcpy(filename, argv[2]);
    int nsteps       = atoi(argv[3]);
    double delta_t   = atof(argv[4]);
    int withGraphics = atoi(argv[5]);
    const double G   = 100.0/N;
    int nData        = N*6;
    double* rawData  = (double*)malloc(nData*sizeof(double));

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

    /*  xPos[i*(nsteps+1)+j] is the x position of planet i at timestep j.
        We use (nsteps+1) to store nsteps data and the initial condition. */
    double* xPos   = (double*)malloc(N*sizeof(double));
    double* yPos   = (double*)malloc(N*sizeof(double));
    double* xSpeed = (double*)malloc(N*sizeof(double));
    double* ySpeed = (double*)malloc(N*sizeof(double));
    double* mass   = (double*)malloc(N*sizeof(double));
    double* bright = (double*)malloc(N*sizeof(double));

    /* inserting initial data into new arrays */
    for (int i=0; i<N; i++){
        xPos[i]   = rawData[i*6+0];
        yPos[i]   = rawData[i*6+1];
        mass[i]   = rawData[i*6+2];
        xSpeed[i] = rawData[i*6+3];
        ySpeed[i] = rawData[i*6+4];
        bright[i] = rawData[i*6+5];
    }

    /* empty arrays to store acceleration information */
    double* xAccelNew = (double*)malloc(N*sizeof(double));
    double* yAccelNew = (double*)malloc(N*sizeof(double));

    /* SIMULATION */
    const float circleRadius = 0.0025, circleColor = 0;
    const int windowWidth = 800;
    float L = 1, W = 1;
    int i = 0;
    if (withGraphics){
        InitializeGraphics(argv[0],windowWidth,windowWidth);
        SetCAxes(0,1);
    }
    
    for (int k=0; k<nsteps; k++){
        calcAccel(N, G, mass, xPos, yPos, xAccelNew, yAccelNew);
        for (int i=0; i<N; i++){

            xSpeed[i] += xAccelNew[i]*delta_t;
            ySpeed[i] += yAccelNew[i]*delta_t;

            xPos[i] += xSpeed[i]*delta_t;
            yPos[i] += ySpeed[i]*delta_t;
        }

        if (withGraphics){
            ClearScreen();
            for (int i=0; i<N; i++){
                DrawCircle(xPos[i], yPos[i], L, W, circleRadius, circleColor);
            }
            Refresh();
            usleep(1000);
        }

        if ((k+1)%((int)(nsteps/15)) == 0){
            printf("Progress: %3.1f %%.\n", ((double)k/nsteps)*100);
        }
    }
    if (withGraphics){
        FlushDisplay();
        CloseDisplay();
    }
    printf("Progress: COMPLETE\n");

    /* preparing data for saving */
    double* output = (double*)malloc(nData*sizeof(double));
    for (int i=0; i<N; i++){
        output[i*6+0] = xPos[i];
        output[i*6+1] = yPos[i];
        output[i*6+2] = mass[i];
        output[i*6+3] = xSpeed[i];
        output[i*6+4] = ySpeed[i];
        output[i*6+5] = bright[i];
    }

    /* SAVING */
    int isWritten = write_doubles_to_file(nData, output, "output.gal");
    if (isWritten ==0){
        printf("Output successfully written to \"output.gal\"\n.");
    }
    else {
        printf("Writing failed\n.");
    }

    /* freeing memory */
    free(rawData);
    free(xPos);
    free(yPos);
    free(xSpeed);
    free(ySpeed);
    free(mass);
    free(bright);
    free(xAccelNew);
    free(yAccelNew);
    free(output);

    return 0;
}
