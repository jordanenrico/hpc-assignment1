#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "graphics.h"
#include "file_operations.h"

double distance(double xMajor, double yMajor, double xMinor, double yMinor){
    return sqrt(pow(xMajor-xMinor, 2)+pow(yMajor-yMinor,2));
}

int animate(char* name, double* xPos, double* yPos, int N, int nsteps){
    
    const float circleRadius=0.0025, circleColor=0;
    const int windowWidth=800;
    
    float L=1, W=1;
    int i = 0;
    
    InitializeGraphics(name,windowWidth,windowWidth);
    SetCAxes(0,1);
    printf("Hit q to quit.\n");

    while (!CheckForQuit()){
        i++;
        i %= nsteps+1;

        ClearScreen();
        int j;
        for (j=0; j<N; j++){
            DrawCircle(xPos[j*(nsteps+1)+i], yPos[j*(nsteps+1)+i], L, W, circleRadius, circleColor);
        }
        Refresh();

        usleep(2000);
    }
    FlushDisplay();
    CloseDisplay();
    return 0;
}

int main(int argc, char *argv[]){
    if(argc!=6){
        printf("ArgumentError: 5 input arguments are expected but %i received.\n", argc-1);
        return -1;
    }

    int i, j, k; // iterators

    int N = atoi(argv[1]);
    char filename[50];
    strcpy(filename, argv[2]);
    int nsteps = atoi(argv[3]);
    double delta_t = atof(argv[4]);
    int withGraphics = atoi(argv[5]);

    const double G = 100/N;
    const double epsilon = 0.001;

    int nData = N*6;
    double rawData[nData];
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
    double* xPos   = (double*)malloc(N*(nsteps+1)*sizeof(double));
    double* yPos   = (double*)malloc(N*(nsteps+1)*sizeof(double));
    double* xSpeed = (double*)malloc(N*(nsteps+1)*sizeof(double));
    double* ySpeed = (double*)malloc(N*(nsteps+1)*sizeof(double));
    double* mass   = (double*)malloc(N*sizeof(double));
    double* bright = (double*)malloc(N*sizeof(double));

    /* inserting initial data into arrays */
    for (i=0; i<N; i++){
        xPos[i*(nsteps+1)+0]   = rawData[i*6+0];
        yPos[i*(nsteps+1)+0]   = rawData[i*6+1];
        xSpeed[i*(nsteps+1)+0] = rawData[i*6+3];
        ySpeed[i*(nsteps+1)+0] = rawData[i*6+4];

        mass[i]   = rawData[i*6+2];
        bright[i] = rawData[i*6+5];
    }

    double xMassPerArea; //stores the x component of the sum of m/r^2
    double yMassPerArea; //stores the y component of the sum of m/r^2
    double majorMass;    //the mass of an object on which the gravity acts
    double majorXPos;    //the x position of majorMass
    double majorYPos;    //the y position of majorMass
    double minorMass;    //the mass of an object which affects majorMass
    double minorXPos;    //the x position of minorMass
    double minorYPos;    //the y position of minorMass
    double accelX;       //the x component of majorMass acceleration caused by gravity
    double accelY;       //the y component of majorMass acceleration caused by gravity
    double newXSpeed, newYSpeed, newXPos, newYPos;

    /* SIMULATION */
    for (k=0; k<nsteps; k++){
        for (i=0; i<N; i++){
            xMassPerArea = 0;
            yMassPerArea = 0;
            majorMass = mass[i];
            majorXPos = xPos[i*(nsteps+1)+k];
            majorYPos = yPos[i*(nsteps+1)+k];

            for (j=0; j<N; j++){
                if (i!=j){
                    minorMass = mass[j];
                    minorXPos = xPos[j*(nsteps+1)+k];
                    minorYPos = yPos[j*(nsteps+1)+k];

                    xMassPerArea += minorMass/pow(distance(majorXPos, majorYPos, minorXPos, minorYPos) + epsilon, 3)*(majorXPos-minorXPos);
                    yMassPerArea += minorMass/pow(distance(majorXPos, majorYPos, minorXPos, minorYPos) + epsilon, 3)*(majorYPos-minorYPos);
                }
            }

            accelX = -G*majorMass*xMassPerArea;
            accelY = -G*majorMass*yMassPerArea;

            newXSpeed = xSpeed[i*(nsteps+1)+k]+accelX*delta_t; //v'=v+at
            newYSpeed = ySpeed[i*(nsteps+1)+k]+accelY*delta_t;
            xSpeed[i*(nsteps+1)+k+1] = newXSpeed;
            ySpeed[i*(nsteps+1)+k+1] = newYSpeed;

            newXPos = xPos[i*(nsteps+1)+k]+newXSpeed*delta_t; //s'=s+vt
            newYPos = yPos[i*(nsteps+1)+k]+newYSpeed*delta_t;
            xPos[i*(nsteps+1)+k+1] = newXPos;
            yPos[i*(nsteps+1)+k+1] = newYPos;
        }
    }

    /* ANIMATION */
    if (withGraphics){
        animate(argv[0], xPos, yPos, N, nsteps);
    }

    /* PREPARING DATA TO BE SAVED */
    double* output = (double*)malloc(N*6*sizeof(double));
    for (i=0; i<N; i++){
        output[i*6+0] = xPos[i*(nsteps+1)+nsteps];
        output[i*6+1] = yPos[i*(nsteps+1)+nsteps];
        output[i*6+2] = mass[i];
        output[i*6+3] = xSpeed[i*(nsteps+1)+nsteps];
        output[i*6+4] = ySpeed[i*(nsteps+1)+nsteps];
        output[i*6+5] = bright[i];
    }

    /* SAVING */
    int isWritten = write_doubles_to_file(nData, output, "output.gal");
    if (isWritten ==0){
        printf("Output successfully written to \"output.gal\"\n.");
    }
    else{
        printf("Writing failed\n.");
    }

    free(xPos);
    free(yPos);
    free(xSpeed);
    free(ySpeed);
    free(mass);
    free(bright);
    free(output);

    return 0;
}